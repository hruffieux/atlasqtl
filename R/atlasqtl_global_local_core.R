# This file is part of the `atlasqtl` R package:
#     https://github.com/hruffieux/atlasqtl
#
# Internal core function to call the variational algorithm for global-local 
# hotspot propensity modelling. 
# See help of `atlasqtl` function for details.
#
atlasqtl_global_local_core_ <- function(Y, X, shr_fac_inv, anneal, df, tol, 
                                        maxit, verbose, list_hyper, list_init, n_cores,
                                        checkpoint_path = NULL, 
                                        trace_path = NULL, full_output = FALSE, 
                                        thinned_elbo_eval = TRUE,   
                                        debug = FALSE, batch = "y") {
  
  n <- nrow(Y)
  p <- ncol(X)
  q <- ncol(Y)
  
  if (any(is.na(Y))) {
    
    mis_pat <- ifelse(is.na(Y), 0, 1)
    Y[is.na(Y)] <- 0
    X_norm_sq <- crossprod(X^2, mis_pat)
    
    cp_X_rm <- lapply(1:q, function(k) {
      if (any(mis_pat[,k] == 0)) {
        ind <- which(mis_pat[,k] == 0)
        crossprod(X[ind,, drop = FALSE])
      } else {
        matrix(0, nrow = p, ncol = p)
      }
    })
    
  } else {
    
    mis_pat <- X_norm_sq <- cp_X_rm <- NULL
    
  }
  
  Y_norm_sq <- colSums(Y^2) # must be after the if as some entries of y set to 0 when missing values
  cp_X <- crossprod(X)
  cp_Y_X <- crossprod(Y, X)
  
  if(verbose){
    cat(n_cores, " cores are used. \n")
  }

  if(n_cores!=1){
    RcppParallel::setThreadOptions(numThreads = n_cores)
  }
  
  # Gathering initial variational parameters. Do it explicitly.
  # with() function not used here as the objects will be modified later.
  #
  gam_vb <- list_init$gam_vb
  mu_beta_vb <- list_init$mu_beta_vb # mu_beta_vb[s, t] = variational posterior mean of beta_st | gamma_st = 1 
                                     # (beta_vb[s, t] = varational posterior mean of beta_st)
  sig02_inv_vb <- list_init$sig02_inv_vb  # horseshoe global precision parameter
  sig2_beta_vb <- list_init$sig2_beta_vb  # sig2_beta_vb[t]  = variational variance of beta_st | gamma_st = 1 
                                          # (same for all s as X is standardized)
  sig2_theta_vb <- list_init$sig2_theta_vb 
  tau_vb <-list_init$tau_vb
  theta_vb <- list_init$theta_vb
  zeta_vb <- list_init$zeta_vb
  
  rm(list_init)
  
  theta_plus_zeta_vb <- sweep(tcrossprod(theta_vb, rep(1, q)), 2, zeta_vb, `+`)
  log_Phi_theta_plus_zeta <- pnorm(theta_plus_zeta_vb, log.p = TRUE)
  log_1_min_Phi_theta_plus_zeta <- pnorm(theta_plus_zeta_vb, log.p = TRUE, lower.tail = FALSE)
  
  # Preparing trace saving if any
  #
  trace_ind_max <- trace_var_max <- NULL
  
  # Preparing annealing if any
  #
  anneal_scale <- TRUE
  
  if (is.null(anneal)) {
    annealing <- FALSE
    c <- c_s <- 1 # c_s for scale parameters
    it_init <- 1 # first non-annealed iteration 
  } else {
    annealing <- TRUE
    ladder <- get_annealing_ladder_(anneal, verbose)
    c <- ladder[1]
    c_s <- ifelse(anneal_scale, c, 1)
    it_init <- anneal[3] # first non-annealed iteration 
  }
  
  eps <- .Machine$double.eps^0.5
  
  if (thinned_elbo_eval) {
    times_conv_sched <- c(1, 5, 10, 50) 
    batch_conv_sched <- c(1, 10, 25, 50) 
  } else {
    times_conv_sched <- 1
    batch_conv_sched <- 1
  }

  
  ind_batch_conv <- length(batch_conv_sched) + 1 # so that, the first time, it enters in the loop below 
  batch_conv <- 1 
  
  with(list_hyper, { # list_init not used with the with() function to avoid
                     # copy-on-write for large objects
    
    # Response-specific parameters: objects derived from t02
    #
    t02_inv <- 1 / t02
    sig2_zeta_vb <- update_sig2_c0_vb_(p, t02, c = c) # stands for a diagonal matrix of size d with this value on the (constant) diagonal
    
    vec_sum_log_det_zeta <- - q * (log(t02) + log(p + t02_inv))
    
    
    # Stored/precomputed objects
    #
    beta_vb <- update_beta_vb_(gam_vb, mu_beta_vb)
    m2_beta <- update_m2_beta_(gam_vb, mu_beta_vb, sig2_beta_vb, sweep = TRUE) # first time keep sweep = TRUE even when missing data, since uses the initial parameter sig2_beta_vb which is a vector.
    
    cp_X_Xbeta <- update_cp_X_Xbeta_(cp_X, beta_vb, cp_X_rm)
    
    # Fixed VB parameter
    #
    nu_xi_inv_vb <- 1 # no change with annealing 
    
    converged <- FALSE
    lb_new <- -Inf
    it <- 0
    
    while ((!converged) & (it < maxit)) {
      
      lb_old <- lb_new
      it <- it + 1
      
      if (verbose != 0 &  (it == 1 | it %% max(5, batch_conv) == 0)) 
        cat(paste0("Iteration ", format(it), "... \n"))
      
      # % #
      nu_vb <- update_nu_vb_(nu, sum(gam_vb), c = c)
      rho_vb <- update_rho_vb_(rho, m2_beta, tau_vb, c = c)
      
      sig2_inv_vb <- nu_vb / rho_vb
      # % #
      
      # % #
      eta_vb <- update_eta_vb_(n, eta, gam_vb, mis_pat, c = c)
      kappa_vb <- update_kappa_vb_(n, Y_norm_sq, cp_Y_X, cp_X_Xbeta, kappa, 
                                   beta_vb, m2_beta, sig2_inv_vb, X_norm_sq, c = c)
      
      tau_vb <- eta_vb / kappa_vb
      
      sig2_beta_vb <- update_sig2_beta_vb_(n, sig2_inv_vb, tau_vb, X_norm_sq, c = c)
      
      log_tau_vb <- update_log_tau_vb_(eta_vb, kappa_vb)
      log_sig2_inv_vb <- update_log_sig2_inv_vb_(nu_vb, rho_vb)
      
      # different possible batch-coordinate ascent schemes:
      #
      if (batch == "y") { # optimal scheme
        
        # C++ Eigen call for expensive updates 
        #
        # Shuffle updates order at each iteration
        #
        # shuffled_ind <- as.integer(sample(0:(p-1))) # Zero-based index in C++
        # sample_q <- as.integer(sample(0:(q-1))) # Zero-based index in C++
        shuffled_ind <- as.integer(0:(p-1)) # Zero-based index in C++
        sample_q <- as.integer(0:(q-1)) # Zero-based index in C++
        
        if(n_cores == 1){
          if (is.null(mis_pat)) {
            coreDualLoop(cp_X, cp_Y_X, gam_vb, log_Phi_theta_plus_zeta,
                         log_1_min_Phi_theta_plus_zeta, log_sig2_inv_vb, log_tau_vb,
                         beta_vb, cp_X_Xbeta, mu_beta_vb, sig2_beta_vb, tau_vb,
                         shuffled_ind, sample_q = sample_q, c = c)
          } else {
            coreDualMisLoop(cp_X, cp_X_rm, cp_Y_X, gam_vb, log_Phi_theta_plus_zeta, 
                            log_1_min_Phi_theta_plus_zeta, log_sig2_inv_vb, log_tau_vb, 
                            beta_vb, cp_X_Xbeta, mu_beta_vb, sig2_beta_vb, tau_vb, 
                            shuffled_ind, sample_q = sample_q, c = c)
          }
        }else{
          if (is.null(mis_pat)) {
            coreDualLoopParResp(cp_X, cp_Y_X, gam_vb, log_Phi_theta_plus_zeta,
                                log_1_min_Phi_theta_plus_zeta, log_sig2_inv_vb, log_tau_vb,
                                beta_vb, cp_X_Xbeta, mu_beta_vb, sig2_beta_vb, tau_vb,
                                shuffled_ind, sample_q = sample_q, c = c)
          }else{
            # not working 
            coreDualMisLoopParResp(cp_X, cp_X_rm, cp_Y_X, gam_vb, log_Phi_theta_plus_zeta, 
                                   log_1_min_Phi_theta_plus_zeta, log_sig2_inv_vb, log_tau_vb, 
                                   beta_vb, cp_X_Xbeta, mu_beta_vb, sig2_beta_vb, tau_vb, 
                                   shuffled_ind, sample_q = sample_q, c = c)
          }
        }
        
        
        
      } else if (batch == "0"){ # no batch, used only internally (slower)
        # schemes "x" of "x-y" are not batch concave
        # hence not implemented as they may diverge
        
        for (k in sample(1:q)) {
          
          
          if (is.null(mis_pat)) {
            
            for (j in sample(1:p)) {
              
              cp_X_Xbeta[,k] <- cp_X_Xbeta[,k] - beta_vb[j, k] * cp_X[j, ]
              
              mu_beta_vb[j, k] <- c * sig2_beta_vb[k] * tau_vb[k] * (cp_Y_X[k, j] - cp_X_Xbeta[j, k])
              
              gam_vb[j, k] <- exp(-log_one_plus_exp_(c * (pnorm(theta_vb[j] + zeta_vb[k], lower.tail = FALSE, log.p = TRUE) -
                                                            pnorm(theta_vb[j] + zeta_vb[k], log.p = TRUE) -
                                                            log_tau_vb[k] / 2 - log_sig2_inv_vb / 2 -
                                                            mu_beta_vb[j, k] ^ 2 / (2 * sig2_beta_vb[k]) -
                                                            log(sig2_beta_vb[k]) / 2)))
              
              beta_vb[j, k] <- gam_vb[j, k] * mu_beta_vb[j, k]
              
              cp_X_Xbeta[,k] <- cp_X_Xbeta[,k] + beta_vb[j, k] * cp_X[j, ]
              
            }
            
          } else {
            
            for (j in sample(1:p)) {

              cp_X_Xbeta[,k] <- cp_X_Xbeta[,k] - beta_vb[j, k] * (cp_X[j, ] - cp_X_rm[[k]][j, ])

              mu_beta_vb[j, k] <- c * sig2_beta_vb[j, k] * tau_vb[k] * (cp_Y_X[k, j] - cp_X_Xbeta[j, k])

              gam_vb[j, k] <- exp(-log_one_plus_exp_(c * (pnorm(theta_vb[j] + zeta_vb[k], lower.tail = FALSE, log.p = TRUE) -
                                                            pnorm(theta_vb[j] + zeta_vb[k], log.p = TRUE) -
                                                            log_tau_vb[k] / 2 - log_sig2_inv_vb / 2 -
                                                            mu_beta_vb[j, k] ^ 2 / (2 * sig2_beta_vb[j, k]) -
                                                            log(sig2_beta_vb[j, k]) / 2)))

              beta_vb[j, k] <- gam_vb[j, k] * mu_beta_vb[j, k]

              cp_X_Xbeta[,k] <- cp_X_Xbeta[,k] + beta_vb[j, k] * (cp_X[j, ] - cp_X_rm[[k]][j, ])

            }
            
          }
        }
        
      } else {
        
        stop ("Batch scheme not defined. Exit.")
        
      }
      
      m2_beta <- update_m2_beta_(gam_vb, mu_beta_vb, sig2_beta_vb, mis_pat = mis_pat)
      
      Z <- update_Z_(gam_vb, theta_plus_zeta_vb, log_1_min_Phi_theta_plus_zeta, log_Phi_theta_plus_zeta, c = c) 
      
      # keep this order!
      #  
      L_vb <- c_s * sig02_inv_vb * shr_fac_inv * (theta_vb^2 + sig2_theta_vb - 2 * theta_vb * m0 + m0^2) / 2 / df 
      rho_xi_inv_vb <- c_s * (A2_inv + sig02_inv_vb) 
      
      if (annealing & anneal_scale) {
        
        lam2_inv_vb <- update_annealed_lam2_inv_vb_(L_vb, c_s, df)
        
      } else {
        
        Q_app <- Q_approx_vec(L_vb)
        
        if (df == 1) {
          
          lam2_inv_vb <- 1 / (Q_app * L_vb) - 1
          
        } else if (df == 3) {
          
          lam2_inv_vb <- exp(-log(3) - log(L_vb) + log(1 - L_vb * Q_app) - log(Q_app * (1 + L_vb) - 1)) - 1 / 3
          
        } else {
          # also works for df = 3 but might be slightly less efficient than the above
          
          exponent <- (df + 1) / 2
          
          lam2_inv_vb <- sapply(1:p, function(j) {
            
            exp(log(compute_integral_hs_(df, L_vb[j] * df, m = exponent, n = exponent, Q_ab = Q_app[j])) -
                  log(compute_integral_hs_(df, L_vb[j] * df, m = exponent, n = exponent - 1, Q_ab = Q_app[j])))
            
          })
          
        }
        
      }
      
      xi_inv_vb <- nu_xi_inv_vb / rho_xi_inv_vb
      
      sig2_theta_vb <- update_sig2_c0_vb_(q, 1 / (sig02_inv_vb * lam2_inv_vb * shr_fac_inv), c = c)
      
      theta_vb <- update_theta_vb_(Z, m0, sig02_inv_vb * lam2_inv_vb * shr_fac_inv, sig2_theta_vb,
                                   vec_fac_st = NULL, zeta_vb, is_mat = FALSE, c = c)
      
      nu_s0_vb <- update_nu_vb_(1 / 2, p, c = c_s)
      
      rho_s0_vb <- c_s * (xi_inv_vb + 
                            sum(lam2_inv_vb * shr_fac_inv * (theta_vb^2 + sig2_theta_vb - 2 * theta_vb * m0 + m0^2)) / 2) 
      
      sig02_inv_vb <- as.numeric(nu_s0_vb / rho_s0_vb)
      
      zeta_vb <- update_zeta_vb_(Z, theta_vb, n0, sig2_zeta_vb, t02_inv,
                                 is_mat = FALSE, c = c) 
      
      theta_plus_zeta_vb <- sweep(tcrossprod(theta_vb, rep(1, q)), 2, zeta_vb, `+`)
      log_Phi_theta_plus_zeta <- pnorm(theta_plus_zeta_vb, log.p = TRUE)
      log_1_min_Phi_theta_plus_zeta <- pnorm(theta_plus_zeta_vb, log.p = TRUE, lower.tail = FALSE)  
      
      if (verbose == 2 && (it == 1 | it %% max(5, batch_conv) == 0)) {
        
        cat(paste0("Variational hotspot propensity global scale: ", 
                   format(sqrt(rho_s0_vb / (nu_s0_vb - 1) / shr_fac_inv), digits = 3), ".\n"))
        cat("Approximate variational hotspot propensity local scale: \n")
        print(summary(sqrt(1 / lam2_inv_vb)))
        cat("\n")
        
      }
      
      if (!is.null(trace_path) && (it == 1 | it %% 25 == 0)) {
        
        list_traces <- plot_trace_var_hs_(lam2_inv_vb, sig02_inv_vb, q, it, 
                                          trace_ind_max, trace_var_max, trace_path)
        
        trace_ind_max <- list_traces$trace_ind_max
        trace_var_max <- list_traces$trace_var_max
        
      }
      
      
      if (annealing) {
        
        if (verbose != 0 & (it == 1 | it %% 5 == 0))
          cat(paste0("Temperature = ", format(1 / c, digits = 4), "\n\n"))
        
        sig2_zeta_vb <- c * sig2_zeta_vb
        
        c <- ifelse(it < length(ladder), ladder[it + 1], 1)
        c_s <- ifelse(anneal_scale, c, 1)
        
        sig2_zeta_vb <- sig2_zeta_vb / c
        
        if (isTRUE(all.equal(c, 1))) {
          
          annealing <- FALSE
          
          if (verbose != 0)
            cat("** Exiting annealing mode. **\n\n")
        }
        
        
      } else {
        
        
        if (it <= it_init + 1 | it %% batch_conv == 0 | it %% batch_conv == 1) { 
          # it <= it_init + 1 evaluate the ELBO for the first two non-annealed iterations
          # to (also) evaluate convergence between two consecutive iterations
          
          lb_new <- elbo_global_local_(Y, A2_inv, beta_vb, df, eta, eta_vb, gam_vb,
                                       kappa, kappa_vb, L_vb, lam2_inv_vb, 
                                       log_1_min_Phi_theta_plus_zeta, log_Phi_theta_plus_zeta, 
                                       m0, m2_beta, n0, nu, nu_s0_vb, nu_vb, nu_xi_inv_vb,
                                       Q_app, rho, rho_s0_vb, rho_vb, rho_xi_inv_vb,
                                       shr_fac_inv, sig02_inv_vb, sig2_beta_vb,
                                       sig2_inv_vb,  sig2_theta_vb, sig2_zeta_vb,
                                       t02_inv, tau_vb, theta_vb, vec_sum_log_det_zeta,
                                       xi_inv_vb, zeta_vb, X_norm_sq, Y_norm_sq, cp_Y_X, cp_X_Xbeta, mis_pat)
          
          if (verbose != 0 & (it == it_init | it %% max(5, batch_conv) == 0))
            cat(paste0("ELBO = ", format(lb_new), "\n\n"))
          
          if (debug && lb_new + eps < lb_old)
            stop("ELBO not increasing monotonically. Exit. ")
          
          diff_lb <- abs(lb_new - lb_old)
          
          sum_exceed <- sum(diff_lb > (times_conv_sched * tol))
          
          if (sum_exceed == 0) {
            
            converged <- TRUE
            
          } else if (ind_batch_conv > sum_exceed) {
            
            ind_batch_conv <- sum_exceed
            batch_conv <- batch_conv_sched[ind_batch_conv]
            
          }
          
        }
        
        checkpoint_(it, checkpoint_path, beta_vb, gam_vb, theta_vb, zeta_vb, 
                    converged, lb_new, lb_old,
                    lam2_inv_vb = lam2_inv_vb, sig02_inv_vb = sig02_inv_vb,
                    names_x = colnames(X), names_y = colnames(Y))
        
      }
      
    }
    
    checkpoint_clean_up_(checkpoint_path)
    
    
    if (verbose != 0) {
      if (converged) {
        cat(paste0("Convergence obtained after ", format(it), " iterations. \n",
                   "Optimal marginal log-likelihood variational lower bound ",
                   "(ELBO) = ", format(lb_new), ". \n\n"))
      } else {
        warning("Maximal number of iterations reached before convergence. Exit.")
      }
    }
    
    lb_opt <- lb_new
    s02_vb <- rho_s0_vb / (nu_s0_vb - 1) / shr_fac_inv # inverse gamma mean, with embedded shrinkage
    
    if (full_output) { # for internal use only
      
      create_named_list_(beta_vb, eta_vb, gam_vb, kappa_vb, lam2_inv_vb, 
                         nu_s0_vb, nu_vb, nu_xi_inv_vb, rho_s0_vb, rho_vb, 
                         rho_xi_inv_vb, shr_fac_inv, sig02_inv_vb, sig2_beta_vb, 
                         sig2_inv_vb,  sig2_theta_vb, sig2_zeta_vb, tau_vb, 
                         theta_vb, cp_Y_X, cp_X, cp_X_Xbeta, xi_inv_vb, zeta_vb)
      
    } else {
      
      names_x <- colnames(X)
      names_y <- colnames(Y)
      names_n <- rownames(Y)
      
      rownames(gam_vb) <- rownames(beta_vb) <- names_x
      colnames(gam_vb) <- colnames(beta_vb) <- names_y
      names(theta_vb) <- names_x
      names(zeta_vb) <- names_y
      names(lam2_inv_vb) <- names_x
      
      diff_lb <- abs(lb_opt - lb_old)
      
      create_named_list_(beta_vb, gam_vb, theta_vb, zeta_vb, 
                         n, p, q, anneal, converged, it, maxit, tol, lb_opt, 
                         diff_lb)
      
    }
  })
  
}



# Internal function which implements the marginal log-likelihood variational
# lower bound (ELBO) corresponding to the `atlasqtl_struct_core` algorithm.
#
elbo_global_local_ <- function(Y, A2_inv, beta_vb, df, eta, eta_vb, gam_vb, 
                               kappa, kappa_vb, L_vb, lam2_inv_vb, 
                               log_1_min_Phi_theta_plus_zeta, log_Phi_theta_plus_zeta, m0, m2_beta, 
                               n0, nu, nu_s0_vb, nu_vb, nu_xi_inv_vb, 
                               Q_app, rho, rho_s0_vb, rho_vb, rho_xi_inv_vb, 
                               shr_fac_inv, sig02_inv_vb, sig2_beta_vb, 
                               sig2_inv_vb,  sig2_theta_vb, sig2_zeta_vb, 
                               t02_inv, tau_vb, theta_vb, vec_sum_log_det_zeta, 
                               xi_inv_vb, zeta_vb, X_norm_sq, Y_norm_sq, cp_Y_X, 
                               cp_X_Xbeta, mis_pat) {
  
  n <- nrow(Y)
  p <- length(L_vb)
  
  # needed for monotonically increasing elbo.
  #
  eta_vb <- update_eta_vb_(n, eta, gam_vb, mis_pat)
  kappa_vb <- update_kappa_vb_(n, Y_norm_sq, cp_Y_X, cp_X_Xbeta, kappa, beta_vb, 
                               m2_beta, sig2_inv_vb, X_norm_sq)
  
  nu_vb <- update_nu_vb_(nu, sum(gam_vb))
  rho_vb <- update_rho_vb_(rho, m2_beta, tau_vb)
  
  log_tau_vb <- update_log_tau_vb_(eta_vb, kappa_vb)
  log_sig2_inv_vb <- update_log_sig2_inv_vb_(nu_vb, rho_vb)
  
  log_sig02_inv_vb <- update_log_sig2_inv_vb_(nu_s0_vb, rho_s0_vb)
  log_xi_inv_vb <- update_log_sig2_inv_vb_(nu_xi_inv_vb, rho_xi_inv_vb)
  
  
  elbo_A <- e_y_(n, kappa, kappa_vb, log_tau_vb, m2_beta, sig2_inv_vb, tau_vb, mis_pat)
  
  elbo_B <- e_beta_gamma_(gam_vb, log_1_min_Phi_theta_plus_zeta, log_Phi_theta_plus_zeta, log_sig2_inv_vb, 
                          log_tau_vb, zeta_vb, 
                          theta_vb, m2_beta, sig2_beta_vb, sig2_zeta_vb,
                          sig2_theta_vb, sig2_inv_vb, tau_vb)
  
  elbo_C <- e_theta_hs_(lam2_inv_vb, L_vb, log_sig02_inv_vb + log(shr_fac_inv), 
                        m0, theta_vb, Q_app, sig02_inv_vb * shr_fac_inv, 
                        sig2_theta_vb, df)
  
  elbo_D <- e_zeta_(zeta_vb, n0, sig2_zeta_vb, t02_inv, vec_sum_log_det_zeta)
  
  elbo_E <- e_tau_(eta, eta_vb, kappa, kappa_vb, log_tau_vb, tau_vb)
  
  elbo_F <- e_sig2_inv_hs_(xi_inv_vb, nu_s0_vb, log_xi_inv_vb, log_sig02_inv_vb, 
                           rho_s0_vb, sig02_inv_vb) # sig02_inv_vb
  
  elbo_G <- e_sig2_inv_(1 / 2, nu_xi_inv_vb, log_xi_inv_vb, A2_inv, 
                        rho_xi_inv_vb, xi_inv_vb) # xi_inv_vb
  
  elbo_H <- e_sig2_inv_(nu, nu_vb, log_sig2_inv_vb, rho, rho_vb, sig2_inv_vb)
  
  as.numeric(elbo_A + elbo_B + elbo_C + elbo_D + elbo_E + elbo_F + elbo_G + elbo_H)
  
}

