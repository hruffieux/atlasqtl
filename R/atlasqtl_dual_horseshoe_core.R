# This file is part of the `atlasqtl` R package:
#     https://github.com/hruffieux/atlasqtl
#
# Internal core function to call the variational algorithm for hotspot 
# propensity control. Sparse regression with identity link, no fixed covariates.
# See help of `atlasqtl` function for details.
#
atlasqtl_horseshoe_core_ <- function(Y, X, anneal, df, tol, maxit, verbose, 
                                     list_hyper, list_init, 
                                     checkpoint_path = NULL, trace_path = NULL,
                                     full_output = FALSE, debug = FALSE, 
                                     batch = "y") {

  # Y must have been centered, and X standardized.
  #
  n <- nrow(Y)
  p <- ncol(X)
  q <- ncol(Y)
  
  
  # Gathering initial variational parameters
  #
  gam_vb <- list_init$gam_vb
  mu_beta_vb <- list_init$mu_beta_vb # mu_beta_vb[s, t] = variational posterior mean of beta_st | gamma_st = 1 
                                     # (beta_vb[s, t] = varational posterior mean of beta_st)
  sig2_beta_vb <- list_init$sig2_beta_vb  # sig2_beta_vb[t]  = variational variance of beta_st | gamma_st = 1 
                                          # (same for all s as X is standardized)
  tau_vb <-list_init$tau_vb
  
  rm(list_init)
  
  
  # Preparing trace saving if any
  #
  trace_ind_max <- trace_var_max <- NULL
  
  with(list_hyper, { # list_init not used with the with() function to avoid
                     # copy-on-write for large objects
    
    shr_fac_inv <- q # = 1 / shrinkage_factor for global variance
    
    # Preparing annealing if any
    #
    anneal_scale <- TRUE
    
    if (is.null(anneal)) {
      annealing <- FALSE
      c <- c_s <- 1 # c_s for scale parameters
    } else {
      annealing <- TRUE
      ladder <- get_annealing_ladder_(anneal, verbose)
      c <- ladder[1]
      c_s <- ifelse(anneal_scale, c, 1)
    }
    
    eps <- .Machine$double.eps^0.5
    
    # Variance initialization
    #
    sig02_inv_vb <- rgamma(1, shape = max(p, q), rate = 1) 
    
    # Some hyperparameters
    #
    A2_inv <- 1 #  hyperparameter 
    
    # Choose m0 so that, `a priori' (i.e. before optimization), E_p_gam is as specified by the user. 
    # In fact, we assume that the variance of theta (s0^2 in the hyperparameter doc) 
    # is very small so that the shift is negligeable: we set m0 to 0.
    #
    m0 <- rep(0, p)
    
    # Parameter initialization here for the top level 
    #
    theta_vb <- rnorm(p, sd = 1 / sqrt(sig02_inv_vb * shr_fac_inv)) 
    sig2_theta_vb <- 1 / (q + rgamma(p, shape = sig02_inv_vb * shr_fac_inv, rate = 1)) # initial guess assuming lam2_inv_vb = 1
    
    zeta_vb <- rnorm(q, mean = n0, sd = sqrt(t02))
    
    
    # Response-specific parameters: objects derived from t02
    #
    T0_inv <- 1 / t02
    sig2_zeta_vb <- update_sig2_c0_vb_(p, t02, c = c) # stands for a diagonal matrix of size d with this value on the (constant) diagonal
    
    vec_sum_log_det_zeta <- - q * (log(t02) + log(p + T0_inv))
    
    
    # Stored/precomputed objects
    #
    beta_vb <- update_beta_vb_(gam_vb, mu_beta_vb)
    m2_beta <- update_m2_beta_(gam_vb, mu_beta_vb, sig2_beta_vb, sweep = TRUE)
    
    mat_x_m1 <- update_mat_x_m1_(X, beta_vb)
    
    # Fixed VB parameter
    #
    nu_xi_inv_vb <- 1 # no change with annealing 
    
    converged <- FALSE
    lb_new <- -Inf
    it <- 0
    
    while ((!converged) & (it < maxit)) {
      
      lb_old <- lb_new
      it <- it + 1
      
      if (verbose & (it == 1 | it %% 5 == 0))
        cat(paste("Iteration ", format(it), "... \n", sep = ""))
      
      # % #
      nu_vb <- update_nu_vb_(nu, sum(gam_vb), c = c)
      rho_vb <- update_rho_vb_(rho, m2_beta, tau_vb, c = c)
      
      sig2_inv_vb <- nu_vb / rho_vb
      # % #
      
      # % #
      eta_vb <- update_eta_vb_(n, eta, gam_vb, c = c)
      kappa_vb <- update_kappa_vb_(Y, kappa, mat_x_m1, beta_vb, m2_beta, sig2_inv_vb, c = c)
      
      tau_vb <- eta_vb / kappa_vb
      # % #
      
      sig2_beta_vb <- update_sig2_beta_vb_(n, sig2_inv_vb, tau_vb, c = c)
      
      log_tau_vb <- update_log_tau_vb_(eta_vb, kappa_vb)
      log_sig2_inv_vb <- update_log_sig2_inv_vb_(nu_vb, rho_vb)
      
      
      # different possible batch-coordinate ascent schemes:
      #
      if (batch == "y") { # optimal scheme
        
        log_Phi_theta_plus_zeta <- sapply(zeta_vb, function(zeta_k) {
          pnorm(theta_vb + zeta_k, log.p = TRUE)})
        
        log_1_min_Phi_theta_plus_zeta <- sapply(zeta_vb, function(zeta_k) {
          pnorm(theta_vb + zeta_k, lower.tail = FALSE, log.p = TRUE)})
        
        # C++ Eigen call for expensive updates
        shuffled_ind <- as.numeric(sample(0:(p-1))) # Zero-based index in C++
        
        coreDualLoop(X, Y, gam_vb, log_Phi_theta_plus_zeta,
                     log_1_min_Phi_theta_plus_zeta, log_sig2_inv_vb,
                     log_tau_vb, beta_vb, mat_x_m1, mu_beta_vb,
                     sig2_beta_vb, tau_vb, shuffled_ind, c = c)
        
      } else if (batch == "0"){ # no batch, used only internally
        # schemes "x" of "x-y" are not batch concave
        # hence not implemented as they may diverge
        
        for (k in sample(1:q)) {
          
          for (j in sample(1:p)) {
            
            mat_x_m1[, k] <- mat_x_m1[, k] - X[, j] * beta_vb[j, k]
            
            mu_beta_vb[j, k] <- c * sig2_beta_vb[k] * tau_vb[k] * crossprod(Y[, k] - mat_x_m1[, k], X[, j])
            
            gam_vb[j, k] <- exp(-log_one_plus_exp_(c * (pnorm(theta_vb[j] + zeta_vb[k], lower.tail = FALSE, log.p = TRUE) -
                                                          pnorm(theta_vb[j] + zeta_vb[k], log.p = TRUE) -
                                                          log_tau_vb[k] / 2 - log_sig2_inv_vb / 2 -
                                                          mu_beta_vb[j, k] ^ 2 / (2 * sig2_beta_vb[k]) -
                                                          log(sig2_beta_vb[k]) / 2)))
            
            beta_vb[j, k] <- gam_vb[j, k] * mu_beta_vb[j, k]
            
            mat_x_m1[, k] <- mat_x_m1[, k] + X[, j] * beta_vb[j, k]
            
          }
        }
        
      } else {
        
        stop ("Batch scheme not defined. Exit.")
        
      }
      
      
      m2_beta <- update_m2_beta_(gam_vb, mu_beta_vb, sig2_beta_vb, sweep = TRUE)
      
      Z <- update_Z_(gam_vb, sweep(tcrossprod(theta_vb, rep(1, q)), 2,
                                        zeta_vb, `+`), c = c) # we use info_ so that the second argument is a matrix
      
      # keep this order!
      #  
      L_vb <- c_s * sig02_inv_vb * shr_fac_inv * (theta_vb^2 + sig2_theta_vb - 2 * theta_vb * m0 + m0^2) / 2 / df 
      rho_xi_inv_vb <- c_s * (A2_inv + sig02_inv_vb) 
        
      
      if (annealing & anneal_scale) {
        
        lam2_inv_vb <- update_annealed_lam2_inv_vb_(L_vb, c_s, df)
        
      } else {
        
        if (df == 1) {
          
          Q_app <- sapply(L_vb, function(L_vb_s) Q_approx(L_vb_s))  # TODO implement a Q_approx for vectors
          
          lam2_inv_vb <- 1 / (Q_app * L_vb) - 1
          
        } else if (df == 3) {
          
          Q_app <- sapply(L_vb, function(L_vb_s) Q_approx(L_vb_s))
          
          lam2_inv_vb <- exp(-log(3) - log(L_vb) + log(1 - L_vb * Q_app) - log(Q_app * (1 + L_vb) - 1)) - 1 / 3
          
        } else {
          # also works for df = 3 but might be slightly less efficient than the above
          
          Q_app <- sapply(L_vb, function(L_vb_s) Q_approx(L_vb_s))
          
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
      
      zeta_vb <- update_zeta_vb_(Z, theta_vb, n0, sig2_zeta_vb, T0_inv,
                                     is_mat = FALSE, c = c) 
      
      
      if (verbose && (it == 1 | it %% 5 == 0)) {
        
        cat(paste0("Updated global variance: ", format(rho_s0_vb / (nu_s0_vb - 1) / shr_fac_inv, digits = 4), ".\n"))
        cat("Updated local variational parameter 1 / lam2_inv_vb for local variances: \n")
        print(summary(1 / lam2_inv_vb))
        cat("\n")
      
      }
      
      if (!is.null(trace_path) && (it == 1 | it %% 25 == 0)) {
        
        list_traces <- plot_trace_var_hs_(lam2_inv_vb, sig02_inv_vb, q, it, trace_ind_max, trace_var_max, trace_path)
        trace_ind_max <- list_traces$trace_ind_max
        trace_var_max <- list_traces$trace_var_max
        
      }
      
      
      if (annealing) {
        
        if (verbose & (it == 1 | it %% 5 == 0))
          cat(paste("Temperature = ", format(1 / c, digits = 4), "\n\n", sep = ""))
        
        sig2_zeta_vb <- c * sig2_zeta_vb
        
        c <- ifelse(it < length(ladder), ladder[it + 1], 1)
        c_s <- ifelse(anneal_scale, c, 1)
        
        sig2_zeta_vb <- sig2_zeta_vb / c
        
        if (isTRUE(all.equal(c, 1))) {
          
          annealing <- FALSE
          
          if (verbose)
            cat("** Exiting annealing mode. **\n\n")
        }
        
        
      } else {
        
        lb_new <- elbo_horseshoe_(Y, xi_inv_vb, A2_inv, lam2_inv_vb, eta, eta_vb, L_vb, gam_vb, kappa, kappa_vb, nu,
                                       nu_vb, nu_xi_inv_vb, nu_s0_vb, m0, n0, zeta_vb,
                                       theta_vb, rho, rho_vb, rho_xi_inv_vb, rho_s0_vb, Q_app, sig2_beta_vb,
                                       sig02_inv_vb, sig2_theta_vb, sig2_inv_vb, sig2_zeta_vb,
                                       T0_inv, tau_vb, beta_vb, m2_beta, mat_x_m1,
                                       vec_sum_log_det_zeta, df, shr_fac_inv)
        
        if (verbose & (it == 1 | it %% 5 == 0))
          cat(paste("ELBO = ", format(lb_new), "\n\n", sep = ""))
        
        if (debug && lb_new + eps < lb_old)
          stop("ELBO not increasing monotonically. Exit. ")
        
        converged <- (abs(lb_new - lb_old) < tol)
        
        
        checkpoint_(it, checkpoint_path, gam_vb, converged, lb_new, lb_old, 
                    lam2_inv_vb = lam2_inv_vb, zeta_vb = zeta_vb, theta_vb = theta_vb, 
                    sig02_inv_vb = sig02_inv_vb)
      }
      
    }
    
    checkpoint_clean_up_(checkpoint_path)
    
    
    if (verbose) {
      if (converged) {
        cat(paste("Convergence obtained after ", format(it), " iterations. \n",
                  "Optimal marginal log-likelihood variational lower bound ",
                  "(ELBO) = ", format(lb_new), ". \n\n", sep = ""))
      } else {
        warning("Maximal number of iterations reached before convergence. Exit.")
      }
    }
    
    lb_opt <- lb_new
    s02_vb <- rho_s0_vb / (nu_s0_vb - 1) / shr_fac_inv # inverse gamma mean, with embedded shrinkage
    
    if (full_output) { # for internal use only
      
      create_named_list_(xi_inv_vb, A2_inv, lam2_inv_vb, eta, eta_vb, L_vb, gam_vb, kappa, kappa_vb, nu,
                         nu_vb, nu_xi_inv_vb, nu_s0_vb, m0, n0, zeta_vb,
                         theta_vb, rho, rho_vb, rho_xi_inv_vb, rho_s0_vb, Q_app, sig2_beta_vb,
                         sig02_inv_vb, s02_vb, sig2_theta_vb, sig2_inv_vb, sig2_zeta_vb,
                         T0_inv, tau_vb, beta_vb, m2_beta, mat_x_m1,
                         vec_sum_log_det_zeta, df, shr_fac_inv)
      
    } else {
      
      names_x <- colnames(X)
      names_y <- colnames(Y)
      
      rownames(gam_vb) <- names_x
      colnames(gam_vb) <- names_y
      names(theta_vb) <- names_x
      names(zeta_vb) <- names_y
      names(lam2_inv_vb) <- names_x
      
      diff_lb <- abs(lb_opt - lb_old)
      
      create_named_list_(beta_vb, gam_vb, theta_vb, zeta_vb, 
                         converged, it, lb_opt, diff_lb)
      
    }
  })
  
}



# Internal function which implements the marginal log-likelihood variational
# lower bound (ELBO) corresponding to the `atlasqtl_struct_core` algorithm.
#
elbo_horseshoe_ <- function(Y, xi_inv_vb, A2_inv, lam2_inv_vb, eta, eta_vb, L_vb, 
                                 gam_vb, kappa, kappa_vb, nu, nu_vb, 
                                 nu_xi_inv_vb, nu_s0_vb, m0, n0, zeta_vb,
                                 theta_vb, rho, rho_vb, rho_xi_inv_vb, rho_s0_vb, 
                                 Q_app, sig2_beta_vb, sig02_inv_vb, sig2_theta_vb, 
                                 sig2_inv_vb, sig2_zeta_vb, T0_inv, tau_vb, 
                                 beta_vb, m2_beta, mat_x_m1, vec_sum_log_det_zeta, 
                                 df, shr_fac_inv) {
  
  n <- nrow(Y)
  p <- length(L_vb)
  
  # needed for monotonically increasing elbo.
  #
  eta_vb <- update_eta_vb_(n, eta, gam_vb)
  kappa_vb <- update_kappa_vb_(Y, kappa, mat_x_m1, beta_vb, m2_beta, sig2_inv_vb)
  
  nu_vb <- update_nu_vb_(nu, sum(gam_vb))
  rho_vb <- update_rho_vb_(rho, m2_beta, tau_vb)
  
  log_tau_vb <- update_log_tau_vb_(eta_vb, kappa_vb)
  log_sig2_inv_vb <- update_log_sig2_inv_vb_(nu_vb, rho_vb)
  
  log_sig02_inv_vb <- update_log_sig2_inv_vb_(nu_s0_vb, rho_s0_vb)
  log_xi_inv_vb <- update_log_sig2_inv_vb_(nu_xi_inv_vb, rho_xi_inv_vb)
  
  
  elbo_A <- e_y_(n, kappa, kappa_vb, log_tau_vb, m2_beta, sig2_inv_vb, tau_vb)
  
  
  elbo_B <- e_beta_gamma_(gam_vb, log_sig2_inv_vb, log_tau_vb,
                               zeta_vb, theta_vb, m2_beta,
                               sig2_beta_vb, sig2_zeta_vb,
                               sig2_theta_vb, sig2_inv_vb, tau_vb)
  
  elbo_C <- e_theta_hs_(lam2_inv_vb, L_vb, log_sig02_inv_vb + log(shr_fac_inv), m0, theta_vb, 
                        Q_app, sig02_inv_vb * shr_fac_inv, sig2_theta_vb, df)
  
  elbo_D <- e_zeta_(zeta_vb, n0, sig2_zeta_vb, T0_inv, vec_sum_log_det_zeta)
  
  elbo_E <- e_tau_(eta, eta_vb, kappa, kappa_vb, log_tau_vb, tau_vb)
  
  elbo_F <- sum(e_sig2_inv_hs_(xi_inv_vb, nu_s0_vb, log_xi_inv_vb, log_sig02_inv_vb, rho_s0_vb, sig02_inv_vb)) # sig02_inv_vb
  
  elbo_G <- sum(e_sig2_inv_(1 / 2, nu_xi_inv_vb, log_xi_inv_vb, A2_inv, rho_xi_inv_vb, xi_inv_vb)) # xi_inv_vb
  
  elbo_H <- e_sig2_inv_(nu, nu_vb, log_sig2_inv_vb, rho, rho_vb, sig2_inv_vb)
  
  as.numeric(elbo_A + elbo_B + elbo_C + elbo_D + elbo_E + elbo_F + elbo_G + elbo_H)
  
}

