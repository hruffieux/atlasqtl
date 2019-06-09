# This file is part of the `atlasqtl` R package:
#     https://github.com/hruffieux/atlasqtl
#
# Internal core function to call the variational algorithm for global hotspot 
# propensity modelling. 
# See help of `atlasqtl` function for details.
#
atlasqtl_global_core_ <- function(Y, X, shr_fac_inv, anneal, df, tol, maxit, 
                                  verbose, list_hyper, list_init, 
                                  checkpoint_path = NULL, full_output = FALSE, 
                                  debug = FALSE, batch = "y") {
  
  
  q <- ncol(Y)
  n <- nrow(Y)
  p <- ncol(X)
  
  
  # Gathering initial variational parameters
  #
  gam_vb <- list_init$gam_vb
  mu_beta_vb <- list_init$mu_beta_vb
  sig02_inv_vb <- list_init$sig02_inv_vb
  sig2_beta_vb <- list_init$sig2_beta_vb
  tau_vb <-list_init$tau_vb
  
  theta_vb <- list_init$theta_vb
  zeta_vb <- list_init$zeta_vb
  
  rm(list_init)
  
  theta_plus_zeta_vb <- sweep(tcrossprod(theta_vb, rep(1, q)), 2, zeta_vb, `+`)
  log_Phi_theta_plus_zeta <- pnorm(theta_plus_zeta_vb, log.p = TRUE)
  log_1_min_Phi_theta_plus_zeta <- pnorm(theta_plus_zeta_vb, log.p = TRUE, lower.tail = FALSE) 
  
  # Preparing annealing if any
  #
  anneal_scale <- TRUE # if TRUE, scale parameters s02 and lam2_inv_vb also annealed.
  
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
  
  nu_s0 <- rho_s0 <- 1 / 2 # gives rise to a Cauchy prior for theta if = 1/2, otherwise, Student t if rho_s0 = 1 / (2*q)
  
  
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
    m2_beta <- update_m2_beta_(gam_vb, mu_beta_vb, sig2_beta_vb, sweep = TRUE)
    
    mat_x_m1 <- update_mat_x_m1_(X, beta_vb)
    
    
    converged <- FALSE
    lb_new <- -Inf
    it <- 0
    
    
    while ((!converged) & (it < maxit)) {
      
      lb_old <- lb_new
      it <- it + 1
      
      if (verbose != 0 & (it == 1 | it %% 5 == 0))
        cat(paste0("Iteration ", format(it), "... \n"))
      
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
      
      Z <- update_Z_(gam_vb, theta_plus_zeta_vb, c = c) # we use info_ so that the second argument is a matrix
      
      # keep this order!
      #
      sig2_theta_vb <- update_sig2_c0_vb_(q, 1 / sig02_inv_vb / shr_fac_inv, c = c)
      
      theta_vb <- update_theta_vb_(Z, m0, sig02_inv_vb * shr_fac_inv, sig2_theta_vb,
                                   vec_fac_st = NULL, zeta_vb, is_mat = FALSE, c = c)
      
      zeta_vb <- update_zeta_vb_(Z, theta_vb, n0, sig2_zeta_vb, t02_inv,
                                 is_mat = FALSE, c = c) # update_zeta_vb_(Z, theta_vb, sig2_zeta_vb)
      
      theta_plus_zeta_vb <- sweep(tcrossprod(theta_vb, rep(1, q)), 2, zeta_vb, `+`)
      log_Phi_theta_plus_zeta <- pnorm(theta_plus_zeta_vb, log.p = TRUE)
      log_1_min_Phi_theta_plus_zeta <- pnorm(theta_plus_zeta_vb, log.p = TRUE, lower.tail = FALSE)  
      
      nu_s0_vb <- c_s * (nu_s0 + p / 2) - c_s + 1 # implement annealing
      rho_s0_vb <- c_s * (rho_s0 + sum(sig2_theta_vb + theta_vb^2 - 2 * theta_vb * m0 + m0^2) / 2)
      
      sig02_inv_vb <- as.numeric(nu_s0_vb / rho_s0_vb)
      
      if (verbose == 2 & (it == 1 | it %% 5 == 0)) {
        
        cat(paste0("Variational hotspot propensity global scale: ", 
                   format(sqrt(rho_s0_vb / (nu_s0_vb - 1) / shr_fac_inv), digits = 3), ".\n"))
        
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
        
        lb_new <- elbo_global_(Y, beta_vb, eta, eta_vb, gam_vb, kappa, kappa_vb, 
                               log_1_min_Phi_theta_plus_zeta, log_Phi_theta_plus_zeta, 
                               m0, m2_beta, mat_x_m1, n0, nu, nu_s0, nu_s0_vb, 
                               nu_vb, rho, rho_s0, rho_s0_vb, rho_vb, 
                               shr_fac_inv, sig02_inv_vb, sig2_beta_vb, 
                               sig2_inv_vb, sig2_theta_vb, sig2_zeta_vb, 
                               t02_inv, tau_vb, theta_vb, vec_sum_log_det_zeta, 
                               zeta_vb)
        
        if (verbose != 0 & (it == 1 | it %% 5 == 0)) 
          cat(paste0("ELBO = ", format(lb_new), "\n\n"))
        
        
        if (debug && lb_new + eps < lb_old)
          stop("ELBO not increasing monotonically. Exit. ")
        
        converged <- (abs(lb_new - lb_old) < tol)
        
        checkpoint_(it, checkpoint_path, gam_vb, converged, lb_new, lb_old, 
                    zeta_vb = zeta_vb, theta_vb = theta_vb, 
                    sig02_inv_vb = sig02_inv_vb)
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
      
      create_named_list_(beta_vb, eta_vb, gam_vb, kappa_vb, nu_s0_vb, nu_vb,  
                         rho_s0_vb, rho_vb, shr_fac_inv, sig02_inv_vb, 
                         sig2_beta_vb, sig2_inv_vb, sig2_theta_vb, sig2_zeta_vb, 
                         tau_vb, theta_vb, zeta_vb)
      
    } else {
      
      names_x <- colnames(X)
      names_y <- colnames(Y)
      
      rownames(gam_vb) <- names_x
      colnames(gam_vb) <- names_y
      names(theta_vb) <- names_x
      names(zeta_vb) <- names_y
      
      diff_lb <- abs(lb_opt - lb_old)
      
      create_named_list_(beta_vb, gam_vb, theta_vb, zeta_vb, converged, it, 
                         lb_opt, diff_lb)
      
    }
  })
  
}



# Internal function which implements the marginal log-likelihood variational
# lower bound (ELBO) corresponding to the `atlasqtl_struct_core` algorithm.
#
elbo_global_ <- function(Y, beta_vb, eta, eta_vb, gam_vb, kappa, kappa_vb, 
                         log_1_min_Phi_theta_plus_zeta, log_Phi_theta_plus_zeta, m0, 
                         m2_beta, mat_x_m1, n0, nu, nu_s0, nu_s0_vb, nu_vb, rho, 
                         rho_s0, rho_s0_vb, rho_vb, shr_fac_inv, sig02_inv_vb, 
                         sig2_beta_vb, sig2_inv_vb, sig2_theta_vb, sig2_zeta_vb, 
                         t02_inv, tau_vb, theta_vb, vec_sum_log_det_zeta, 
                         zeta_vb) {
  
  n <- nrow(Y)
  p <- length(theta_vb)
  
  # needed for monotonically increasing elbo.
  #
  eta_vb <- update_eta_vb_(n, eta, gam_vb)
  kappa_vb <- update_kappa_vb_(Y, kappa, mat_x_m1, beta_vb, m2_beta, sig2_inv_vb)
  
  nu_vb <- update_nu_vb_(nu, sum(gam_vb))
  rho_vb <- update_rho_vb_(rho, m2_beta, tau_vb)
  
  log_tau_vb <- update_log_tau_vb_(eta_vb, kappa_vb)
  log_sig2_inv_vb <- update_log_sig2_inv_vb_(nu_vb, rho_vb)
  
  log_sig02_inv_vb <- update_log_sig2_inv_vb_(nu_s0_vb, rho_s0_vb)
  
  vec_sum_log_det_theta <- p * (log_sig02_inv_vb + log(shr_fac_inv) + log(sig2_theta_vb)) # E(log(det(sig02_inv))) + log(det(sig2_theta_vb_bl))
  
  elbo_A <- e_y_(n, kappa, kappa_vb, log_tau_vb, m2_beta, sig2_inv_vb, tau_vb)
  
  elbo_B <- e_beta_gamma_(gam_vb, log_1_min_Phi_theta_plus_zeta, log_Phi_theta_plus_zeta, log_sig2_inv_vb, 
                          log_tau_vb, zeta_vb, 
                          theta_vb, m2_beta, sig2_beta_vb, sig2_zeta_vb,
                          sig2_theta_vb, sig2_inv_vb, tau_vb)
  
  elbo_C <- e_theta_(m0, theta_vb, shr_fac_inv * sig02_inv_vb, sig2_theta_vb, 
                     vec_sum_log_det_theta)
  
  elbo_D <- e_zeta_(zeta_vb, n0, sig2_zeta_vb, t02_inv, vec_sum_log_det_zeta)
  
  elbo_E <- e_tau_(eta, eta_vb, kappa, kappa_vb, log_tau_vb, tau_vb)
  
  elbo_F <- e_sig2_inv_(nu, nu_vb, log_sig2_inv_vb, rho, rho_vb, sig2_inv_vb)
  
  elbo_G <- e_sig2_inv_(nu_s0, nu_s0_vb, log_sig02_inv_vb, rho_s0, rho_s0_vb, 
                        sig02_inv_vb)
  
  
  as.numeric(elbo_A + elbo_B + elbo_C + elbo_D + elbo_E + elbo_F + elbo_G)
  
}

