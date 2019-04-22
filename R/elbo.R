# This file is part of the `atlasqtl` R package:
#     https://github.com/hruffieux/atlasqtl
#
# Internal functions gathering the ELBO terms common to core algorithms.
#
########################################################
## E log p(beta, gamma | rest) - E log q(beta, gamma) ##
########################################################

e_beta_gamma_dual_ <- function(gam_vb, log_sig2_inv_vb, log_tau_vb,
                               mu_rho_vb, mu_theta_vb, m2_beta,
                               sig2_beta_vb, sig2_rho_vb,
                               sig2_theta_vb, sig2_inv_vb, tau_vb) {
  
  eps <- .Machine$double.eps^0.75 # to control the argument of the log when gamma is very small
  
  q <- length(tau_vb)
  
  mat_struct <- sweep(tcrossprod(mu_theta_vb, rep(1, q)), 2, mu_rho_vb, `+`)
  
  sum(log_sig2_inv_vb * gam_vb / 2 +
        sweep(gam_vb, 2, log_tau_vb, `*`) / 2 -
        sweep(m2_beta, 2, tau_vb, `*`) * sig2_inv_vb / 2 +
        gam_vb * pnorm(mat_struct, log.p = TRUE) +
        (1 - gam_vb) * pnorm(mat_struct, lower.tail = FALSE, log.p = TRUE) -
        sig2_rho_vb / 2 + 1 / 2 * sweep(gam_vb, 2, log(sig2_beta_vb) + 1, `*`) -
        gam_vb * log(gam_vb + eps) - (1 - gam_vb) * log(1 - gam_vb + eps) -
        sig2_theta_vb / 2)
  
}


########################################
## E log p(rho | rest) - E log q(rho) ##
########################################

e_rho_ <- function(mu_rho_vb, n0, sig2_rho_vb, T0_inv, vec_sum_log_det_rho) {
  
  q <- length(mu_rho_vb)
  
  (vec_sum_log_det_rho - # vec_sum_log_det_rho = log(det(T0_inv)) + log(det(sig2_rho_vb))
    T0_inv * crossprod(mu_rho_vb - n0) -
    q * T0_inv * sig2_rho_vb + q) / 2 # trace of a product
  
}


##################################################
## E log p(sig2_inv | rest) - E log q(sig2_inv) ##
##################################################

e_sig2_inv_ <- function(lambda, lambda_vb, log_sig2_inv_vb, nu, nu_vb, sig2_inv_vb) {
  
  (lambda - lambda_vb) * log_sig2_inv_vb - (nu - nu_vb) * sig2_inv_vb +
    lambda * log(nu) - lambda_vb * log(nu_vb) - lgamma(lambda) + lgamma(lambda_vb)
  
}


e_sig2_inv_hs_ <- function(a_inv_vb, lambda_s0_vb, log_a_inv_vb, log_S0_inv_vb, 
                           nu_s0_vb, S0_inv_vb) {
  
  - 1/2 * log_S0_inv_vb - a_inv_vb * S0_inv_vb + log_a_inv_vb / 2 - lgamma(1 / 2) -
    (lambda_s0_vb - 1) * log_S0_inv_vb + nu_s0_vb * S0_inv_vb - 
    lambda_s0_vb * log(nu_s0_vb) + lgamma(lambda_s0_vb)
  
}


########################################
## E log p(tau | rest) - E log q(tau) ##
########################################

e_tau_ <- function(eta, eta_vb, kappa, kappa_vb, log_tau_vb, tau_vb) {
  
  sum((eta - eta_vb) * log_tau_vb - (kappa - kappa_vb) * tau_vb +
        eta * log(kappa) - eta_vb * log(kappa_vb) - lgamma(eta) + lgamma(eta_vb))
  
}


############################################
## E log p(theta | rest) - E log q(theta) ##
############################################

e_theta_ <- function(m0, mu_theta_vb, S0_inv, sig2_theta_vb, vec_sum_log_det) {
  
  p <- length(mu_theta_vb)
  
  sum(vec_sum_log_det - S0_inv * crossprod(mu_theta_vb - m0) -
        p * S0_inv * sig2_theta_vb + p) / 2 
  
}


e_theta_hs_ <- function(b_vb, G_vb, log_S0_inv_vb, m0, mu_theta_vb, Q_app, 
                        S0_inv_vb, sig2_theta_vb, df) {
  
  if (df == 1) {
    
    sum(log_S0_inv_vb / 2 - S0_inv_vb * b_vb *
          (mu_theta_vb^2 + sig2_theta_vb - 2 * m0 * mu_theta_vb + m0^2) / 2 +
          (log(sig2_theta_vb) + 1) / 2 - log(pi) + G_vb * b_vb + log(Q_app))
    
    
  } else if (df == 3) {
    
    # G_vb is tilde G_vb, i.e., G_vb / df
    
    log_B <- log(9) - log(Q_app * (1 + G_vb) - 1)
    
    sum(log(6) + log(3) / 2 - log(pi) - log_B + df * G_vb * b_vb +
          log_S0_inv_vb / 2  - S0_inv_vb * b_vb *
          (mu_theta_vb^2 + sig2_theta_vb - 2 * m0 * mu_theta_vb + m0^2) / 2 +
          (log(sig2_theta_vb) + 1) / 2)
    
  } else {
    
    # valid for df = odd number, so should also be valid for df = 1 and 3, 
    # but the above is slightly more efficient
    #
    p <- length(b_vb)
    
    exponent <- (df + 1) / 2
    
    # G_vb is tilde G_vb, i.e., G_vb / df 
    log_B <- - log(sapply(1:p, function(j) {
      compute_integral_hs_(df, G_vb[j] * df, m = exponent, n = exponent - 1, Q_ab = Q_app[j])}))
    
    sum(-log(pi) / 2 - lgamma(df / 2) + df * log(df) / 2 + lfactorial((df - 1)/2) - 
          log_B + df * G_vb * b_vb  +
          log_S0_inv_vb / 2 - S0_inv_vb * b_vb *
          (mu_theta_vb^2 + sig2_theta_vb - 2 * m0 * mu_theta_vb + m0^2) / 2 +
          (log(sig2_theta_vb) + 1) / 2)
    
  }
  
  
}

#######################
## E log p(y | rest) ##
#######################

e_y_ <- function(n, kappa, kappa_vb, log_tau_vb, m2_beta, sig2_inv_vb, tau_vb) {
  
  sum(-n / 2 * log(2 * pi) + n / 2 * log_tau_vb - tau_vb *
        (kappa_vb - colSums(m2_beta) * sig2_inv_vb / 2 - kappa))
  
}
