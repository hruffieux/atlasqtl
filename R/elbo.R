# This file is part of the `atlasqtl` R package:
#     https://github.com/hruffieux/atlasqtl
#
# Internal functions gathering the ELBO terms common to core algorithms.
#
########################################################
## E log p(beta, gamma | rest) - E log q(beta, gamma) ##
########################################################

e_beta_gamma_ <- function(gam_vb, log_1_pnorm, log_pnorm, log_sig2_inv_vb, log_tau_vb,
                          zeta_vb, theta_vb, m2_beta,
                          sig2_beta_vb, sig2_zeta_vb,
                          sig2_theta_vb, sig2_inv_vb, tau_vb) {
  
  eps <- .Machine$double.eps^0.75 # to control the argument of the log when gamma is very small
  
  q <- length(tau_vb)
  
  arg <- log_sig2_inv_vb * gam_vb / 2 +
        sweep(gam_vb, 2, log_tau_vb, `*`) / 2 -
        sweep(m2_beta, 2, tau_vb, `*`) * sig2_inv_vb / 2 +
        gam_vb * log_pnorm +
        (1 - gam_vb) * log_1_pnorm - #log(1 - exp(log_pnorm)) is unstable
        sig2_zeta_vb / 2 -  gam_vb * log(gam_vb + eps) - 
        (1 - gam_vb) * log(1 - gam_vb + eps) - sig2_theta_vb / 2
  
  sig2_beta_vb <- as.matrix(sig2_beta_vb)
  if (ncol(sig2_beta_vb) > 1) {
    sum(arg + 1 / 2 * gam_vb * (log(sig2_beta_vb) + 1))
  } else {
    sum(arg + 1 / 2 * sweep(gam_vb, 2, log(sig2_beta_vb) + 1, `*`))
  }

}


##################################################
## E log p(sig2_inv | rest) - E log q(sig2_inv) ##
##################################################

e_sig2_inv_ <- function(nu, nu_vb, log_sig2_inv_vb, rho, rho_vb, sig2_inv_vb) {
  
  (nu - nu_vb) * log_sig2_inv_vb - (rho - rho_vb) * sig2_inv_vb +
    nu * log(rho) - nu_vb * log(rho_vb) - lgamma(nu) + lgamma(nu_vb)
  
}


e_sig2_inv_hs_ <- function(xi_inv_vb, nu_s0_vb, log_xi_inv_vb, log_sig02_inv_vb, 
                           rho_s0_vb, sig02_inv_vb) {
  
  - 1/2 * log_sig02_inv_vb - xi_inv_vb * sig02_inv_vb + log_xi_inv_vb / 2 - lgamma(1 / 2) -
    (nu_s0_vb - 1) * log_sig02_inv_vb + rho_s0_vb * sig02_inv_vb - 
    nu_s0_vb * log(rho_s0_vb) + lgamma(nu_s0_vb)
  
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

e_theta_ <- function(m0, theta_vb, sig02_inv, sig2_theta_vb, vec_sum_log_det) {
  
  p <- length(theta_vb)
  
  sum(vec_sum_log_det - sig02_inv * crossprod(theta_vb - m0) -
        p * sig02_inv * sig2_theta_vb + p) / 2 
  
}


e_theta_hs_ <- function(lam2_inv_vb, L_vb, log_sig02_inv_vb, m0, theta_vb, Q_app, 
                        sig02_inv_vb, sig2_theta_vb, df) {
  
  if (df == 1) {
    
    sum(log_sig02_inv_vb / 2 - sig02_inv_vb * lam2_inv_vb *
          (theta_vb^2 + sig2_theta_vb - 2 * m0 * theta_vb + m0^2) / 2 +
          (log(sig2_theta_vb) + 1) / 2 - log(pi) + L_vb * lam2_inv_vb + log(Q_app))
    
    
  } else if (df == 3) {
    
    # L_vb is tilde L_vb, i.e., L_vb / df
    
    log_B <- log(9) - log(Q_app * (1 + L_vb) - 1)
    
    sum(log(6) + log(3) / 2 - log(pi) - log_B + df * L_vb * lam2_inv_vb +
          log_sig02_inv_vb / 2  - sig02_inv_vb * lam2_inv_vb *
          (theta_vb^2 + sig2_theta_vb - 2 * m0 * theta_vb + m0^2) / 2 +
          (log(sig2_theta_vb) + 1) / 2)
    
  } else {
    
    # valid for df = odd number, so should also be valid for df = 1 and 3, 
    # but the above is slightly more efficient
    #
    p <- length(lam2_inv_vb)
    
    exponent <- (df + 1) / 2
    
    # L_vb is tilde L_vb, i.e., L_vb / df 
    log_B <- - log(sapply(1:p, function(j) {
      compute_integral_hs_(df, L_vb[j] * df, m = exponent, n = exponent - 1, Q_ab = Q_app[j])}))
    
    sum(-log(pi) / 2 - lgamma(df / 2) + df * log(df) / 2 + lfactorial((df - 1)/2) - 
          log_B + df * L_vb * lam2_inv_vb  +
          log_sig02_inv_vb / 2 - sig02_inv_vb * lam2_inv_vb *
          (theta_vb^2 + sig2_theta_vb - 2 * m0 * theta_vb + m0^2) / 2 +
          (log(sig2_theta_vb) + 1) / 2)
    
  }
  
  
}


#######################
## E log p(y | rest) ##
#######################

e_y_ <- function(n, kappa, kappa_vb, log_tau_vb, m2_beta, sig2_inv_vb, tau_vb, 
                 cs_mis_pat = NULL) {
  
  if (is.null(cs_mis_pat)) {
    arg <- -n / 2 * log(2 * pi) + n / 2 * log_tau_vb
  } else {
    arg <- cs_mis_pat * (log_tau_vb - log(2 * pi)) / 2
  }
  
  sum(arg - tau_vb * (kappa_vb - colSums(m2_beta) * sig2_inv_vb / 2 - kappa))
  
}


##########################################
## E log p(zeta | rest) - E log q(zeta) ##
##########################################

e_zeta_ <- function(zeta_vb, n0, sig2_zeta_vb, t02_inv, vec_sum_log_det_zeta) {
  
  q <- length(zeta_vb)
  
  (vec_sum_log_det_zeta - # vec_sum_log_det_zeta = log(det(t02_inv)) + log(det(sig2_zeta_vb))
    t02_inv * crossprod(zeta_vb - n0) -
    q * t02_inv * sig2_zeta_vb + q) / 2 # trace of a product
  
}

