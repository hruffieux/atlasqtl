# This file is part of the `atlasqtl` R package:
#     https://github.com/hruffieux/atlasqtl
#
# Internal functions gathering the variational updates for the core algorithms.
# Besides improving code readability via modular programming, the main purpose
# is to avoid copy-and-paste programming, as most of these updates (or slightly
# modified versions) are used more than once in the different core algorithms.
# For this reason, we choose to create functions for most variational updates,
# even for those consisting in very basic operations.
# Note that we don't modularize the body of the core for loops for performance
# reasons.

####################
## beta's updates ##
####################

update_beta_vb_ <- function(gam_vb, mu_beta_vb) gam_vb * mu_beta_vb

update_m2_beta_ <- function(gam_vb, mu_beta_vb, sig2_beta_vb, sweep = FALSE) {
  
  if(sweep) {
    
    sweep(mu_beta_vb ^ 2, 2, sig2_beta_vb, `+`) * gam_vb
    
  } else {
    
    (mu_beta_vb ^ 2 + sig2_beta_vb) * gam_vb
    
  }
  
}

update_sig2_beta_vb_ <- function(n, sig2_inv_vb, tau_vb = NULL, c = 1) {
  
  if(is.null(tau_vb)) {
    
    1 / (c * (n - 1 + sig2_inv_vb))
    
  } else {
    
    1 / (c * (n - 1 + sig2_inv_vb) * tau_vb)
    
  }
}

update_mat_x_m1_ <- function(X, beta_vb) X %*% beta_vb


####################
## b's updates ##
####################

update_annealed_b_vb_ <- function(G_vb, c, df) { # here G_vb <- c * G_vb / df
  
  if (df == 1) {
    
    gsl::gamma_inc(- c + 2, G_vb) / (gsl::gamma_inc(- c + 1, G_vb) * G_vb) - 1
    
  } else { # also works for df = 1, but slightly less efficient
    
    (gamma(c * (df - 1) / 2 + 2) * gamma(c) * gsl::hyperg_1F1(c * (df - 1) / 2 + 2, 3 - c, G_vb) / (c - 1) / (c - 2) / gamma(c * (df + 1) / 2) +
       gamma(2 - c) * G_vb^(c - 2) * gsl::hyperg_1F1(c * (df + 1) / 2, c - 1, G_vb) ) /
      (gamma(c * (df - 1) / 2 + 1) * gamma(c) * gsl::hyperg_1F1(c * (df - 1) / 2 + 1, 2 - c, G_vb) / (c - 1) / gamma(c * (df + 1) / 2) +
         gamma(1 - c) * G_vb^(c - 1) * gsl::hyperg_1F1(c * (df + 1) / 2, c, G_vb) ) / df
    
  }
  
}


########################
## c0 and c's updates ##
########################

update_sig2_c0_vb_ <- function(d, s02, c = 1) 1 / (c * (d + (1/s02)))


#####################
## zeta's updates ##
#####################

update_zeta_vb_ <- function(W, mat_add, n0, sig2_zeta_vb, T0_inv, is_mat = FALSE, c = 1) {
  
  
  if (is_mat) {
    as.vector(c * sig2_zeta_vb * (colSums(W) + T0_inv * n0 - colSums(mat_add))) # mat_add <- sweep(mat_v_mu, 1, zeta_vb, `-`)
  } else {
    # as.vector(sig2_zeta_vb %*% (colSums(W) + T0_inv %*% n0 - sum(theta_vb)))
    # sig2_zeta_vb and T0_inv is stored as a scalar which represents the value on the diagonal of the corresponding diagonal matrix
    as.vector(c * sig2_zeta_vb * (colSums(W) + T0_inv * n0 - sum(mat_add))) # mat_add = theta_vb
  }
  
}

#####################
## sigma's updates ##
#####################

update_lambda_vb_ <- function(lambda, sum_gam, c = 1) c * (lambda + sum_gam / 2) - c + 1

update_nu_vb_ <- function(nu, m2_beta, tau_vb, c = 1) c * as.numeric(nu + crossprod(tau_vb, colSums(m2_beta)) / 2)

update_log_sig2_inv_vb_ <- function(lambda_vb, nu_vb) digamma(lambda_vb) - log(nu_vb)


###################
## tau's updates ##
###################

update_eta_vb_ <- function(n, eta, gam_vb, c = 1) c * (eta + n / 2 + colSums(gam_vb) / 2) - c + 1

update_kappa_vb_ <- function(Y, kappa, mat_x_m1, beta_vb, m2_beta, sig2_inv_vb, c = 1) {
  
  n <- nrow(Y)
  
  c * (kappa + (colSums(Y^2) - 2 * colSums(Y * mat_x_m1)  +
                  (n - 1 + sig2_inv_vb) * colSums(m2_beta) +
                  colSums(mat_x_m1^2) - (n - 1) * colSums(beta_vb^2))/ 2)
  
}

update_log_tau_vb_ <- function(eta_vb, kappa_vb) digamma(eta_vb) - log(kappa_vb)


#####################
## theta's updates ##
#####################

update_theta_vb_ <- function(W, m0, S0_inv, sig2_theta_vb, vec_fac_st,
                                mat_add = 0, is_mat = FALSE, c = 1) {
  
  if (is.null(vec_fac_st)) {
    
    # S0_inv and sig2_zeta_vb are stored as scalars which represent the values on the diagonal of the corresponding diagonal matrix
    
    if (is_mat) {
      
      theta_vb <- c * sig2_theta_vb * (rowSums(W) + S0_inv * m0 - rowSums(mat_add)) # mat_add = sweep(mat_v_mu, 1, theta_vb, `-`)
      
    } else {
      
      theta_vb <- c * sig2_theta_vb * (rowSums(W) + S0_inv * m0 - sum(mat_add)) # mat_add = zeta_vb
      
    }
    
    
  } else {
    
    if (c != 1)
      stop("Annealing not implemented when Sigma_0 is not the identity matrix.")
    
    bl_ids <- unique(vec_fac_st)
    n_bl <- length(bl_ids)
    
    if (is_mat) {
      
      theta_vb <- unlist(lapply(1:n_bl, function(bl) {
        sig2_theta_vb[[bl]] %*% (rowSums(W[vec_fac_st == bl_ids[bl], , drop = FALSE]) +
                                   S0_inv[[bl]] %*% m0[vec_fac_st == bl_ids[bl]] -
                                   rowSums(mat_add[vec_fac_st == bl_ids[bl], , drop = FALSE]))  # mat_add = sweep(mat_v_mu, 1, theta_vb, `-`)
      }))
    } else {
      
      theta_vb <- unlist(lapply(1:n_bl, function(bl) {
        sig2_theta_vb[[bl]] %*% (rowSums(W[vec_fac_st == bl_ids[bl], , drop = FALSE]) +
                                   S0_inv[[bl]] %*% m0[vec_fac_st == bl_ids[bl]] -
                                   sum(mat_add)) # mat_add = zeta_vb
      }))
    }
    
  }
  
}


#################
## W's updates ##
#################

update_W_info_ <- function(gam_vb, mat_v_mu, c = 1) {
  
  sqrt_c <- sqrt(c)
  
  (gam_vb * (inv_mills_ratio_(1, sqrt_c * mat_v_mu) - inv_mills_ratio_(0, sqrt_c * mat_v_mu)) +
    inv_mills_ratio_(0, sqrt_c * mat_v_mu)) / sqrt_c + mat_v_mu
  
}
