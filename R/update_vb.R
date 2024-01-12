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

update_m2_beta_ <- function(gam_vb, mu_beta_vb, sig2_beta_vb, sweep = FALSE, mis_pat = NULL) {
  
  if(sweep | is.null(mis_pat)) {
    
    sweep(mu_beta_vb ^ 2, 2, sig2_beta_vb, `+`) * gam_vb
    
  } else {
    
    (mu_beta_vb ^ 2 + sig2_beta_vb) * gam_vb
    
  }
  
}

update_sig2_beta_vb_ <- function(n, sig2_inv_vb, tau_vb = NULL, X_norm_sq = NULL, c = 1) {
  
  if(is.null(tau_vb)) {
    
    if (is.null(X_norm_sq))
      1 / (c * (n - 1 + sig2_inv_vb))
    else
      1 / (c * (X_norm_sq + sig2_inv_vb))
    
  } else {
    
    if (is.null(X_norm_sq))
      1 / (c * (n - 1 + sig2_inv_vb) * tau_vb)
    else
      1 / (c * sweep(X_norm_sq + sig2_inv_vb, 2, tau_vb, `*`))
    
  }
}

update_X_beta_vb_ <- function(X, beta_vb) X %*% beta_vb

update_cp_betaX_X_ <- function(cp_X, beta_vb) crossprod(beta_vb, cp_X)


####################
## b's updates ##
####################

update_annealed_lam2_inv_vb_ <- function(L_vb, c, df) { # here L_vb <- c * L_vb / df
  
  if (df == 1) {
    
    gsl::gamma_inc(- c + 2, L_vb) / (gsl::gamma_inc(- c + 1, L_vb) * L_vb) - 1
    
  } else { # also works for df = 1, but slightly less efficient
    
    (gamma(c * (df - 1) / 2 + 2) * gamma(c) * gsl::hyperg_1F1(c * (df - 1) / 2 + 2, 3 - c, L_vb) / (c - 1) / (c - 2) / gamma(c * (df + 1) / 2) +
       gamma(2 - c) * L_vb^(c - 2) * gsl::hyperg_1F1(c * (df + 1) / 2, c - 1, L_vb) ) /
      (gamma(c * (df - 1) / 2 + 1) * gamma(c) * gsl::hyperg_1F1(c * (df - 1) / 2 + 1, 2 - c, L_vb) / (c - 1) / gamma(c * (df + 1) / 2) +
         gamma(1 - c) * L_vb^(c - 1) * gsl::hyperg_1F1(c * (df + 1) / 2, c, L_vb) ) / df
    
  }
  
}


########################
## c0 and c's updates ##
########################

update_sig2_c0_vb_ <- function(d, s02, c = 1) 1 / (c * (d + (1/s02)))


#####################
## zeta's updates ##
#####################

update_zeta_vb_ <- function(Z, mat_add, n0, sig2_zeta_vb, t02_inv, is_mat = FALSE, c = 1) {
  
  
  if (is_mat) {
    as.vector(c * sig2_zeta_vb * (colSums(Z) + t02_inv * n0 - colSums(mat_add))) # mat_add <- sweep(mat_v_mu, 1, zeta_vb, `-`)
  } else {
    # as.vector(sig2_zeta_vb %*% (colSums(Z) + t02_inv %*% n0 - sum(theta_vb)))
    # sig2_zeta_vb and t02_inv is stored as a scalar which represents the value on the diagonal of the corresponding diagonal matrix
    as.vector(c * sig2_zeta_vb * (colSums(Z) + t02_inv * n0 - sum(mat_add))) # mat_add = theta_vb
  }
  
}

#####################
## sigma's updates ##
#####################

update_nu_vb_ <- function(nu, sum_gam, c = 1) c * (nu + sum_gam / 2) - c + 1

update_rho_vb_ <- function(rho, m2_beta, tau_vb, c = 1) c * as.numeric(rho + crossprod(tau_vb, colSums(m2_beta)) / 2)

update_log_sig2_inv_vb_ <- function(nu_vb, rho_vb) digamma(nu_vb) - log(rho_vb)


###################
## tau's updates ##
###################

update_eta_vb_ <- function(n, eta, gam_vb, mis_pat = NULL, c = 1) {

  if (is.null(mis_pat))
    c * (eta + n / 2 + colSums(gam_vb) / 2) - c + 1
  else
    c * (eta + colSums(mis_pat) / 2 + colSums(gam_vb) / 2) - c + 1
  
}

update_kappa_vb_ <- function(n, Y_norm_sq, cp_Y_X, cp_X, kappa, X_beta_vb, # remove X_beta_vb from arguments once version for missingness implemented
                             beta_vb, m2_beta, sig2_inv_vb, 
                             X_norm_sq = NULL, mis_pat = NULL, c = 1) {
  
  stopifnot(!xor(is.null(X_norm_sq), is.null(mis_pat)))
  
  if (is.null(mis_pat)) {
    
    # c * (kappa + (Y_norm_sq - 2 * colSums(Y * X_beta_vb)  +
    #                 (n - 1 + sig2_inv_vb) * colSums(m2_beta) +
    #                 colSums(X_beta_vb^2) - (n - 1) * colSums(beta_vb^2))/ 2)

    c * (kappa + (Y_norm_sq - 2 * colSums(beta_vb * t(cp_Y_X))  +
                    (n - 1 + sig2_inv_vb) * colSums(m2_beta) +
                    diag(crossprod(beta_vb, cp_X %*% beta_vb)) - (n - 1) * colSums(beta_vb^2))/ 2)
    
  } else {
    
    c * (kappa + (Y_norm_sq - 2 * colSums(Y * X_beta_vb)  +
                    sig2_inv_vb * colSums(m2_beta) + colSums(X_norm_sq * m2_beta) +
                    colSums(X_beta_vb^2 * mis_pat) - colSums(X_norm_sq * beta_vb^2))/ 2)
    
  }
  
}

update_log_tau_vb_ <- function(eta_vb, kappa_vb) digamma(eta_vb) - log(kappa_vb)


#####################
## theta's updates ##
#####################

update_theta_vb_ <- function(Z, m0, sig02_inv, sig2_theta_vb, vec_fac_st,
                                mat_add = 0, is_mat = FALSE, c = 1) {
  
  if (is.null(vec_fac_st)) {
    
    # sig02_inv and sig2_zeta_vb are stored as scalars which represent the values on the diagonal of the corresponding diagonal matrix
    
    if (is_mat) {
      
      theta_vb <- c * sig2_theta_vb * (rowSums(Z) + sig02_inv * m0 - rowSums(mat_add)) # mat_add = sweep(mat_v_mu, 1, theta_vb, `-`)
      
    } else {
      
      theta_vb <- c * sig2_theta_vb * (rowSums(Z) + sig02_inv * m0 - sum(mat_add)) # mat_add = zeta_vb
      
    }
    
    
  } else {
    
    if (c != 1)
      stop("Annealing not implemented when Sigma_0 is not the identity matrix.")
    
    bl_ids <- unique(vec_fac_st)
    n_bl <- length(bl_ids)
    
    if (is_mat) {
      
      theta_vb <- unlist(lapply(1:n_bl, function(bl) {
        sig2_theta_vb[[bl]] %*% (rowSums(Z[vec_fac_st == bl_ids[bl], , drop = FALSE]) +
                                   sig02_inv[[bl]] %*% m0[vec_fac_st == bl_ids[bl]] -
                                   rowSums(mat_add[vec_fac_st == bl_ids[bl], , drop = FALSE]))  # mat_add = sweep(mat_v_mu, 1, theta_vb, `-`)
      }))
    } else {
      
      theta_vb <- unlist(lapply(1:n_bl, function(bl) {
        sig2_theta_vb[[bl]] %*% (rowSums(Z[vec_fac_st == bl_ids[bl], , drop = FALSE]) +
                                   sig02_inv[[bl]] %*% m0[vec_fac_st == bl_ids[bl]] -
                                   sum(mat_add)) # mat_add = zeta_vb
      }))
    }
    
  }
  
}


#################
## Z's updates ##
#################

update_Z_ <- function(gam_vb, mat_v_mu, log_1_pnorm, log_pnorm, c = 1) {
  
  if (!isTRUE(all.equal(c, 1))) {
    
    sqrt_c <- sqrt(c)
    
    log_pnorm <- pnorm(sqrt_c * mat_v_mu, log.p = TRUE)
    log_1_pnorm <- pnorm(sqrt_c * mat_v_mu, log.p = TRUE, lower.tail = FALSE)
    
  } else {
    
    sqrt_c <- 1
  }

  imr0 <- inv_mills_ratio_(0, sqrt_c * mat_v_mu, log_1_pnorm, log_pnorm)
  (gam_vb * (inv_mills_ratio_(1, sqrt_c * mat_v_mu, log_1_pnorm, log_pnorm) - imr0) + imr0) / sqrt_c + mat_v_mu

}
