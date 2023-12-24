# This file is part of the `atlasqtl` R package:
#     https://github.com/hruffieux/atlasqtl
#

# Diverse utility functions implementing sanity checks, basic preprocessing,
# and ticks to prevent overflow/underflow.
#


check_natural_ <- function(x, eps = .Machine$double.eps^0.75){
  if (any(x < eps | abs(x - round(x)) > eps)) {
    stop(paste0(deparse(substitute(x)),
                " must be natural."))
  }
}

check_positive_ <- function(x, eps = .Machine$double.eps^0.75){
  if (any(x < eps)) {
    err_mess <- paste0(deparse(substitute(x)), " must be positive, greater than ",
                       format(eps, digits = 3), ".")
    if (length(x) > 1) err_mess <- paste0("All entries of ", err_mess)
    stop(err_mess)
  }
}

check_zero_one_ <- function(x){
  if (any(x < 0) | any(x > 1)) {
    err_mess <- paste0(deparse(substitute(x)), " must lie between 0 and 1.")
    if (length(x) > 1) err_mess <- paste0("All entries of ", err_mess)
    stop(err_mess)
  }
}

check_structure_ <- function(x, struct, type, size = NULL,
                             null_ok = FALSE,  inf_ok = FALSE, na_ok = FALSE) {
  if (type == "double") {
    bool_type <-  is.double(x)
    type_mess <- "a double-precision "
  } else if (type == "integer") {
    bool_type <- is.integer(x)
    type_mess <- "an integer "
  } else if (type == "numeric") {
    bool_type <- is.numeric(x)
    type_mess <- "a numeric "
  } else if (type == "logical") {
    bool_type <- is.logical(x)
    type_mess <- "a boolean "
  } else if (type == "string") {
    bool_type <- is.character(x)
    type_mess <- "string "
  }

  bool_size <- TRUE # for case size = NULL (no assertion on the size/dimension)
  size_mess <- ""
  if (struct == "vector") {
    bool_struct <- is.vector(x) & (length(x) > 0) # not an empty vector
    if (!is.null(size)) {
      bool_size <- length(x) %in% size
      size_mess <- paste0(" of length ", paste0(size, collapse=" or "))
    }
  } else if (struct == "matrix") {
    bool_struct <- is.matrix(x) & (length(x) > 0) # not an empty matrix
    if (!is.null(size)) {
      bool_size <- all(dim(x) == size)
      size_mess <- paste0(" of dimension ", size[1], " x ", size[2])
    }
  }

  correct_obj <- bool_struct & bool_type & bool_size

  bool_null <- is.null(x)

  if (!is.list(x) & type != "string") {
    na_mess <- ""
    if (!na_ok) {
      if (!bool_null) correct_obj <- correct_obj & !any(is.na(x))
      na_mess <- " without missing value"
    }

    inf_mess <- ""
    if (!inf_ok) {
      if (!bool_null) correct_obj <- correct_obj & all(is.finite(x[!is.na(x)]))
      inf_mess <- ", finite"
    }
  } else {
    na_mess <- ""
    inf_mess <- ""
  }

  null_mess <- ""
  if (null_ok) {
    correct_obj <- correct_obj | bool_null
    null_mess <- " or must be NULL"
  }

  if(!(correct_obj)) {
    stop(paste0(deparse(substitute(x)), " must be a non-empty ", type_mess, struct,
               size_mess, inf_mess, na_mess, null_mess, "."))
  }
}


create_named_list_ <- function(...) {
  setNames(list(...), as.character(match.call()[-1]))
}


get_annealing_ladder_ <- function(anneal, verbose) {

  # ladder set following:
  # Importance Tempering, Robert B. Gramacy & Richard J. Samworth, pp.9-10, arxiv v4

  k_m <- 1 / anneal[2]
  m <- anneal[3]

  if(anneal[1] == 1) {

    type <- "geometric"

    delta_k <- k_m^(1 / (1 - m)) - 1

    ladder <- (1 + delta_k)^(1 - m:1)

  } else if (anneal[1] == 2) { # harmonic spacing

    type <- "harmonic"

    delta_k <- ( 1 / k_m - 1) / (m - 1)

    ladder <- 1 / (1 + delta_k * (m:1 - 1))

  } else { # linear spacing

    type <- "linear"

    delta_k <- (1 - k_m) / (m - 1)

    ladder <- k_m + delta_k * (1:m - 1)
  }

  if (verbose != 0)
    cat(paste0("** Annealing with ", type," spacing ** \n\n"))

  ladder

}


log_one_plus_exp_ <- function(x) { # computes log(1 + exp(x)) avoiding
                                   # numerical overflow
  m <- x
  m[x < 0] <- 0

  log(exp(x - m) + exp(- m)) + m
}

log_det <- function(list_mat) {

  if (is.list(list_mat)) {
    sapply(list_mat, function(mat) {
      log_det <- determinant(mat, logarithm = TRUE)
      log_det$modulus * log_det$sign
    })
  } else {
    log_det <- determinant(list_mat, logarithm = TRUE)
    log_det$modulus * log_det$sign
  }

}


inv_mills_ratio_ <- function(y, U, log_1_pnorm_U, log_pnorm_U) {
  
  stopifnot(y %in% c(0, 1))
  
  # writing explicitly the formula for pnorm(, log = TRUE) is faster...
  if (y == 1) {
    
    m <- exp(-U^2/2 - log(sqrt(2*pi)) - log_pnorm_U)
    m[m < -U] <- -U
    
  } else {
    
    m <- - exp(-U^2/2 - log(sqrt(2*pi)) - log_1_pnorm_U)
    m[m > -U] <- -U
    
  }
  
  m
  
}


log_sum_exp_ <- function(x) {
  # Computes log(sum(exp(x))

  if ( max(abs(x)) > max(x) )
    offset <- min(x)
  else
    offset <- max(x)
  log(sum(exp(x - offset))) + offset

}


# entropy_ <- function(Y, U) {
#
#   log((2 * pi * exp(1))^(1/2) *
#            exp(Y * pnorm(U, log.p = TRUE) +
#                  (1-Y) * pnorm(U, lower.tail = FALSE, log.p = TRUE))) -
#   U * inv_mills_ratio_(Y, U) / 2
#
# }


# Functions for hyperparameter settings
#
E_Phi_X <- function(mu, s2, lower_tail = TRUE) {

  pnorm(mu / sqrt(1 + s2), lower.tail = lower_tail)

}

E_Phi_X_2 <- function(mu, s2) {

  pnorm(mu / sqrt(1 + s2)) -
    2 * PowerTOST::OwensT(mu / sqrt(1 + s2), 1 / sqrt(1 + 2 * s2))

}

get_V_p_t <- function(mu, s2, p) {
  p * (p - 1) * E_Phi_X_2(mu, s2) -
    p^2 * E_Phi_X(mu, s2)^2 +
    p * E_Phi_X(mu, s2)
}


get_mu <- function(E_p_t, s2, p) {

  sqrt(1 + s2) * qnorm(E_p_t / p)

}


get_n0_t02 <- function(q, p, p_star) {
  
  E_p_t <- p_star[1]
  V_p_t <- min(p_star[2], floor(2 * p / 3))

  dn <- 1e-6
  up <- 1e5
  
  # Get n0 and t02 similarly as for a_omega_t and b_omega_t in HESS 
  # (specify expectation and variance of number of active predictors per response)
  #
  # Look at : gam_st | theta_s = 0
  #
  tryCatch(t02 <- uniroot(function(x)
    get_V_p_t(get_mu(E_p_t, x, p), x, p) - V_p_t,
    interval = c(dn, up))$root,
    error = function(e) {
      stop(paste0("No hyperparameter values matching the expectation and variance ",
                  "of the number of active predictors per responses supplied in p0. ",
                  "Please change p0."))
    })
  
  # n0 sets the level of sparsity.
  n0 <- get_mu(E_p_t, t02, p)
  n0 <- rep(n0, q)
  
  create_named_list_(n0, t02)
}



rm_constant_ <- function(mat, verbose) {

  bool_cst <- is.nan(colSums(mat))

  if (any(bool_cst)) {

    rmvd_cst <- colnames(mat)[bool_cst]

    if (verbose != 0) {
      if (sum(bool_cst) < 50) {
        cat(paste0("Variable(s) ", paste0(rmvd_cst, collapse=", "),
                   " constant across subjects. \n",
                   "Removing corresponding column(s) and saving its/their id(s) ",
                   "in the function output ... \n\n"))
      } else {
        cat(paste0(sum(bool_cst), " variables constant across subjects. \n",
                   "Removing corresponding column(s) and saving their ids ",
                   "in the function output ... \n\n"))
      }
    }

    mat <- mat[, !bool_cst, drop = FALSE]
  } else {
    rmvd_cst <- NULL
  }

  create_named_list_(mat, bool_cst, rmvd_cst)
}

rm_collinear_ <- function(mat, verbose) {

  bool_coll <- duplicated(mat, MARGIN = 2)

  if (any(bool_coll)) {

    mat_coll <- mat[, bool_coll, drop = FALSE]
    rmvd_coll <- colnames(mat_coll)

    if (verbose != 0) {
      if (length(rmvd_coll) < 50) {
        cat(paste0("Presence of collinear variable(s). ",
                   paste0(rmvd_coll, collapse=", "), " redundant. \n",
                   "Removing corresponding column(s) and saving its/their id(s) ",
                   "in the function output ... \n"))
      } else {
        cat(paste0("Presence of collinear variables. ", length(rmvd_coll),
                   " redundant.\n", "Removing corresponding columns and saving ",
                   "their ids in the function output ... \n"))
      }
    }

    # associate to each removed replicate the name of the covariate with which
    # it is duplicated and that is kept in the dataset
    bool_with_coll <- duplicated(mat, MARGIN = 2, fromLast = TRUE) & !bool_coll
    mat_with_coll <- mat[, bool_with_coll, drop = FALSE]

    assoc_coll <- colnames(mat_with_coll)[match(data.frame(mat_coll),
                                                data.frame(mat_with_coll))]
    names(rmvd_coll) <- assoc_coll

    mat <- mat[, !bool_coll, drop = FALSE]

  } else {
    rmvd_coll <- NULL
  }

  create_named_list_(mat, bool_coll, rmvd_coll)
}


Q_approx <- function(x, eps1 = 1e-30, eps2 = 1e-7) {
  
  if(x <= 1) {
    
    gsl::expint_E1(x) * exp(x)
    
  } else {
    
    f_p <- eps1
    C_p <- eps1
    D_p <- 0
    Delta <- 2 + eps2
    j <- 1
    
    while( abs(Delta-1) >= eps2 ) {
      
      j <- j+1
      
      D_c <- x + 2*j - 1 - ((j-1)^{2}) * D_p
      C_c <- x + 2*j - 1 - ((j-1)^{2}) / C_p
      D_c <- 1 / D_c
      
      Delta <- C_c * D_c
      f_c <- f_p * Delta
      f_p <- f_c
      C_p <- C_c
      D_p <- D_c
    }
    
    1/(x + 1 + f_c)
  }
}


Q_approx_vec <- function(x, eps1 = 1e-30, eps2 = 1e-7) {
  
  qapprox <- rep(NA, length(x))
  
  x_lower <- x[x <= 1]
  
  if (length(x_lower) > 0) {
    qapprox[x <= 1] <- gsl::expint_E1(x_lower) * exp(x_lower)
  }
  
  
  x_upper <- x[x > 1]
  
  if (length(x_upper) > 0) {
    
    f_p <- eps1
    C_p <- eps1
    D_p <- 0
    Delta <- 2 + eps2
    j <- 1
    
    
    while( max(abs(Delta-1)) >= eps2 ) {
      
      j <- j+1
      
      D_c <- x_upper + 2*j - 1 - ((j-1)^{2}) * D_p
      C_c <- x_upper + 2*j - 1 - ((j-1)^{2}) / C_p
      D_c <- 1 / D_c
      
      Delta <- C_c * D_c
      f_c <- f_p * Delta
      f_p <- f_c
      C_p <- C_c
      D_p <- D_c
    }
    
    qapprox[x > 1] <- 1/(x_upper + 1 + f_c)
    
  }
  
  qapprox
  
}

compute_integral_hs_ <- function(alpha, beta, m, n, Q_ab) {
  
  # computes int_0^infty x^n (1 + alpha * x)^(-m) * exp(- beta * x) dx
  # for m = n or m = n + 1, n natural n > 0, beta > 0 (if n = 0, then = Q_ab)
  
  # Q_ab = Q_approx(alpha / beta) # precomputed to avoid computing it too many times
  # = exp(beta/alpha) * E_1(beta/alpha)
  
  
  if (m == n) {
    
    out <- alpha^(-n) * beta^(-1) 
    
    if (n == 1) {
      
      out <- out - alpha^(-2) * Q_ab
      
    } else if (n == 2) {
      
      out <- out - 2 * alpha^(-3) * Q_ab + alpha^(-3) - alpha^(-4) * beta * Q_ab
      
    } else if (n == 3) {
      
      v1 <- c(-3 * log(alpha) - log(beta),
              log(3) - 4 * log(alpha),
              -5 * log(alpha) - log(2) + log(beta))
      
      
      v2 <- c(log(3) - 4 * log(alpha) + log(Q_ab),
              log(3) - 5 * log(alpha) + log(beta) + log(Q_ab),
              -4 * log(alpha) - log(2),
              -6 * log(alpha) - log(2) + 2 * log(beta) + log(Q_ab))
      
      out <- exp(log_sum_exp_(v1)) - exp(log_sum_exp_(v2))
      
    } else if (n == 4){
      
      v1 <- c(-4 * log(alpha) - log(beta),
              log(4) - 5 * log(alpha),
              log(2) - 5 * log(alpha),
              log(2) - 7 * log(alpha) + 2 * log(beta) + log(Q_ab),
              -7 * log(alpha) - log(6) + 2 * log(beta),
              -5 * log(alpha) - log(3))
      
      v2 <- c(log(4) - 5 * log(alpha) + log(Q_ab),
              log(4) - 6 * log(alpha) + log(beta) + log(Q_ab),
              log(2) - 7 * log(alpha) + log(beta),
              -6 * log(alpha) - log(6) + log(beta))
      
      out <- exp(log_sum_exp_(v1)) - exp(log_sum_exp_(v2))
      
    } else {
      
      v1 <- c(-n * log(alpha) - log(beta), 
              -n * log(alpha) + log(n) + unlist(sapply(2:(n-1), function(k) {
                -k * log(alpha) - lfactorial(k-1) +
                  sapply(seq(1, k-1, by = 2), function(j) {
                    lfactorial(j-1) + (k-j-1) * log(beta) + j * log(alpha)}) })),
              -2*n *log(alpha) - lfactorial(n-1) +
                sapply(seq(1, n-1, by = 2), function(j) {
                  lfactorial(j-1) + (n-j-1) * log(beta) + j * log(alpha)})
      )
      
      v2 <- c(log(n) - (n + 1) * log(alpha) + log(Q_ab), 
              -n * log(alpha) + log(n) + unlist(sapply(3:(n-1), function(k) { # k = 2 doesn't contribute
                -k * log(alpha) - lfactorial(k-1) +
                  sapply(seq(2, k-1, by = 2), function(j) {
                    lfactorial(j-1) * (k-j-1) * log(beta) + j * log(alpha)}) })),
              -n * log(alpha) + log(n) + sapply(2:(n-1), function(k) {
                -k * log(alpha) - lfactorial(k-1) + (k-1) * log(beta) + log(Q_ab)}),
              - 2*n * log(alpha) - lfactorial(n-1) +
                sapply(seq(2, n-1, by = 2), function(j) {
                  lfactorial(j-1) + (n-j-1) * log(beta) + j * log(alpha)}), 
              -2*n * log(alpha) - lfactorial(n-1) + (n-1) * log(beta) + log(Q_ab)
      )
      
      out <- exp(log_sum_exp_(v1)) - exp(log_sum_exp_(v2))
      
    }
    
    
  } else if (m == n + 1) { # not stable for n >= 4, i.e., can't be used for df = 9 and higher.
    
    if (n == 1) {
      
      out <- alpha^(-2) * Q_ab - alpha^(-2) +  alpha^(-3) * beta * Q_ab
      
    } else if (n == 2){
      
      v1 <- c(-3 * log(alpha) + log(Q_ab),
              -3 * log(alpha) - log(2),
              -5 * log(alpha) - log(2) + 2 * log(beta) + log(Q_ab),
              -4 * log(alpha) + log(2) + log(beta) + log(Q_ab)
      )
      
      v2 <- c(-4 * log(alpha) - log(2) + log(beta),
              -3 * log(alpha) + log(2)
      )
      
      
      out <- exp(log_sum_exp_(v1)) - exp(log_sum_exp_(v2))
      
    } else {
      
      
      v1 <- c(-(n+1) * log(alpha) + log(Q_ab),
              -(2*n+1) * log(alpha) - lfactorial(n) +
                sapply(seq(2, n, by = 2), function(j) {
                  lfactorial(j-1) + (n-j) * log(beta) + j * log(alpha)}),
              -(2*n+1) * log(alpha) - lfactorial(n) + n * log(beta) + log(Q_ab),
              -n * log(alpha) + log(n) + unlist(sapply(2:(n-1), function(k) {
                -(1+k) * log(alpha) - lfactorial(k) +
                  sapply(seq(2, k, by = 2), function(j) {
                    lfactorial(j-1) + (k-j) * log(beta) + j * log(alpha) })
              })),
              -n * log(alpha) + log(n) + sapply(1:(n-1), function(k) {
                -(1+k) * log(alpha) - lfactorial(k) + k * log(beta) + log(Q_ab)
              })
      )
      
      
      v2 <- c(-(2*n+1) * log(alpha) - lfactorial(n) +
                sapply(seq(1, n, by = 2), function(j) {
                  lfactorial(j-1) + (n-j) * log(beta) + j * log(alpha)}),
              -n * log(alpha) + log(n) + unlist(sapply(1:(n-1), function(k) {
                -(1+k) * log(alpha) - lfactorial(k) +
                  sapply(seq(1, k, by = 2), function(j) {
                    lfactorial(j-1) + (k-j) * log(beta) + j * log(alpha) })
              }))
              
      )
      
      out <- exp(log_sum_exp_(v1)) - exp(log_sum_exp_(v2))
      
    }
    
  } else {
    
    stop("Invalid value of m, must be n or n + 1.")
    
  }
  
  out
}


checkpoint_ <- function(it, checkpoint_path, 
                        beta_vb, gam_vb, theta_vb, zeta_vb, 
                        converged, lb_new, lb_old, 
                        lam2_inv_vb = NULL, sig02_inv_vb = NULL, rate = 100,
                        names_x = NULL, names_y = NULL){
  
  if (!is.null(checkpoint_path) && it %% rate == 0) {
    
    diff_lb <- abs(lb_new - lb_old)
    
    if (!is.null(names_x)) {
      
      rownames(gam_vb) <- rownames(beta_vb) <- names(theta_vb) <- names_x

      if (!is.null(lam2_inv_vb)) {
        names(lam2_inv_vb) <- names_x
      }

    }
    
    if (!is.null(names_y)) {
      
      colnames(gam_vb) <- colnames(beta_vb) <- names(zeta_vb) <- names_y
      
    }
    
    tmp_vb <- create_named_list_(beta_vb, gam_vb, theta_vb, zeta_vb, converged, it, lb_new, diff_lb, 
                                 lam2_inv_vb, sig02_inv_vb)

    file_save <- paste0(checkpoint_path, "tmp_output_it_", it, ".RData")
    
    save(tmp_vb, file = file_save)
    
    old_file_clean_up <- paste0(checkpoint_path, "tmp_output_it_", it - 2 * rate, ".RData") # keep only the last two for comparison
    
    if (file.exists(old_file_clean_up)) 
      file.remove(old_file_clean_up)

  }
    
}


checkpoint_clean_up_ <- function(checkpoint_path) {
  
  if (!is.null(checkpoint_path)) {
    
    old_files_clean_up <- list.files(path = checkpoint_path, pattern = "tmp_output_it_")
    
    sapply(old_files_clean_up, function(ff) {
      if (file.exists(file.path(checkpoint_path, ff))) 
        file.remove(file.path(checkpoint_path, ff))
    })
    
  } 
  
}
 

plot_trace_var_hs_ <- function(lam2_inv_vb, sig02_inv_vb, shr_fac_inv, it, trace_ind_max, trace_var_max, path_trace) {
  
  x <- 1 / sig02_inv_vb * 1 / lam2_inv_vb / shr_fac_inv
  
  vec_ind_max <- which(x == max(x))
  
  nb_max <- 4
  for (i in 2:nb_max) {
    
    vec_ind_max <- c(vec_ind_max, which(x == max(x[-vec_ind_max])))
    
  }
  
  trace_ind_max <- rbind(trace_ind_max, vec_ind_max)
  trace_var_max <- rbind(trace_var_max, x[vec_ind_max])
  rownames(trace_var_max)[nrow(trace_var_max)] <- rownames(trace_ind_max)[nrow(trace_ind_max)] <- it
  
  # display
  list_changepoints <- lapply(1:nb_max, function(i) {
    1 + which(diff(trace_ind_max[, i]) != 0)
  })
  
  png(file.path(path_trace, "traces_top_local_x_global_parameters.png"), width = 5000, height = 4000,
      res = 600, type="cairo")
  
  matplot(rownames(trace_var_max), trace_var_max, type = "o", col = "black", bg = 2:5, pch = 16,
          main = "Trace 1 / mu_s0_inv_vb x 1 / lam2_inv_vb x shr_factor", xlab = "Iteration",
          ylab = paste0("Traces for the ", nb_max, " largest variataional parameters related to the hotspot variances"))
  for (i in 1:nb_max) {
    points(rownames(trace_var_max)[list_changepoints[[i]]], trace_var_max[list_changepoints[[i]], i], col = "blue", pch = 16)
  }
  abline(h = 5, col = "red", lty = 3)
  legend("topleft", legend = "Current predictor different from that of previous iterations on the path.", col = "blue", pch = 16, bty = "n", cex = 0.8)
  
  dev.off()
  
  create_named_list_(trace_ind_max, trace_var_max)
  
}


add_collinear_back_ <- function(beta_vb, gam_vb, theta_vb, initial_colnames_X, rmvd_coll_x, verbose) {
  
  stopifnot(length(rmvd_coll_x) > 0)
  
  p_all <- length(initial_colnames_X) # all the (non-constant) variables, including the collinear variables
                                      # (we don't add any posterior summary for the constant variables...)
  
  if (verbose) {
    cat("Adding the posterior summary to the final output for the collinear X variables: \n")
    print(rmvd_coll_x)
  }

  gam_vb_wo_coll <- gam_vb
  beta_vb_wo_coll <- beta_vb
  theta_vb_wo_coll <- theta_vb
  
  # matrices / vectors with the results for the collinear X variables added back
  #
  gam_vb <- beta_vb <- matrix(NA, nrow = p_all, ncol = ncol(gam_vb_wo_coll))
  theta_vb <- rep(NA, p_all)

  # boolean vector indicating the position of the collinear X variables
  #
  bool_rmvd <- initial_colnames_X %in% rmvd_coll_x
  
  # original results (no result for the redundant X variables)
  #
  gam_vb[!bool_rmvd,] <- gam_vb_wo_coll
  rownames(gam_vb)[!bool_rmvd] <- rownames(gam_vb_wo_coll)
  
  beta_vb[!bool_rmvd,] <- beta_vb_wo_coll
  rownames(beta_vb)[!bool_rmvd] <- rownames(beta_vb_wo_coll)
  
  theta_vb[!bool_rmvd] <- theta_vb_wo_coll
  names(theta_vb)[!bool_rmvd] <- names(theta_vb_wo_coll)
  
  # results for redundant X variables rebuilt from those X variables to which they are collinear
  #
  gam_vb_coll_var <- gam_vb[match(names(rmvd_coll_x), initial_colnames_X),, drop = FALSE]
  rownames(gam_vb_coll_var) <- rmvd_coll_x
  
  beta_vb_coll_var <- beta_vb[match(names(rmvd_coll_x), initial_colnames_X),, drop = FALSE]
  rownames(beta_vb_coll_var) <- rmvd_coll_x
  
  theta_vb_coll_var <- theta_vb[match(names(rmvd_coll_x), initial_colnames_X)]
  names(theta_vb_coll_var) <- rmvd_coll_x
  
  # results for redundant X variables added back to the matrices / vectors
  #
  gam_vb[match(rmvd_coll_x, initial_colnames_X),] <- gam_vb_coll_var
  rownames(gam_vb)[match(rmvd_coll_x, initial_colnames_X)] <- rownames(gam_vb_coll_var)
  colnames(gam_vb) <- colnames(gam_vb_wo_coll)
  
  beta_vb[match(rmvd_coll_x, initial_colnames_X),] <- beta_vb_coll_var
  rownames(beta_vb)[match(rmvd_coll_x, initial_colnames_X)] <- rownames(beta_vb_coll_var)
  colnames(beta_vb) <- colnames(beta_vb_wo_coll)
  
  theta_vb[match(rmvd_coll_x, initial_colnames_X)] <- theta_vb_coll_var
  names(theta_vb)[match(rmvd_coll_x, initial_colnames_X)] <- names(theta_vb_coll_var)
  
  create_named_list_(beta_vb, gam_vb, theta_vb)
  
}


