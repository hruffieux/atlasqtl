# This file is part of the `atlasqtl` R package:
#     https://github.com/hruffieux/atlasqtl
#

#' Gather model hyperparameters provided by the user.
#'
#' This function is used to provide hyperparameter values for 
#' \code{\link{atlasqtl}}.
#'
#' The \code{\link{atlasqtl}} function can also be used with default 
#' hyperparameter values (i.e., without using \code{\link{set_hyper}}) by 
#' setting its argument \code{list_hyper} to \code{NULL}.
#'
#' @param q Number of responses.
#' @param p Number of candidate predictors.
#' @param eta Vector of length 1 or q. Provides the shape hyperparameter 
#'   \eqn{\eta} for the gamma prior distribution of the response residual 
#'   precision, \eqn{\tau}. If of length 1, the provided value is repeated q 
#'   times. 
#' @param kappa Vector of length 1 or q. Provides the rate hyperparameter
#'   \eqn{\kappa} for the gamma prior distribution of the response residual
#'   precision, \eqn{\tau}. If of length 1, the provided value is repeated q 
#'   times. 
#' @param n0 Vector of length 1 or q, prior mean for the response effects.
#' @param nu Vector of length 1 providing the shape hyperparameter \eqn{\nu} for 
#'   the gamma prior distribution of \eqn{\sigma^{-2}}. \eqn{\sigma} represents 
#'   the typical size of nonzero effects.
#' @param rho Vector of length 1 providing the rate hyperparameter \eqn{\rho}
#'   for the prior distribution of \eqn{\sigma^{-2}}. \eqn{\sigma} represents
#'   the typical size of nonzero effects.
#' @param t02 Prior variance for the response effects.
#'
#' @return An object of class "\code{hyper}" preparing user hyperparameter in a
#'   form that can be passed to the \code{\link{atlasqtl}} function.
#'
#' @examples
#' seed <- 123; set.seed(seed)
#'
#' ###################
#' ## Simulate data ##
#' ###################
#'
#' ## Examples using small problem sizes:
#' ##
#' n <- 200; p <- 50; p_act <- 10; q <- 100; q_act <- 50
#' 
#' # Candidate predictors (subject to selection)
#' #
#' # Here example with common genetic variants under Hardy-Weinberg equilibrium
#' #
#' X_act <- matrix(rbinom(n * p_act, size = 2, p = 0.25), nrow = n)
#' X_inact <- matrix(rbinom(n * (p - p_act), size = 2, p = 0.25), nrow = n)
#'
#' # shuffle indices 
#' shuff_x_ind <- sample(p)
#' shuff_y_ind <- sample(q)
#' 
#' X <- cbind(X_act, X_inact)[, shuff_x_ind]
#'
#' # Association pattern and effect sizes
#' #
#' pat <- matrix(FALSE, ncol = q, nrow = p)
#' bool_x <- shuff_x_ind <= p_act
#' bool_y <- shuff_y_ind <= q_act
#' 
#' pat_act <- beta_act <- matrix(0, nrow = p_act, ncol = q_act)
#' pat_act[sample(p_act * q_act, floor(p_act * q_act / 5))] <- 1
#' beta_act[as.logical(pat_act)] <-  rnorm(sum(pat_act))
#' 
#' pat[bool_x, bool_y] <- pat_act
#'
#' # Gaussian responses
#' #
#' Y_act <- matrix(rnorm(n * q_act, mean = X_act %*% beta_act), nrow = n)
#' Y_inact <- matrix(rnorm(n * (q - q_act)), nrow = n)
#'
#' Y <- cbind(Y_act, Y_inact)[, shuff_y_ind]
#' 
#' #############################
#' ## Specify hyperparameters ##
#' #############################
#' 
#' list_hyper <- set_hyper(q, p, eta = 1, kappa = 1, n0 = -2, nu = 1, rho = 1, 
#'                         t02 = 0.1)
#'                         
#' ########################
#' ## Infer associations ##
#' ########################
#' 
#' p0 <- c(mean(colSums(pat)), 10)
#'
#' res_atlas <- atlasqtl(Y, X, p0, list_hyper = list_hyper, user_seed = seed)
#'
#' @seealso  \code{\link{set_init}}, \code{\link{atlasqtl}}
#'
#' @export
#'
set_hyper <- function(q, p, eta, kappa, n0, nu, rho, t02) {

  check_structure_(q, "vector", "numeric", 1)
  check_natural_(q)

  check_structure_(p, "vector", "numeric", 1)
  check_natural_(p)

  check_structure_(n0, "vector", "double", c(1, q))
  if (length(n0) == 1) n0 <- rep(n0, q)

  check_structure_(t02, "vector", "double", 1)
  check_positive_(t02)

  check_structure_(nu, "vector", "double", 1)
  check_positive_(nu)

  check_structure_(rho, "vector", "double", 1)
  check_positive_(rho)

  check_structure_(eta, "vector", "double", c(1, q))
  check_positive_(eta)
  if (length(eta) == 1) eta <- rep(eta, q)

  check_structure_(kappa, "vector", "double", c(1, q))
  check_positive_(kappa)
  if (length(kappa) == 1) kappa <- rep(kappa, q)

  m0 <- 0 # prior mean of horseshoe distribution is 0
  
  A2_inv <- 1 # horseshoe global scale hyperprior 
  
  q_hyper <- q
  p_hyper <- p

  list_hyper <- create_named_list_(q_hyper, p_hyper, A2_inv, eta, kappa, m0, n0, 
                                   nu, rho, t02)

  class(list_hyper) <- "hyper"

  list_hyper

}


# Internal function setting default model hyperparameters when not provided by
# the user.
#
auto_set_hyper_ <- function(Y, p, p0) {

  q <- ncol(Y)

  nu <- 1e-2
  rho <- 1

  eta <- 1 / median(apply(Y, 2, var, na.rm = TRUE)) 
  if (!is.finite(eta)) eta <- 1e3
  eta <- rep(eta, q)
  kappa <- rep(1, q)

  E_p_t <- p0[1]
  V_p_t <- p0[2]
  
  dn <- 1e-6
  up <- 1e5

  # Get n0 and t02 by specifying a prior expectation and variance for number of 
  # predictors associated with each response
  #
  tryCatch(t02 <- uniroot(function(x)
    get_V_p_t(get_mu(E_p_t, x, p), x, p) - V_p_t,
    interval = c(dn, up))$root,
    error = function(e) {
      stop(paste0("No hyperparameter values matching the expectation and variance ",
                  "of the number of active predictors per responses supplied in p0.",
                  "Please change p0."))
    })

  # n0 sets the level of sparsity.
  #
  n0 <- get_mu(E_p_t, t02, p)
  n0 <- rep(n0, q)

  check_positive_(t02)
  
  m0 <- 0     # prior mean of horseshoe distribution is 0
  A2_inv <- 1 # horseshoe global scale hyperprior 
  

  q_hyper <- q
  p_hyper <- p
  
  list_hyper <- create_named_list_(q_hyper, p_hyper, A2_inv, eta, kappa, m0, n0, 
                                   nu, rho, t02)

  class(list_hyper) <- "out_hyper"

  list_hyper

}

#' Gather initial variational parameters provided by the user.
#'
#' This function must be used to provide initial values for the variational
#' parameters by \code{\link{atlasqtl}}.
#'
#' The \code{\link{atlasqtl}} function can also be used with default initial
#' parameter choices (without using \code{\link{set_init}}) by setting
#' its argument \code{list_init} to \code{NULL}.
#'
#' @param q Number of responses.
#' @param p Number of candidate predictors.
#' @param gam_vb Matrix of size p x q with initial values for the variational
#'   parameter yielding posterior probabilities of inclusion.
#' @param mu_beta_vb Matrix of size p x q with initial values for the
#'   variational parameter yielding the posterior mean of regression estimates 
#'   for predictor-response pairs included in the model.
#' @param sig02_inv_vb Initial parameter for the hotspot propensity global 
#'   precision.
#' @param sig2_beta_vb Vector of length q. Initial values for the variational 
#' parameter yielding the posterior variance of regression estimates for
#'   predictor-response pairs included in the model.
#' @param sig2_theta_vb Vector of length p. Initial values for the variational 
#'   parameter yielding the posterior variance of the hotspot propensities.
#' @param tau_vb Vector of length q with initial values for the variational 
#'   parameter yielding estimates for the response residual precisions. 
#' @param theta_vb Vector of length p. Initial values for the variational 
#'   parameter yielding the posterior mean of the hotspot propensities.
#' @param zeta_vb Vector of length q. Initial values for the variational 
#'   parameter yielding the posterior mean of the response propensities.
#'
#' @return An object of class "\code{init}" preparing user initial values for
#'   the variational parameters in a form that can be passed to the
#'   \code{\link{atlasqtl}} function.
#'
#' @examples
#' seed <- 123; set.seed(seed)
#'
#' ###################
#' ## Simulate data ##
#' ###################
#'
#' # Example with small problem sizes:
#' #
#' n <- 200; p <- 50; p_act <- 10; q <- 100; q_act <- 50
#'
#' # Candidate predictors (subject to selection)
#' #
#' # Here example with common genetic variants under Hardy-Weinberg equilibrium
#' #
#' X_act <- matrix(rbinom(n * p_act, size = 2, p = 0.25), nrow = n)
#' X_inact <- matrix(rbinom(n * (p - p_act), size = 2, p = 0.25), nrow = n)
#'
#' # shuffle indices 
#' shuff_x_ind <- sample(p)
#' shuff_y_ind <- sample(q)
#' 
#' X <- cbind(X_act, X_inact)[, shuff_x_ind]
#'
#' # Association pattern and effect sizes
#' #
#' pat <- matrix(FALSE, ncol = q, nrow = p)
#' bool_x <- shuff_x_ind <= p_act
#' bool_y <- shuff_y_ind <= q_act
#' 
#' pat_act <- beta_act <- matrix(0, nrow = p_act, ncol = q_act)
#' pat_act[sample(p_act * q_act, floor(p_act * q_act / 5))] <- 1
#' beta_act[as.logical(pat_act)] <-  rnorm(sum(pat_act))
#' 
#' pat[bool_x, bool_y] <- pat_act
#'
#' # Gaussian responses
#' #
#' Y_act <- matrix(rnorm(n * q_act, mean = X_act %*% beta_act), nrow = n)
#' Y_inact <- matrix(rnorm(n * (q - q_act)), nrow = n)
#'
#' Y <- cbind(Y_act, Y_inact)[, shuff_y_ind]
#'
#'
#' ################################
#' ## Specify initial parameters ##
#' ################################
#' 
#' tau_vb <- rep(1, q)
#' 
#' gam_vb <- matrix(rbeta(p * q, shape1 = 1, shape2 = 4 * q - 1), nrow = p)
#' 
#' mu_beta_vb <- matrix(rnorm(p * q), nrow = p)
#' sig2_beta_vb <- 1 / rgamma(q, shape = 2, rate = 1)
#' 
#' sig02_inv_vb <- rgamma(1, shape = max(p, q), rate = 1)
#' 
#' theta_vb <- rnorm(p, sd = 1 / sqrt(sig02_inv_vb * q))
#' sig2_theta_vb <- 1 / (q + rgamma(p, shape = sig02_inv_vb * q, rate = 1))
#' 
#' zeta_vb <- rnorm(q, mean = -2, sd = 0.1)
#' 
#' list_init <- set_init(q, p, gam_vb, mu_beta_vb, sig02_inv_vb, sig2_beta_vb, 
#'                       sig2_theta_vb, tau_vb, theta_vb, zeta_vb)
#'                       
#'                       
#' ########################
#' ## Infer associations ##
#' ########################
#' 
#' p0 <- c(mean(colSums(pat)), 10)
#' 
#' res_atlas <- atlasqtl(Y, X, p0, list_init = list_init)
#'
#' @seealso  \code{\link{set_hyper}}, \code{\link{atlasqtl}}
#'
#' @export
#'
set_init <- function(q, p, gam_vb, mu_beta_vb, sig02_inv_vb, sig2_beta_vb, 
                     sig2_theta_vb, tau_vb, theta_vb, zeta_vb) {

  check_structure_(q, "vector", "numeric", 1)
  check_natural_(q)

  check_structure_(p, "vector", "numeric", 1)
  check_natural_(p)

  check_structure_(gam_vb, "matrix", "double", c(p, q))
  check_zero_one_(gam_vb)

  check_structure_(mu_beta_vb, "matrix", "double", c(p, q))

  check_structure_(sig02_inv_vb, "vector", "numeric", 1)
  check_positive_(sig02_inv_vb)
  
  check_structure_(sig2_beta_vb, "vector", "double", q)
  check_positive_(sig2_beta_vb)
  
  check_structure_(sig2_theta_vb, "vector", "double", p)
  check_positive_(sig2_theta_vb)
  
  check_structure_(tau_vb, "vector", "double", q)
  check_positive_(tau_vb)
  
  check_structure_(theta_vb, "vector", "double", p)
  
  check_structure_(zeta_vb, "vector", "double", q)
  
  q_init <- q
  p_init <- p
  
  list_init <- create_named_list_(q_init, p_init, gam_vb, mu_beta_vb,
                                  sig02_inv_vb, sig2_beta_vb, sig2_theta_vb, 
                                  tau_vb, theta_vb, zeta_vb)
  
  class(list_init) <- "init"

  list_init
}


# Internal function setting default starting values when not provided by the user.
#
auto_set_init_ <- function(Y, p, p0, shr_fac_inv, user_seed) {

  q <- ncol(Y)

  if (!is.null(user_seed)) set.seed(user_seed)

  E_p_t <- p0[1]
  V_p_t <- p0[2]

  dn <- 1e-6
  up <- 1e5

  # Get n0 and t02 by specifying a prior expectation and variance for number of 
  # predictors associated with each response
  #
  tryCatch(t02 <- uniroot(function(x)
    get_V_p_t(get_mu(E_p_t, x, p), x, p) - V_p_t,
    interval = c(dn, up))$root,
    error = function(e) {
      stop(paste0("No hyperparameter values matching the expectation and variance ",
                  "of the number of active predictors per responses supplied in p0.",
                  "Please change p0."))
    })

  # n0 sets the level of sparsity.
  n0 <- get_mu(E_p_t, t02, p)

  # Look at : gam_st
  #
  s02 <- 1e-4
  check_positive_(t02)

  gam_vb <- matrix(pnorm(rnorm(p * q, mean = n0, sd = s02 + t02)), nrow = p)                                            

  mu_beta_vb <- matrix(rnorm(p * q), nrow = p)

  sig2_inv_vb <- 1e-2

  tau_vb <- 1 / median(apply(Y, 2, var, na.rm = TRUE))
  if (!is.finite(tau_vb)) tau_vb <- 1e3
  tau_vb <- rep(tau_vb, q)

  sig2_beta_vb <- 1 / rgamma(q, shape = 2, rate = 1 / (sig2_inv_vb * tau_vb))
  
  sig02_inv_vb <- rgamma(1, shape = max(p, q), rate = 1)

  theta_vb <- rnorm(p, sd = 1 / sqrt(sig02_inv_vb * shr_fac_inv))
  sig2_theta_vb <- 1 / (q + rgamma(p, shape = sig02_inv_vb * shr_fac_inv, rate = 1)) # initial guess assuming lam2_inv_vb = 1
  
  zeta_vb <- rnorm(q, mean = n0, sd = sqrt(t02))
  
  q_init <- q
  p_init <- p
 
  list_init <- create_named_list_(q_init, p_init, gam_vb, mu_beta_vb,
                                  sig02_inv_vb, sig2_beta_vb, sig2_theta_vb, 
                                  tau_vb, theta_vb, zeta_vb)

  class(list_init) <- "out_init"

  list_init
  
}



#' Evaluate the approximation in the hyperprior elicitation for the number 
#' of predictors associated with each response
#'
#' This function is used to guide hyperprior elicitation before running
#' \code{\link{atlasqtl}}. The errors arising from the approximation in the 
#' hyper prior elicitation are estimated by Monte-Carlo simulation.
#'
#' @param p0 Vector of size 2 whose entries are the prior expectation and 
#'   variance of the number of predictors associated with each response.
#' @param p Number of candidate predictors.
#' @param q Number of responses.
#' @param n_draws Number of draws used for Monte-Carlo simulation (default 1e5).
#' @param n_cpus Number of cores used (default serial).
#'
#' @return A list of errors for the prior mean and prior standard deviation of 
#'   the number of predictors associated with each response. 
#'   
#' @examples
#' seed <- 123; set.seed(seed)
#
#' n <- 200; p <- 100; q <- 20000
#' p0 <- c(1, 10)
#' n_draws <- 1e5 # must be large for accurate estimation
#' 
#' map_hyperprior_elicitation(p0, p, q, n_draws, n_cpus = 1)
#'
#' @seealso \code{\link{atlasqtl}}
#'
#' @export
#'
map_hyperprior_elicitation <- function(p0, p, q, n_draws = 1e5, n_cpus = 1) {
  
  check_structure_(p0, "vector", "numeric", 2)
  E_p <- p0[1]
  V_p <- p0[2]
  check_positive_(E_p)
  check_positive_(V_p)
  
  check_natural_(p)
  check_natural_(q)
  
  check_natural_(n_draws)
  
  if (n_draws < 1e3) {
    warning("The number of draws may be too small for accurate Monte Carlo estimation.")
  }
  
  if (E_p > p) {
    stop(paste0("The prior mean number of predictors associated with each response, ", 
                "E_p, must be smaller than the total number of candidate predictors, p."))
  }
  
  # Only used for throwing a message if E_p and V_p are incompatible
  #
  list_n0_t02 <- get_n0_t02(1, p, c(E_p, V_p))
  n0 <- list_n0_t02$n0
  t02 <- list_n0_t02$t02
  
  if (n_cpus > parallel::detectCores()) {
    stop("The number of CPUs must be smaller than the number of CPUs on the machine.")
  }
  
  lambda <- LaplacesDemon::rhalfcauchy(n_draws, scale = 1)
  sigma0 <- LaplacesDemon::rhalfcauchy(n_draws, scale = 1 / sqrt(q))
  
  E_Phi_X_hs <- pnorm(n0 / sqrt(1 + t02 + lambda^2 * sigma0^2))
  E_Phi_X_2_hs <- E_Phi_X_hs - 2 * 
    unlist(parallel::mclapply(1:n_draws, function(ii) { 
      PowerTOST::OwensT(n0 / sqrt(1 + t02 + lambda[ii]^2 * sigma0[ii]^2), 1 / sqrt(1 + 2 * t02 + 2 * lambda[ii]^2 * sigma0[ii]^2))
    }, mc.cores = n_cpus))
  
  E_p_hs <- mean(p * pnorm(n0 / sqrt(1 + t02 + lambda^2 * sigma0^2)))
  V_p_hs <- mean(p * (p - 1) * E_Phi_X_2_hs - p^2 * E_Phi_X_hs^2 + p * E_Phi_X_hs)
  
  error_E_p <- abs(E_p_hs - E_p)
  error_sd_p <- abs(sqrt(V_p_hs) - sqrt(V_p))
  
  create_named_list_(error_E_p, error_sd_p)
  
}
