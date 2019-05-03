# This file is part of the `atlasqtl` R package:
#     https://github.com/hruffieux/atlasqtl
#

#' Gather model hyperparameters provided by the user.
#'
#' This function must be used to provide hyperparameter values for the model
#' used in \code{\link{atlasqtl}}.
#'
#' The \code{\link{atlasqtl}} function can also be used with default hyperparameter
#' choices (without using \code{\link{set_hyper}}) by setting its argument
#' \code{list_hyper} to \code{NULL}.
#'
#' @param q Number of responses.
#' @param p Number of candidate predictors.
#' @param lambda Vector of length 1 providing the values of hyperparameter
#'   \eqn{\lambda} for the prior distribution of \eqn{\sigma^{-2}}. \eqn{\sigma}
#'   represents the typical size of nonzero effects.
#' @param nu Vector of length 1 providing the values of hyperparameter \eqn{\nu}
#'   for the prior distribution of \eqn{\sigma^{-2}}. \eqn{\sigma} represents
#'   the typical size of nonzero effects.
#' @param eta Vector of length 1 or q. Provides the values of
#'   hyperparameter \eqn{\eta} for the prior distributions of the response 
#'   residual precisions, \eqn{\tau}. If of length 1, the provided
#'   value is repeated q times. 
#' @param kappa Vector of length 1 or q. Provides the values of hyperparameter
#'   \eqn{\kappa} for the prior distributions of the response residual
#'   precisions, \eqn{\tau}. If of length 1, the provided value is repeated q 
#'   times. 
#' @param n0 Vector of length 1 or q, prior mean for the response effects.
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
#' n <- 200; p <- 200; p_act <- 20; q <- 20; q_act <- 15
#'
#' ## Candidate predictors (subject to selection)
#' ##
#' # Here we simulate common genetic variants (but any type of candidate
#' # predictors can be supplied).
#' # 0 = homozygous, major allele, 1 = heterozygous, 2 = homozygous, minor allele
#'
#' X_act <- matrix(rbinom(n * p_act, size = 2, p = 0.25), nrow = n)
#' X_inact <- matrix(rbinom(n * (p - p_act), size = 2, p = 0.25), nrow = n)
#'
#' shuff_x_ind <- sample(p)
#' X <- cbind(X_act, X_inact)[, shuff_x_ind]
#'
#' bool_x_act <- shuff_x_ind <= p_act
#'
#' pat_act <- beta <- matrix(0, nrow = p_act, ncol = q_act)
#' pat_act[sample(p_act*q_act, floor(p_act*q_act/5))] <- 1
#' beta[as.logical(pat_act)] <-  rnorm(sum(pat_act))
#'
#' ## Gaussian responses
#' ##
#' Y_act <- matrix(rnorm(n * q_act, mean = X_act %*% beta, sd = 0.5), nrow = n)
#' Y_inact <- matrix(rnorm(n * (q - q_act), sd = 0.5), nrow = n)
#' shuff_y_ind <- sample(q)
#' Y <- cbind(Y_act, Y_inact)[, shuff_y_ind]
#'
#' ########################
#' ## Infer associations ##
#' ########################
#' 
#' list_hyper <- set_hyper(q, p, lambda = 1, nu = 1, eta = 1, kappa = 1, 
#'                         n0 = -2.5, t02 = 0.1)
#'
#' vb <- atlasqtl(Y = Y, X = X, p0 = c(5, 25), list_hyper = list_hyper, 
#'                user_seed = seed)
#'
#' @seealso  \code{\link{set_init}}, \code{\link{atlasqtl}}
#'
#' @export
#'
set_hyper <- function(q, p, lambda, nu, eta, kappa, n0, t02) {

  check_structure_(q, "vector", "numeric", 1)
  check_natural_(q)

  check_structure_(p, "vector", "numeric", 1)
  check_natural_(p)

  check_structure_(n0, "vector", "double", c(1, q))
  if (length(n0) == 1) n0 <- rep(n0, q)

  check_structure_(t02, "vector", "double", 1)
  check_positive_(t02)

  check_structure_(lambda, "vector", "double", 1)
  check_positive_(lambda)

  check_structure_(nu, "vector", "double", 1)
  check_positive_(nu)

  check_structure_(eta, "vector", "double", c(1, q))
  check_positive_(eta)
  if (length(eta) == 1) eta <- rep(eta, q)

  check_structure_(kappa, "vector", "double", c(1, q))
  check_positive_(kappa)
  if (length(kappa) == 1) kappa <- rep(kappa, q)

  q_hyper <- q
  p_hyper <- p

  list_hyper <- create_named_list_(q_hyper, p_hyper, eta, kappa, lambda, nu, n0, t02)

  class(list_hyper) <- "hyper"

  list_hyper

}


# Internal function setting default model hyperparameters when not provided by
# the user.
#
auto_set_hyper_ <- function(Y, p, p0) {

  q <- ncol(Y)

  lambda <- 1e-2
  nu <- 1

  # hyperparameter set using the data Y
  eta <- 1 / median(apply(Y, 2, var)) #median to be consistent when doing permutations
  if (!is.finite(eta)) eta <- 1e3
  eta <- rep(eta, q)
  kappa <- rep(1, q)

  E_p_t <- p0[1]
  V_p_t <- p0[2]
  
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
                  "of the number of active predictors per responses supplied in p0.",
                  "Please change p0."))
    })

  # n0 sets the level of sparsity.
  n0 <- - get_mu(E_p_t, t02, p)
  n0 <- rep(n0, q)

  check_positive_(t02)

  q_hyper <- q
  p_hyper <- p
  
  list_hyper <- create_named_list_(q_hyper, p_hyper, eta, kappa, lambda, nu, n0, t02)

  class(list_hyper) <- "out_hyper"

  list_hyper

}

#' Gather initial variational parameters provided by the user.
#'
#' This function must be used to provide initial values for the variational
#' parameters used in \code{\link{atlasqtl}}.
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
#'   variational parameter yielding regression coefficient estimates for
#'   predictor-response pairs included in the model.
#' @param sig2_beta_vb Vector of length q. For these values are the same
#'   for all the predictors (as a result of the predictor variables being
#'   standardized before running the variational algorithm). 
#' @param tau_vb Vector of length q with initial values for the variational 
#'   parameter yielding estimates for the continuous response residual 
#'   precisions. 
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
#' ## Examples using small problem sizes:
#' ##
#' n <- 200; p <- 200; p_act <- 20; q <- 20; q_act <- 15
#'
#' ## Candidate predictors (subject to selection)
#' ##
#' # Here we simulate common genetic variants (but any type of candidate
#' # predictors can be supplied).
#' # 0 = homozygous, major allele, 1 = heterozygous, 2 = homozygous, minor allele
#'
#' X_act <- matrix(rbinom(n * p_act, size = 2, p = 0.25), nrow = n)
#' X_inact <- matrix(rbinom(n * (p - p_act), size = 2, p = 0.25), nrow = n)
#'
#' shuff_x_ind <- sample(p)
#' X <- cbind(X_act, X_inact)[, shuff_x_ind]
#'
#' bool_x_act <- shuff_x_ind <= p_act
#'
#' pat_act <- beta <- matrix(0, nrow = p_act, ncol = q_act)
#' pat_act[sample(p_act*q_act, floor(p_act*q_act/5))] <- 1
#' beta[as.logical(pat_act)] <-  rnorm(sum(pat_act))
#'
#' ## Gaussian responses
#' ##
#' Y_act <- matrix(rnorm(n * q_act, mean = X_act %*% beta, sd = 0.5), nrow = n)
#' Y_inact <- matrix(rnorm(n * (q - q_act), sd = 0.5), nrow = n)
#' shuff_y_ind <- sample(q)
#' Y <- cbind(Y_act, Y_inact)[, shuff_y_ind] 
#'
#'
#' ########################
#' ## Infer associations ##
#' ########################
#'
#' ## Continuous responses
#' ##
#'
#' # gam_vb chosen so that the prior mean number of responses associated with
#' # each candidate predictor is 1/4.
#' gam_vb <- matrix(rbeta(p * q, shape1 = 1, shape2 = 4*q-1), nrow = p)
#' mu_beta_vb <- matrix(rnorm(p * q), nrow = p)
#' tau_vb <- 1 / apply(Y, 2, var)
#' sig2_beta_vb <- 1 / rgamma(q, shape = 2, rate = 1 / tau_vb)
#'
#' list_init <- set_init(q, p, gam_vb, mu_beta_vb, sig2_beta_vb, tau_vb)
#' 
#' vb_g <- atlasqtl(Y = Y, X = X, p0 = c(5, 25), list_init = list_init)
#'
#' @seealso  \code{\link{set_hyper}}, \code{\link{atlasqtl}}
#'
#' @export
#'
set_init <- function(q, p, gam_vb, mu_beta_vb, sig2_beta_vb, tau_vb) {

  check_structure_(q, "vector", "numeric", 1)
  check_natural_(q)

  check_structure_(p, "vector", "numeric", 1)
  check_natural_(p)

  check_structure_(gam_vb, "matrix", "double", c(p, q))
  check_zero_one_(gam_vb)

  check_structure_(mu_beta_vb, "matrix", "double", c(p, q))

  check_structure_(sig2_beta_vb, "vector", "double", q)
 
  check_structure_(tau_vb, "vector", "double", q)
  check_positive_(tau_vb)

  check_positive_(sig2_beta_vb)

  q_init <- q
  p_init <- p
  
  list_init <- create_named_list_(q_init, p_init, gam_vb, mu_beta_vb,
                                  sig2_beta_vb, tau_vb)

  class(list_init) <- "init"

  list_init
}


# Internal function setting default starting values when not provided by the user.
#
auto_set_init_ <- function(Y, p, p0, user_seed) {

  # Initialisation not modified for dual = TRUE (should not matter, but maybe change this)

  q <- ncol(Y)

  if (!is.null(user_seed)) set.seed(user_seed)

  E_p_t <- p0[1]
  V_p_t <- p0[2]

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
                  "of the number of active predictors per responses supplied in p0.",
                  "Please change p0."))
    })


  # n0 sets the level of sparsity.
  n0 <- - get_mu(E_p_t, t02, p)

  # Look at : gam_st
  #
  s02 <- 1e-4

  check_positive_(t02)

  gam_vb <- matrix(pnorm(rnorm(p * q, mean = n0, sd = s02 + t02)), nrow = p)                                            

  mu_beta_vb <- matrix(rnorm(p * q), nrow = p)

  sig2_inv_vb <- 1e-2

  tau_vb <- 1 / median(apply(Y, 2, var))
  if (!is.finite(tau_vb)) tau_vb <- 1e3
  tau_vb <- rep(tau_vb, q)

  sig2_beta_vb <- 1 / rgamma(q, shape = 2, rate = 1 / (sig2_inv_vb * tau_vb))

  q_init <- q
  p_init <- p
 
  list_init <- create_named_list_(q_init, p_init, gam_vb, mu_beta_vb,
                                  sig2_inv_vb, sig2_beta_vb, tau_vb)

  class(list_init) <- "out_init"

  list_init
  
}
