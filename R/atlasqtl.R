# This file is part of the `atlasqtl` R package:
#     https://github.com/hruffieux/atlasqtl
#

#' Fit a flexible hierarchical model for hotspot detection using annealed 
#' variational inference. 
#'
#' @param Y Response data matrix of dimension n x q, where n is the number of
#'   samples and q is the number of response variables.
#' @param X Input matrix of dimension n x p, where p is the number of candidate
#'   predictors. \code{X} cannot contain NAs. No intercept must be supplied.
#' @param p0 Vector of size 2 whose arguments are the prior expectation and 
#'   variance of the number of predictors associated with each response.
#'   Must be \code{NULL} if \code{list_init} and \code{list_hyper}
#'   are both non-\code{NULL}.
#' @param list_hyper An object of class "\code{hyper}" containing the model
#'   hyperparameters. Must be constructed using the \code{\link{set_hyper}}
#'   function or set to \code{NULL} for default hyperparameters.
#' @param list_init An object of class "\code{init}" containing the initial
#'   variational parameters. Must be constructed using the 
#'   \code{\link{set_init}} function or set to \code{NULL} for a default 
#'   initialization.
#' @param user_seed Seed set for reproducible default choices of hyperparameters
#'   (if \code{list_hyper} is \code{NULL}) and initial variational parameters
#'   (if \code{list_init} is \code{NULL}). Default is \code{NULL}, no
#'   seed set.
#' @param tol Tolerance for the stopping criterion (default is 0.1).
#' @param maxit Maximum number of iterations allowed (default is 1000).
#' @param anneal Parameters for annealing scheme. Must be a vector whose first
#'   entry is the type of schedule: 1 = geometric spacing (default), 
#'   2 = harmonic spacing or 3 = linear spacing, the second entry is the initial 
#'   temperature (default is 2), and the third entry is the temperature grid 
#'   size (default is 10). If \code{NULL}, no annealing is performed.
#' @param save_hyper If \code{TRUE}, the hyperparameters used for the model are
#'   saved as output.
#' @param save_init If \code{TRUE}, the initial variational parameters used for
#'   the inference are saved as output.
#' @param verbose If \code{TRUE}, messages are displayed during execution.
#' @param checkpoint_path Path where to save temporary checkpoint outputs. 
#'   Default is \code{NULL}, for no checkpointing.
#' @param trace_path Path where to save trace plot for the variance of hotspot
#'   propensities. Default is \code{NULL}, for no trace saved.
#'   
#' @details The columns of \code{Y} are centered before running the variational 
#' algorithm, and the columns of \code{X} are standardized.
#'
#' @return An object of class "\code{vb}" containing the following variational
#'   estimates and settings:
#'  \item{gam_vb}{Posterior inclusion probability matrix of dimension p x q.
#'                Entry (s, t) corresponds to the posterior probability of
#'                association between candidate predictor s and response t.}
#'  \item{mu_theta_vb}{Vector of length p containing the posterior mean of
#'                     theta. Entry s corresponds to the propensity of candidate
#'                     predictor s to be included in the model. \code{NULL} if
#'                     \code{dual} is \code{FALSE}.}
#'  \item{converged}{A boolean indicating whether the algorithm has converged
#'                   before reaching \code{maxit} iterations.}
#'  \item{it}{Final number of iterations.}
#'  \item{lb_opt}{Optimized variational lower bound for the marginal
#'                log-likelihood (ELBO).}
#'  \item{diff_lb}{Difference in ELBO between the last and penultimate
#'                 iterations. Useful to convergence diagnostic when 
#'                 \code{maxit} has been reached.}
#'  \item{rmvd_cst_x}{Vectors containing the indices of constant variables in 
#'                    \code{X} removed prior to the analysis.}
#'  \item{rmvd_coll_x}{Vectors containing the indices of variables in \code{X} 
#'                     removed prior to the analysis because collinear to other
#'                     variables. The entry name indicates the corresponding 
#'                     variable kept in the analysis (i.e., that causing the 
#'                     collinearity for the entry in question).}
#'  \item{list_hyper, list_init}{If \code{save_hyper}, resp. \code{save_init},
#'                               \code{TRUE}, hyperparameters, resp. initial
#'                               variational parameters, used for inference are
#'                               saved as output.}
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
#' n <- 200; p <- 100; p_act <- 10; q <- 100; q_act <- 50
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
#' ########################
#' ## Infer associations ##
#' ########################
#'
#' # Expectation and variance for the prior number of predictors associated with
#' # each response
#' #
#' p0 <- c(mean(colSums(pat)), 10)
#' 
#' res_atlas <- atlasqtl(Y = Y, X = X, p0 = p0, user_seed = seed)
#'
#' @references
#' Helene Ruffieux, Anthony C. Davison, Jorg Hager, Jamie Inshaw, Benjamin P. 
#' Fairfax, Sylvia Richardson, Leonardo Bottolo, A global-local approach for 
#' detecting hotspots in multiple-response regression, arXiv:1811.03334, 2018.
#'
#' @seealso \code{\link{set_hyper}}, \code{\link{set_init}}.
#'
#' @export
#'
atlasqtl <- function(Y, X, p0, list_hyper = NULL, list_init = NULL,
                     user_seed = NULL, tol = 1e-1, maxit = 1000, 
                     anneal = c(1, 2, 10), save_hyper = FALSE, save_init = FALSE, 
                     verbose = TRUE, checkpoint_path = NULL, trace_path = NULL) {
  
  if (verbose) cat("== Preparing the data ... \n")
  
  check_annealing_(anneal)
  
  dat <- prepare_data_(Y, X, user_seed, tol, maxit, verbose, checkpoint_path, 
                       trace_path)
  
  bool_rmvd_x <- dat$bool_rmvd_x
  
  X <- dat$X
  Y <- dat$Y
  
  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(Y)
  
  names_x <- colnames(X)
  names_y <- colnames(Y)
  
  if (verbose) cat("... done. == \n\n")

  if (is.null(list_hyper) | is.null(list_init)) {
    
    check_structure_(p0, "vector", "numeric", 2)
    check_positive_(p0)
    
  } else {
    
    if (!is.null(p0))
      warning(paste("Provided argument p0 not used, as both list_hyper ",
                    "and list_init were provided.", sep = ""))
    
  }
  
  
  
  if (verbose) cat("== Preparing the hyperparameters ... \n\n")
  
  list_hyper <- prepare_list_hyper_(list_hyper, Y, p, p0, bool_rmvd_x, 
                                    names_x, names_y, verbose)
  
  if (verbose) cat("... done. == \n\n")
  
  if (verbose) cat("== Preparing the parameter initialization ... \n\n")
  
  list_init <- prepare_list_init_(list_init, Y, p, p0, bool_rmvd_x, 
                                  user_seed, verbose)
  
  if (verbose) cat("... done. == \n\n")
  
  if (verbose){
    cat(paste("============================================================== \n",
              "== Variational inference for sparse multivariate regression == \n",
              "============================================================== \n\n",
              sep = ""))
  }
  
  
  hs <- TRUE
    if (hs) {
      
      df <- 1
      vb <- atlasqtl_dual_horseshoe_core_(Y, X, list_hyper, list_init$gam_vb,
                                          list_init$mu_beta_vb, 
                                          list_init$sig2_beta_vb,
                                          list_init$tau_vb, df, tol, maxit, 
                                          anneal, verbose,
                                          checkpoint_path = checkpoint_path,
                                          trace_path = trace_path)
      
    } else {
      
      vb <- atlasqtl_dual_prior_core_(Y, X, list_hyper, list_init$gam_vb,
                                      list_init$mu_beta_vb, list_init$sig2_beta_vb,
                                      list_init$tau_vb, tol, maxit, anneal, 
                                      verbose, checkpoint_path = checkpoint_path)
      
    }
         
  vb$p0 <- p0
  
  vb$rmvd_cst_x <- dat$rmvd_cst_x
  vb$rmvd_coll_x <- dat$rmvd_coll_x
 
  if (save_hyper) vb$list_hyper <- list_hyper
  if (save_init) vb$list_init <- list_init
  
  class(vb) <- "atlasqtl"
  
  if (verbose) cat("... done. == \n\n")
  
  vb
  
}
