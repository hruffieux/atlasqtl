# This file is part of the `atlasqtl` R package:
#     https://github.com/hruffieux/atlasqtl
#

# Internal function implementing sanity checks and needed preprocessing before
# the application of the different `atlasqtl_*_core` algorithms.
#
prepare_data_ <- function(Y, X, tol, maxit, user_seed, verbose, checkpoint_path, 
                          trace_path) {

  check_structure_(user_seed, "vector", "numeric", 1, null_ok = TRUE)
  
  check_structure_(tol, "vector", "numeric", 1)
  check_positive_(tol, eps=.Machine$double.eps)

  check_structure_(maxit, "vector", "numeric", 1)
  check_natural_(maxit)

  check_structure_(X, "matrix", "numeric")
  
  if (!is.null(checkpoint_path) && !dir.exists(checkpoint_path))
    stop(paste0("The directory specified in checkpoint_path does not exist. ", 
                "Please make sure to provide a valid path."))
  
  if (!is.null(trace_path) && !dir.exists(trace_path)) 
    stop(paste0("The directory specified in trace_path does not exist. ", 
                "Please make sure to provide a valid path."))
  
  n <- nrow(X)
  p <- ncol(X)

  check_structure_(Y, "matrix", "numeric", na_ok = TRUE)
  
  q <- ncol(Y)
  
  if (nrow(Y) != n) 
    stop("X and Y must have the same number of samples.")
  
  if (sum(!is.na(Y)) / (n * q) < 0.05)
    stop("Too few non-NA values in matrix Y. Exit.")
  
  ind_low <- apply(Y, 2, function(y) sum(!is.na(y)) / n < 0.025)
  if (any(ind_low))
      stop(paste0("Column(s) ", which(ind_low), " of matrix Y have more than ",
                  "97.5% missing values, and should be removed. Exit."))

  if (is.null(rownames(X)) & is.null(rownames(Y)))
    rownames(X) <- rownames(Y) <- paste0("Ind_", 1:n)
  else if (is.null(rownames(X))) rownames(X) <- rownames(Y)
  else if (is.null(rownames(Y))) rownames(Y) <- rownames(X)
  else if (any(rownames(X) != rownames(Y)))
    stop("The provided rownames of X and Y must be the same.")

  if (is.null(colnames(X))) colnames(X) <- paste0("Cov_x_", 1:p)
  if (is.null(colnames(Y))) colnames(Y) <- paste0("Resp_", 1:q)

  X <- scale(X)

  list_X_cst <- rm_constant_(X, verbose)
  X <- list_X_cst$mat
  bool_cst_x <- list_X_cst$bool_cst
  rmvd_cst_x <- list_X_cst$rmvd_cst
  
  initial_colnames_X <- colnames(X) # names (and ordering) of the variables before 
                                    # having removed the collinear variables 
                                    # but after having removed the constant variables

  list_X_coll <- rm_collinear_(X, verbose)
  X <- list_X_coll$mat

  bool_coll_x <- list_X_coll$bool_coll
  rmvd_coll_x <- list_X_coll$rmvd_coll

  bool_rmvd_x <- bool_cst_x
  bool_rmvd_x[!bool_cst_x] <- bool_coll_x

  p <- ncol(X) # get the number of variables after removal of constant or
               # collinear variables 
  
  if (p < 1) stop(paste0("There must be at least 1 non-constant candidate ", 
                         "predictor stored in X."))

  Y <- scale(Y, center = TRUE, scale = FALSE)

  create_named_list_(Y, X, bool_rmvd_x, initial_colnames_X, rmvd_cst_x, rmvd_coll_x) 

}


check_verbose_ <- function(verbose) {
  
  if (!(verbose %in% 0:2))
    stop(paste0("The verbose argument must be set to 0, 1 or 2."))
  
}    

# Internal function implementing sanity checks for the annealing schedule 
# specification.
#
check_annealing_ <- function(anneal) {

  check_structure_(anneal, "vector", "numeric", 3, null_ok = TRUE)

  if (!is.null(anneal)) {

    check_natural_(anneal[c(1, 3)])
    check_positive_(anneal[2])

    if (!(anneal[1] %in% 1:3))
      stop(paste0("The annealing spacing scheme must be set to 1 for geometric ", 
                  "2 for harmonic or 3 for linear spacing."))

    if (anneal[2] < 1.5)
      stop(paste0("Initial annealing temperature very small. May not be large ",
                  "enough for a successful exploration. Please increase it or ", 
                  "select no annealing."))

    if (anneal[3] > 1000)
      stop(paste0("Temperature grid size very large. This may be unnecessarily ",
                  "computationally demanding. Please decrease it."))
    
  }

}


# Internal function implementing sanity checks and needed preprocessing for the
# model hyperparameters before the application of the different `atlasqtl_*_core`
# algorithms.
#
prepare_list_hyper_ <- function(list_hyper, Y, p, p0, bool_rmvd_x, names_x, 
                                names_y, verbose) {

  q <- ncol(Y)
  
  if (is.null(list_hyper)) {

    if (verbose != 0) cat("list_hyper set automatically. \n")

    list_hyper <- auto_set_hyper_(Y, p, p0)

  } else {

    if (!inherits(list_hyper, c("hyper", "out_hyper")))
      stop(paste0("The provided list_hyper must be an object of class ", 
                  "``hyper'' or ``out_hyper''. \n *** you must either use the ", 
                  "function set_hyper to set your own hyperparameters or ", 
                  "list_hyper to NULL for automatic choice. ***"))

    if (inherits(list_hyper, "hyper")) {

      p_hyper_match <- length(bool_rmvd_x)

    } else {
      
      p_hyper_match <- p
    
    }


    if (list_hyper$q_hyper != q)
      stop(paste0("The dimensions (q) of the provided hyperparameters ",
                 "(list_hyper) are not consistent with that of Y.\n"))

    if (list_hyper$p_hyper != p_hyper_match)
      stop(paste0("The dimensions (p) of the provided hyperparameters ",
                 "(list_hyper) are not consistent with that of X.\n"))
    
    if (!is.null(names(list_hyper$eta)) && names(list_hyper$eta) != names_y)
      stop(paste0("Provided names for the entries of eta do not match the ", 
                  "colnames of the continuous variables in Y"))

    if (!is.null(names(list_hyper$kappa)) && names(list_hyper$kappa) != names_y)
      stop(paste0("Provided names for the entries of kappa do not match the ", 
           "colnames of the continuous variables in Y"))

  }

  class(list_hyper) <- "out_hyper"

  list_hyper
}


# Internal function implementing sanity checks and needed preprocessing for the
# starting values before the application of the different `atlasqtl_*_core`
# algorithms.
#
prepare_list_init_ <- function(list_init, Y, p, p0, bool_rmvd_x, shr_fac_inv, 
                               user_seed, verbose) {

  q <- ncol(Y)
  n <- nrow(Y)

  if (is.null(list_init)) {

    if (!is.null(user_seed) & verbose != 0) cat(paste0("Seed set to user_seed ",
                                                 user_seed,". \n"))

    if (verbose != 0) cat(paste0("list_init set automatically. \n"))

    list_init <- auto_set_init_(Y, p, p0, shr_fac_inv, user_seed)

  } else {
    
    if (!is.null(user_seed))
      warning("user_seed not used since a non-NULL list_init was provided. \n")

    if (!inherits(list_init, c("init", "out_init")))
      stop(paste0("The provided list_init must be an object of class ``init'' ", 
                  "or `` out_init''. \n *** you must either use the function ", 
                  "set_init to set your own initialization or set the argument ", 
                  "list_init to NULL for automatic initialization. ***"))

    if (inherits(list_init, "init")) {
      p_init_match <- length(bool_rmvd_x)
    } else {
      p_init_match <- p
    }

    if (inherits(list_init, "init")) {
      p_init_match <- length(bool_rmvd_x)
    } else {
      p_init_match <- p
    }

    if (list_init$q_init != q)
      stop(paste0("The dimensions (q) of the provided initial parameters ",
                 "(list_init) are not consistent with that of Y.\n"))

    if (list_init$p_init != p_init_match)
      stop(paste0("The dimensions (p) of the provided initial parameters ",
                  "(list_init) are not consistent with that of X.\n"))

    if (inherits(list_init, "init")) {

      # drops the rows corresponding to the removed constant and collinear
      # predictors in X (if any)
      list_init$gam_vb <- list_init$gam_vb[!bool_rmvd_x,, drop = FALSE]
      list_init$mu_beta_vb <- list_init$mu_beta_vb[!bool_rmvd_x,, drop = FALSE]

    }

  }

  class(list_init) <- "out_init"

  list_init
}