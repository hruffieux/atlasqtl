#' Print function for S3 class "atlasqtl"
#' 
#' This function provides functionality for printing a summary of the atlasqtl
#' run. 
#' 
#' @param x an object with S3 class \code{"atlasqtl"}
#' @param ... further arguments passed to or from other methods.
#' 
#' @seealso \code{\link{atlasqtl}}
#' 
#' @export 
#' 
print.atlasqtl <- function(x, ...) {

  if (x$converged) {
    cat(paste0("****************************************************** \n", 
               "Successful convergence after ", x$it, " iterations, using a\n", 
               "tolerance of ", x$tol, " on the absolute changes in the ELBO\n", 
               "****************************************************** \n\n"))
    
    if (!is.null(x$anneal)) {
      
      if (x$anneal[1] == 1) {
        anneal_type <- "Geometric"
      } else if (x$anneal[1] == 2) {
        anneal_type <- "Harmonic"
      } else if (x$anneal[1] == 3) {
        anneal_type <- "Linear"
      }
      cat(paste0(anneal_type, " annealing on the inverse temperature was\n", 
                 "applied for the first ", x$anneal[3], " iterations, with initial\n",
                 "temperature of ", x$anneal[2],
                 ifelse(all(x$anneal == c(1, 2, 10)), " (default).\n", ".\n")))
    }
    
  } else {
    cat(paste0("************************************************ \n", 
               "Unsuccessful convergence after ", x$maxit, " iterations. \n",
               "Difference between last two consecutive values\n", 
               "of the ELBO: ", format(x$diff_lb, digits = 3), ".\n\n",
               "Try increasing:\n", 
               "- the maximum number of iterations (maxit) or\n", 
               "- the convergence threshold (tol). \n",
               "************************************************ \n\n"))
  }
}

# do summary and plot functions.