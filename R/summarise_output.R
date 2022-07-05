#' Print function for S3 class "atlasqtl"
#' 
#' This function provides functionality for printing basic information about the 
#' atlasqtl run. 
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
               "tolerance of ", x$tol, " on the absolute changes in the ELBO.\n", 
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
                 ifelse(all(x$anneal == c(1, 2, 10)), " (default).\n\n", ".\n\n")))
    }
    
    cat(paste0("Number of samples: ", x$n, ";\n",
               "Number of (non-redundant) candidate predictors: ", x$p, ";\n",
               "Number of responses: ", x$q, ";\n",
               "Prior expectation for the number of predictors\n", 
               "associated with each response: ", x$p0[1], " (sd: ", 
               format(sqrt(x$p0[2]), digits = 2), ").\n\n"))
    
    cat(paste0("The posterior quantities inferred by ATLASQTL can\n",
               "be accessed as list elements from the `atlasqtl` S3\n", 
               "object, and a summary can obtained using the\n",  
               "`summary` function.\n\n"))
    
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


#' Summary function for S3 class "atlasqtl"
#' 
#' This function provides functionality for printing a basic summary of the 
#' atlasqtl posterior output.
#' 
#' @param object an object with S3 class \code{"atlasqtl"}
#' @param thres threshold to be applied on the posterior probabilities of 
#'              association (PPI) to define hotspots (default is 0.5, for the 
#'              Median Probability Model).
#' @param fdr_adjust if TRUE the threshold is applied on FDR estimates computed
#'                   from the "atlasqtl" posterior probabilities of association,
#'                   otherwise it is directly applied on the posterior 
#'                   probabilities of association (default is FALSE).
#' @param full_summary Whether or not to print a detailed summary. Default is TRUE.
#' @param ... additional arguments affecting the summary produced.
#' 
#' @seealso \code{\link{atlasqtl}}
#' 
#' @export 
#' 
summary.atlasqtl <- function(object, thres, fdr_adjust = FALSE, full_summary = TRUE, ...) {
  
    cat(paste0("****************************************************** \n", 
               "* ATLASQTL: posterior summary for variable selection *\n", 
               "****************************************************** \n\n"))
    
    if (full_summary) {
      cat("Posterior probabilities pairwise association, pr(gamma_st = 1 | y)\n")
      print(summary(as.vector(object$gam_vb)))
      cat("\nPosterior mean of pairwise regression coefficients, E(beta_st | y)\n")
      print(summary(as.vector(object$beta_vb))) 
      cat("\nPosterior mean of hotspot propensities, E(theta_s | y)\n ")
      print(summary(object$theta_vb)) 
      cat("\n\n")
    }
    
    if (fdr_adjust) {
      mat_fdr <- assign_bFDR(object$gam_vb)
      nb_pairwise <- sum(mat_fdr < thres)
      rs_thres <- rowSums(mat_fdr < thres)
    } else {
      nb_pairwise <- sum(object$gam_vb > thres)
      rs_thres <- rowSums(object$gam_vb > thres)
    }
    cat(paste0("Using ", ifelse(fdr_adjust, "an FDR control of ", "an PPI threshold of "), 
               100*thres, "%:\n"))
    cat("---------------------------\n")
    cat(paste0("\nNb of pairwise (predictor-response) associations: ", nb_pairwise), "\n")
    cat(paste0("\nNb of predictors associated with at least one response \n", 
               "(active predictors): ", sum(rs_thres>0)), "\n")
    cat("\nHotspot sizes (nb of responses associated with each \nactive predictors):\n")
    print(summary(rs_thres[rs_thres>0])) 
    
    if (any(rs_thres>0)) {
      sorted_hotspot_sizes <- sort(rs_thres, decreasing = TRUE)
      top <- min(sum(rs_thres>0), 6)
      
      top_half <- ifelse(top > 3, 3, top) 
      
      cat(paste0("\nTop hotspots: \n", paste0(names(sorted_hotspot_sizes)[1:top_half], 
                                              " (size ", sorted_hotspot_sizes[1:top_half], ")", collapse = ", "), 
                 ifelse(top>3, ", ", ". "), "\n"))
      
      if (top > 3) {
        
        cat(paste0(names(sorted_hotspot_sizes)[(top_half+1):top], 
                                                " (size ", sorted_hotspot_sizes[(top_half+1):top], ")", collapse = ", "))
        
      }
    }

    
}




#' Plot function for S3 class "atlasqtl"
#' 
#' This function provides functionality for plotting a so-called Manhattan plot 
#' with the position and size of "hotspot" predictors (i.e., predictors 
#' associated with multiple responses).
#' 
#' @param x an object with S3 class \code{"atlasqtl"}.
#' @param thres threshold to be applied on the posterior probabilities of 
#'              association to define hotspots (default is 0.5, for the 
#'              Median Probability Model).
#' @param fdr_adjust if TRUE the threshold is applied on FDR estimates computed
#'                   from the "atlasqtl" posterior probabilities of association,
#'                   otherwise it is directly applied on the posterior 
#'                   probabilities of association (default is FALSE).
#' @param pch type of points.
#' @param ylim_max upper limit for y-axis (must be set to NULL if \code{add} is TRUE).
#' @param main plot title (must be set to NULL if \code{add} is TRUE).
#' @param xlab x-axis title (must be set to NULL if \code{add} is TRUE).
#' @param ylab y-axis title (must be set to NULL if \code{add} is TRUE).
#' @param add whether the plot should be overlaid on an existing plot.
#' @param ... additional arguments.
#' 
#' @seealso \code{\link{atlasqtl}}
#' 
#' @export 
#' 
plot.atlasqtl <- function(x, thres = 0.5, fdr_adjust = FALSE, pch = 20, 
                          ylim_max = NULL, main = "Hotspot manhattan plot", 
                          xlab = "Predictors", 
                          ylab = "Nb of associated responses",
                          add = FALSE, ...) { 
  
  if (fdr_adjust) {
    mat_fdr <- assign_bFDR(x$gam_vb)
    rs_thres <- rowSums(mat_fdr < thres)
  } else {
    rs_thres <- rowSums(x$gam_vb > thres)
  }
  
  if (!add) {
    plot(rs_thres, pch = pch, 
         ylim = c(0, ifelse(is.null(ylim_max), max(rs_thres), ylim_max)),
         main = main, 
         xlab = xlab,
         ylab = ylab)
  } else { 
    points(rs_thres, pch = pch)
  }
  
}


#' Compute Bayesian FDR estimates for pairwise association between the 
#' predictors and responses
#' 
#' This function computes Bayesian FDR estimates from posterior probabilities of 
#' association.
#' 
#' @param mat_ppi a matrix of posterior probabilities of association (e.g, gam_vb
#'                from the "atlasqtl" object).
#' 
#' @seealso \code{\link{atlasqtl}}
#' 
#' @export 
#' 
assign_bFDR <- function(mat_ppi) {
  
  vec_ppi <- as.vector(mat_ppi)
  ind <- order(vec_ppi, decreasing = TRUE)
  vec_ppi_ord <- vec_ppi[ind]
  
  vec_fdr_ord <- cumsum(1-vec_ppi_ord) / 1:length(vec_ppi)
  
  ind_back <- order(ind)
  vec_fdr <- vec_fdr_ord[ind_back]
  
  mat_fdr <- matrix(vec_fdr, ncol = ncol(mat_ppi))
  rownames(mat_fdr) <- rownames(mat_ppi)
  colnames(mat_fdr) <- colnames(mat_ppi)
  
  mat_fdr
}
