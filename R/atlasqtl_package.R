#' atlasqtl: a package for combined predictor and outcome selection in
#' high-dimensional set-ups using variational inference
#'
#' The atlasqtl package provides an efficient variational algorithm for
#' simultaneous variable selection of predictors and associated outcomes based
#' on a sparse multivariate regression model (H. Ruffieux, A. C. Davison,
#' J. Hager, I. Irincheeva, Efficient inference for genetic association studies
#' with multiple outcomes, Biostatistics, 2017). The methods from this package
#' have been used on large genetic datasets from molecular quantitative trait
#' atlasqtl (QTL) problems with over 200K single nucleotide polymorphisms (SNPs),
#' hundreds of molecular expression levels and hundreds of samples.
#'
#' @section atlasqtl functions: set_hyper, set_init, generate_null, atlasqtl,
#'   set_blocks, set_cv, set_groups, set_struct.
#'
#' @docType package
#' @name atlasqtl-package
#' @useDynLib atlasqtl, .registration = TRUE
#' @import RcppEigen
#' @importFrom Rcpp evalCpp
#' @importFrom stats cor dnorm median pnorm qnorm rbeta rbinom rgamma rnorm setNames uniroot var
#' @importFrom grDevices dev.off png
#' @importFrom graphics abline legend matplot points
NULL