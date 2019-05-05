#' atlasqtl: a package for fAsT gLobal-locAl hotSpot QTL detection
#'
#' Global-local hierachical approach for joint detection of hotspot predictors 
#' in multiple-response regression. Inference uses closed-form variational 
#' updates coupled with a simulated annealing algorithm to enhance exploration 
#' of highly multimodal spaces. This software allows joint inference at 
#' large-scale, e.g., for molecular quantitative trait locus (QTL) studies where 
#' the hotspots are genetic variants regulating several molecular traits 
#' simultaneously. See H. Ruffieux, et al., A global-local approach for 
#' detecting hotspots in multiple-response regression, 2018, arxiv:1811.03334). 
#'
#' @section atlasqtl functions: set_hyper, set_init, atlasqtl.
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
