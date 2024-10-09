#' atlasqtl: an R package for variable selection in sparse regression with hierarchically-related responses
#'
#' Flexible sparse regression for variable selection in large predictor and 
#' response settings, based on a series of hierarchically-related spike-and-slab 
#' submodels. The model is also tailored to the detection of hotspots, namely, 
#' predictors associated with multiple responses, which it represents using a 
#' global-local horseshoe specification. Inference uses closed-form variational 
#' updates coupled with a simulated annealing algorithm to enhance exploration of 
#' highly multimodal spaces. This software allows joint inference at large scale, 
#' e.g., for molecular quantitative trait locus (QTL) studies, where the hotspots 
#' are genetic variants regulating several molecular traits simultaneously. 
#' See H. Ruffieux, A. C. Davison, J. Hager, J. Inshaw, B. Fairfax, S. Richardson, 
#' and L. Bottolo. A global-local approach for detecting hotspots in multiple 
#' response regression. The Annals of Applied Statistics, 14:905-928, 2020.
#'
#' @section atlasqtl functions: atlasqtl, assign_bFDR, map_hyperprior_elicitation, print.atlasqtl, set_hyper, set_init, summary.atlasqtl.
#'
#' "_PACKAGE"
#' @name atlasqtl-package
#' @useDynLib atlasqtl, .registration = TRUE
#' @import RcppEigen
#' @importFrom Rcpp evalCpp
#' @importFrom stats cor dnorm median pnorm qnorm rbeta rbinom rgamma rnorm setNames uniroot var
#' @importFrom grDevices dev.off png
#' @importFrom graphics abline legend matplot points
#' @importFrom RcppParallel RcppParallelLibs
NULL
