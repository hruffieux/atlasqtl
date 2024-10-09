#ifndef ATLASQTL_UTILS_H_
#define ATLASQTL_UTILS_H_

#include <RcppEigen.h>
#include <RcppParallel.h>
#include "atlasqtl_types.h"

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppParallel)]]

using namespace Rcpp;

Arr1D logOnePlusExp(const  Arr1D& x);

#endif // ATLASQTL_UTILS_H_
