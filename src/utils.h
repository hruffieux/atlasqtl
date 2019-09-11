#ifndef ATLASQTL_UTILS_H_
#define ATLASQTL_UTILS_H_

#include <RcppEigen.h>
#include "atlasqtl_types.h"

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;

Arr1D logOnePlusExp(const  Arr1D& x);

#endif // ATLASQTL_UTILS_H_
