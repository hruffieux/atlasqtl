/*
 *
 * This file is part of the `atlasqtl` R package:
 *     https://github.com/hruffieux/atlasqtl
 *
 * Functions for computationally expensive updates in algorithms without
 * external information.
 *
 * These functions use Eigen::Map to pass large matrices by reference from R.
 * Given dimensionalities involved in some applications, copying such matrices
 * would imply a prohibitive RAM overconsumption.
 *
 */

#include "utils.h"

// for atlasqtl_core function
// [[Rcpp::export]]
void coreDualLoop(const MapMat cp_X,
                  const MapMat cp_Y_X,
                  MapArr2D gam_vb,
                  const MapArr2D log_Phi_theta_plus_zeta,
                  const MapArr2D log_1_min_Phi_theta_plus_zeta,
                  const double log_sig2_inv_vb,
                  const MapArr1D log_tau_vb,
                  MapMat m1_beta,
                  MapMat cp_betaX_X,
                  MapArr2D mu_beta_vb,
                  const MapArr1D sig2_beta_vb,
                  const MapArr1D tau_vb,
                  const Eigen::VectorXi shuffled_ind,
                  const double c = 1) {
  
  const Arr1D cst = -(log_tau_vb + log_sig2_inv_vb + log(sig2_beta_vb) )/ 2;
  
  for (int i = 0; i < shuffled_ind.size(); ++i) {
    
    int j = shuffled_ind[i];
    
    // cout << cp_betaX_X.rows();
    
    cp_betaX_X.noalias() -= m1_beta.row(j).transpose() * cp_X.row(j);
    
    mu_beta_vb.row(j) = c * sig2_beta_vb * tau_vb * (cp_Y_X.col(j) - cp_betaX_X.col(j)).array();
    
    gam_vb.row(j) = exp(-logOnePlusExp(c * (log_1_min_Phi_theta_plus_zeta.row(j) -
      log_Phi_theta_plus_zeta.row(j) - mu_beta_vb.row(j).square() / (2 * sig2_beta_vb.transpose()) +
      cst.transpose())));
    
    m1_beta.row(j) = mu_beta_vb.row(j) * gam_vb.row(j);
    
    cp_betaX_X.noalias() += m1_beta.row(j).transpose() * cp_X.row(j);
    
  }
  
}


// for atlasqtl_core function handling missing values in Y
// [[Rcpp::export]]
void coreDualMisLoop(const MapMat cp_X,
                     const List cp_X_rm,
                     const MapMat cp_Y_X,
                     MapArr2D gam_vb,
                     const MapArr2D log_Phi_theta_plus_zeta,
                     const MapArr2D log_1_min_Phi_theta_plus_zeta,
                     const double log_sig2_inv_vb,
                     const MapArr1D log_tau_vb,
                     MapMat m1_beta,
                     MapMat cp_betaX_X, //MapArr2D X_beta_vb,
                     MapArr2D mu_beta_vb,
                     const MapArr2D sig2_beta_vb,
                     const MapArr1D tau_vb,
                     const Eigen::VectorXi shuffled_ind,
                     const double c = 1) {
  
  const Arr1D cst = -(log_tau_vb + log_sig2_inv_vb)/ 2;
  int q = cp_Y_X.rows();
  
  for (int i = 0; i < shuffled_ind.size(); ++i) {
    
    int j = shuffled_ind[i];
    
    cp_betaX_X.noalias() -= m1_beta.row(j).transpose() * cp_X.row(j);
    for (int k = 0; k < q; ++k) {
      // MapMat cp_X_rm_k = as<MapMat>(cp_X_rm[k]);
      cp_betaX_X.row(k) += m1_beta(j, k) * as<MapMat>(cp_X_rm[k]).row(j);
    }
    mu_beta_vb.row(j) = c * sig2_beta_vb.row(j).transpose() * tau_vb *
      (cp_Y_X.col(j) - cp_betaX_X.col(j)).array();
    
    gam_vb.row(j) = exp(-logOnePlusExp(c * (log_1_min_Phi_theta_plus_zeta.row(j) -
      log_Phi_theta_plus_zeta.row(j) - mu_beta_vb.row(j).square() / (2 * sig2_beta_vb.row(j)) -
      log(sig2_beta_vb.row(j)) / 2 + cst.transpose())));
    
    m1_beta.row(j) = mu_beta_vb.row(j) * gam_vb.row(j);
    
    cp_betaX_X.noalias() += m1_beta.row(j).transpose() * cp_X.row(j);
    for (int k = 0; k < q; ++k) {
      // MapMat cp_X_rm_k = as<MapMat>(cp_X_rm[k]);
      cp_betaX_X.row(k) -= m1_beta(j, k) * as<MapMat>(cp_X_rm[k]).row(j);
    }
    
  }
  
}
