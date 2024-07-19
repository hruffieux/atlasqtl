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
/*
 double crossprod(const Eigen::VectorXd &one, const Eigen::VectorXd &two) {
 // R Crossproduct of vectors = sum of element-wise products
 // Assume one.size() == two.size() without checking
 double sum = 0;
 for (int i = 0; i < one.size(); i++) {
 sum += one[i] * two[i];
 }
 return sum;
 }
 */

double logOnePlusExp(double x) {
  double m = x;
  if (x < 0)
    m = 0;
  return log(exp(x-m) + exp(-m)) + m;
}


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
                  const Eigen::VectorXi sample_q,
                  const double c = 1) {
  
  // Rcout << "one\n";
  
  const Arr1D cst = -(log_tau_vb + log_sig2_inv_vb + log(sig2_beta_vb) )/ 2;
  
  for (int a = 0; a < sample_q.size(); a++) {
    int k = sample_q[a];
    
    // Rcout << "two" << a << " " << k <<  "\n";
    // Rcout << sample_q << "\n";
    
    for (int b = 0; b < shuffled_ind.size(); b++) {
      int j = shuffled_ind[b];
      
      //Rcout << "three: " << a << " " << k << " " << b << " " << j << "\n";
      
      double m1_beta_jk = m1_beta(j, k);
      //double cp_betaX_X_kj = cp_betaX_X(k, j) - m1_beta_jk * cp_X(j,j);
      double cp_betaX_X_jk = cp_betaX_X(j, k) - m1_beta_jk * cp_X(j,j);
      
      mu_beta_vb(j, k) = c * sig2_beta_vb[k] * tau_vb[k] * (cp_Y_X(k, j) - cp_betaX_X_jk);
      
      gam_vb(j, k) = exp(-logOnePlusExp(c * (log_1_min_Phi_theta_plus_zeta(j, k) - log_Phi_theta_plus_zeta(j, k)
                                               - mu_beta_vb(j, k)*mu_beta_vb(j, k) / (2 * sig2_beta_vb[k])
                                               + cst[k])));
                                               
                                               m1_beta(j, k) = gam_vb(j, k) * mu_beta_vb(j, k);
                                               
                                               cp_betaX_X.col(k) += (m1_beta(j, k) - m1_beta_jk) * cp_X.col(j);
                                               
                                               
    } 
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
                     const Eigen::VectorXi sample_q,
                     const double c = 1) {
  
  const Arr1D cst = -(log_tau_vb + log_sig2_inv_vb)/ 2;
  // int q = cp_Y_X.rows();
  
  for (int a = 0; a < sample_q.size(); a++) {
    int k = sample_q[a];
    MapMat cp_X_rm_k = as<MapMat>(cp_X_rm[k]);
    
    for (int b = 0; b < shuffled_ind.size(); b++) {
      int j = shuffled_ind[b];
      
      
      double m1_beta_jk = m1_beta(j, k);
      //double cp_betaX_X_kj = cp_betaX_X(k, j) - m1_beta(j, k)*(cp_X(j, j) - cp_X_rm_k(j, j));
      double cp_betaX_X_jk = cp_betaX_X(j, k) - m1_beta_jk*(cp_X(j, j) - cp_X_rm_k(j, j));
      
      // cp_betaX_X.row(k) += -m1_beta(j, k) * (cp_X.row(j) - cp_X_rm_k.row(j));
      
      mu_beta_vb(j, k) = c * sig2_beta_vb(j,k) * tau_vb[k] * (cp_Y_X(k, j) - cp_betaX_X_jk);
      
      gam_vb(j, k) = exp(-logOnePlusExp(c * (log_1_min_Phi_theta_plus_zeta(j, k) -
        log_Phi_theta_plus_zeta(j, k) - mu_beta_vb(j, k)*mu_beta_vb(j, k) / (2 * sig2_beta_vb(j, k)) -
        log(sig2_beta_vb(j, k)) / 2 + cst[k])));
      
      m1_beta(j, k) = gam_vb(j, k) * mu_beta_vb(j, k);
      cp_betaX_X.col(k) += (m1_beta(j, k) - m1_beta_jk) * (cp_X.col(j) - cp_X_rm_k.col(j)); 
      
      
    } 
  }
  
}
