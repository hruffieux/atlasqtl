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
                  MapMat cp_X_Xbeta,
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
      
      
      double m1_beta_jk = m1_beta(j, k);
      //double cp_X_Xbeta_kj = cp_X_Xbeta(k, j) - m1_beta_jk * cp_X(j,j);
      double cp_X_Xbeta_jk = cp_X_Xbeta(j, k) - m1_beta_jk * cp_X(j,j);
      
      mu_beta_vb(j, k) = c * sig2_beta_vb[k] * tau_vb[k] * (cp_Y_X(k, j) - cp_X_Xbeta_jk);
      
      gam_vb(j, k) = exp(-logOnePlusExp(c * (log_1_min_Phi_theta_plus_zeta(j, k) - log_Phi_theta_plus_zeta(j, k)
                                               - mu_beta_vb(j, k)*mu_beta_vb(j, k) / (2 * sig2_beta_vb[k])
                                               + cst[k])));
                                               
       m1_beta(j, k) = gam_vb(j, k) * mu_beta_vb(j, k);
       
       cp_X_Xbeta.col(k) += (m1_beta(j, k) - m1_beta_jk) * cp_X.col(j);
       
                                               
    } 
  }
}



// bool is_in_vector(int k, const std::vector<int>& vec) {
//   return std::find(vec.begin(), vec.end(), k) != vec.end();
// }

// for atlasqtl_core function handling missing values in Y
// [[Rcpp::export]]
void coreDualMisLoop(const MapMat cp_X,
                     const MapMat sq_X,
                     const MapMat cp_Y_X,
                     MapArr2D gam_vb,
                     const MapArr2D log_Phi_theta_plus_zeta,
                     const MapArr2D log_1_min_Phi_theta_plus_zeta,
                     const double log_sig2_inv_vb,
                     const MapArr1D log_tau_vb,
                     MapMat m1_beta,
                     MapMat cp_X_Xbeta, //MapArr2D X_beta_vb,
                     MapArr2D mu_beta_vb,
                     const MapArr2D sig2_beta_vb,
                     const MapArr1D tau_vb,
                     const Eigen::VectorXi shuffled_ind,
                     const Eigen::VectorXi sample_q,
                     const MapMat X,
                     const Eigen::VectorXi k_index_vec, //storing the index of k in vec_ind_k_mis_zero_ind
                     const List& list_ind_mis,
                     const double c = 1) {
  
  
  const Arr1D cst = -(log_tau_vb + log_sig2_inv_vb)/ 2;
  int q = cp_Y_X.rows();
  int p =  X.cols();
  int n = X.rows();
  
  for (int a = 0; a < sample_q.size(); a++) {
    int k = sample_q[a]; // k as column, j as row
    
    if (k_index_vec[k] >= 0) {
      
      int k_mis = k_index_vec[k]; 
      int n_mis = as<Eigen::VectorXi>(list_ind_mis[k_mis]).size(); // number of NA for response k
      // Eigen::MatrixXd X_ind_k_mis(n_mis, p);
      // for (int ind = 0; ind < n_mis; ++ind)  {
      //   // construct this matrix (haven't found an easy way to subset X, i.e., select a subset of rows X[list_ind_mis[[i_list_ind_mis]], ])
      //   int ind_n_mis = as<Eigen::VectorXi>(list_ind_mis[k_mis])(ind)-1; // minus 1 as zero-based
      //   X_ind_k_mis.row(ind) = X.row(ind_n_mis);
      // }
      
      
      Eigen::VectorXi ind_n_mis = Eigen::Map<Eigen::VectorXi>(INTEGER(list_ind_mis[k_mis]), n_mis) - Eigen::VectorXi::Constant(n_mis, 1);
      // Eigen::MatrixXd X_ind_k_mis = X(ind_n_mis, Eigen::all);
      
      for (int b = 0; b < shuffled_ind.size(); b++) {
        int j = shuffled_ind[b];
        
        double m1_beta_jk = m1_beta(j, k);
        // double cp_X_Xbeta_jk = cp_X_Xbeta(j, k) - m1_beta_jk*(cp_X(j, j)-X_ind_k_mis.col(j).array().square().sum());
        // X_ind_k_mis.col(j).array().square().sum() can be replaced by extracting from sq_X directly
        double cp_X_Xbeta_jk = cp_X_Xbeta(j, k) - m1_beta_jk*(cp_X(j, j) - sq_X(ind_n_mis, j).sum());
        
        //calculating X_ind_k_mis.transpose() * X_ind_k_mis.col(j) directly without extracting X_ind_k_mis
        Eigen::MatrixXd cp_X_mis_j = Eigen::VectorXd::Zero(p, 1);
        for (int col=0; col < p; col++){
          for (int row=0; row < n_mis; row++){
            cp_X_mis_j(col, 0) += X(ind_n_mis[row], col) * X(ind_n_mis[row], j);
          }
        }
        
        
        mu_beta_vb(j, k) = c * sig2_beta_vb(j,k) * tau_vb[k] * (cp_Y_X(k, j) - cp_X_Xbeta_jk);
        
        gam_vb(j, k) = exp(-logOnePlusExp(c * (log_1_min_Phi_theta_plus_zeta(j, k) -
          log_Phi_theta_plus_zeta(j, k) - mu_beta_vb(j, k)*mu_beta_vb(j, k) / (2 * sig2_beta_vb(j, k)) -
          log(sig2_beta_vb(j, k)) / 2 + cst[k])));
        
        m1_beta(j, k) = gam_vb(j, k) * mu_beta_vb(j, k);
        
        // cp_X_Xbeta.col(k) += (m1_beta(j, k) - m1_beta_jk) * (cp_X.col(j) - (X_ind_k_mis.transpose() * X_ind_k_mis.col(j)));
        cp_X_Xbeta.col(k) += (m1_beta(j, k) - m1_beta_jk) * (cp_X.col(j) - cp_X_mis_j);
        
      }
      
    }else{
      
      for (int b = 0; b < shuffled_ind.size(); b++) {
        int j = shuffled_ind[b];
        
        double m1_beta_jk = m1_beta(j, k);
        double cp_X_Xbeta_jk = cp_X_Xbeta(j, k) - m1_beta_jk*cp_X(j, j);
        
        mu_beta_vb(j, k) = c * sig2_beta_vb(j,k) * tau_vb[k] * (cp_Y_X(k, j) - cp_X_Xbeta_jk);
        
        gam_vb(j, k) = exp(-logOnePlusExp(c * (log_1_min_Phi_theta_plus_zeta(j, k) -
          log_Phi_theta_plus_zeta(j, k) - mu_beta_vb(j, k)*mu_beta_vb(j, k) / (2 * sig2_beta_vb(j, k)) -
          log(sig2_beta_vb(j, k)) / 2 + cst[k])));
        
        m1_beta(j, k) = gam_vb(j, k) * mu_beta_vb(j, k);
        cp_X_Xbeta.col(k) += (m1_beta(j, k) - m1_beta_jk) * cp_X.col(j);
      }
    }
  }
}




