#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

using namespace Eigen;
using namespace Rcpp;

//' Random Walk with Restart using Eigen Sparse Matrix
//'
//' @param A column-normalized sparse matrix
//' @param v initial restart vector
//' @param restart restart probability (e.g., 0.5)
//' @param threshold convergence threshold
//' @param max_iter maximal number of iterations
//' @return numeric vector of stationary probabilities
//' @noRd
// [[Rcpp::export]]
NumericVector rwr_eigen_cpp(const Eigen::MappedSparseMatrix<double>& A, 
                            const Eigen::Map<Eigen::VectorXd>& v, 
                            double restart, 
                            double threshold, 
                            int max_iter) {
    int n = A.rows();
    VectorXd u = v;
    VectorXd u_old = VectorXd::Zero(n);
    
    // Pre-calculate matrices and vectors for iteration
    SparseMatrix<double> trans_A = A * (1.0 - restart);
    VectorXd restart_v = v * restart;
    
    int iter = 0;
    while(iter < max_iter) {
        u_old = u;
        // Sparse matrix multiplication, Eigen optimizes this in C++
        u = trans_A * u_old + restart_v; 
        
        // Check for convergence (L1 norm of difference)
        if((u - u_old).cwiseAbs().sum() < threshold) {
            break;
        }
        iter++;
    }
    
    return Rcpp::wrap(u);
}
