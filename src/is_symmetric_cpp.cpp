#include <Rcpp.h>
#include <cmath>

using namespace Rcpp;

//' C++ implementation to check if a matrix is symmetric
//' 
//' @param m A matrix to check.
//' @param tol A Numeric scalar >= 0. Differences smaller than tol are ignored.
//' @return Returns TRUE if the matrix is symmetric and FALSE otherwise.
//' @export
//[[Rcpp::export]]
bool is_symmetric_cpp(NumericMatrix m, double tol = 1e-12) {
  if(m.nrow() != m.ncol())
    return 0;
  
  int p = m.nrow();
  
  // Check the off-diagonal.
  for(int j = 0; j < (p - 1); j++) {
    for(int i = j + 1; i < p; i++) {
      if(std::abs((double) (m(i, j) - m(j, i))) > tol) {
        return 0;
      }
    }
  }
  
  return 1;
}