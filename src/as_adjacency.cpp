#include <Rcpp.h>

using namespace Rcpp;

//' C++ implementation to convert a matrix to an adjacency matrix
//' 
//' Sets all nonzero values to 1.
//' @param m A matrix to convert.
//' @return Returns the adjacency matrix.
//' @export
//[[Rcpp::export]]
NumericMatrix overwrite_as_adjacency_cpp(NumericMatrix m, double tol = 10^-13) {
  if(m.nrow() != m.ncol())
    stop("m is not a square matrix.");
    
  int p = m.nrow();
  
  for(int i = 0; i < p; i++) {
    m(i, i) = 0;
  }
  
  // Use values from the lower-diagonal.
  for(int j = 0; j < (p - 1); j++) {
    for(int i = j + 1; i < p; i++) {
      if(abs(m(i, j)) < tol) {
        m(i, j) = 0;
        m(j, i) = 0;
      } else {
        m(i, j) = 1;
        m(j, i) = 1;
      }
    }
  }
  
  return m;
}