#include <Rcpp.h>

using namespace Rcpp;

//' C++ implementation to check if a matrix is an adjacency matrix
//' 
//' @param m A matrix to check.
//' @return Returns TRUE if the matrix is an adjacency matrix and FALSE 
//' otherwise.
//' @export
//[[Rcpp::export]]
bool is_adjacency_cpp(NumericMatrix m) {
  if(m.nrow() != m.ncol())
    return 0;
  
  int p = m.nrow();
  
  // Check the diagonal.
  for(int i = 0; i < p; i++) {
    if(m(i, i) != 0) {
      return 0;
    }
  }
  
  // Check the off-diagonal.
  for(int j = 0; j < (p - 1); j++) {
    for(int i = j + 1; i < p; i++) {
      if(m(i, j) != m(j, i)) {
        return 0;
      }
      if(m(i, j) != 1 && m(i, j) != 0) {
        return 0;
      }
    }
  }
  
  return 1;
}