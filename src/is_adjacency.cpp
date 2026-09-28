#include <Rcpp.h>

using namespace Rcpp;

//' C++ implementation to check if a matrix is an adjacency matrix
//' 
//' @param m A matrix to check.
//' @return Returns 0 if the matrix is an adjacency matrix. If the matrix is 
//' not square, returns 1; if the diagonal entries are not all zero, returns 2; 
//' if the matrix is not symmetric, returns 3; if the matrix contains values
//' other than 0 or 1, returns 4.
//' @export
//[[Rcpp::export]]
int check_adjacency_cpp(NumericMatrix m) {
  if(m.nrow() != m.ncol())
    return 1;
  
  int p = m.nrow();
  
  // Check the diagonal.
  for(int i = 0; i < p; i++) {
    if(m(i, i) != 0) {
      return 2;
    }
  }
  
  // Check the off-diagonal.
  for(int j = 0; j < (p - 1); j++) {
    for(int i = j + 1; i < p; i++) {
      if(m(i, j) != m(j, i)) {
        return 3;
      }
      if(m(i, j) != 1 && m(i, j) != 0) {
        return 4;
      }
    }
  }
  
  return 0;
}
