#include <Rcpp.h>

using namespace Rcpp;

//' C++ implementation for creating a ring lattice
//' 
//' @param p The number of nodes in the lattice.
//' @param neig_size The neighborhood side within which nodes are connected.
//' @return Returns the adjacency matrix for the ring lattice.
//' @export
//[[Rcpp::export]]
NumericMatrix ring_lattice_cpp(int p, int neig_size) {
  int size = neig_size; // Make a copy.
  if(size > (p / 2.0)) {
    size = (int) p / 2.0;
  }
  
  NumericMatrix adj(p, p);
  
  for(int i = 0; i < p; i++) {
    for(int j = 1; j <= size; j++) {
      int neig = i + j;
      if(neig >= p) {
        neig = neig - p;
      }
      
      adj(i, neig) += 1.0;
      adj(neig, i) += 1.0;
    }
  }
  
  return adj;
}