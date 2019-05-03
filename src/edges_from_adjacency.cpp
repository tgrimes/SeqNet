#include <Rcpp.h>

using namespace Rcpp;

//' C++ implementation for obtaining an edge list from adjacency matrix
//' 
//' @param adj An adjacency matrix.
//' @return Returns a matrix with 2 columns containing the indicies in the
//' lower-triangle of the matrix that are nonzero.
//' @export
//[[Rcpp::export]]
NumericMatrix edges_from_adjacency_cpp(NumericMatrix adj) {
  int p = adj.nrow();
  int n_edges = 0;
  for(int j = 0; j < (p - 1); j++) {
    for(int i = j + 1; i < p; i++) {
      if(adj(i, j) != 0) {
        n_edges++;
      }
    }
  }
  NumericMatrix edges(n_edges, 2);
  int index = 0;
  
  if(n_edges > 0) {
    for(int j = 0; j < (p - 1); j++) {
      for(int i = j + 1; i < p; i++) {
        if(adj(i, j) != 0) {
          edges(index, 0) = j + 1;
          edges(index, 1) = i + 1;
          index++;
        }
      }
    }
  } 
  
  return edges;
}