#include <Rcpp.h>
#include <stack>

using namespace Rcpp;

//' C++ implementation to obtain connected components in a graph.
//' 
//' @param adj An adjacency matrix.
//' @return Returns a matrix with 2 columns containing the indicies in the
//' lower-triangle of the matrix that are nonzero.
//' @export
//[[Rcpp::export]]
IntegerVector components_in_adjacency(NumericMatrix adj) {
  int p = adj.nrow();
  LogicalVector visited(p); 
  IntegerVector group(p); 
  
  for(int i = 0; i < p; i++) {
    visited[i] = FALSE;
    group[i] = 0;
  }
  
  std::stack <int> node;
  int component = 0;
  int k;
  for(int i = 0; i < p; i++) {
    // If node has not been visted, push the node and investigate its neighbors.
    if(!visited[i]) {
      node.push(i);
      component++;
    }
    while(!node.empty()) {
      k = node.top();
      node.pop();
      // If node has not been visited, then push all of its neighbors.
      if(!visited[k]) {
        // Set node to current component.
        group[k] = component;
        for(int j = 0; j < p; j++) {
          // Push neighbors (to be visited) and set component.
          if(adj(k, j) != 0) {
            node.push(j);
          }
        }
        // The node has now been visited.
        visited[k] = 1;
      }
    }
    // After exploring all neighbors, go to next component.
  }
  
  return group;
}