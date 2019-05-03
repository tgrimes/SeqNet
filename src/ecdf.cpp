#include <Rcpp.h>

using namespace Rcpp;

//' C++ implementation of empirical CDF
//' 
//' Constructs the empirical CDF, F, for a set of observations, x, and 
//' returns F(x).
//' @param x The observation to construct the empirical CDF from.
//' @return Returns the values for F(x).
//' @export
//[[Rcpp::export]]
NumericVector ecdf_cpp(NumericVector x) {
  int n = x.size();
  NumericVector obs = x + 0.0;
  NumericVector ans(n);
  std::sort(obs.begin(), obs.end());
  for (int i = 0; i < n; ++i){
    ans[i] = (std::upper_bound(obs.begin(), obs.end(), x[i]) - obs.begin());
  }
  return ans/((double) n);
}