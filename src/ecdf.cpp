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
  int nobs = x.size();
  NumericVector ans(nobs);
  std::sort(x.begin(), x.end());
  for (int i = 0; i < nobs; ++i){
    ans[i] = (std::upper_bound(x.begin(), x.end(), x[i]) - x.begin());
  }
  return ans/((double) nobs);
}