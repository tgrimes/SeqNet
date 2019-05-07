#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

//' C++ implementation of qzinb
//' 
//' @param p A vector of probabilities
//' @param mu The dispersion paramater used in dnbinom.
//' @param size The distribution mean.
//' @param rho The zero-inflation parameter.
//' @param lower Logical; if TRUE, then probabilities are P(X <= x). 
//' Otherwise, P(X > x).
//' @param log Logical; if TRUE, then exp(p) is used.
//' @export
// [[Rcpp::export]]
NumericMatrix overwrite_with_qzinb_cpp(NumericMatrix p, NumericVector size, NumericVector mu,
                          NumericVector rho,
                          bool lower = true, bool log = false ) {
  if(!lower) {
    p = 1 - p;
  }
  
  if(log) {
    for(int i = 0; i < p.ncol(); i++) {
      for(int j = 0; j < p.nrow(); j++) {
        p(i, j) = exp(p(i, j));
      }
    }
  }
    
  double x;
  for(int i = 0; i < p.ncol(); i++) {
    if(rho[i] == 0) {
      p(_, i) = Rcpp::qnbinom_mu(p(_, i), 
        size[i], mu[i], lower, log);
    } else {
      p(_, i) = Rcpp::qnbinom_mu(Rcpp::pmax(0, (p(_, i) - rho[i]) / (1 - rho[i])), 
        size[i], mu[i], lower, log);
    }
  }
  return p;
}
