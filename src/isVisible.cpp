#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
LogicalVector isVisible(NumericVector x) {
  int n = x.size();
  LogicalVector out(n);
  out[0] = true;
  
  double max_tangent = -9999;
  
  for(int i = 1; i < n; ++i) {
    double this_tangent = x[i];
    
    if (this_tangent > max_tangent) {
      max_tangent = this_tangent;
      out[i] = true;
    }
  }
  return out;
}