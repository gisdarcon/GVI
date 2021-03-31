#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector tangents(int x1, int y1, double height0, IntegerVector xy2, NumericVector dsm_profile) {
  NumericVector out(dsm_profile.size());
  
  for(int i = 0; i < xy2.size(); ++i) {
    //# Distance traveled so far
    double distance_traveled = sqrt((x1 - xy2[i])*(x1 - xy2[i]) + (y1 - xy2[i+1])*(y1 - xy2[i+1]));
    //# Calculate tangent
    out[i/2] = (dsm_profile[i/2] - height0) / (distance_traveled);
    i+=1;
  }
  return out;
}