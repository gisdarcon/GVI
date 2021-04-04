#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List countGroups(IntegerMatrix xyVisible, IntegerVector uniqueXY) {
  IntegerVector dxy = seq(1, max(uniqueXY));
  IntegerVector visibleTotal(dxy.size());
  IntegerVector visibleGreen(dxy.size());
  
  for(int i = 0; i < xyVisible.nrow(); ++i) {
    int this_dxy = xyVisible(i,0);
    visibleTotal[this_dxy-1] += 1;
    visibleGreen[this_dxy-1] += xyVisible(i,1);
  }
  
  List ret;
  ret["dxy"] = dxy;
  ret["visibleTotal"] = visibleTotal;
  ret["visibleGreen"] = visibleGreen;
  return ret;
}