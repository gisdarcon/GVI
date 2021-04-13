#include <Rcpp.h>
using namespace Rcpp;

int Sign(const int dxy)
{
  if(dxy<0) return -1; 
  else if(dxy>0) return 1; 
  else return 0;
}

IntegerMatrix bresenham(const int x1, const int y1, const int x2, const int y2)
{
  int Dx = x2 - x1;
  int Dy = y2 - y1;
  
  //# Increments
  int Sx = Sign(Dx); 
  int Sy = Sign(Dy);
  
  //# Segment length
  Dx = abs(Dx); 
  Dy = abs(Dy); 
  
  int D;
  if(Dx>=Dy) D = Dx;
  else D = Dy;
  
  
  //# Initial remainder
  double R = D / 2;
  
  int X = x1;
  int Y = y1;
  
  IntegerMatrix out((D+1),2);
  
  if(Dx > Dy) {   
    //# Main loop
    for(int i=0; i<=D; i++) {   
      out(i,0) = Y;
      out(i,1) = X;
      
      //# Update (X, Y) and R
      X+= Sx; R+= Dy; //# Lateral move
      if (R >= Dx) {
        Y+= Sy; 
        R-= Dx; //# Diagonal move
      }
    }
  } else {   
    //# Main loop
    for(int i=0; i<=D; i++) {     
      out(i,0) = Y;
      out(i,1) = X;
      
      //# Update (X, Y) and R
      Y+= Sy; R+= Dx; //# Lateral move
      if(R >= Dy) {    
        X+= Sx; 
        R-= Dy; //# Diagonal move
      }
    }
  }
  return out;
}

NumericVector tangents(const int x1, const int y1, const double height0, const IntegerMatrix xy2, const NumericVector dsm_profile) {
  NumericVector out(xy2.nrow());
  
  for(int i = 0; i < xy2.nrow(); ++i) {
    //# Distance traveled so far
    const double distance_traveled = sqrt((x1 - xy2(i,1))*(x1 - xy2(i,1)) + (y1 - xy2(i,0))*(y1 - xy2(i,0)));
    
    //# Calculate tangent
    out[i] = (dsm_profile[i] - height0) / (distance_traveled);
  }
  return out;
}

LogicalVector isVisible(const NumericVector x) {
  const int n = x.size();
  LogicalVector out(n);
  out[0] = true;
  
  double max_tangent = -9999;
  double this_tangent;
  
  for(int i = 1; i < n; ++i) {
    this_tangent = x[i];
    
    if (this_tangent > max_tangent) {
      max_tangent = this_tangent;
      out[i] = true;
    }
  }
  return out;
}

//' @title Line of Sight
//' @description Computes a single Line of Sight (LoS) using a Digital Surface Model (DSM) matrix. Uses a C++ implementation of Bresenham's line algorithm.
//'
//' @param rc1 integer matrix; Row and column of the targets location in the DSM matrix.
//' @param r0 integer; Specifying the row of the starting location in the DSM matrix.
//' @param c0 integer; Specifying the column of the starting location in the DSM matrix.
//' @param observerHeight numeric > 0; Height of the observer (e.g. 1.7).
//' @param dsm_mat numeric matrix; Matrix derived from a DSM.
//'
//' @return Cell numbers of all visible cells of the DSM along the line. 
//' 
//' @export
// [[Rcpp::export]]
std::vector<int> LoSCpp(const IntegerMatrix rc1, const int r0, const int c0, const double observerHeight, const NumericMatrix dsm_mat) {
  std::vector<int> visibleCells;
  
  for(int i = 0; i < rc1.nrow(); ++i) {
    //# End points
    int r1 = rc1(i,0);
    int c1 = rc1(i,1);
    
    //# Get XY coordinates of all points in the matrix along a line from r0,c0 to r1,c1
    IntegerMatrix XYmat = bresenham(c0, r0, c1, r1);
    
    //# Pull DSM height values from matrix
    NumericVector dsm_profile(XYmat.nrow());
    for(int j = 0; j < XYmat.nrow(); ++j) {
      dsm_profile[j] = dsm_mat(XYmat(j,0),XYmat(j,1));
    }
    
    //# Is visible? Current tangent must be greater than max. tangent
    NumericVector tangent = tangents(c0, r0, observerHeight, XYmat, dsm_profile);
    LogicalVector visibile = isVisible(tangent);
      
    //# Get visible cell numbers
    for(int j = 0; j < XYmat.nrow(); ++j) {
      if(visibile[j] == TRUE) {
        int thisCel = (XYmat(j,1)-1)*dsm_mat.nrow() + XYmat(j,0);
        visibleCells.push_back(thisCel);
      }
    }
  }
  
  return visibleCells;
}

