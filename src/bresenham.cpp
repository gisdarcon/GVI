#include <Rcpp.h>
using namespace Rcpp;

int Sign(int dxy)
{
  if(dxy<0) return -1; 
  else if(dxy>0) return 1; 
  else return 0;
}

// [[Rcpp::export]]
IntegerVector bresenham(int x1, int y1, int x2, int y2)
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
  
  IntegerVector out(2*(D+1));
  
  if(Dx > Dy) {   
    //# Main loop
    for(int I=0; I<=D; I++) {   
      out[2*I] = X;
      out[2*I+1] = Y;
      
      //# Update (X, Y) and R
      X+= Sx; R+= Dy; //# Lateral move
      if (R >= Dx) {
        Y+= Sy; 
        R-= Dx; //# Diagonal move
      }
    }
  } else {   
    //# Main loop
    for(int I=0; I<=D; I++) {     
      out[2*I] = X;
      out[2*I+1] = Y;
      
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