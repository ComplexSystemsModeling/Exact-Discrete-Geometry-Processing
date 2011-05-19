#include <itkMatrix.h>

#include "vnl/vnl_det.h" 
#include "vnl/vnl_matrix_fixed.h" 

#include <iostream> 
#include <cstdlib>

int main( int argc, char** argv )
{
  if( argc < 3 )
    {
    std::cout << "usage: <exe> epsilon_x epsilon_y";
    std::cout << std::endl;
    }

  typedef itk::Matrix< double, 4,4 > MatrixType;
  MatrixType M;

  char **stringPtr;
  double epsilon_x = std::strtod( argv[1], stringPtr );
  double epsilon_y = std::strtod( argv[2], stringPtr );
  double ax = 0.0;
  double ay = 0.0;
  double bx = 1.0;
  double by = 0.0;
  double cx = 1.0;
  double cy = 1.0;
  double dx = 0.0 + epsilon_x;
  double dy = 1.0 + epsilon_y;

  M(0,0) = ax;
  M(0,1) = ay;
  M(0,2) = ax*ax+ay*ay;
  M(0,3) = 1.0;
  M(1,0) = bx;
  M(1,1) = by;
  M(1,2) = bx*bx+by*by;
  M(1,3) = 1.0;
  M(2,0) = cx;
  M(2,1) = cy;
  M(2,2) = cx*cx+cy*cy;
  M(2,3) = 1.0;
  M(3,0) = dx;
  M(3,1) = dy;
  M(3,2) = dx*dx+dy*dy;
  M(3,3) = 1.0;
 
  std::cout << "M: " << std::endl;
  std::cout << M << std::endl;

  double det = vnl_det( M.GetVnlMatrix() );

  std::cout << "Det(M): " << det << std::endl;

  return(det>0?0:1);
}
