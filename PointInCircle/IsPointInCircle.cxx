#include <itkMatrix.h>

#include "vnl/vnl_det.h" 
#include "vnl/vnl_matrix_fixed.h" 

#include <iostream> 
#include <cstdlib>

extern "C"
{
double incircle(double* pa, double* pb, double* pc, double* pd);
}

int main( int argc, char** argv )
{
  if( argc < 4 )
    {
    std::cout << "usage: <exe> epsilon_x epsilon_y IsExact";
    std::cout << std::endl;
    }

  char **stringPtr;
  double epsilon_x = std::strtod( argv[1], stringPtr );
  double epsilon_y = std::strtod( argv[2], stringPtr );
  std::cout << "Epsilon_x: " << epsilon_x << std::endl;
  std::cout << "Epsilon_y: " << epsilon_y << std::endl;

  double ax = 1.0;
  double ay = 0.0;
  double bx = 1.0;
  double by = 1.0;
  double cx = 0.0;
  double cy = 1.0;
  double dx = epsilon_x;
  double dy = epsilon_y;
  std::cout << "dx: " << dx << std::endl;
  std::cout << "dy: " << dy << std::endl;

  double det;

  if( atoi( argv[3] ) == 0 )
    { 

    typedef itk::Matrix< double, 4,4 > MatrixType;
    MatrixType M;

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

    det = vnl_det( M.GetVnlMatrix() );
    std::cout << "Det(M): " << det << std::endl;

    }
  else
    {
    double * pa = new double[2];
    double * pb = new double[2];
    double * pc = new double[2];
    double * pd = new double[2]; 
    pa[0] = ax; pa[1] = ay;
    pb[0] = bx; pb[1] = by;
    pc[0] = cx; pc[1] = cy;
    pd[0] = dx; pd[1] = dy;

    det = incircle( pa, pb, pc, pd );
    std::cout << "Det(M): " << det << std::endl;
    }

  return(det>=0?0:1);

}
