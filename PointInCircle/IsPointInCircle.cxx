#include <itkMatrix.h>
#include <itkPoint.h>

#include "vnl/vnl_det.h" 
#include "vnl/vnl_matrix_fixed.h" 

#include <iostream> 
#include <cstdlib>

//------------------------------------------------------------------------------
// Skewchuck code
//
extern "C"
{
double incircle(double* pa, double* pb, double* pc, double* pd);
double orient2d(double* pa, double* pb, double* pc);
}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
// The test functor to map skewchuck code to ITK's API
//
template< typename PointType >
bool
IsInside(
  PointType TrianglePoint1,
  PointType TrianglePoint2,
  PointType TrianglePoint3,
  PointType PointToTest )
{

  double * pa = new double[2];
  double * pb = new double[2];
  double * pc = new double[2];
  double * pd = new double[2]; 
  pa[0] = TrianglePoint1[0];
  pa[1] = TrianglePoint1[1];
  pb[0] = TrianglePoint2[0];
  pb[1] = TrianglePoint2[1];
  pc[0] = TrianglePoint3[0];
  pc[1] = TrianglePoint3[1];
  pd[0] = PointToTest[0];
  pd[1] = PointToTest[1];

  // orientation test
  double orientation = orient2d( pa, pb, pc );

  // incircle test - the result is multiplied by the orientation test result
  double det = incircle( pa, pb, pc, pd ) * orientation;
  
  // NOTE ALEX: replace by debug macro
  std::cout << "Det(M): " << det << std::endl;
 
  // NOTE ALEX: zero, which means the point is ON the circle is considered IN.
  return( det>=0?0:1 );
}
//------------------------------------------------------------------------------



//------------------------------------------------------------------------------
// test code
//
int main( int argc, char** argv )
{

  if( argc < 4 )
    {
    std::cout << "usage: <exe> epsilon_x epsilon_y IsExact";
    std::cout << std::endl;
    }

  //--------------
  // input parsing
  //--------------

  // stringPtr init
  char value = 0;
  char* ptr;
  ptr = &value;
  char** stringPtr;
  stringPtr = &ptr;

  double epsilon_x = std::strtod( argv[1], stringPtr );
  double epsilon_y = std::strtod( argv[2], stringPtr );
  std::cout << "Epsilon_x: " << epsilon_x << std::endl;
  std::cout << "Epsilon_y: " << epsilon_y << std::endl;

  //--------------
  // test data
  //--------------


  // NOTE ALEX: define mesh, qemesh types
  // then both pointtypes
  // with dim values of 2,3 and more
  // then with pixeltype of int, float, double
  // and test
  // use templated test function

  // for now, best case scenario, life is good
  typedef itk::Point< float, 2 > PointType;
  PointType mpa, mpb, mpc, mpd;
  mpa[0] = mpc[0] = mpc[1] = mpb[1] = 1.0;
  mpa[1] = mpb[0] =                   0.0;
  mpd[0] = epsilon_x;
  mpd[1] = epsilon_y;

  // if the triangle is defined clock wise,
  // and you do not take into accont orientation,
  // then the test does not work
  double ax = 1.0;
  double ay = 0.0;
  double cx = 1.0;
  double cy = 1.0;
  double bx = 0.0;
  double by = 1.0;
  double dx = epsilon_x;
  double dy = epsilon_y;
  std::cout << "dx: " << dx << std::endl;
  std::cout << "dy: " << dy << std::endl;


  //----------------
  // non exact test
  //---------------
  if( atoi( argv[3] ) == 0 )
    { 
    double det;

    // orientation test - determination of the sign of ad-bc
    double a = ax-cx;
    double b = ay-cy;
    double c = bx-cx;
    double d = by-cy;
    double orientation = a*d-b*c;


    // NOTE ALEX: this should be extracted/inferred from the mesh/point types
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

    // determinant computation - the result is multiplied by 'orientation'
    det = vnl_det( M.GetVnlMatrix() ) * orientation;
    std::cout << "Det(M): " << det << std::endl;

    return( det>=0?0:1);
    }

  //----------------
  // exact test
  //---------------
  else
    {
    return IsInside< PointType >( mpa, mpb, mpc, mpd );
    }

}
