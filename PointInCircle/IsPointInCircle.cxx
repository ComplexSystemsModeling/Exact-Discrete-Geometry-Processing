#include <itkMatrix.h>
#include <itkPoint.h>
#include <itkMesh.h>
#include <itkQuadEdgeMesh.h>
#include <itkTriangleCell.h>
#include <itkCellInterface.h>

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
  const PointType& TrianglePoint1,
  const PointType& TrianglePoint2,
  const PointType& TrianglePoint3,
  const PointType& PointToTest )
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
//  Non exact, ITK native version
//
template< typename RealType >
bool
IsInsideNotExact(
  const RealType& ax,
  const RealType& ay,
  const RealType& bx,
  const RealType& by,
  const RealType& cx,
  const RealType& cy,
  const RealType& dx,
  const RealType& dy )
{
  double det;

  // orientation test - determination of the sign of ad-bc
  double a = ax-cx;
  double b = ay-cy;
  double c = bx-cx;
  double d = by-cy;
  double orientation = a*d-b*c;

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

  return( det>=0?0:1 );
}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
// Convenience Functor
//
// fetch triangle from the container,
// fetch pointsIds from triangleCell
// fetch points from point container
// and test
//
template<
  typename MeshType, 
  typename MeshPointerType, 
  typename CellIdentifier, 
  typename PointType >
bool
TestPointInTriangleInMesh(
  MeshPointerType myMesh,
  const CellIdentifier& TriangleId,
  const PointType& PointToTest,
  bool ExactTest
  )
{
  typedef typename MeshType::CellType   CellType; // abstract
  typedef typename CellType::CellAutoPointer CellAutoPointer; // abstract

  typedef typename MeshType::PixelType  RealType;
  typedef typename MeshType::CellTraits CellTraits;

  typedef itk::CellInterface< RealType, CellTraits > CellInterfaceType;
  typedef itk::TriangleCell< CellInterfaceType >     TriangleCellType;
  typedef typename TriangleCellType::PointIdIterator PointIdIterator;

  PointType mpa, mpb, mpc;

  CellAutoPointer cellIterator;
  if( myMesh->GetCell( TriangleId, cellIterator ) )
    {
    if( cellIterator->GetType() == 2 )
      { 
      // we have a triangle, let s get the points
      PointIdIterator pointIdIterator = cellIterator->PointIdsBegin();
      // should check the return value
      myMesh->GetPoint( *pointIdIterator, &mpa );
      pointIdIterator++;
      myMesh->GetPoint( *pointIdIterator, &mpb );
      pointIdIterator++;
      myMesh->GetPoint( *pointIdIterator, &mpc );
      if( !ExactTest )
        {
        return IsInsideNotExact( //< RealType >(
          mpa[0], mpa[1],
          mpb[0], mpb[1],
          mpc[0], mpc[1],
          PointToTest[0], PointToTest[1] );
        }
      else
        {
        return IsInside( mpa, mpb, mpc, PointToTest);
        // result_2 = IsInside< PointType >( mpa, mpb, mpc, mpd );
        }
      }
    else
      {
      // oops - this is not a triangle
      return( false );
      }
    }
  else
   {
   // oops, the cell ID was not found in the container
   return( false );
   }
}
//------------------------------------------------------------------------------



//------------------------------------------------------------------------------
// test function
//
// main subroutine of the test
// for a given meshtype,
// test the given point against the triangle, either given as 3 points
// or as a triangle
// illustrate the code using the itk APIs
//
template< typename MeshType >
bool
test(
  const double& eps_x,
  const double& eps_y,
  const bool& ExactTest
  )
{
  // sanity check
  // check dimension, get out if <1 and issue a warning if >2

  // infer types from template
  typedef typename MeshType::PixelType  RealType;
  typedef typename MeshType::PointType  PointType;
  typedef typename MeshType::CellType   CellType; // abstract
  typedef typename MeshType::CellTraits CellTraits;

  typedef typename MeshType::CellIdentifier  CellIdentifier;
  typedef typename MeshType::PointIdentifier PointIdentifier;

  typedef typename CellType::CellAutoPointer CellAutoPointer; // abstract

  typedef itk::CellInterface< RealType, CellTraits > CellInterfaceType;
  typedef itk::TriangleCell< CellInterfaceType >     TriangleCellType;
  typedef typename TriangleCellType::PointIdIterator PointIdIterator;

  //--------------------------------
  //  direct method, just use points
  //--------------------------------

  bool result_1 = false;

  PointType mpa, mpb, mpc, mpd;
  mpa[0] = mpc[0] = mpc[1] = mpb[1] = 1.0;
  mpa[1] = mpb[0] =                   0.0;
  mpd[0] = eps_x;
  mpd[1] = eps_y;

  if( !ExactTest )
    {
    result_1 = IsInsideNotExact(
      mpa[0], mpa[1],
      mpb[0], mpb[1],
      mpc[0], mpc[1],
      mpd[0], mpd[1] );
    }
  else
    {
    result_1 = IsInside( mpa, mpb, mpc, mpd );
    }

  //----------------------------------
  // closer to real case, use triangle
  //----------------------------------

  bool result_2 = false;

  // create an empty mesh
  typedef typename MeshType::Pointer MeshPointerType;
  MeshPointerType mesh = MeshType::New();

  // add points
  mesh->SetPoint( 0, mpa );
  mesh->SetPoint( 1, mpb );
  mesh->SetPoint( 2, mpc );

  // add cell
  CellAutoPointer dummyAbstractCell;
  TriangleCellType * dummyCell = new TriangleCellType();
  dummyAbstractCell.TakeOwnership( dummyCell ); // polymorphism
  PointIdentifier dummyCellPoints[3] = {0,1,2};
  dummyAbstractCell->SetPointIds( dummyCellPoints );
  mesh->SetCell( 0, dummyAbstractCell ); // invalidate the cell

  // now this code below should be close to what the user will do
  CellIdentifier myTriangle = 0;
  result_2 = TestPointInTriangleInMesh< MeshType >( mesh, myTriangle, mpd, ExactTest );

  return( result_1 | result_2 );

}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
// main
//
// Test for the given epsilon values, exact or not depending on user
// for a collection of different types
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

  bool TestExact = atoi(argv[3]);
  std::cout << "Exact Testing? " << TestExact << std::endl;

  //--------------
  // test data
  //--------------

  return(
      test< itk::Mesh< float,  2 > >( epsilon_x, epsilon_y, TestExact )
    | test< itk::Mesh< float,  3 > >( epsilon_x, epsilon_y, TestExact )
    | test< itk::Mesh< float,  4 > >( epsilon_x, epsilon_y, TestExact )
    | test< itk::Mesh< double, 2 > >( epsilon_x, epsilon_y, TestExact )
    | test< itk::Mesh< double, 3 > >( epsilon_x, epsilon_y, TestExact )
    | test< itk::Mesh< double, 4 > >( epsilon_x, epsilon_y, TestExact )
    | test< itk::QuadEdgeMesh< float,  2 > >( epsilon_x, epsilon_y, TestExact )
    | test< itk::QuadEdgeMesh< float,  3 > >( epsilon_x, epsilon_y, TestExact )
    | test< itk::QuadEdgeMesh< float,  4 > >( epsilon_x, epsilon_y, TestExact )
    | test< itk::QuadEdgeMesh< double, 2 > >( epsilon_x, epsilon_y, TestExact )
    | test< itk::QuadEdgeMesh< double, 3 > >( epsilon_x, epsilon_y, TestExact )
    | test< itk::QuadEdgeMesh< double, 4 > >( epsilon_x, epsilon_y, TestExact )
    );
}

