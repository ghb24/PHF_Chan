#include "grid_3d_class.h"
#include <iostream>

using namespace std;

int main( int argc, char* argv[] )
{
  int nx = 4;
  int ny = 4; 
  int nz = 4;
  double lx = -1.0e0;
  double ux = 1.0e0;
  double ly = -1.0e0;
  double uy = 1.0e0;
  double lz = -1.0e0;
  double uz = 1.0e0;
  phf::coulumb_grid::grid_3d func( nx, lx, ux, 
		                  ny, ly, uy,
				  nz, lz, uz);

  // the way to get the coordinate ix = 2;
  int ix = 2;
  cout << func.get_x(ix) << endl;


  // I set all the grids with initial value 1.5
  cout << func(0, 0, 0) + func( 0,1,2) << endl; // so it should be 1.5+1.5 = 3

  // I change this f(0,1,2) to 2.0
  func(0,1,2) = 2.0e0;
  cout << func(0,1,2) << endl;
  // so it's now 2.0
  // CONCLUSION: To access this ixth, iyth and izth value of the func, just use func( ix, iy, iz)

}
