#include <iostream>
#include "grid_3d_class.h"

using namespace std;

namespace phf { namespace coulomb_grid {

void grid_3d::set_ngrid( int nx, int ny, int nz )
{
  this-> num_x = nx;
  this-> num_y = ny;
  this-> num_z = nz;
  this-> num_grid = num_x * num_y * num_z;
}

void grid_3d::set_grid_x( int nx, double lx, double ux)
{
  this->set_nx(nx);
  this->set_lx(lx);
  this->set_ux(ux);
  this-> incr_x = ( get_ux() - get_lx() )/ get_nx();
  double value = get_lx();
  for( int i = 0; i < get_nx(); i++ ){
    this->grid_x.push_back( value );
    value += incr_x;
  }

}

void grid_3d::set_grid_y( int ny, double ly, double uy)
{
  this->set_ny(ny);
  this->set_ly(ly);
  this->set_uy(uy);
  this-> incr_y = ( get_uy() - get_ly() )/ get_ny();
  double value = get_ly();
  for( int i = 0; i < get_ny(); i++ ){
    this->grid_y.push_back( value );
    value += incr_y;
  }

}

void grid_3d::set_grid_z( int nz, double lz, double uz)
{
  this->set_nz(nz);
  this->set_lz(lz);
  this->set_uz(uz);
  this-> incr_z = ( get_uz() - get_lz() )/ get_nz();
  double value = get_lz();
  for( int i = 0; i < get_nz(); i++ ){
    this->grid_z.push_back( value );
    value += incr_z;
  }
}

void grid_3d::init_gridv( int nx, int ny, int nz )
{
  this->set_ngrid( nx, ny, nz );
  double value = 1.5e0; 
  for( int i = 0; i < get_ngrid(); i++ ){
   grid_value.push_back( value );
  }
//  this->grid_value.resize( this->get_ngrid() );
}

void coord_grid()
{
 // To be done
}

double& grid_3d::get_gridv( int ix, int iy, int iz ){
 int pos = ix * get_ny() * get_nz() + iy * get_nz() + iz;
 return this-> grid_value.at( pos );
}

double& grid_3d::operator() ( int ix, int iy, int iz )
{
 return this->get_gridv( ix, iy, iz );
}

grid_3d::grid_3d( int n )
{
 // to be done
}

grid_3d::grid_3d( int n, double c)
{
  // to be done
}

grid_3d::grid_3d( int nx, int ny, int nz)
{
 // to be done
}

grid_3d::grid_3d( int nx, int ny, int nz, double c )
{ 
  // to be done
} 

grid_3d::grid_3d( int nx, int ny, int nz, double cx, double cy, double cz )
{
 // to be done
}

grid_3d::grid_3d( int nx, double lx, double ux,
                  int ny, double ly, double uy,
                  int nz, double lz, double uz )
{

  cout << " constructor grid_3d::grid_3d( nx, lx, ux, ny, ly, uy, nz, lz, yz )" << endl;
  this->set_grid_x( nx, lx, ux );
  this->set_grid_y( ny, ly, uy );
  this->set_grid_z( nx, lz, uz );

  cout.precision(10);
  cout << " check : " << endl;
  cout << " nx = " << get_nx() << " lower_x = " << get_lx() << " upper_x = " << get_ux() << " increment x = " << get_incrx() << endl;
  cout << " ny = " << get_ny() << " lower_y = " << get_ly() << " upper_y = " << get_uy() << " increment y = " << get_incry() << endl;
  cout << " nz = " << get_nz() << " lower_z = " << get_lz() << " upper_z = " << get_uz() << " increment z = " << get_incrz() << endl;

  this->init_gridv( get_nx(), get_ny(), get_nz() );

  cout << " n_grid = " << get_ngrid() << endl;

}

} } // end of phf::coulomb_grid
