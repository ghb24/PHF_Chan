#include <vector>

namespace phf { namespace coulumb_grid { 

class grid_3d
{
// Constructors and destructors
public:
 grid_3d(){};
 grid_3d( int n );
 grid_3d( int n, double c);
 grid_3d( int nx, int ny, int nz );
 grid_3d( int nx, int ny, int nz, double c );
 grid_3d( int nx, int ny, int nz, double cx, double cy, double cz );
 grid_3d( int nx, double lx, double ux,
          int ny, double ly, double uy,
	  int nz, double lz, double uz); // main contructor to use, feed it with the n, upper limit and lower limit 

 ~grid_3d(){};

public:
 int num_x, num_y, num_z;  // number of the points in the grid for x, y, z
 double incr_x, incr_y, incr_z; // the increment in the grid of the x, y, z directions
 
 double lower_x, upper_x; // upper and lower limits of x grid
 double lower_y, upper_y; // upper and lower limits of y grid
 double lower_z, upper_z; // upper and lower limits of z grid

public:
 std::vector<double> grid_x; // the actual grid of x, size = num_x
 std::vector<double> grid_y; // the actual grid of y, size = num_y
 std::vector<double> grid_z; // the actual grid of z, size = num_z

public:
 int num_grid;
 std::vector<double> grid_value; // stores the function values f(x,y,z), size = num_x * num_y * num_z

 std::vector<std::vector<double> > coord;
// functions
public:
 // set and return the size of the grid
 void set_nx( int nx ){ this-> num_x = nx; }
 void set_ny( int ny ){ this-> num_y = ny; }
 void set_nz( int nz ){ this-> num_z = nz; }

 int get_nx(){ return this-> num_x; }
 int get_ny(){ return this-> num_y; }
 int get_nz(){ return this-> num_z; }

 // set and return the increment of grid
 void set_incrx( double incr ){ this-> incr_x = incr; }
 void set_incry( double incr ){ this-> incr_y = incr; }
 void set_incrz( double incr ){ this-> incr_z = incr; }

 double get_incrx(){ return this-> incr_x; }
 double get_incry(){ return this-> incr_y; }
 double get_incrz(){ return this-> incr_z; }

// set and return the upper and lower limits of the grid
 void set_lx( double lx ){ this-> lower_x = lx; }
 void set_ux( double ux ){ this-> upper_x = ux; }
 void set_ly( double ly ){ this-> lower_y = ly; }
 void set_uy( double uy ){ this-> upper_y = uy; }
 void set_lz( double lz ){ this-> lower_z = lz; }
 void set_uz( double uz ){ this-> upper_z = uz; }
 
 double get_lx(){ return this-> lower_x; }
 double get_ux(){ return this-> upper_x; }
 double get_ly(){ return this-> lower_y; }
 double get_uy(){ return this-> upper_y; }
 double get_lz(){ return this-> lower_z; }
 double get_uz(){ return this-> upper_z; }

// get the coordinate
 double get_x( int ix ){ return grid_x.at(ix); }
 double get_y( int iy ){ return grid_y.at(iy); }
 double get_z( int iz ){ return grid_z.at(iz); }

// set num of the grid, this thing will also be done in the constructor
 void set_ngrid( int nx, int ny, int nz );
 int get_ngrid(){ return this->num_grid; }

// set the actual x, y, z grid which are stored in the std::vector<double>
 void set_grid_x( int nx, double lx, double ux);
 void set_grid_y( int ny, double ly, double uy);
 void set_grid_z( int nz, double lz, double uz);

// set the size of the f(x, y, z)
 void init_gridv( int nx, int ny, int nz);   

 void coord_grid();

// set the value of the f(x, y, z)
 double& get_gridv( int ix, int iy, int iz );
 double& operator() (int ix, int iy, int iz);

 
};


} } // end of phf::coulumb_grid
