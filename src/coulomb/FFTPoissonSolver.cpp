#include <fftw3.h>
#include <iostream>
#include <cmath>
#include "grid_3d_class.h"

extern "C" {
    void dscal_(const int *n, const double *da, double *dx, const int *incx);
}

namespace cls { // Coulomb Lattice Summation
static double k_dot_k(double k[3])
{
    return k[0] * k[0] + k[1] * k[1] + k[2] * k[2];
}

/*
 **************************************************
 * for Poisson equation FFT solver
 * size of array kcoords >= size_of [grid3d.num_grid] * 3
 */
static void kspace_grids_generate(const phf::coulomb_grid::grid_3d& grid3d,
                                  const FLattice& lattice,
                                  double *kcoords)
{
    int x, y, z;
    for (int ix = 0; ix < grid3d.num_x; ix++) {
        if (ix < grid3d.num_x/2) {
            x = ix;
        } else {
            x = ix - grid3d.num_x;
        }
        for (int iy = 0; iy < grid3d.num_y; iy++) {
            if (iy < grid3d.num_y/2) {
                y = iy;
            } else {
                y = iy - grid3d.num_y;
            }
            for (int iz = 0; iz < grid3d.num_z/2+1; iz++) {
                if (iz < grid3d.num_z/2) {
                    z = iz;
                } else {
                    z = iz - grid3d.num_z;
                }
                kcoords[0] = lattice.K[0][0] * x
                           + lattice.K[1][0] * y
                           + lattice.K[2][0] * z;
                kcoords[1] = lattice.K[0][1] * x
                           + lattice.K[1][1] * y
                           + lattice.K[2][1] * z;
                kcoords[2] = lattice.K[0][2] * x
                           + lattice.K[1][2] * y
                           + lattice.K[2][2] * z;
                kcoords += 3;
            }
        }
    }
}


int poisson_kspace_solver(double* density,double* potential,
                          const phf::coulomb_grid::grid_3d& grid3d,
                          const FLattice& lattice)
{
    int nx = grid3d.num_x;
    int ny = grid3d.num_y;
    int nz = grid3d.num_z;
    fftw_complex *k_den = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * (nx * ny * (nz/2+1)));
    fftw_plan p = fftw_plan_dft_r2c_3d(nx, ny, nz, density, k_den, FFTW_ESTIMATE);
    fftw_execute(p);
    fftw_destroy_plan(p);

    double k2_inv;
    std::vector<double> k_grids(grid3d.num_grid * 3);
    kspace_grids_generate(grid3d, lattice, &k_grids[0]);
    for (int k_id = 1; k_id < nx * ny * (nz/2+1); k_id++) {
        k2_inv = 4 * M_PI / k_dot_k(&k_grids[k_id * 3]);
        k_den[k_id][0] = k_den[k_id][0] * k2_inv;
        k_den[k_id][1] = k_den[k_id][1] * k2_inv;
    }
    k_den[0][0] = 0;
    k_den[0][1] = 0;

    p = fftw_plan_dft_c2r_3d(nx, ny, nz, k_den, potential, FFTW_ESTIMATE);
    fftw_execute(p);
    fftw_destroy_plan(p);
    int INC1 = 1;
    double n_inv = 1.l / grid3d.num_grid;
    dscal_(&grid3d.num_grid, &n_inv, potential, &INC1);

    fftw_free(k_den);
    return 1;
}

} // end namespace cls
