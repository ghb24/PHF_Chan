/*
 * File: CouLatSum.h
 * Author: Qiming Sun <osirpt.sun@gmail.com>
 */

#include "../PhfSolidDef.h"
#include "../PhfBasisSet.h"
#include "../PhfOperator.h"
#include "grid_3d_class.h"

namespace cls { // Coulomb Lattice Summation

// value of density on each grid in a unit cell,
// double *density is an array with the length of the number of grids in
// a unit cell
void density_unit_cell(double *density);
void init_env(const FSolidModel& solid, const FOpMatrix& den_mat);
void del_env();


int poisson_kspace_solver(double* density,double* potential,
                          const phf::coulomb_grid::grid_3d& grid3d,
                          const FLattice& lattice);

void set_num_grids_in_unit_cell(const int n_a1, const int n_a2, const int n_a3);

double coul_matrix(const FSolidModel& solid, const FOpMatrix& den_mat,
                   FOpMatrix& coul_mat);
}
