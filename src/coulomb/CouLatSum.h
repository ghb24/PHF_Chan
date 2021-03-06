/*
 * File: CouLatSum.h
 * Author: Qiming Sun <osirpt.sun@gmail.com>
 */

#include "../PhfSolidDef.h"
#include "../PhfBasisSet.h"
#include "../PhfOperator.h"

namespace cls { // Coulomb Lattice Summation

typedef void (*FCellIter_t)(const int grid_id, const int *cell_id,
                            double *v, double *mat, double *ao);
typedef void (*FGridIter_t)(const int grid_id,
                            double *v, double *mat, double *ao);
typedef int (*FCellTruncate_t)(int *);
// value of density on the given grid, grid is indexed by grid_id
double density_at_grid(const int grid_id);
// value of density on each grid in a unit cell,
// double *density is an array with the length of the number of grids in
// a unit cell
void density_unit_cell(double *density);
void init_env(const FSolidModel& solid, const FOpMatrix& den_mat);
void del_env();



void set_num_grids_in_unit_cell(const int n_a1, const int n_a2, const int n_a3);

double coul_matrix(const FSolidModel& solid, const FOpMatrix& den_mat,
                   FOpMatrix& coul_mat);
}
