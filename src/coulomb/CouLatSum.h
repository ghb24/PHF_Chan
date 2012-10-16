/*
 * File: CouLatSum.h
 * Author: Qiming Sun <osirpt.sun@gmail.com>
 */

#include "../PhfSolidDef.h"
#include "../PhfBasisSet.h"
#include "../PhfOperator.h"

namespace cls { // Coulomb Lattice Summation
typedef void (*FCellIter_t)(const int *grid_id, const int *cell_id,
                            double *v, double *mat, double *ao);
typedef void (*FGridIter_t)(const int *grid_id,
                            double *v, double *mat, double *ao);
typedef int (*FuncTerm_t)(const int *);


int init_coul_fft(const FSolidModel& solid, const FOpMatrix& den_mat);

int del_coul_fft();

int set_den_cutoff(const double den_cutoff);

int set_background_charge(const double chg);

double density_at_grid(const int grid_id[3]);

void density_unit_cell(double *density);

double coul_matrix(const FSolidModel& solid, const FOpMatrix& den_mat,
                   FOpMatrix& coul_mat);

}
