/*
 * File: GridDen.cc
 * Author: Qiming Sun <osirpt.sun@gmail.com>
 */

#include <iostream>
#include <assert.h>
#include "GridDen.h"
#include "grid_3d_class.h"
#include "fblas.h"

//extern "C" {
//#include <cblas.h>
//#include <stdio.h>
//#include <stdlib.h>
//}


namespace cls { // Coulomb Lattice Summation

static struct MAYBE_ENV_T { 
    den_cutoff;
    background_chg;
    double cell_parameter[9]; // three vectors in lattice
    phf::coulomb::grid_3d *ptr_grid3d;
    FBasisSet *ptr_orb_basis;
    FUnitCell *ptr_unit_cell;
    FLattice *ptr_lattice;
};

static struct MAYBE_ENV_T _env;

/*
 **************************************************
 * local funcs
 */
// sum mu[n_mu] * mat[n_mu,n_nu] * nu[n_nu]
static double dtrace_vmv(int n_mu, int n_nu, double *mu, double *nu, double *mat)
{
    const int INC1 = 1;
    const double D0 = 0;
    const double D1 = 1;
    const char TRANS_T = 'T';
    mv = new double[n_nu];
    double s;
    dgemv_(&TRANS_T, &n_mu, &n_nu, &D1, mat, &n_mu, mu, &INC1, &D0, mv, &INC1);
    s = ddot_(&n_nu, mv, &INC1, nv, &INC1);
    delete mv[];
    return s;
}

// loop over cells, and trace over mat to get scalar or from scalar to mat,
static void cells_accumulator(FuncIter_t func, FuncTerm_t fterm,
                              const int *grid_id,
                              double *v, double *mat, double *mu, double *nu)
{
    int *cell_id = 0;
    while (!(*fterm)(cell_id)) {
        (*func)(grid_id, cell_id, v, mat, mu, nu);
        cell_id++;
    }
}

static void grids_accumulator(FuncIter_t func,
                              double *v, double *mat, double *mu, double *nu)
{
    for (int *grid_id = 0; grid_id < tot_grids; grid_id++) {
        (*func)(gird_id, NONE v, mat, mu, nu);
    }
}

/************************************************/
static double ao_at_grid(const int ao_id, const int *grid_id, const int *cell_id)
{
    double ao = 0;
    // call grid class
    return ao;
}
static int aos_on_cell(const int *grid_id, const int *cell_id, double *aos)
{
    for (int mu = 0; mu < num_ao; mu++) {
        aos[mu] = ao_at_grid(grid_id, cell_id);
    }
    return 0;
}

static int number_of_aos_in_one_cell()
{
    return 0;
}


/************************************************/
static int is_cell_far(const int *cell_id)
{
    if (cell_id > SUPER_CELL_INFO) {
        return 1; // remote cell
    } else {
        return 0;
    }
}

static void density_grid_iter(const int *grid_id, const int *cell_id,
                              double *v, double *dmat, double *mu, double *nu)
{
    //TODO: call grids_accumulator here to make grids inner-loop
    aos_on_cell(grid_id, cell_id, ao_nu);
    p_celldmat = density_matrix_pannel(cell_id);
    *v += dtrace_vmv(num_ao, num_ao, ao_mu, ao_nu, p_celldmat);
}
static double density_at_grid(const int *grid_id, double *dmat,
                              double *v, double *dmat, double *mu, double *nu)
{
    int num_ao = number_of_aos_in_one_cell();
    double *ao_mu = new double[num_ao];
    double *ao_nu = new double[num_ao];
    double den;
    aos_on_cell(grid_id, cell_id, ao_mu);
    cells_accumulator(density_grid_iter, is_cell_far,
                      grid_id, &den, dmat, ao_mu, ao_nu);
    delete ao_mu[];
    delete ao_nu[];
    return den;
}

static int coul_matrix_grid_iter(const int *grid_id, cell_id,
                                 double *v, double *dmat, double *mu, double *nu)
{
    const int INC1 = 1;
    const double D0 = 0;
    const double D1 = 1;
    const char TRANS_T = 'T';
    int num_ao = number_of_aos_in_one_cell();
    double *ao_mu = new double[num_ao];
    double *ao_nu = new double[num_ao];
    dger_(&num_ao, &num_ao, v, mu, &INC1, nu, &INC1, dmat, ao_mu);
    delete ao_mu[];
    delete ao_nu[];
    return den;
}
static int coul_matrix_cell_iter(const int *grid_id, const int *cell_id,
                                 double *v, double *mat)
{
    int num_ao = number_of_aos_in_one_cell();
    double *ao_mu = new double[num_ao];
    double *ao_nu = new double[num_ao];
    double den;
    grids_accumulator(coul_matrix_grid_iter, v, mat, ao_mu, ao_nu);
    delete ao_mu[];
    delete ao_nu[];
    return den;
}

static int coul_matrix(double *v, double *mat)
{
    int num_ao = number_of_aos_in_one_cell();
    double *ao_mu = new double[num_ao];
    double *ao_nu = new double[num_ao];
    aos_on_cell(grid_id, cell_id, ao_mu);
    cells_accumulator(coul_matrix_cell_iter, is_cell_far,
                      grid_id, &den, dmat, ao_mu, ao_nu);
    delete ao_mu[];
    delete ao_nu[];
    return 0;
}



/*
 **************************************************
 * global functions
 */
int init_coul_fft()
{
    _env.ptr_grid3d = new phf::coulomb::grid_3d();
    nao = pOrbBasis->nBasisFn(),
}

int del_coul_fft()
{
    delete _env.ptr_grid3d;
}

int set_den_cutoff(const double den_cutoff)
{
    cls::_env.den_cutoff = den_cutoff;
    return 0
}

int set_background_charge(const double chg)
{
    cls::_env.background_chg = chg;
    return 0
}

// grid_id of one cell
double density_at_grid(const int *grid_id)
{
    return cls::density_at_grid(grid_id)
        - cls::_env.background_chg;
}

} // end namespace cls
