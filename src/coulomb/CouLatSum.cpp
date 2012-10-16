/*
 * File: GridDen.cc
 * Author: Qiming Sun <osirpt.sun@gmail.com>
 */

#include <iostream>
#include <assert.h>
#include "CouLatSum.h"
#include "grid_3d_class.h"
//#include "../lib/CxAlgebra.h"
#include "fblas.h"

//extern "C" {
//#include <cblas.h>
//#include <stdio.h>
//#include <stdlib.h>
//}


namespace cls { // Coulomb Lattice Summation

struct maybe_env_t {
    const FSolidModel *p_solid;
    const FSuperCell *p_super_cell;
    const FLattice *p_lattice;
    const FUnitCell *p_unit_cell;
    const FBasisSet *p_bas;
    const FOpMatrix *p_den_mat;
    phf::coulomb_grid::grid_3d *p_grid3d;

    int nao_unit_cell;
    double den_cutoff;
    double background_chg;
    double grid_weight;
};

static struct maybe_env_t *_env;

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
    double *mv = new double[n_nu];
    double s;
    dgemv_(&TRANS_T, &n_mu, &n_nu, &D1, mat, &n_mu, mu, &INC1, &D0, mv, &INC1);
    s = ddot_(&n_nu, mv, &INC1, nu, &INC1);
    delete[] mv;
    return s;
}

/*
 * data flow:
 * -> cells_accumulator -> grids_accumulator ->
 */
static int cell_offset_by_cell_id(const int cell_id[3])
{
    int nx = _env->p_super_cell->Size[0];
    int ny = _env->p_super_cell->Size[1];
    int nz = _env->p_super_cell->Size[2];
    return cell_id[2] * ny * nz + cell_id[1] * nz + cell_id[0];
}
// loop over cells, and trace over mat to get scalar or from scalar to mat,
static void cells_accumulator(FCellIter_t func, FuncTerm_t fterm,
                              const int grid_id[3],
                              double *v, double *mat, double *ao)
{
    int cell_id[3];
    // move the 0-cell to the center of the supercell
    int nc_x0 = -_env->p_super_cell->Size[0]/2;
    int nc_y0 = -_env->p_super_cell->Size[1]/2;
    int nc_z0 = -_env->p_super_cell->Size[2]/2;
    int nc_x1 = (_env->p_super_cell->Size[0]-1)/2;
    int nc_y1 = (_env->p_super_cell->Size[1]-1)/2;
    int nc_z1 = (_env->p_super_cell->Size[2]-1)/2;
    for (int ix = nc_x0; ix < nc_x1; ix++)
        for (int iy = nc_y0; iy < nc_y1; iy++)
            for (int iz = nc_z0; iz < nc_z1; iz++) {
                cell_id[0] = ix;
                cell_id[1] = iy;
                cell_id[2] = iz;
                if (!(*fterm)(cell_id)) {
                    (*func)(grid_id, cell_id, v, mat, ao);
                }
            }
}

static void grids_accumulator(FGridIter_t func,
                              double *v, double *mat, double *ao)
{
    int grid_id[3];
    for (int ix = 0; ix < _env->p_grid3d->num_x; ix++)
        for (int iy = 0; iy < _env->p_grid3d->num_y; iy++)
            for (int iz = 0; iz < _env->p_grid3d->num_z; iz++) {
                grid_id[0] = ix;
                grid_id[1] = iy;
                grid_id[2] = iz;
                (*func)(grid_id, v, mat, ao);
                v++;
            }
}

/************************************************/
static void cell_coord_by_cell_id(const int cell_id[3], double cell_coord[3])
{
    cell_coord[0] = cell_id[0] * _env->p_lattice->T[0][0]
                  + cell_id[1] * _env->p_lattice->T[1][0]
                  + cell_id[2] * _env->p_lattice->T[2][0];
    cell_coord[1] = cell_id[0] * _env->p_lattice->T[0][1]
                  + cell_id[1] * _env->p_lattice->T[1][1]
                  + cell_id[2] * _env->p_lattice->T[2][1];
    cell_coord[2] = cell_id[0] * _env->p_lattice->T[0][2]
                  + cell_id[1] * _env->p_lattice->T[1][2]
                  + cell_id[2] * _env->p_lattice->T[2][2];
}
static void grid_coord_by_grid_id(const int grid_id[3], double grid_coord[3])
{
    grid_coord[0] = _env->p_grid3d->grid_x[grid_id[0]];
    grid_coord[1] = _env->p_grid3d->grid_y[grid_id[1]];
    grid_coord[2] = _env->p_grid3d->grid_z[grid_id[2]];
}
//FIXME ao_at_grid and ao_unit_cell to generate value of AOs on each grid
// the grid_coord should be constrained in 0-th cell
static double ao_at_grid(const int ao_id, const double grid_coord[3])
{
    double ao = 0;
    // call grid class
//    _env->p_bas->FD(eval_basis_fn_on_grid)();
//        (double *pOut, FORTINT const &nCompSt,
//   FORTINT *pCentersOut, FORTINT *pMap, FORTINT &nMap,
//   FBasisSet const &Basis, double  (*pGridPt)[3], FORTINT const &nGridPt,
//   FORTINT const &DerivOrder, FORTINT &iContext)
    return ao;
}

static void ao_unit_cell(const double grid_coord[3], const int cell_id[3],
                         double *ao)
{
    for (int i = 0; i < _env->nao_unit_cell; i++) {
        //ao[i] = ao_at_grid(grid_coord, cell_id);
    }
}
// for the grid in 0-th cell, value of all AOs in the super_cell
static void ao_super_cell(const double grid_coord[3], double *ao)
{
    int cell_id[3];
    for (int ix = 0; ix < _env->p_super_cell->Size[0]; ix++)
        for (int iy = 0; iy < _env->p_super_cell->Size[1]; iy++)
            for (int iz = 0; iz < _env->p_super_cell->Size[2]; iz++) {
                cell_id[0] = ix;
                cell_id[1] = iy;
                cell_id[2] = iz;
                //grid_cell_coord_by_cell_id(grid_coord, cell_id, coord);
                //cell_coord_by_cell_id(cell_id, coord);
                ao_unit_cell(grid_coord, cell_id, ao);
                ao += _env->nao_unit_cell;
    }
}


/************************************************/
static int is_cell_far(const int cell_id[3])
{
    return 0;
    // TODO:
    //if (cell_id > SUPER_CELL_INFO) {
    //    return 1; // remote cell
    //} else {
    //    return 0;
    //}
}

// for a given grid_id in 0-th cell,
// iterate this fn in cells_accumulator to run over all the image grids
static void density_cell0(const int grid_id[3], const int cell_id[3],
                          double *v, double *dmat, double *ao)
{
    //TODO: call grids_accumulator here to make grids inner-loop

    double grid_coord[3];
    grid_coord_by_grid_id(grid_id, grid_coord);
    double cell_coord[3];
    cell_coord_by_cell_id(cell_id, cell_coord);
    double coord[3];
    coord[0] = grid_coord[0] - cell_coord[0];
    coord[1] = grid_coord[1] - cell_coord[1];
    coord[2] = grid_coord[2] - cell_coord[2];
    ao_super_cell(coord, ao);

    int n_mu = _env->p_den_mat->nRows;
    int n_nu = _env->p_den_mat->nCols;
    *v += dtrace_vmv(n_mu, n_nu, ao, ao, dmat);
}
// constrain grid_id in 0-th cell, v is the value of density at the
// given grid
static void density_grid_iter(const int grid_id[3],
                                double *v, double *dmat, double *ao)
{
    cells_accumulator(density_cell0, is_cell_far,
                      grid_id, v, (double *)&_env->p_den_mat->Data[0], ao);
    *v += -_env->background_chg;
}

// for a given grid_id in 0-th cell,
// iterate this fn in cells_accumulator to run over all the image grids
static void coul_matrix_cell_iter(const int grid_id[3], const int cell_id[3],
                                  double *v, double *dmat, double *ao)
{
    double grid_coord[3];
    grid_coord_by_grid_id(grid_id, grid_coord);
    double cell_coord[3];
    cell_coord_by_cell_id(cell_id, cell_coord);
    double coord[3];
    coord[0] = grid_coord[0] + cell_coord[0];
    coord[1] = grid_coord[1] + cell_coord[1];
    coord[2] = grid_coord[2] + cell_coord[2];
    ao_super_cell(coord, ao);

    dmat += cell_offset_by_cell_id(cell_id)
        * _env->nao_unit_cell * _env->nao_unit_cell;
    int n_mu = _env->p_den_mat->nRows;
    int n_nu = _env->p_den_mat->nCols;
    const int INC1 = 1;
    double vw = *v * _env->grid_weight;
    dger_(&n_mu, &n_nu, &vw, ao, &INC1, ao, &INC1, dmat, &n_mu);
}

static void coul_matrix_grid_iter(const int grid_id[3],
                                  double *v, double *dmat, double *ao)
{
    // for any grids in 0-th cell, accumulate the contribution of image grids
    cells_accumulator(coul_matrix_cell_iter, is_cell_far,
                      grid_id, v, dmat, ao);
}

static void coul_matrix_acc(double *real_space_pot, double *vmat)
{
    double *ao = new double[_env->p_den_mat->nCols];
    // accumulate all the grids in 0-th cell
    grids_accumulator(coul_matrix_grid_iter, real_space_pot, vmat, ao);
    delete[] ao;
}



/*
 **************************************************
 * global functions
 */
int init_coul_fft(const FSolidModel& solid, const FOpMatrix& den_mat)
{
    _env->p_solid = &solid;
    _env->p_super_cell = &solid.SuperCell;
    _env->p_lattice = &solid.Lattice;
    _env->p_unit_cell = &solid.UnitCell;
    //_env->p_bas = &solid.UnitCell.OrbBasis;
    _env->p_den_mat = &den_mat;

    //_env->nao_unit_cell = solid.UnitCell.OrbBasis.nFn;
    _env->nao_unit_cell = den_mat.nCols;

    double nele = 0; //FIXME = total_number_of_electron_in_one_cell
    double vol = solid.UnitCell.Volume;
    set_background_charge(nele/vol);

    _env->p_grid3d = new phf::coulomb_grid::grid_3d();
    _env->grid_weight = vol / (_env->p_grid3d->num_x
                             * _env->p_grid3d->num_y
                             * _env->p_grid3d->num_z);
    return 0;
}

int del_coul_fft()
{
    delete _env->p_grid3d;
    return 0;
}

int set_den_cutoff(const double den_cutoff)
{
    cls::_env->den_cutoff = den_cutoff;
    return 0;
}

int set_background_charge(const double chg)
{
    cls::_env->background_chg = chg;
    return 0;
}

// grid_id of one cell
double density_at_grid(const int grid_id[3])
{
    double *ao = new double[_env->p_den_mat->nCols];
    double den;
    density_grid_iter(grid_id, &den, (double *)&_env->p_den_mat->Data[0], ao);
    delete[] ao;
    return den;
}

void density_unit_cell(double *density)
{
    double *ao = new double[_env->p_den_mat->nCols];
    cls::grids_accumulator(cls::density_grid_iter,
                           density, (double *)&_env->p_den_mat->Data[0], ao);
    delete[] ao;
}

// return the coulomb energy
double coul_matrix(const FSolidModel& solid, const FOpMatrix& den_mat,
                   FOpMatrix& coul_mat)
{
    init_coul_fft(solid, den_mat);

    int num_grid_in_unit_cell = _env->p_grid3d->num_x
                              * _env->p_grid3d->num_y
                              * _env->p_grid3d->num_z;
    double *real_space_density = new double[num_grid_in_unit_cell];
    double *real_space_pot = new double[num_grid_in_unit_cell];
    // calculate the density
    density_unit_cell(real_space_density);
    // FFT Poisson solver
    //TODO: poisson_solver(real_space_density, real_space_pot);
    // Coulomb matrix
    coul_matrix_acc(real_space_pot, &coul_mat[0]);
    delete[] real_space_density;
    delete[] real_space_pot;

    del_coul_fft();
    const int INC1 = 1;
    int n = coul_mat.size();
    int ncells = _env->p_super_cell->Size[0]
               * _env->p_super_cell->Size[1]
               * _env->p_super_cell->Size[2];
    return ddot_(&n, &coul_mat[0], &INC1, &den_mat[0], &INC1) / ncells;
}

} // end namespace cls
