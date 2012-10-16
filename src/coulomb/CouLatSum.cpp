/*
 * File: CouLatSum.cpp
 * Author: Qiming Sun <osirpt.sun@gmail.com>
 */

#include <iostream>
#include <assert.h>
#include "../lib/CtIo.h"
#include "CouLatSum.h"
#include "nfft.h"
//#include "grid_3d_class.h"
//#include "../lib/CxAlgebra.h"
#include "fblas.h"

#define DEBUG

#if defined DEBUG
#include <math.h>
//#define LOGGER(X,Y)    xout << "DEBUG:" << X << std::endl
#define LOGGER(X,Y)     std::cout << #X << ":" << Y << std::endl
#else
#define LOGGER(X,Y)
#endif


namespace cls { // Coulomb Lattice Summation

struct maybe_env_t {
    FSolidModel *p1;
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

static int cell_offset_by_cell_id(const int cell_id[3])
{
    //int nx = _env->p_super_cell->Size[0];
    int ny = _env->p_super_cell->Size[1];
    int nz = _env->p_super_cell->Size[2];
    return cell_id[2] * ny * nz + cell_id[1] * nz + cell_id[0];
}
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

/*
 * data flow:
 * -> cells_accumulator -> grids_accumulator ->
 */
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
//FIXME ao_at_grid and ao_unit_cell to generate value of AOs on each grid
// the grid_coord should be constrained in 0-th cell
static double ao_at_grid(const int ao_id, const double grid_coord[3])
{
    double ao = 0.01;
    // call grid class
//    _env->p_bas->FD(eval_basis_fn_on_grid)();
//        (double *pOut, FORTINT const &nCompSt,
//   FORTINT *pCentersOut, FORTINT *pMap, FORTINT &nMap,
//   FBasisSet const &Basis, double  (*pGridPt)[3], FORTINT const &nGridPt,
//   FORTINT const &DerivOrder, FORTINT &iContext)
    return ao;
}

// for the grid in 0-th cell, value of all AOs in the unit cell
static void ao_unit_cell(const double grid_coord[3], double *ao)
{
    for (int i = 0; i < _env->nao_unit_cell; i++) {
        ao[i] = ao_at_grid(i, grid_coord);
    }
}

// for the grid in 0-th cell, value of all AOs in the super cell
static void ao_super_cell(const double grid_coord[3], double *ao)
{
    int cell_id[3];
    double coord[3];
    for (int ix = 0; ix < _env->p_super_cell->Size[0]; ix++)
        for (int iy = 0; iy < _env->p_super_cell->Size[1]; iy++)
            for (int iz = 0; iz < _env->p_super_cell->Size[2]; iz++) {
                cell_id[0] = ix;
                cell_id[1] = iy;
                cell_id[2] = iz;
                cell_coord_by_cell_id(cell_id, coord);
                coord[0] += grid_coord[0];
                coord[1] += grid_coord[1];
                coord[2] += grid_coord[2];
                ao_unit_cell(grid_coord, ao);
                ao += _env->nao_unit_cell;
    }
}


/************************************************/
static int is_cell_far(const int cell_id[3])
{
    return 0;
    // TODO:
    //if (SUPER_CELL_INFO or cutoff_info) {
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

    double cell_coord[3];
    double coord[3];
    grid_coord_by_grid_id(grid_id, coord);
    cell_coord_by_cell_id(cell_id, cell_coord);
    coord[0] += -cell_coord[0];
    coord[1] += -cell_coord[1];
    coord[2] += -cell_coord[2];
    ao_super_cell(coord, ao);

    int n_mu = _env->p_den_mat->nRows;
    int n_nu = _env->p_den_mat->nCols;
    *v += dtrace_vmv(n_mu, n_nu, ao, ao, dmat);
}

// grid_id is constrained in 0-th cell,
// v is the value of density corresponding to the given grid_id
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
    double cell_coord[3];
    double coord[3];
    grid_coord_by_grid_id(grid_id, coord);
    cell_coord_by_cell_id(cell_id, cell_coord);
    coord[0] += cell_coord[0];
    coord[1] += cell_coord[1];
    coord[2] += cell_coord[2];
    ao_super_cell(coord, ao);

    dmat += cell_offset_by_cell_id(cell_id)
        * _env->nao_unit_cell * _env->nao_unit_cell;
    int n_mu = _env->p_den_mat->nRows;
    int n_nu = _env->p_den_mat->nCols;
    const int INC1 = 1;
    double vw = *v * _env->grid_weight;
    dger_(&n_mu, &n_nu, &vw, ao, &INC1, ao, &INC1, dmat, &n_mu);
}

// iterate this fn in grids_accumulator to run over all grids in the unit cell
static void coul_matrix_grid_iter(const int grid_id[3],
                                  double *v, double *dmat, double *ao)
{
    // for any grids in 0-th cell, accumulate the contribution of image grids
    cells_accumulator(coul_matrix_cell_iter, is_cell_far,
                      grid_id, v, dmat, ao);
}

static void coul_matrix_acc(double *real_space_pot, double *vmat)
{
    LOGGER(DEBUG, "coul_matrix_acc starts");
    double *ao = new double[_env->p_den_mat->nCols];
    // accumulate all the grids in 0-th cell
    grids_accumulator(coul_matrix_grid_iter, real_space_pot, vmat, ao);
    delete[] ao;
}



/*
 **************************************************
 * global functions
 */
void init_env(const FSolidModel& solid, const FOpMatrix& den_mat)
{
    LOGGER(DEBUG, "Init cls_env");
    _env = new maybe_env_t;

    _env->p_solid = &solid;
    _env->p_super_cell = &solid.SuperCell;
    _env->p_lattice = &solid.Lattice;
    _env->p_unit_cell = &solid.UnitCell;
    _env->p_bas = &solid.UnitCell.OrbBasis;
    _env->p_den_mat = &den_mat;

    //_env->nao_unit_cell = solid.UnitCell.OrbBasis.nFn;
    _env->nao_unit_cell = den_mat.nRows;
    LOGGER(DEBUG, "nao_unit_cell = " << _env->nao_unit_cell);

    double nele = 0;
    for (TArray<FORTINT>::const_iterator n = solid.UnitCell.Elements.begin();
         n != solid.UnitCell.Elements.end(); n++) {
        nele += *n;
    }
    LOGGER(DEBUG, "nele in unit cell = " << nele);

    double vol = solid.UnitCell.Volume;
    set_background_charge(nele/vol);
    LOGGER(DEBUG, "unit cell vol = " << vol);
    LOGGER(DEBUG, "unit cell average charge = " << _env->background_chg);

    // set up mesh grids
    _env->p_grid3d = new phf::coulomb_grid::grid_3d();
    _env->p_grid3d->set_ngrid(2, 3, 4);
    // FIXME
    _env->p_grid3d->set_grid_x(_env->p_grid3d->num_x, 0, _env->p_lattice->T[0][0]);
    _env->p_grid3d->set_grid_y(_env->p_grid3d->num_y, 0, _env->p_lattice->T[1][1]);
    _env->p_grid3d->set_grid_z(_env->p_grid3d->num_z, 0, _env->p_lattice->T[2][2]);
    LOGGER(DEBUG, "unit cell grid numbers (x,y,z,total) = ("
           << _env->p_grid3d->num_x << ","
           << _env->p_grid3d->num_y << ","
           << _env->p_grid3d->num_z << ","
           << _env->p_grid3d->num_grid << ")");

    _env->grid_weight = vol / _env->p_grid3d->num_grid;
    LOGGER(DEBUG, "unit cell grid weight = " << _env->grid_weight);
    LOGGER(DEBUG, "init_env ends");
}

void del_env()
{
    delete _env->p_grid3d;
    delete _env;
}

// TODO: density is obtained by looping over all the unit cells in the
// super cell, cutoff may help to reduce the number of unit cells to be
// looped
void set_den_cutoff(const double den_cutoff)
{
    cls::_env->den_cutoff = den_cutoff;
}

void set_background_charge(const double chg)
{
    cls::_env->background_chg = chg;
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

//void density_unit_cell(double *density)
//{
//    double *ao = new double[_env->p_den_mat->nCols];
//    cls::grids_accumulator(cls::density_grid_iter,
//                           density, (double *)&_env->p_den_mat->Data[0], ao);
//    delete[] ao;
//}

double density_unit_cell(double *density)
{
    LOGGER(DEBUG, "density on each grid starts");
    double nele = 0;
    int grid_id[3];
    for (int ix = 0; ix < _env->p_grid3d->num_x; ix++)
        for (int iy = 0; iy < _env->p_grid3d->num_y; iy++)
            for (int iz = 0; iz < _env->p_grid3d->num_z; iz++) {
                grid_id[0] = ix;
                grid_id[1] = iy;
                grid_id[2] = iz;
                *density = density_at_grid(grid_id);
                nele += _env->grid_weight * *density;
                density++;
            }
    LOGGER(INFO, "sum(weight * density) => number of nele = " << nele);
#if defined DEBUG
    LOGGER(DEBUG, "density_unit_cell ends");
#else
    if (abs(nele) > 1e-5) {
        throw "Incorrect number of electron";
    }
#endif
    return nele;
}

// return the coulomb energy
double coul_matrix(const FSolidModel& solid, const FOpMatrix& den_mat,
                   FOpMatrix& coul_mat)
{
    LOGGER(DEBUG, "coul_matrix starts");
#if defined DEBUG
    FOpMatrix dm_tmp = den_mat;
    for (unsigned int i =0; i<dm_tmp.size(); i++) {
        dm_tmp[i] = cos(i) * .1;
    }
    init_env(solid, dm_tmp);
#else
    init_env(solid, den_mat);
#endif

    double *real_space_density = new double[_env->p_grid3d->num_grid];
    double *real_space_pot = new double[_env->p_grid3d->num_grid];
    // calculate the density
    density_unit_cell(real_space_density);
    // FFT Poisson solver
    LOGGER(DEBUG, "call poisson solver");
    poisson_kspace_solver(real_space_density, real_space_pot, _env->p_grid3d);
    // Coulomb matrix
    coul_matrix_acc(real_space_pot, &coul_mat[0]);
    delete[] real_space_density;
    delete[] real_space_pot;

    del_env();
    const int INC1 = 1;
    int n = coul_mat.size();
    int ncells = _env->p_super_cell->Size[0]
               * _env->p_super_cell->Size[1]
               * _env->p_super_cell->Size[2];
    double e = ddot_(&n, &coul_mat[0], &INC1, &den_mat[0], &INC1) / ncells;
    LOGGER(INFO, "e_coul = " << e);
    LOGGER(DEBUG, "coul_mat ends");
    return e;
}

} // end namespace cls
