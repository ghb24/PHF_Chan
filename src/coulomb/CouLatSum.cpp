/*
 * File: CouLatSum.cpp
 * Author: Qiming Sun <osirpt.sun@gmail.com>
 */

#include <iostream>
#include <vector>
#include <sstream>
#include <assert.h>
#include "../lib/CtIo.h"
#include "CouLatSum.h"
#include "nfft.h"
#include "../PhfBasisSet.h"
//#include "grid_3d_class.h"
//#include "../lib/CxAlgebra.h"
#include "fblas.h"

#define DEBUG

#if defined DEBUG
#include <math.h>
//#define LOGGER(X,Y)    xout << "DEBUG:" << X << std::endl
#define LOGGER(X,Y)     std::cout << #X << ": " << Y << std::endl
#else
#define LOGGER(X,Y)
#endif


namespace cls { // Coulomb Lattice Summation

struct maybe_env_t {
    const FSolidModel *p_solid;
    const FSuperCell *p_super_cell;
    const FLattice *p_lattice;
    const FUnitCell *p_unit_cell;
    const FOpMatrix *p_den_mat;
    phf::coulomb_grid::grid_3d *p_grid3d;
     
    int nao_unit_cell;
    int n_unit_cell_in_super_cell;
    double den_cutoff;
    double background_chg;
    double grid_weight;
    double *grids_coord_in_unit_cell;
    //TODO:cells_coord_in_super_cell;

    double* pOrbVal;
    FORTINT nComp;
    FORTINT *pCenterIndices;
    FORTINT *pMap;
    FORTINT nMap;
    const FBasisSet *p_bas;
    double (*GridPt)[3];
    FORTINT nGridPt;
    FORTINT nDiffBf;
    FORTINT iContext;
    FORTINT ic;
};

static struct maybe_env_t *_env;
static double _default_num_grids_in_unit_cell[3] = {20, 20, 20};

/*
 **************************************************
 * local funcs
 ***************** accumulator ********************
 * data flow:
 * -> cells_accumulator -> grids_accumulator ->
 */
static void cells_accumulator(FCellIter_t func, FCellTruncate_t ftrunc,
                              const int grid_id,
                              double *v, double *mat, double *ao)
{
    int cell_id[3];
    for (int ix = 0; ix < _env->p_super_cell->Size[0]; ix++)
        for (int iy = 0; iy < _env->p_super_cell->Size[1]; iy++)
            for (int iz = 0; iz < _env->p_super_cell->Size[2]; iz++) {
                cell_id[0] = ix;
                cell_id[1] = iy;
                cell_id[2] = iz;
                if (!(*ftrunc)(cell_id)) {
                    (*func)(grid_id, cell_id, v, mat, ao);
                }
            }
}

static void grids_accumulator(FGridIter_t func,
                              double *v, double *mat, double *ao)
{
    int grid_id = 0;
    for (int ix = 0; ix < _env->p_grid3d->num_x; ix++)
        for (int iy = 0; iy < _env->p_grid3d->num_y; iy++)
            for (int iz = 0; iz < _env->p_grid3d->num_z; iz++) {
                (*func)(grid_id, v, mat, ao);
                grid_id++;
                v++;
            }
}
/***************** accumulator end ****************/

/*
 ********** FCellTruncate_t functions *************
 * return 1 to truncate, drop the cell associated with cell_id
 * return 0 to truncate, drop the cell associated with cell_id
 * see the function cells_accumulator
 */
// TODO:
//static int is_cell_far(int cell_id[3])
//{
//    if (SUPER_CELL_INFO or cutoff_info) {
//        return 1; // remote cell
//    } else {
//        return 0;
//    }
//}

static int cover_whole_super_cell(int cell_id[3])
{
    return 0;
}
/********** FCellTruncate_t functions end *********/

/************************************************/
static void dset0(const int n, double *a)
{
    for (int i = 0; i < n; i++) {
        a[i] = 0;
    }
}
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

//ABORTstatic int cell_offset_by_cell_id(const int cell_id[3])
//ABORT{
//ABORT    //int nx = _env->p_super_cell->Size[0];
//ABORT    int ny = _env->p_super_cell->Size[1];
//ABORT    int nz = _env->p_super_cell->Size[2];
//ABORT    return cell_id[2] * ny * nz + cell_id[1] * nz + cell_id[0];
//ABORT}
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
static void grid_coord_by_grid_id(const int grid_id, double grid_coord[3])
{
    grid_coord[0] = _env->grids_coord_in_unit_cell[grid_id * 3 + 0];
    grid_coord[1] = _env->grids_coord_in_unit_cell[grid_id * 3 + 1];
    grid_coord[2] = _env->grids_coord_in_unit_cell[grid_id * 3 + 2];
}

// for the grid in 0-th cell, value of all AOs in the super cell
// TODO: save the output of FD(eval_basis_fn_on_grid) once, and load it here
static void ao_super_cell(const double grid_coord[3], double *ao)
{
    int cell_id[3];
    double *const coords = new double[_env->n_unit_cell_in_super_cell * 3];
    double *ao_not0 = new double[_env->n_unit_cell_in_super_cell * _env->nao_unit_cell];
    FORTINT *p_map = new FORTINT[_env->nao_unit_cell];
    FORTINT *pCenterIndices = new FORTINT[_env->nao_unit_cell];
    FORTINT nmap;
    double *pcoord = coords;
    for (int ix = 0; ix < _env->p_super_cell->Size[0]; ix++)
        for (int iy = 0; iy < _env->p_super_cell->Size[1]; iy++)
            for (int iz = 0; iz < _env->p_super_cell->Size[2]; iz++) {
                cell_id[0] = ix;
                cell_id[1] = iy;
                cell_id[2] = iz;
                cell_coord_by_cell_id(cell_id, pcoord);
                pcoord[0] = grid_coord[0] - pcoord[0];
                pcoord[1] = grid_coord[1] - pcoord[1];
                pcoord[2] = grid_coord[2] - pcoord[2];
                pcoord += 3;
            }
    dset0(_env->n_unit_cell_in_super_cell * _env->nao_unit_cell, ao);
    FORTINT deriv_order = 0;
    FORTINT stride_for_deriv_order = 1;
    FD(eval_basis_fn_on_grid)(ao_not0, stride_for_deriv_order,
                              pCenterIndices, p_map, nmap, *(_env->p_bas),
                              (double (*)[3])coords,
                              (FORTINT)(_env->n_unit_cell_in_super_cell), 
                              deriv_order, _env->ic);
    // ao_not0 is saved as Fortran_array(grid,bas)
    double *pao_not0;
    for (int i = 0; i < _env->n_unit_cell_in_super_cell; i++) {
        pao_not0 = ao_not0 + i;
        for (int imap = 0; imap < nmap; imap++) {
            ao[p_map[imap]-1] = *pao_not0; // p_map is in fortran index
            pao_not0 += _env->n_unit_cell_in_super_cell; // next ao
        }
        ao += _env->nao_unit_cell;
    }
    delete[] coords;
    delete[] ao_not0;
    delete[] p_map;
    delete[] pCenterIndices;
}


// for a given grid_id in 0-th cell,
// iterate this fn in cells_accumulator to run over all the image grids
static void density_cell0(const int grid_id, const int cell_id[3],
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
static void density_grid_iter(const int grid_id,
                              double *v, double *dmat, double *ao)
{
    *v = 0;
    cells_accumulator(density_cell0, cover_whole_super_cell,
                      grid_id, v, (double *)&_env->p_den_mat->Data[0], ao);
    *v += -_env->background_chg;
}

// grid_id of one cell
double density_at_grid(const int grid_id)
{
    std::vector<double> ao(_env->p_den_mat->nCols);
    double den;
    density_grid_iter(grid_id,
                      &den,
                      const_cast<double *>(&(*_env->p_den_mat)[0]), &ao[0]);
    return den;
}

void density_unit_cell(double *density)
{
    std::vector<double> ao(_env->p_den_mat->nCols);
    grids_accumulator(density_grid_iter,
                      density,
                      const_cast<double *>(&(*_env->p_den_mat)[0]), &ao[0]);
}
//double density_unit_cell(double *density)
//{
//    LOGGER(DEBUG, "density on each grid starts");
//    double nele = 0;
//    int grid_id = 0;
//    for (int ix = 0; ix < _env->p_grid3d->num_x; ix++)
//        for (int iy = 0; iy < _env->p_grid3d->num_y; iy++)
//            for (int iz = 0; iz < _env->p_grid3d->num_z; iz++) {
//                density[grid_id] = density_at_grid(grid_id);
//                nele += _env->grid_weight * density[grid_id];
//                grid_id++;
//            }
//    return nele;
//}

// for a given grid_id in 0-th cell,
// iterate this fn in cells_accumulator to run over all the image grids
static void coul_matrix_cell_iter(const int grid_id, const int cell_id[3],
                                  double *v, double *dmat, double *ao)
{
    double cell_coord[3];
    double coord[3];
    int cell_id_shift[3];
    // FIXME: shift the 0th-cell to the center of the super cell?
    // division should be rounded down
    cell_id_shift[0] = cell_id[0] - _env->p_super_cell->Size[0]/2;
    cell_id_shift[1] = cell_id[1] - _env->p_super_cell->Size[1]/2;
    cell_id_shift[2] = cell_id[2] - _env->p_super_cell->Size[2]/2;
    grid_coord_by_grid_id(grid_id, coord);
    cell_coord_by_cell_id(cell_id_shift, cell_coord);
    coord[0] += cell_coord[0];
    coord[1] += cell_coord[1];
    coord[2] += cell_coord[2];
    ao_super_cell(coord, ao);

    int n_mu = _env->p_den_mat->nRows;
    int n_nu = _env->p_den_mat->nCols;
    const int INC1 = 1;
    double vw = *v * _env->grid_weight;
    dger_(&n_mu, &n_nu, &vw, ao, &INC1, ao, &INC1, dmat, &n_mu);
}

// iterate this fn in grids_accumulator to run over all grids in the unit cell
static void coul_matrix_grid_iter(const int grid_id,
                                  double *v, double *dmat, double *ao)
{
    cells_accumulator(coul_matrix_cell_iter, cover_whole_super_cell,
                      grid_id, v, dmat, ao);
}

static void coul_matrix_acc(double *real_space_pot, double *vmat)
{
    LOGGER(DEBUG, "coul_matrix_acc starts");
    std::vector<double> ao(_env->p_den_mat->nCols);
    dset0((int)_env->p_den_mat->size(), vmat);
    // accumulate all the grids in 0-th cell
    grids_accumulator(coul_matrix_grid_iter, real_space_pot, vmat, &ao[0]);
}

// init after construct obj grid3d
void init_grids_coord_in_unit_cell()
{
    _env->grids_coord_in_unit_cell = new double[_env->p_grid3d->num_grid * 3];
    double *p = _env->grids_coord_in_unit_cell;
    for (int ix = 0; ix < _env->p_grid3d->num_x; ix++)
        for (int iy = 0; iy < _env->p_grid3d->num_y; iy++)
            for (int iz = 0; iz < _env->p_grid3d->num_z; iz++) {
                p[0] = _env->p_lattice->T[0][0] * ix / _env->p_grid3d->num_x
                     + _env->p_lattice->T[1][0] * iy / _env->p_grid3d->num_y
                     + _env->p_lattice->T[2][0] * iz / _env->p_grid3d->num_z;
                p[1] = _env->p_lattice->T[0][1] * ix / _env->p_grid3d->num_x
                     + _env->p_lattice->T[1][1] * iy / _env->p_grid3d->num_y
                     + _env->p_lattice->T[2][1] * iz / _env->p_grid3d->num_z;
                p[2] = _env->p_lattice->T[0][2] * ix / _env->p_grid3d->num_x
                     + _env->p_lattice->T[1][2] * iy / _env->p_grid3d->num_y
                     + _env->p_lattice->T[2][2] * iz / _env->p_grid3d->num_z;
                p += 3;
            }
}

void init_env(const FSolidModel& solid, const FOpMatrix& den_mat)
{
    LOGGER(DEBUG, "Init cls_env");
    _env = new maybe_env_t;

    _env->p_solid = &solid;
    _env->p_super_cell = &solid.SuperCell;

    _env->p_lattice = &solid.Lattice;
    LOGGER(DEBUG, "lattice vector 1 = ("
           << solid.Lattice.T[0][0] << ","
           << solid.Lattice.T[0][1] << ","
           << solid.Lattice.T[0][2] << ")");
    LOGGER(DEBUG, "lattice vector 2 = ("
           << solid.Lattice.T[1][0] << ","
           << solid.Lattice.T[1][1] << ","
           << solid.Lattice.T[1][2] << ")");
    LOGGER(DEBUG, "lattice vector 3 = ("
           << solid.Lattice.T[2][0] << ","
           << solid.Lattice.T[2][1] << ","
           << solid.Lattice.T[2][2] << ")");

    _env->p_unit_cell = &solid.UnitCell;
    LOGGER(DEBUG, "atoms in unit cell");
    for (unsigned int i = 0; i <_env->p_unit_cell->Elements.size(); i++) {
        LOGGER(DEBUG, "atom " << i << " nuc = "
               << _env->p_unit_cell->Elements[i] << " @ ("
               << _env->p_unit_cell->Coords[i][0] << ","
               << _env->p_unit_cell->Coords[i][1] << ","
               << _env->p_unit_cell->Coords[i][2] << ")");
    }

    _env->p_bas = &solid.UnitCell.OrbBasis;
    for (unsigned int i = 0; i < _env->p_bas->Centers.size(); i++) {
        LOGGER(DEBUG, "bas " << i
               << ", nFn = " << _env->p_bas->Groups[i].nFn
               << ", nExp = " << _env->p_bas->Groups[i].nExp
               << ", nCo =  " << _env->p_bas->Groups[i].nCo
               << ", l =  " << _env->p_bas->Groups[i].l
               << ", Type = " << _env->p_bas->Groups[i].Type
               << ", iCen = " << _env->p_bas->Groups[i].iCen
               << ", Range = " << _env->p_bas->Data[_env->p_bas->Groups[i].iCen]);
        LOGGER(DEBUG, "     "
               << " [iExp..] = " << _env->p_bas->Data[_env->p_bas->Groups[i].iExp]
               << ",[iCo..] = " << _env->p_bas->Data[_env->p_bas->Groups[i].iCo]);
    }

    _env->p_den_mat = &den_mat;
    LOGGER(DEBUG, "dim of den_mat = (" << den_mat.nRows
           << "," << den_mat.nCols << ")");

    //_env->nao_unit_cell = solid.UnitCell.OrbBasis.nFn;
    _env->nao_unit_cell = den_mat.nRows;
    LOGGER(DEBUG, "nao_unit_cell = " << _env->nao_unit_cell);

    _env->n_unit_cell_in_super_cell = _env->p_super_cell->Size[0]
                                    * _env->p_super_cell->Size[1]
                                    * _env->p_super_cell->Size[2];
    LOGGER(DEBUG, "num unit cell (x,y,z,total) = ("
           << _env->p_super_cell->Size[0] << ","
           << _env->p_super_cell->Size[1] << ","
           << _env->p_super_cell->Size[2] << ","
           << _env->n_unit_cell_in_super_cell << ")");

    double nele = 0;
    for (TArray<FORTINT>::const_iterator n = solid.UnitCell.Elements.begin();
         n != solid.UnitCell.Elements.end(); n++) {
        nele += *n;
    }
    LOGGER(DEBUG, "nele in unit cell = " << nele);
    double vol = solid.UnitCell.Volume;
    LOGGER(DEBUG, "unit cell vol = " << vol);
    _env->background_chg = nele / vol;
    LOGGER(DEBUG, "set background charge to unit cell average charge = "
           << _env->background_chg);

    // set up mesh grids
    _env->p_grid3d =
        new phf::coulomb_grid::grid_3d(const_cast<FSolidModel&>(solid),
                                       _default_num_grids_in_unit_cell[0],
                                       _default_num_grids_in_unit_cell[1],
                                       _default_num_grids_in_unit_cell[2]);
#if defined DEBUG
    _env->p_grid3d->set_ngrid(1, 1, 2);
#endif
    LOGGER(DEBUG, "unit cell grid numbers (x,y,z,total) = ("
           << _env->p_grid3d->num_x << ","
           << _env->p_grid3d->num_y << ","
           << _env->p_grid3d->num_z << ","
           << _env->p_grid3d->num_grid << ")");
    init_grids_coord_in_unit_cell();

    _env->grid_weight = vol / _env->p_grid3d->num_grid;
    LOGGER(DEBUG, "unit cell grid weight = " << _env->grid_weight);

    LOGGER(DEBUG, "create_integral_context");
    _env->ic = FD(create_integral_context)(0,0, 1e-10);

    LOGGER(DEBUG, "init_env ends");
}

void del_env()
{
    LOGGER(DEBUG, "destroy_integral_context");
    FD(destroy_integral_context)(_env->ic);

    delete[] _env->grids_coord_in_unit_cell;
    delete _env->p_grid3d;
    delete _env;
}

/*
 **************************************************
 * for Poisson equation FFT solver
 * size of array kcoords >= size_of [_env->p_grid3d->num_grid][3]
 */
static void kspace_grids_generator(double *kcoords)
{
    for (int ix = 0; ix < _env->p_grid3d->num_x; ix++)
        for (int iy = 0; iy < _env->p_grid3d->num_y; iy++)
            for (int iz = 0; iz < _env->p_grid3d->num_z; iz++) {
                kcoords[0] = _env->p_lattice->K[0][0] * ix / _env->p_grid3d->num_x
                           + _env->p_lattice->K[1][0] * iy / _env->p_grid3d->num_y
                           + _env->p_lattice->K[2][0] * iz / _env->p_grid3d->num_z;
                kcoords[1] = _env->p_lattice->K[0][1] * ix / _env->p_grid3d->num_x
                           + _env->p_lattice->K[1][1] * iy / _env->p_grid3d->num_y
                           + _env->p_lattice->K[2][1] * iz / _env->p_grid3d->num_z;
                kcoords[2] = _env->p_lattice->K[0][2] * ix / _env->p_grid3d->num_x
                           + _env->p_lattice->K[1][2] * iy / _env->p_grid3d->num_y
                           + _env->p_lattice->K[2][2] * iz / _env->p_grid3d->num_z;
                kcoords += 3;
            }
}

/*
 **************************************************
 * global functions
 */

// TODO: density is obtained by looping over all the unit cells in the
// super cell, cutoff may help to reduce the number of unit cells to be
// looped
// void set_den_cutoff(const double den_cutoff)
// {
//     cls::_env->den_cutoff = den_cutoff;
// }

void set_num_grids_in_unit_cell(const int n_a1, const int n_a2, const int n_a3)
{
    _default_num_grids_in_unit_cell[0] = n_a1;
    _default_num_grids_in_unit_cell[1] = n_a2;
    _default_num_grids_in_unit_cell[2] = n_a3;
}

// return the coulomb energy
double coul_matrix(const FSolidModel& solid, const FOpMatrix& den_mat,
                   FOpMatrix& coul_mat)
{
    LOGGER(DEBUG, "coul_matrix starts");
#if defined DEBUG
    FOpMatrix dm_tmp = den_mat;
    for (unsigned int i =0; i<dm_tmp.size(); i++) {
        //dm_tmp[i] = cos(i);
        dm_tmp[i] = i*.1;
    }
    init_env(solid, dm_tmp);
#else
    init_env(solid, den_mat);
#endif

    std::vector<double> real_space_density(_env->p_grid3d->num_grid);
    std::vector<double> real_space_pot(_env->p_grid3d->num_grid);

    //TODO:int ngrid3d = (_env->p_grid3d->num_x) * _env->n_unit_cell_in_super_cell;
    //TODO:_env->nComp = 1;
    //TODO:FORTINT nbasis = _env->p_bas->nFn;
    //TODO:_env->pCenterIndices = new FORTINT [nbasis];
    //TODO:_env->pMap = new FORTINT [nbasis];
    //TODO:_env->nGridPt = _env->n_unit_cell_in_super_cell * _env->p_grid3d->num_x;
    //TODO:_env->pOrbVal = new double [ _env->nGridPt * _env->nComp * nbasis ];
    //TODO://for( int i = 0; i < 3; i++ ){
    //TODO: _env->GridPt = (double (*)[3])(new double [ ngrid3d*3 ]);
    //TODO://}
    //TODO:_env->nDiffBf = 0;
    //TODO:// calculate the basis function values in the super cell grid
    //TODO:FD(eval_basis_fn_on_grid)( _env->pOrbVal, _env->nComp, _env->pCenterIndices, _env->pMap, _env->nMap,
    //TODO:  *(_env->p_bas), _env->GridPt, _env->nGridPt, _env->nDiffBf, _env->ic );

    // calculate the density
    density_unit_cell(&real_space_density[0]);
    double nele = 0;
    for (unsigned int i = 0; i < real_space_density.size(); i++)
        nele += real_space_density[i];
    nele *= _env->grid_weight;
    LOGGER(INFO, "sum(weight * density) => number of nele = " << nele);
#if defined DEBUG
    LOGGER(DEBUG, "density_unit_cell ends");
#else
    if (std::abs(nele) > 1e-5) {
        std::stringstream msg;
        msg << "Incorrect number of electron! " << nele;
        std::cerr << msg.str() << std::endl;
        throw msg.str();
    }
#endif
    // FFT Poisson solver
    LOGGER(DEBUG, "call poisson solver");
    std::vector<double> kspace_grids(_env->p_grid3d->num_grid * 3);
    kspace_grids_generator(&kspace_grids[0]);
    poisson_kspace_solver(&real_space_density[0], &real_space_pot[0],
                          _env->p_grid3d, &kspace_grids[0]);
    // Coulomb matrix
    coul_matrix_acc(&real_space_pot[0], &coul_mat[0]);

    del_env();

    // Coulomb energy
    const int INC1 = 1;
    int n = coul_mat.size();
#if defined DEBUG
    double e = ddot_(&n, &coul_mat[0], &INC1, &dm_tmp[0], &INC1)
            / _env->n_unit_cell_in_super_cell;
#else
    double e = ddot_(&n, &coul_mat[0], &INC1, &den_mat[0], &INC1)
            / _env->n_unit_cell_in_super_cell;
#endif
    LOGGER(DEBUG, "coul_mat ends");
    return e;
}

} // end namespace cls
