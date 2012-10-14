/*
 * File: GridDen.cc
 * Author: Qiming Sun <osirpt.sun@gmail.com>
 */

#include <iostream>
#include <assert.h>
#include "GridDen.h"

//#include "grid.h" // wf's class

extern "C" {
#include <cblas.h>
#include <stdio.h>
#include <stdlib.h>
}

namespace fftcoul {
/*
 * local funcs
 */

double
accumlator(double (NumInt::*fgrid)(const int, xc_func_type&, xc_func_type&,
                                   double&, double *, const double *, const int) const,
                    xc_func_type& func_x, xc_func_type& func_c,
                    double& e, double *v, const double *dm) const
{
    int id, nn;
    double nelec = 0;
    double nelec_priv, e_priv;
    double *v_priv;
    if (_relativity == XC_RELATIVISTIC) {
        nn = _num_cgto * _num_cgto * 4; // v, spvsp
    } else {
        nn = _num_cgto * _num_cgto;
    }
    e = 0;
    dset0(nn, v);
//#if defined OMP
//#pragma omp parallel default(none) \
//    shared(nn, nelec, fgrid, func_x, func_c, e, v, dm) \
//    private(id, nelec_priv, e_priv, v_priv)
//#endif
    {
        nelec_priv = 0;
        e_priv = 0;
        v_priv = new double[nn];
        dset0(nn, v_priv);
//#if defined OMP
//#pragma omp for nowait schedule(dynamic, 20)
//#endif
        for (id = 0; id < _tot_num_grid; id++) {
            nelec_priv += _weights[id]
                * (this->*fgrid)(id, func_x, func_c, e_priv, v_priv, dm, nn);
        }
//#if defined OMP
//#pragma omp critical
//#endif
        {
            nelec += nelec_priv;
            e += e_priv;
            cblas_daxpy(nn, 1, v_priv, 1, v, 1);
        }
        delete[] v_priv;
    }
    return nelec;
}

static int cells_looper()
{
}

static int grids_looper()
{
    for (int cell_id = 0; cell_id < ..; cell_id++) {
    }
}

static int density_at_grid(const int tot_cell, const int tot_grid)
{
    cells_looper(grids_looper);
}


/*
 * class defs
 */


GridDen::GridDen()
{
    _den_cutoff = 1e-10;
    _bg_chg = 0;
}

GridDen::~GridDen()
{}

int
GridDen::set_den_cutoff(const double den_cutoff)
{
    _den_cutoff = den_cutoff;
    return 0
}

int
GridDen::set_background_charge(const double chg)
{
    _bg_chg = chg;
    return 0
}

double
GridDen::density_at_grid(const int grid_id)
{
}

int
GridDen::density_matrix(const int flags, CELLINFO, GRIDINFO, double *mat)
{
    return fftcoul::density_matrix()
}




} // end namespace fftcoul
