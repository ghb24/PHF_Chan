/*
 * File: GridDen.h
 * Author: Qiming Sun <osirpt.sun@gmail.com>
 */

#include <vector>
#include "basis.h"
extern "C" {
#include <xc.h>
}

#define NONE 0

typedef void (*FuncIter_t)(const int *grid_id, const int *cell_id,
                           double *v, double *mat, double *mu, double *nu);
typedef void (*FuncTerm_t)(const int *);



int set_den_cutoff(const double);

int set_background_charge(const double);

double dentsity_at_grid(const int *grid_id);
//TODO:int density_matrix(const int)
