/*
 * File: GridDen.h
 * Author: Qiming Sun <osirpt.sun@gmail.com>
 */

#include <vector>
#include "basis.h"
extern "C" {
#include <xc.h>
}

using namespace std;

class GridDen {
public:
    GridDen();
    ~GridDen();

    int set_den_cutoff(const double);
    int set_background_charge(const double);
    double dentsity_at_grid(const int grid_id)
    int density_matrix(const int)
private:
    double _bg_chg;
    double _den_cutoff;
};

