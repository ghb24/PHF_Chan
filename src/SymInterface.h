/*
 * Oct. 13, 2012, Written by Naoki Nakatani
 * ----------------------------------------------------------------------------------------------------
 * Skeleton interfaces to symmetry operation for unit-cell, super-cell, and orbitals
 * They can be called from both Fortran and C/C++ but will be changed to fit global data structure
 */
#ifndef SYM_ITERFACE_H
#define SYM_ITERFACE_H

/*
 * Initializing Lattice
 * + latVec : real-space lattice vector (Real(3), input)
 * + nRept  : repetition vector (Integer(3), input)
 * + spcVec : real-space super-cell vector (Real(3), output)
 * + recVec : reciprocal lattice vector (Real(3), output)
 * + volume : volume of super-cell (Real, output)
 */
void InitLattice(double* latVec, int* nRept, double* spcVec, double* recVec, double* volume);

/*
 * Compute Symmetry-Unique Atoms
 * + x      : cart.-coord. x[iAtom] (Real(nAtoms), input)
 * + y      : cart.-coord. y[iAtom] (Real(nAtoms), input)
 * + z      : cart.-coord. z[iAtom] (Real(nAtoms), input)
 * + nAtoms : # of atoms in unit-cell
 * + iSymAtom
 *          : list of symm-unique iAtom (Integer(nSymAtoms), output)
 * + nSymAtoms
 *          : # of symm-unique atoms (Integer, output)
 */
void GetSymUniqueAtoms(double* x, double* y, double* z, int* nAtoms, int* iSymAtom, int* nSymAtoms);

/*
 * Compute Translational Symmetry
 * + x      : cart.-coord. x[iAtom] (Real(nAtoms), input)
 * + y      : cart.-coord. y[iAtom] (Real(nAtoms), input)
 * + z      : cart.-coord. z[iAtom] (Real(nAtoms), input)
 * + nAtoms : # of atoms in unit-cell
 * + nx     : list of translation nx[iTrAtom] (Real(nTrAtoms), output)
 * + ny     : list of translation ny[iTrAtom] (Real(nTrAtoms), output)
 * + nz     : list of translation nz[iTrAtom] (Real(nTrAtoms), output)
 * + weight : weight of iTrAtom (Real(nTrAtoms), output)
 * + nTrAtoms
 *          : # of atoms in super-cell
 */
void GetTransSymAtoms(double* x, double* y, double* z, int* nAtoms,
                      double* nx, double* ny, double* nz, double* weight, int* nTrAtoms);

/*
 * Transformation between AO and SO
 * + pAO    : data matrix spanned by AO (Real(nAO*nAO), in/out)
 * + nAO    : # of AO in unit-cell (Integer, in/out)
 * + pSO    : data matrix spanned by SO (Real(nSO*nSO), in/out)
 * + nSO    : # of SO in unit-cell (Integer, in/out)
 * + iDir   : direction of transformation (0: AO->SO, 1:SO->AO)
 * (more input parameters might be required)
 */
void SymTrans(double* pAO, int* nAO, double* pSO, int* nSO, int* iDir);

#endif
