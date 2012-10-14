#ifndef PHF_OPERATOR_H
#define PHF_OPERATOR_H

#include "PhfTypes.h"

/// represents the matrix elements of a translationally symmetric 1e operator,
/// in either a raw basis, or a symmetry adapted basis (same storage format).
///
/// Notes:
///   - This operator effectively represents a matrix of format
///         nAoSuperCell x nAoSuperCell
///     where nAoSuperCell is the number of basis functions in the entire
///     super-cell (nAo x Lattice.Ts.size()).
///
///   - Due to symmetry, we can bind the first index to the first unit cell
///     without loss of information: If T is a lattice translation (or other
///     spatial symmetry group operation), then
///         <\mu|op|\nu> = < T\mu|op|T\nu>.
///     Thus all matrix elements can be recovered from <\mu|op|(all T)\nu>
///     with \mu,\nu in the first unit-cell. (For this we need the action of a
///     symmetry element T on the basis functions.)
///
///
/// For a start, the following representation is used:
///   - The type is stored as 3-dimensional array, with format
///        nAo x nAo x SuperCell.Ts.size()
///     where nAo is the number of basis functions in the unit-cell and
///     the last index iterates through translations in the super-cell
///     between the first unit-cell and the second unit cell.
///
///   - That is: The first index refers to UnitCell.OrbBasis, while the second
///     refers to SuperCell.OrbBasis.
///
///   - Later we would like to restrict the first index to symmetry-unique
///     basis functions (full space group) additionally; and to support the
///     transformation to the fully symmetry adapted basis (in which the operator
///     is block-diagonal). For this reason the storage format might change.
template<class FScalar>
struct TOpMatrix
{



};
typedef TOpMatrix<double>
   FOpMatrix;
typedef TOpMatrix<complex_double>
   FOpMatrixC;


#endif // PHF_OPERATOR_H
