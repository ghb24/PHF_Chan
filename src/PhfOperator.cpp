#include <ostream>
#include "lib/CtMatrix.h" // for PrintMatrixGen
#include "PhfOperator.h"

template<class FScalar>
TOpMatrix<FScalar>::TOpMatrix(FSolidModel const &M)
{
   nRows = M.UnitCell.OrbBasis.nFn;
   nCols = M.SuperCell.OrbBasis.nFn;

   Data.resize(nRows * nCols);
};

template<class FScalar>
void TOpMatrix<FScalar>::Print(std::ostream &xout, std::string const &Caption)
{
   ct::PrintMatrixGen(xout, &Data[0], nRows, 1, nCols, nRows, Caption);
};


// explicitly instantiate the templates for the types we will use.
// This binds their implementation to this file.
template class TOpMatrix<double>;
// template class TOpMatrix<complex_double>;
