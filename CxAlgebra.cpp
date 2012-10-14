/* Copyright (c) 2012  Gerald Knizia
 * 
 * This file is part of the bfint program
 * (See http://www.theochem.uni-stuttgart.de/~knizia)
 * 
 * bfint is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3.
 * 
 * bfint is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with bfint (LICENSE). If not, see http://www.gnu.org/licenses/
 */

#include <stdexcept>
#include <stdlib.h>

#include "CxAlgebra.h"
#include "CxDefs.h" // for assert.

namespace ct {

// Out = f * A * B
void Mxm(double *pOut, ptrdiff_t iRowStO, ptrdiff_t iColStO,
         double const *pA, ptrdiff_t iRowStA, ptrdiff_t iColStA,
         double const *pB, ptrdiff_t iRowStB, ptrdiff_t iColStB,
         size_t nRows, size_t nLink, size_t nCols, bool AddToDest, double fFactor )
{
    assert( iRowStO == 1 || iColStO == 1 );
    assert( iRowStA == 1 || iColStA == 1 );
    assert( iRowStB == 1 || iColStB == 1 );
    // ^- otherwise dgemm directly not applicable. Would need local copy
    // of matrix/matrices with compressed strides.

//     if ( nRows == 1 || nLink == 1 || nCols == 1 ) {
//         if ( !AddToDest )
//             for ( uint ic = 0; ic < nCols; ++ ic )
//                 for ( uint ir = 0; ir < nRows; ++ ir )
//                     pOut[ir*iRowStO + ic*iColStO] = 0;
//
//         for ( uint ic = 0; ic < nCols; ++ ic )
//             for ( uint ir = 0; ir < nRows; ++ ir )
//                 for ( uint il = 0; il < nLink; ++ il )
//                     pOut[ir*iRowStO + ic*iColStO] += fFactor * pA[ir*iRowStA + il*iColStA] * pB[il*iRowStB + ic*iColStB];
//         return;
//     }

    double
        Beta = AddToDest? 1.0 : 0.0;
    char
        TransA, TransB;
    FORTINT
        lda, ldb,
        ldc = (iRowStO == 1)? iColStO : iRowStO;

    if ( iRowStA == 1 ) {
        TransA = 'N'; lda = iColStA;
    } else {
        TransA = 'T'; lda = iRowStA;
    }
    if ( iRowStB == 1 ) {
        TransB = 'N'; ldb = iColStB;
    } else {
        TransB = 'T'; ldb = iRowStB;
    }

    DGEMM( TransA, TransB, nRows, nCols, nLink,
        fFactor, pA, lda, pB, ldb, Beta, pOut, ldc );
}

// note: both H and S are overwritten. Eigenvectors go into H.
void DiagonalizeGen(double *pEw, double *pH, uint ldH, double *pS, uint ldS, uint N)
{
    FORTINT info = 0, nWork = 128*N;
    double *pWork = (double*)::malloc(sizeof(double)*nWork);
    DSYGV(1, 'V', 'L', N, pH, ldH, pS, ldS, pEw, pWork, nWork, info );
    ::free(pWork);
    if ( info != 0 ) throw std::runtime_error("dsygv failed.");
};

void Diagonalize(double *pEw, double *pH, uint ldH, uint N)
{
    FORTINT info = 0, nWork = 128*N;
    double *pWork = (double*)::malloc(sizeof(double)*nWork);
    DSYEV('V', 'L', N, pH, ldH, pEw, pWork, nWork, info );
    ::free(pWork);
    if ( info != 0 ) throw std::runtime_error("dsyev failed.");
};

}

// kate: indent-width 4
