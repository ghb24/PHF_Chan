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

#ifndef CT8k_MATRIX_H
#define CT8k_MATRIX_H

#include "AicCommon.h" // FStackMatrix sits on FMemoryStack
#include "CxAlgebra.h"
#include "CtCommon.h"

#include <string>
#include <sstream>

namespace ct {

typedef double
    FScalar;

using std::size_t;

/// Describes the layout of a matrix somewhere in memory, which this object
/// DOES NOT own.
struct FMatrixView
{
    FScalar
        *pData; ///< start of data of this matrix.
    uint
        nRows, ///< number of rows (first index, fast)
        nCols, ///< number of columns (second index, slow)
        nRowSt, ///< number of words between adjacent rows. In default-alignment, nRowSt == 1
        nColSt;  ///< number of words between adjacent cols. In default-alignment, nColSt == nRows
    // Note: One of both strides should be 1.

    /// element retrieval operator: m(nRow,nCol) -> element stored at that
    /// part of the matrix. Note: indices are 0-based!.
    inline FScalar& operator () ( uint iRow, uint iCol ){
        assert_fast( iRow <= nRows );
        assert_fast( iCol <= nCols );
        return pData[nRowSt * iRow + nColSt * iCol];
    };
    inline FScalar const& operator () ( uint iRow, uint iCol ) const {
        return const_cast<FMatrixView *const>(this)->operator()( iRow, iCol );
    };

    FScalar &operator [] (uint iEntry) { return pData[iEntry]; }
    FScalar const &operator []  (uint iEntry) const { return pData[iEntry]; }

    FMatrixView()
        : pData(0)
    {}

    FMatrixView( double *pData_, uint nRows_, uint nCols_, uint nRowSt_ = 1, uint nColSt_ = 0 )
        : pData(pData_), nRows(nRows_), nCols(nCols_),
          nRowSt(nRowSt_), nColSt( (nColSt_ != 0)? nColSt_ : nRows )
    {};

    inline uint GetStridedSize() const;

    void Print( std::ostream &out, std::string const &Name = "" ) const;
    /// sets matrix to zero
    void Clear();
    /// sets matrix to identity
    void SetIdentity();
    /// returns copy of diagonal elements of the matrix. Only valid for square matrices.
    std::vector<FScalar> GetDiagonal() const;

    bool IsSquare() const { return nRows == nCols; };
    bool IsSymmetric(FScalar Thresh) const;

    /// if nRowSt==1, sets space between nRow and nColSt to zero (in order to
    /// not have possible NaNs/Denorms in spaces in the matrices which are only
    /// there for alignment purposes; e.g., ntb vs nt strides.
    void ClearEmptySpace();

    /// returns size of this matrix would require in triangular storage (e.g.,
    /// storing only the lower triangular part). If Sign == -1, size is
    /// calculated for anti-symmetric matrix.
    size_t TriangularStorageSize(int Sign = 1) const;
    size_t TriangularStorageSize1() const { return TriangularStorageSize(1); };
    void TriangularExpand(int Sign = 1, FScalar *pTriangData = 0);
    void TriangularReduce(int Sign = 1, FScalar *pTriangData = 0);
    enum FTriangularTransformFlags {
        TRIANG_Expand,
        TRIANG_Reduce
    };
    void TriangularExpandOrReduce(int Sign, FScalar *pTriangData, FTriangularTransformFlags Direction);
};

uint FMatrixView::GetStridedSize() const {
   if ( nRowSt == 1 )
       return nCols * nColSt;
   else {
       if ( nColSt == 1 )
           return nRows * nRowSt;
       else
           assert(0); // not supported.
           return 0;
   }
};



// export version (can also be used in other files than CtMatrix.cpp)
void AssertCompatible1( FMatrixView const &A, FMatrixView const &B );

uint const MXM_AddToDest = 0x1u;

/// Creates a view representing the transpose of the input matrix.
/// The returned view still references the original data (i.e., In and Out are aliased)
FMatrixView Transpose( FMatrixView const &In );
/// Out += f * In
void Add( FMatrixView &Out, FMatrixView const &In, FScalar fFactor = 1.0 );
/// Out = f * In
void Move( FMatrixView &Out, FMatrixView const &In, FScalar fFactor = 1.0 );
/// Out = A * B
void Mxm( FMatrixView &Out, FMatrixView const &A, FMatrixView const &B, uint Flags = 0 );
/// Out = f * A * B
void Mxm( FMatrixView &Out, FMatrixView const &A, FMatrixView const &B, FScalar fFactor, uint Flags = 0 );
/// Out = f * A * A^T for symmetric matrix Out [Note: Out will be symmetrized!].
void Syrk( FMatrixView &Out, FMatrixView const &A, FScalar fFactor = 1.0, uint Flags = 0 );
/// Out = f * A^T * A for symmetric matrix Out [Note: Out will be symmetrized!].
void SyrkT( FMatrixView &Out, FMatrixView const &A, FScalar fFactor = 1.0, uint Flags = 0 );

// Out = A * B, Out,B being vectors
void Mxva( FScalar RESTRICT *pOut, FMatrixView const &A, FScalar const RESTRICT *pIn );
// Out += A * B, Out,B being vectors
void Mxvb( FScalar RESTRICT *pOut, FMatrixView const &A, FScalar const RESTRICT *pIn );
/// trace of a matrix
FScalar Trace( FMatrixView const &In );
/// dot-product of two matrices. dot(A^T B) = Tr(A * B).
FScalar Dot( FMatrixView const &A, FMatrixView const &B );
/// scale matrix by a factor in place
void Scale( FMatrixView &Out, FScalar fFactor );
/// diagonalize a matrix in place and store the eigenvalues at pEigenValues.
/// Input matrix must be symmetric and pEigenValues must hold room for nRows==nCols
/// entries.
void Diagonalize( FMatrixView &InOut, FScalar *pEigenValues, FMemoryStack &Mem );
/// As Diagonalize(), but eigenvalues are returned in descending order.
void Diagonalize_LargeEwFirst( FMatrixView &InOut, FScalar *pEigenValues, FMemoryStack &Mem );


// calculate Cholesky factorization M = L * L^T
void CalcCholeskyFactors(FMatrixView M);
// invert lower triangular matrix L. Upper part ignored.
void InvertTriangularMatrix(FMatrixView L);
// Solve L X = B where L is a lower triangular matrix.
void TriangularSolve(FMatrixView RhsSol, FMatrixView const &L);
// set A := L * A where L is a lower triangular matrix (Side == 'L').
// or  A := A * L where L is a lower triangular matrix (Side == 'R')
void TriangularMxm(FMatrixView A, FMatrixView const &L, char Side, double fFactor = 1.0);
// solve linear least squares ||M * x - b||_2 -> min subject to ||x||_2 ->min
// using singular value decomposition of A. Singular values below fThrRel are
// treated as zero. A negative value indicates the usage of machine precision.
// Both M and RhsSol are overwritten.
void LeastSquaresSolve(FMatrixView M, FMatrixView RhsSol, double fThrRel, FMemoryStack &Mem);


// // Solve X L = B where L is a lower triangular matrix.
// void TriangularSolveRight(FMatrixView RhsSol, FMatrixView const &L);
// Solve A X = B, where A = L*L^T and L is a lower triangular matrix.
void CholeskySolve(FMatrixView RhsAndSolution, FMatrixView const &L);
void CholeskySolve(double *pRhsAndSolution, FMatrixView const &L);

void Symmetrize(FMatrixView M, double Phase = 1);

/// creates a view representing the sub-matrix of size nRow x nCol
/// starting at (iRow,iCol).
/// The returned view still references the original data (i.e., In and Out are aliased)
FMatrixView Select( FMatrixView const &In, uint iRow, uint iCol, uint nRows, uint nCols );

// exactly same as above, but allowing for temporaries as first argument (C++ ftw...).
inline void Add0( FMatrixView Out, FMatrixView const &In, FScalar fFactor = 1.0 ) {
    return Add(Out,In,fFactor);
}
inline void Move0( FMatrixView Out, FMatrixView const &In, FScalar fFactor = 1.0 ) {
    return Move(Out, In, fFactor);
}
inline void Mxm0( FMatrixView Out, FMatrixView const &A, FMatrixView const &B, uint Flags = 0 ) {
    return Mxm(Out, A, B, Flags);
}
inline void Mxm0( FMatrixView Out, FMatrixView const &A, FMatrixView const &B, FScalar fFactor, uint Flags = 0 ) {
    return Mxm(Out, A, B, fFactor, Flags);
}
inline void Scale0( FMatrixView Out, FScalar fFactor ) {
    return Scale(Out, fFactor);
}



// An interface to FMatrixView which owns its data.
// Will adjust to input dimensions of operations as long as the actual
// matrix size does not exceed nMaxSize.
struct FStackMatrix : public FMatrixView
{
    FStackMatrix( uint nMaxSize_, FMemoryStack *pMemoryStack_ )
        : nMaxSize(nMaxSize_), pMemory(pMemoryStack_)
    {
        pMemory->Align(8);
        pMemory->Alloc(pData, nMaxSize_);
        Reshape(0,0);
    };

    FStackMatrix( uint nRows_, uint nCols_, FMemoryStack *pMemoryStack_ )
        : nMaxSize(nRows_*nCols_), pMemory(pMemoryStack_)
    {
        pMemory->Align(8);
        pMemory->Alloc(pData, nMaxSize);
        Reshape(nRows_, nCols_);
    };

    ~FStackMatrix(){ Free(); }
    void Free() { if (pData) pMemory->Free(pData); pData = 0; }
    void Reshape( uint nRows_, uint nCols_ ){
        assert( nRows_ * nCols_ <= nMaxSize );
        nRows = nRows_; nCols = nCols_; nRowSt = 1; nColSt = nRows_;
    };
private:
    uint
        nMaxSize;
    FMemoryStack
        *pMemory; // a
};

/// Out = f * In
void Move( FStackMatrix &Out, FMatrixView const &In, FScalar fFactor = 1.0 );
/// Out = A * B
void Mxm( FStackMatrix &Out, FMatrixView const &A, FMatrixView const &B, uint Flags = 0 );
/// Out = f * A * B
void Mxm( FStackMatrix &Out, FMatrixView const &A, FMatrixView const &B, FScalar fFactor, uint Flags = 0 );

// Note: the other operations from FMatrixView should work unchanged
// on FStackMatrix, because they do not require a change in the shape.

void PrintMatrixGen( std::ostream &out, double const *pData,
        uint nRows, uint nRowStride, uint nCols, uint nColStride,
        std::string const &Name );


struct FSmhOptions{
    double
        // ignore eigenvectors with eigenvalue < std::max(ThrAbs, MaxEw * ThrRel)
        ThrAbs, ThrRel;
    char const
        *pDelMsg; // may be 0.
    std::ostream
        *pXout; // may be 0 if pDelMsg is 0.
    FSmhOptions( double ThrAbs_ = 0, double ThrRel_ = 0, char const *pDelMsg_ = 0,
            std::ostream *pXout_ = &xout )
        : ThrAbs(ThrAbs_), ThrRel(ThrRel_), pDelMsg(pDelMsg_), pXout(pXout_)
    {}
};

// calculates M^{-1/2} in place.
void CalcSmhMatrix( FMatrixView M, FMemoryStack &Mem, FSmhOptions const &Opt );
// orthogonalizes C in place, where S is C's overlap matrix.
// both C and S are overwritten.
void OrthSchmidt1(FMatrixView S, FMatrixView C, FMemoryStack &Mem);

// permute rows of matrix such that
//    NewM[iRow,iCol] = OldM[P[iRow], iCol]
void PermuteRows(FMatrixView M, uint const *P, FMemoryStack &Mem);
void Permute(FMatrixView &M, uint const *pRowPerm, uint const *pColPerm, FMemoryStack &Mem);
void WriteMatrixToFile(std::string const &FileName, FMatrixView const &M, std::string const &MatrixName, uint *pRowPerm = 0, uint *pColPerm = 0);

} // namespace ct

#endif // CT8k_MATRIX_H

// kate: space-indent on; tab-indent on; backspace-indent on; tab-width 4; indent-width 4; mixedindent off; indent-mode normal;
