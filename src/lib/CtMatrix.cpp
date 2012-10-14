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
#include <boost/format.hpp>
#include <cmath>
#include <fstream> // only for matrix import/export
using boost::format;

#include "fmt.h"
#include "CtCommon.h"
#include "CtMatrix.h"

#include "CxAlgebra.h" // fortran BLAS/LAPACK externals.


namespace ct {

inline void AssertCompatible( FMatrixView const &A, FMatrixView const &B )
{
    assert( A.nRows == B.nRows );
    assert( A.nCols == B.nCols );
};

void AssertCompatible1( FMatrixView const &A, FMatrixView const &B )
{
    AssertCompatible(A,B);
};


FMatrixView Transpose( FMatrixView const &In )
{
    return FMatrixView( In.pData,
        In.nCols, In.nRows, In.nColSt, In.nRowSt );
};

FMatrixView Select( FMatrixView const &In, uint iRow, uint iCol, uint nRows, uint nCols )
{
    assert(iRow + nRows <= In.nRows && iCol + nCols <= In.nCols);
    return FMatrixView(In.pData + iRow*In.nRowSt + iCol*In.nColSt,
        nRows, nCols, In.nRowSt, In.nColSt);
}

bool FMatrixView::IsSymmetric(FScalar Thresh) const
{
    if ( nRows != nCols )
        return false;
    for ( uint iRow = 0; iRow < nRows; ++ iRow )
        for ( uint iCol = 0; iCol < iRow; ++ iCol )
            if ( std::abs((*this)(iRow, iCol) - (*this)(iCol, iRow)) > Thresh )
                return false;
    return true;
};

// Out += f * In
void Add( FMatrixView &Out, FMatrixView const &In, FScalar fFactor )
{
    AssertCompatible(Out,In);
    for ( uint nCol = 0; nCol < Out.nCols; ++ nCol )
        for ( uint nRow = 0; nRow < Out.nRows; ++ nRow )
            Out(nRow,nCol) += fFactor * In(nRow,nCol);
};

// Out = f * In
void Move( FMatrixView &Out, FMatrixView const &In, FScalar fFactor )
{
    AssertCompatible(Out,In);
    for ( uint nCol = 0; nCol < Out.nCols; ++ nCol )
        for ( uint nRow = 0; nRow < Out.nRows; ++ nRow )
            Out(nRow,nCol) = fFactor * In(nRow,nCol);
};

/// Out = f * A * B
void Mxm( FMatrixView &Out, FMatrixView const &A, FMatrixView const &B,
    FScalar fFactor, uint Flags )
{
    assert( Out.nRows == A.nRows );
    assert( Out.nCols == B.nCols );
    assert( A.nCols == B.nRows );

    assert( A.nRowSt == 1 || A.nColSt == 1 );
    assert( B.nRowSt == 1 || B.nColSt == 1 );
    assert( Out.nRowSt == 1 || Out.nColSt == 1 );
    // ^- otherwise dgemm directly not applicable. Would need local copy
    // of matrix/matrices with compressed strides.

    double
        Beta = (0 != (Flags & MXM_AddToDest))? 1.0 : 0.0;
    char
        TransA = (A.nRowSt == 1? 'N' : 'T'),
        TransB = (B.nRowSt == 1? 'N' : 'T');
    FORTINT
        lda,
        ldb,
//         lda = (A.nRowSt == 1)? A.nColSt : A.nRowSt,
//         ldb = (B.nRowSt == 1)? B.nColSt : B.nRowSt,
        ldc = (Out.nRowSt == 1)? Out.nColSt : Out.nRowSt;
    if ( A.nRowSt == 1 )
        // this is to bypass the screwball checks in MKL for
        // lda >= nrows, even if the stride is 1 because there is
        // only one column.
        lda = std::max(A.nColSt, A.nRows);
    else
        lda = A.nRowSt;
    if ( B.nRowSt == 1 )
        ldb = std::max(B.nColSt, B.nRows);
    else
        ldb = B.nRowSt;

    DGEMM( TransA, TransB, Out.nRows, Out.nCols, A.nCols,
        fFactor, A.pData, lda, B.pData, ldb, Beta, Out.pData, ldc );
}

// Out = A * B
void Mxm( FMatrixView &Out, FMatrixView const &A, FMatrixView const &B, uint Flags )
{
    return Mxm(Out, A, B, 1.0, Flags);
};

/// Out = f * A^T * A for symmetric matrix Out.
void SyrkT( FMatrixView &Out, FMatrixView const &A, FScalar fFactor, uint Flags )
{
//     return Mxm( Out, Transpose(A), A, fFactor, AddToDest );
    assert( Out.nRows == Out.nCols );
    assert( A.nCols == Out.nCols );

    assert( A.nRowSt == 1 || A.nColSt == 1 );
    assert( Out.nRowSt == 1 || Out.nColSt == 1 );
    // ^- otherwise dsyrk directly not applicable. Would need local copy
    // of matrix/matrices with compressed strides.

    double
        Beta = (0 != (Flags & MXM_AddToDest))? 1.0 : 0.0;
    char
        TransA = (A.nRowSt == 1? 'T' : 'N');
    FINTARG
        lda = (A.nRowSt == 1)? A.nColSt : A.nRowSt,
        ldc = (Out.nRowSt == 1)? Out.nColSt : Out.nRowSt;

    DSYRK( 'L', TransA, Out.nRows, A.nRows, fFactor,
        A.pData, lda, Beta, Out.pData, ldc );
    // fix up the upper triangle of the matrix.
    for ( uint iCol = 0; iCol < Out.nCols; ++ iCol )
        for ( uint iRow = 0; iRow < iCol; ++ iRow )
            Out(iRow,iCol) = Out(iCol,iRow);
}

/// Out = f * A * A^T for symmetric matrix Out.
void Syrk( FMatrixView &Out, FMatrixView const &A, FScalar fFactor, uint Flags )
{
    return SyrkT(Out, Transpose(A), fFactor, Flags);
}



// // Out = f * A * B
// void Mxm( FMatrixView &Out, FMatrixView const &A, FMatrixView const &B,
//     FScalar fFactor, bool AddToDest )
// {
//     assert( Out.nRows == A.nRows );
//     assert( Out.nCols == B.nCols );
//     assert( A.nCols == B.nRows );
//
//     if ( !AddToDest )
//         Out.Clear();
//
//     // [pack B]
//     for ( uint nRow = 0; nRow < Out.nRows; ++ nRow ){
//         // [pack A]
//         for ( uint nCol = 0; nCol < Out.nCols; ++ nCol ){
//             for ( uint nLink = 0; nLink < A.nCols; ++ nLink ){
//                 Out(nRow,nCol) += fFactor * A(nRow,nLink) * B(nLink,nCol);
//             }
//         }
//     }
// }


// diagonalize a matrix in place and store the eigenvalues at pEigenValues.
// Input matrix must be symmetric and pEigenValues must hold room for nRows==nCols
// entries.
void Diagonalize( FMatrixView &InOut, FScalar *pEigenValues, FMemoryStack &Mem )
{
    assert_rt(InOut.nRowSt == 1);
    assert(InOut.nRows == InOut.nCols);

//     double
//         *pWork;
//     FORTINT
//         info = 0,
//         nb = 256,
//         // ^- no idea?! actually calling ILAENV to get nb could cause
//         // considerable portability problems due to string passing.
//         lWork = (nb+2) * InOut.nRows;
//
//     Mem.Alloc(pWork, lWork);
//     xout << format("diagonalize: %i x %i [ld=%i]  lwrk: %i") % InOut.nRows % InOut.nRows % InOut.nColSt % lWork << std::endl;
//     DSYEV('V', 'L', InOut.nRows, InOut.pData, InOut.nColSt,
//         pEigenValues, pWork, lWork, info);
//     Mem.Free(pWork);

    FORTINT N = InOut.nRows, nWork, info = 0;
    FORTINT ldH = InOut.nColSt;
    nWork = 128*N;
    double *pWork = (double*)::malloc(sizeof(double)*nWork);
    DSYEV('V', 'L', N, InOut.pData, ldH, pEigenValues, pWork, nWork, info );
    ::free(pWork);
    if ( info != 0 ) throw std::runtime_error("dsyev failed.");


    if ( info != 0 ) {
        std::stringstream str;
        str << "Something went wrong when trying to diagonalize a " << InOut.nRows << "x" << InOut.nCols << " matrix. "
            << "DSYEV returned error code " << info << ".";
        throw std::runtime_error(str.str());
    }
};

// As Diagonalize(), but eigenvalues are returned in descending order.
void Diagonalize_LargeEwFirst( FMatrixView &InOut, FScalar *pEigenValues, FMemoryStack &Mem )
{
    Scale(InOut, -1.0);
    Diagonalize(InOut, pEigenValues, Mem);
    Scale(pEigenValues, -1.0, InOut.nCols);
}


// Out = A * B, Out,B being vectors
void Mxva( FScalar RESTRICT *pOut, FMatrixView const &A, FScalar const RESTRICT *pIn )
{
    // TODO: inline this and call mxva/mxvb, or make a better implementation.
    // These routines are actually used in the coupling coefficient calculations
    // and should be reasonably fast.
    for ( uint nRow = 0; nRow < A.nRows; ++ nRow ){
        FScalar
            d = 0;
        for ( uint nCol = 0; nCol < A.nCols; ++ nCol )
            d += A(nRow,nCol) * pIn[nCol];
        pOut[nRow] = d;
    }
}

// Out += A * B, Out,B being vectors
void Mxvb( FScalar RESTRICT *pOut, FMatrixView const &A, FScalar const RESTRICT *pIn )
{
    FMatrixView Out(pOut, A.nRows,1);
    return Mxm( Out, A, FMatrixView(const_cast<double*>(pIn),A.nCols,1), 1.0, true );

    // TODO: inline this and call mxva/mxvb, or make a better implementation.
    for ( uint nRow = 0; nRow < A.nRows; ++ nRow ){
        FScalar
            d = 0;
        for ( uint nCol = 0; nCol < A.nCols; ++ nCol )
            d += A(nRow,nCol) * pIn[nCol];
        pOut[nRow] += d;
    }
};

// trace of a matrix
FScalar Trace( FMatrixView const &In )
{
    assert(In.nRows == In.nCols);
    FScalar
        res = 0.0;
    for ( uint n = 0; n < In.nRows; ++ n )
        res += In(n,n);
    return res;
};

// dot-product of two matrices. dot(A^T B) = Tr(A * B).
FScalar Dot( FMatrixView const &A, FMatrixView const &B )
{
    AssertCompatible(A,B);
    FScalar
        res = 0.0;
    for ( uint nCol = 0; nCol < A.nCols; ++ nCol )
        for ( uint nRow = 0; nRow < A.nRows; ++ nRow )
            res += A(nRow,nCol) * B(nRow,nCol);
    return res;
};

// scale matrix by a factor in place
void Scale( FMatrixView &Out, FScalar fFactor )
{
    for ( uint nCol = 0; nCol < Out.nCols; ++ nCol )
        for ( uint nRow = 0; nRow < Out.nRows; ++ nRow )
            Out(nRow,nCol) *= fFactor;
};


void FMatrixView::Clear()
{
    for ( uint nCol = 0; nCol < nCols; ++ nCol )
        for ( uint nRow = 0; nRow < nRows; ++ nRow )
            (*this)(nRow,nCol) = 0.0f;
};

void FMatrixView::SetIdentity()
{
    assert( nRows == nCols );
    Clear();
    for ( uint i = 0; i < nRows; ++ i )
       (*this)(i,i) = 1.0f;
};

std::vector<FScalar> FMatrixView::GetDiagonal() const
{
    assert( nRows == nCols );
    std::vector<FScalar>
        r(nRows);
    for ( uint i = 0; i < nRows; ++ i )
       r[i] = (*this)(i,i);
    return r;
};


void FMatrixView::Print( std::ostream &out, std::string const &Name ) const
{
    PrintMatrixGen( out, pData, nRows, nRowSt, nCols, nColSt, Name );
};


void Move( FStackMatrix &Out, FMatrixView const &In, FScalar fFactor ) {
    Out.Reshape(In.nRows, In.nCols);
    Move( *static_cast<FMatrixView*>(&Out), In, fFactor );
};

void Mxm( FStackMatrix &Out, FMatrixView const &A, FMatrixView const &B, uint Flags )
{
    assert( (0==(Flags & MXM_AddToDest)) || (Out.nRows == A.nRows && Out.nCols == B.nCols ) );
    Out.Reshape(A.nRows, B.nCols);
    Mxm( *static_cast<FMatrixView*>(&Out), A, B, Flags );
};

void Mxm( FStackMatrix &Out, FMatrixView const &A, FMatrixView const &B, FScalar fFactor, uint Flags )
{
    assert( (0==(Flags & MXM_AddToDest)) || (Out.nRows == A.nRows && Out.nCols == B.nCols ) );
    Out.Reshape(A.nRows, B.nCols);
    Mxm( *static_cast<FMatrixView*>(&Out), A, B, fFactor, Flags );
};


/// prints a general rectangular matrix specified by memory layout.
/// Element m(r,c) is given by pData[ r * nRowStride + c * nColStride ].
void PrintMatrixGen( std::ostream &out, double const *pData,
        uint nRows, uint nRowStride, uint nCols, uint nColStride,
        std::string const &Name )
{
    using namespace fmt;
    uint
        nFloatWidth = 14, //25, //20, //14,
        nFloatPrec  = 8;  //16; //15; //8;
    if ( Name != "" ){
        out << "  Matrix " << Name << ", " << nRows << "x" << nCols << "." << std::endl;
        if ( nRows * nCols == 0 )
            return;
    }
    out << "           ";
    for ( uint nCol = 0; nCol < nCols; ++ nCol )
        out << fi(nCol,nFloatWidth-3) << "   ";
    out << "\n";
    for ( uint nRow = 0; nRow < nRows; ++ nRow ) {
        out << "    " << fi(nRow,4) << "   ";
        for ( uint nCol = 0; nCol < nCols; ++ nCol )
        {
            double const
                &f = pData[ nRow * nRowStride + nCol * nColStride ];
//             if ( *(const uint32_t*)&f != 0xffffffff )
                out << ff(f, nFloatWidth, nFloatPrec);
//                 out << fe(f, nFloatWidth, nFloatPrec);
/*            else
                out << fc("`", nFloatWidth,1);*/
        };
        out << "\n";
    }
    out.flush();
};


std::size_t FMatrixView::TriangularStorageSize(int Sign) const
{
    assert( nRows == nCols );
    assert( Sign == +1 || Sign == -1 );
//     if ( Sign == 1 )
//         return (nRows*(nRows+1))/2;
//     else
//         return (nRows*(nRows-1))/2;
    // ^- Molpro does not actually omit storing the diagonal of antisymmetric
    //    matrices. By doing it also this way we can keep binary compatibility.
    return (nRows*(nRows+1))/2;
}

void FMatrixView::TriangularExpand(int Sign, FScalar *pTriangData)
{
    assert( nRowSt == 1 && nColSt >= nRows ); // <- required?
    assert( nRows == nCols );
    if ( pTriangData == 0 )
        pTriangData = pData; // in-place.
    std::size_t
        iOff = nRows*(nRows+1)/2;
    uint
        nRowSt_ = nRowSt,
        nColSt_ = nColSt;
    assert( nRowSt == 1 || nColSt == 1 );
    if ( nRowSt_ > nColSt_ ) std::swap(nRowSt_,nColSt_);
    for ( uint iCol_ = nCols; iCol_ != 0; -- iCol_ )
        for ( uint iRow_ = iCol_; iRow_ != 0; -- iRow_ ) {
            iOff -= 1;
            pData[(iRow_-1) * nRowSt_ + (iCol_-1) * nColSt_] = pTriangData[iOff];
        }
    assert(iOff == 0);
    for ( uint iCol_ = 0; iCol_ != nCols; ++ iCol_ )
        for ( uint iRow_ = 0; iRow_ != iCol_; ++ iRow_ )
            pData[iCol_ * nRowSt_ + iRow_ * nColSt_] = Sign * pData[iCol_ * nColSt_ + iRow_ * nRowSt_];
    ClearEmptySpace();
}

void FMatrixView::TriangularReduce(int Sign, FScalar *pTriangData)
{
    assert( nRowSt == 1 && nColSt >= nRows ); // <- required?
    assert( nRows == nCols );
    if ( pTriangData == 0 )
        pTriangData = pData; // in-place.
    std::size_t
        iOff = 0;
    uint
        nRowSt_ = nRowSt,
        nColSt_ = nColSt;
    assert( nRowSt == 1 || nColSt == 1 );
    if ( nRowSt_ > nColSt_ ) std::swap(nRowSt_,nColSt_);
#ifdef _DEBUG
    for ( uint iCol = 0; iCol < nCols; ++ iCol )
        for ( uint iRow = 0; iRow <= iCol; ++ iRow ) {
            assert( std::abs(pData[iRow * nRowSt_ + iCol * nColSt_] -
              (double)Sign * pData[iCol * nRowSt_ + iRow * nColSt_]) < 1e-10 );
        }
    iOff = 0;
#endif
    for ( uint iCol = 0; iCol < nCols; ++ iCol )
        for ( uint iRow = 0; iRow <= iCol; ++ iRow ) {
//             xout << boost::format("pack: T[%3i] <- Q[%3i]  iRow=%2i  iCol=%2i  F=%f") % iOff % (iRow * nRowSt_ + iCol * nColSt_) % iRow % iCol % pData[iRow * nRowSt_ + iCol * nColSt_] << std::endl;
            pTriangData[iOff] = pData[iRow * nRowSt_ + iCol * nColSt_];
            iOff += 1;
        }
//     PrintMatrixGen(xout, pTriangData, 1, 1, iOff, 1, "TRIANGULAR DATA AFTER PACKING:");
}

void FMatrixView::TriangularExpandOrReduce(int Sign, FScalar *pTriangData, FTriangularTransformFlags Direction)
{
    assert( Direction == TRIANG_Expand || Direction == TRIANG_Reduce);
    if ( Direction == TRIANG_Expand )
        return TriangularExpand(Sign, pTriangData);
    else
        return TriangularReduce(Sign, pTriangData);
}

void FMatrixView::ClearEmptySpace()
{
    if ( nRowSt == 1 ) {
        if ( nColSt == nRows )
            return;
        for ( uint iCol = 0; iCol < nCols; ++ iCol )
            for ( uint iEmptyRow = nRows; iEmptyRow != nColSt; ++ iEmptyRow )
                pData[iEmptyRow + iCol*nColSt] = 0;

    } else {
        assert( nColSt == 1 );
        if ( nRowSt == nCols )
            return;
        for ( uint iRow = 0; iRow < nRows; ++ iRow )
            for ( uint iEmptyCol = nCols; iEmptyCol != nRowSt; ++ iEmptyCol )
                pData[iEmptyCol + iRow*nRowSt] = 0;
    }
}



// maybe move this to Matrix.cpp?
void CalcSmhMatrix( FMatrixView M, FMemoryStack &Mem, FSmhOptions const &Opt )
{
    assert_rt(M.nRows == M.nCols);
    uint
        N = M.nRows;
    FStackMatrix
        T1(N, N, &Mem),
        T2(N, N, &Mem);
    Move(T1, M); // should get rid of possible row strides of M (diag2 doesn't like them).
//     Move(T1, M, -1.); // should get rid of possible row strides of M (diag2 doesn't like them).
    double
        *pEw; // eigenvalues
    Mem.Alloc(pEw, N);

//     T1.Print(xout, "(A|J|B)");
    Diagonalize(T1, pEw, Mem);
//     Scale(pEw, -1., N);

/*    xout << "Spectrum for " << Opt.pDelMsg << "\n    ";
    for ( uint i = 0; i < N; ++ i )
        xout << format("  %10.3e") % pEw[i];
    xout << std::endl;*/
    if ( Opt.pXout && Opt.pDelMsg ) {
        std::ostream &xout = *Opt.pXout;
        xout << "" << Opt.pDelMsg << ": ";
        uint
            nLen = 0;
        for ( char const *p = Opt.pDelMsg; *p; ++p ) nLen += 1;
        for ( uint i = 27; i > nLen; -- i ) xout << " ";
        xout << format("  EwMin =%10.3e  EwMax=%10.3e  Ratio=%10.3e")
                % pEw[0] % pEw[N-1] % (pEw[0]/pEw[N-1])
             << std::endl;
    }
    uint
        nDel = 0;
    for ( uint iCol = 0; iCol < N; ++ iCol ){
        double
            Ew = pEw[iCol],
            f = 0;
        if ( Ew < Opt.ThrAbs || Ew < Opt.ThrRel * pEw[N-1] ) {
            f = 0;
            nDel += 1;
        }
        else
            f = 1.0/std::sqrt(std::sqrt(Ew));
        for ( uint iRow = 0; iRow < N; ++ iRow )
            T2(iRow,N-iCol-1) = f * T1(iRow,iCol); // <- invert order: smallest coeffs first.
    }

    if ( nDel != 0 && Opt.pDelMsg ) {
        std::ostream
            &out = *(Opt.pXout? Opt.pXout : &xout);
        out << "Warning: Deleted " << nDel
            << " vectors while constructing " << Opt.pDelMsg << "." << std::endl;
    }

//     Mxm(M, T2, Transpose(T2));
    Syrk(M, T2);
    Mem.Free(pEw);
};


void LinearSolveGenSvd( double *pOut, double *pMat, uint nStr, double *pIn, uint nDim, double Thr, FMemoryStack &Mem )
{
    FStackMatrix
        Binv(nDim, nDim, &Mem),
        Ew(nDim,1, &Mem);
    Move(Binv, FMatrixView(pMat,nDim,nDim,1,nStr));
//     Binv.Print(xout, "SolveSvd: Input");
/*    FMatrixView(pIn,1,nDim).Print(xout,"LINSOLVE: Rhs");
    Binv.Print(xout,"LINSOLVE: B");*/
    Diagonalize( Binv, Ew.pData, Mem );
//     Binv.Print(xout,"LINSOLVE: Bev");
    for ( uint i = 0; i < nDim; ++ i )
        pOut[i] = 0;

    for ( uint iEw_ = nDim; iEw_ != 0; -- iEw_ ){
        uint   iEw = iEw_ - 1;
        double ew = Ew.pData[iEw];
//         _xout0("ew(" << iEw << ") = " << ew);
        if ( std::abs(ew) < Thr )
            continue;
        double
            *pEv = Binv.pData + nDim * iEw,
            d = Dot(pEv, pIn, nDim);
//         _xout0("dot(" << iEw << ") = " << d << "  ew: " << ew);
        Add(pOut, pEv, d/ew, nDim);
//         Add(pOut, pEv, d/(std::exp(ew)-1.), nDim);
    }
//     FMatrixView(pOut,1,nDim).Print(xout,"LINSOLVE: Lhs");
};


/*void OrthSchmidt1(FMatrixView S, FMatrixView C)
{
    assert(S.nRows == S.nCols && C.nRows == C.nCols & S.nRows == C.nRows);



}*/


// solve linear least squares ||M * x - b||_2 -> min subject to ||x||_2 ->min
// using singular value decomposition of A. Singular values below fThrRel are
// treated as zero. A negative value indicates the usage of machine precision.
// Both M and RhsSol are overwritten.
void LeastSquaresSolve(FMatrixView M, FMatrixView RhsSol, double fThrRel, FMemoryStack &Mem)
{
    using std::min;
    using std::max;
    assert(M.nRowSt == 1);
    assert(RhsSol.nRowSt == 1);
    assert(RhsSol.nRows == M.nCols);

    FORTINT
        nRank = -1,
        info = 0,
        // minimal work space required according to docs.
        nWork = 3*min(M.nRows,M.nCols) + max(max(2*min(M.nRows, M.nCols), max(M.nRows, M.nCols)), RhsSol.nCols);
    double
        fWork = 0,
        *pWork,
        *pSig; // singular values (decreasing order).
    Mem.Alloc(pSig, min(M.nRows, M.nCols));
    DGELSS(M.nRows, M.nCols, RhsSol.nCols, M.pData, M.nColSt,
        RhsSol.pData, RhsSol.nColSt, pSig, fThrRel, &nRank,
        &fWork, -1, &info );
    if ( info != 0 ) throw std::runtime_error("dgelss failed in workspace query.");
    nWork = static_cast<FORTINT>(fWork + 0.5);
    Mem.Alloc(pWork, nWork);
    DGELSS(M.nRows, M.nCols, RhsSol.nCols, M.pData, M.nColSt,
        RhsSol.pData, RhsSol.nColSt, pSig, fThrRel, &nRank,
        pWork, nWork, &info );
    if ( info != 0 ) throw std::runtime_error("dgelss failed.");
    Mem.Free(pWork);
    Mem.Free(pSig);
};








// // transpose square matrix in place.
// static void cd_transp_1(double *c, uint nb, uint n)
// {
//     for ( uint i = 0; i < n; ++ i ){
//         for ( uint j = 0; j < i; ++ j ) {
//             double t = c[i + nb*j];
//             c[i + nb*j] = c[j + nb*i];
//             c[j + nb*i] = t;
//         }
//     }
// }

// given the n x n MO overlap matrix s = SMO = C^T SAO C, orthogonalize
// c via Schmidt orthogonalization. s is overwritten.
// nbs and nbc are the lead dimensions of s and c, respectively.
static void OrthSchmidt1(double *S, uint nbs, uint n, double *C, uint nbc, uint nc, FMemoryStack &Mem)
{
    // Let C' := C U with C being the original orbitals, C' being the new,
    // orthogonalized orbitals, and U being an upper triangular matrix.
    // We want U^T S U to be identity.
    // The upper diagonality makes sure that the i'th row of C' is only composed
    // of linear combinations of rows j of C with j <= i. This is equivalent
    // to the Schmidt incremental orthogonalization.
    FORTINT
        info = 0;
    // factorize S into S = L * L^T.
    DPOTRF('L', n, S, nbs, &info);
    if (info != 0) throw std::runtime_error("orths1: dpotrf failed.");
    // solve L^T U = id  (note: L^T U = id implies U^T S U = (U^T L) (L^T S) = id)
    //       -> U = (L^T)^{-1}
    //       -> C' = C * U = C * (L^T)^{-1}
    DTRSM('R','L','T','N',n,nc,1.0,S,nbs,C,nbc,&info);
    if (info != 0) throw std::runtime_error("orths1: dtrtrs failed.");
}

void OrthSchmidt1(FMatrixView S, FMatrixView C, FMemoryStack &Mem)
{
    assert(S.nRows == S.nCols && S.nCols == C.nCols);
    assert_rt(S.nRowSt == 1 && C.nRowSt == 1); // sorry. otherwise need to copy stuff around.

    OrthSchmidt1(S.pData, S.nColSt, S.nRows, C.pData, C.nColSt, C.nCols, Mem);
}


// calculate Cholesky factorization M = L * L^T
void CalcCholeskyFactors(FMatrixView M)
{
    assert(M.nRowSt == 1 && M.nRows == M.nCols);
    FORTINT
        info = 0;
    // factorize M into M = L * L^T.
    DPOTRF('L', M.nRows, M.pData, M.nColSt, &info);
    if (info != 0) throw std::runtime_error("CalcCholeskyMatrix: dpotrf failed.");
    for ( uint i = 0; i < M.nRows; ++ i )
        for ( uint j = i+1; j < M.nRows; ++j )
            M(i,j) = 0;
};

// invert lower triangular matrix L. Upper part ignored.
void InvertTriangularMatrix(FMatrixView L)
{
    assert(L.nRowSt == 1 && L.nRows == L.nCols);
    FORTINT
        info = 0;
    DTRTRI('L','N', L.nRows, L.pData, L.nColSt, &info);
    if (info != 0) throw std::runtime_error("InvertTriangularMatrix: dtrtri failed.");
};

// Solve L X = B where L is a lower triangular matrix.
void TriangularSolve(FMatrixView RhsSol, FMatrixView const &L)
{
    assert(L.nRowSt == 1 && L.nRows == L.nCols);
    assert(RhsSol.nRowSt == 1 || RhsSol.nColSt == 1);
    assert(L.nRows == RhsSol.nRows);
    FORTINT
        info = 0;
    if ( RhsSol.nRowSt == 1 )
        DTRSM('L', 'L', 'N', 'N', RhsSol.nRows, RhsSol.nCols, 1.0, L.pData, L.nColSt, RhsSol.pData, RhsSol.nColSt, &info);
    else
        // A X^T == B^T
        // <-> X A^T == B
        DTRSM('R', 'L', 'T', 'N', RhsSol.nCols, RhsSol.nRows, 1.0, L.pData, L.nColSt, RhsSol.pData, RhsSol.nRowSt, &info);
    if (info != 0) throw std::runtime_error("TriangularSolve: dtrtrs failed.");
};

// // Solve X L = B where L is a lower triangular matrix.
// void TriangularSolveRight(FMatrixView RhsSol, FMatrixView const &L)
// {
//     assert(L.nRowSt == 1 && L.nRows == L.nCols);
//     assert(RhsSol.nRowSt == 1);
//     assert(L.nRows == RhsSol.nCols);
//     FORTINT
//         info = 0;
//     DTRSM('R', 'L', 'N', 'N', RhsSol.nRows, RhsSol.nCols, L.pData, L.nColSt, RhsSol.pData, RhsSol.nColSt, &info);
// };

// set A := L * A where L is a lower triangular matrix (LR = 'L').
// or  A := A * L where L is a lower triangular matrix (LR = 'R')
void TriangularMxm(FMatrixView RhsSol, FMatrixView const &L, char Side, double fFactor)
{
    assert(L.nRowSt == 1 && L.nRows == L.nCols);
    assert(Side == 'L' || Side == 'R');
    if ( Side == 'L' ) {
        assert(RhsSol.nRowSt == 1 || RhsSol.nColSt == 1);
        assert(L.nRows == RhsSol.nRows);
        if ( RhsSol.nRowSt == 1 )
            DTRMM('L', 'L', 'N', 'N', RhsSol.nRows, RhsSol.nCols, fFactor, L.pData, L.nColSt, RhsSol.pData, RhsSol.nColSt);
        else
            // A X^T == B^T
            // <-> X A^T == B
            DTRMM('R', 'L', 'T', 'N', RhsSol.nCols, RhsSol.nRows, fFactor, L.pData, L.nColSt, RhsSol.pData, RhsSol.nRowSt);
    } else {
        assert(RhsSol.nRowSt == 1 || RhsSol.nColSt == 1);
        assert(L.nRows == RhsSol.nCols);

        if ( RhsSol.nRowSt == 1 )
            DTRMM('R', 'L', 'N', 'N', RhsSol.nRows, RhsSol.nCols, fFactor, L.pData, L.nColSt, RhsSol.pData, RhsSol.nColSt);
        else
            // A X^T == B^T
            // <-> X A^T == B
            DTRMM('L', 'L', 'T', 'N', RhsSol.nCols, RhsSol.nRows, fFactor, L.pData, L.nColSt, RhsSol.pData, RhsSol.nRowSt);
    }
};


// Solve A X = B, where A = L*L^T and L is a lower triangular matrix.
void CholeskySolve(FMatrixView RhsSol, FMatrixView const &L)
{
    assert(L.nRowSt == 1 && L.nRows == L.nCols);
    assert(RhsSol.nRowSt == 1);
    assert(L.nRows == RhsSol.nRows);
    FORTINT
        info = 0;
    DPOTRS('L', L.nRows, RhsSol.nCols, L.pData, L.nColSt, RhsSol.pData, RhsSol.nColSt, &info);
    if (info != 0) throw std::runtime_error("CholeskySolve: dpotrs failed.");
};

void CholeskySolve(double *pRhsAndSolution, FMatrixView const &L)
{
    CholeskySolve(FMatrixView(pRhsAndSolution,L.nRows,1), L);
};


void Symmetrize(FMatrixView M, double Phase)
{
   assert(M.IsSquare());
   for ( uint i = 0; i < M.nRows; ++ i )
      for ( uint j = 0; j <= i; ++ j ) {
         FScalar f = .5*(M(i,j) + Phase * M(j,i));
         M(i,j) = f;
         M(j,i) = Phase * f;
      }
}

// permute rows of matrix such that
//    NewM[iRow,iCol] = OldM[P[iRow], iCol]
void PermuteRows(FMatrixView M, uint const *P, FMemoryStack &Mem)
{
   for ( uint iCol = 0; iCol < M.nCols; ++ iCol ) {
      double *p;
      Mem.Alloc(p, M.nRows);
      for ( uint iRow = 0; iRow < M.nRows; ++ iRow )
         p[iRow] = M(P[iRow], iCol);
      for ( uint iRow = 0; iRow < M.nRows; ++ iRow )
         M(iRow, iCol) = p[iRow];
      Mem.Free(p);
   };
};

void Permute(FMatrixView &M, uint const *pRowPerm, uint const *pColPerm, FMemoryStack &Mem)
{
   PermuteRows(M, pRowPerm, Mem);
   PermuteRows(Transpose(M), pColPerm, Mem);
};

void ReadMatrixFromFile(FMatrixView &M, std::string const &FileName, FMemoryStack &Mem)
{
   // format:
   //    <matrix name>
   //    <nRows> x <nCols> <SYMMETRIC/RECTANGULAR>
   //    <irow> <icol> <m[irow,icol]>
   std::string
      MatrixName, Dummy, Type;
   std::ifstream
      File(FileName.c_str());
   if ( !File.good() )
      throw std::runtime_error("failed to open: " + FileName + ".");
   File >> MatrixName;
   uint
      nRows, nCols;
   std::getline(File, Dummy); // matrix name
   File >> nRows >> Dummy >> nCols >> Type;
   if ( Dummy != "x" || (Type != "RECTANGULAR" && Type != "SYMMETRIC") )
      throw std::runtime_error("unexpected format: input matrix " + FileName + ".");
   bool
      Symmetric = false;
   Symmetric = Type == "SYMMETRIC";
   if ( Symmetric && nRows != nCols )
      throw std::runtime_error("symmetric matrices must have equal number of columns and rows: input matrix " + FileName + ".");

   M = FMatrixView(0, nRows, nCols);
   Mem.ClearAlloc(M.pData, M.GetStridedSize());
   while ( true ) {
      int
         iRow, iCol;
      double
         f;
      File >> iRow >> iCol >> f;
      if ( !File.good() )
         break;
      if ( iRow < 1 || (unsigned)iRow > nRows || iCol < 1 || (unsigned)iCol > nCols ) {
         std::stringstream str;
         str << format("index out of bounds: iRow = %i  iCol = %i while reading %i x %i matrix %s.")
            % iRow % iCol % nRows % nCols % FileName;
         throw std::runtime_error(str.str());
      }
      M(iRow-1, iCol-1) = f;
      if ( Symmetric )
         M(iCol-1, iRow-1) = f;
   };
};

void WriteMatrixToFile(std::string const &FileName, FMatrixView const &M, std::string const &MatrixName, uint *pRowPerm, uint *pColPerm)
{
   // format:
   //    <matrix name>
   //    <nRows> x <nCols> <SYMMETRIC/RECTANGULAR>
   //    <irow> <icol> <m[irow,icol]>
   std::ofstream
      File(FileName.c_str());
   if ( !File.good() )
      throw std::runtime_error("failed to open: " + FileName + ".");
   bool
      Symmetric = M.IsSymmetric(1e-10);
   File << format("%s\n%i x %i %s\n")
      % MatrixName
      % M.nRows % M.nCols % (Symmetric? "SYMMETRIC" : "RECTANGULAR");

   for ( uint iCol = 0; iCol < M.nCols; ++ iCol )
      for ( uint iRow = 0; iRow < M.nRows; ++ iRow ) {
         if ( Symmetric && iCol < iRow )
            continue;
         uint
            iRow_ = pRowPerm? pRowPerm[iRow] : iRow,
            iCol_ = pColPerm? pColPerm[iCol] : iCol;
         FScalar
            f = M(iRow_,iCol_);
         if ( Symmetric )
            f = .5 * (f + M(iCol_,iRow_));
         if ( std::abs(f) > 1e-25 )
            File << format("%5i %5i %24.16e\n") % (iRow+1) % (iCol+1) % f;
      }
};



} // namespace ct

// kate: space-indent on; tab-indent on; backspace-indent on; tab-width 4; indent-width 4; mixedindent off; indent-mode normal;
