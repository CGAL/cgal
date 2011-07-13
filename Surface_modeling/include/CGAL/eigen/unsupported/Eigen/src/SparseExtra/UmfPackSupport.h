// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2008-2009 Gael Guennebaud <gael.guennebaud@inria.fr>
//
// Eigen is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 3 of the License, or (at your option) any later version.
//
// Alternatively, you can redistribute it and/or
// modify it under the terms of the GNU General Public License as
// published by the Free Software Foundation; either version 2 of
// the License, or (at your option) any later version.
//
// Eigen is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License or the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License and a copy of the GNU General Public License along with
// Eigen. If not, see <http://www.gnu.org/licenses/>.

#ifndef EIGEN_UMFPACKSUPPORT_H
#define EIGEN_UMFPACKSUPPORT_H

/* TODO extract L, extract U, compute det, etc... */

// generic double/complex<double> wrapper functions:

inline void umfpack_free_numeric(void **Numeric, double)
{ umfpack_di_free_numeric(Numeric); }

inline void umfpack_free_numeric(void **Numeric, std::complex<double>)
{ umfpack_zi_free_numeric(Numeric); }

inline void umfpack_free_symbolic(void **Symbolic, double)
{ umfpack_di_free_symbolic(Symbolic); }

inline void umfpack_free_symbolic(void **Symbolic, std::complex<double>)
{ umfpack_zi_free_symbolic(Symbolic); }

inline int umfpack_symbolic(int n_row,int n_col,
                            const int Ap[], const int Ai[], const double Ax[], void **Symbolic,
                            const double Control [UMFPACK_CONTROL], double Info [UMFPACK_INFO])
{
  return umfpack_di_symbolic(n_row,n_col,Ap,Ai,Ax,Symbolic,Control,Info);
}

inline int umfpack_symbolic(int n_row,int n_col,
                            const int Ap[], const int Ai[], const std::complex<double> Ax[], void **Symbolic,
                            const double Control [UMFPACK_CONTROL], double Info [UMFPACK_INFO])
{
  return umfpack_zi_symbolic(n_row,n_col,Ap,Ai,&internal::real_ref(Ax[0]),0,Symbolic,Control,Info);
}

inline int umfpack_numeric( const int Ap[], const int Ai[], const double Ax[],
                            void *Symbolic, void **Numeric,
                            const double Control[UMFPACK_CONTROL],double Info [UMFPACK_INFO])
{
  return umfpack_di_numeric(Ap,Ai,Ax,Symbolic,Numeric,Control,Info);
}

inline int umfpack_numeric( const int Ap[], const int Ai[], const std::complex<double> Ax[],
                            void *Symbolic, void **Numeric,
                            const double Control[UMFPACK_CONTROL],double Info [UMFPACK_INFO])
{
  return umfpack_zi_numeric(Ap,Ai,&internal::real_ref(Ax[0]),0,Symbolic,Numeric,Control,Info);
}

inline int umfpack_solve( int sys, const int Ap[], const int Ai[], const double Ax[],
                          double X[], const double B[], void *Numeric,
                          const double Control[UMFPACK_CONTROL], double Info[UMFPACK_INFO])
{
  return umfpack_di_solve(sys,Ap,Ai,Ax,X,B,Numeric,Control,Info);
}

inline int umfpack_solve( int sys, const int Ap[], const int Ai[], const std::complex<double> Ax[],
                          std::complex<double> X[], const std::complex<double> B[], void *Numeric,
                          const double Control[UMFPACK_CONTROL], double Info[UMFPACK_INFO])
{
  return umfpack_zi_solve(sys,Ap,Ai,&internal::real_ref(Ax[0]),0,&internal::real_ref(X[0]),0,&internal::real_ref(B[0]),0,Numeric,Control,Info);
}

inline int umfpack_get_lunz(int *lnz, int *unz, int *n_row, int *n_col, int *nz_udiag, void *Numeric, double)
{
  return umfpack_di_get_lunz(lnz,unz,n_row,n_col,nz_udiag,Numeric);
}

inline int umfpack_get_lunz(int *lnz, int *unz, int *n_row, int *n_col, int *nz_udiag, void *Numeric, std::complex<double>)
{
  return umfpack_zi_get_lunz(lnz,unz,n_row,n_col,nz_udiag,Numeric);
}

inline int umfpack_get_numeric(int Lp[], int Lj[], double Lx[], int Up[], int Ui[], double Ux[],
                               int P[], int Q[], double Dx[], int *do_recip, double Rs[], void *Numeric)
{
  return umfpack_di_get_numeric(Lp,Lj,Lx,Up,Ui,Ux,P,Q,Dx,do_recip,Rs,Numeric);
}

inline int umfpack_get_numeric(int Lp[], int Lj[], std::complex<double> Lx[], int Up[], int Ui[], std::complex<double> Ux[],
                               int P[], int Q[], std::complex<double> Dx[], int *do_recip, double Rs[], void *Numeric)
{
  double& lx0_real = internal::real_ref(Lx[0]);
  double& ux0_real = internal::real_ref(Ux[0]);
  double& dx0_real = internal::real_ref(Dx[0]);
  return umfpack_zi_get_numeric(Lp,Lj,Lx?&lx0_real:0,0,Up,Ui,Ux?&ux0_real:0,0,P,Q,
                                Dx?&dx0_real:0,0,do_recip,Rs,Numeric);
}

inline int umfpack_get_determinant(double *Mx, double *Ex, void *NumericHandle, double User_Info [UMFPACK_INFO])
{
  return umfpack_di_get_determinant(Mx,Ex,NumericHandle,User_Info);
}

inline int umfpack_get_determinant(std::complex<double> *Mx, double *Ex, void *NumericHandle, double User_Info [UMFPACK_INFO])
{
  double& mx_real = internal::real_ref(*Mx);
  return umfpack_zi_get_determinant(&mx_real,0,Ex,NumericHandle,User_Info);
}


template<typename _MatrixType>
class SparseLU<_MatrixType,UmfPack> : public SparseLU<_MatrixType>
{
  protected:
    typedef SparseLU<_MatrixType> Base;
    typedef typename Base::Scalar Scalar;
    typedef typename Base::RealScalar RealScalar;
    typedef Matrix<Scalar,Dynamic,1> Vector;
    typedef Matrix<int, 1, _MatrixType::ColsAtCompileTime> IntRowVectorType;
    typedef Matrix<int, _MatrixType::RowsAtCompileTime, 1> IntColVectorType;
    typedef SparseMatrix<Scalar,Lower|UnitDiag> LMatrixType;
    typedef SparseMatrix<Scalar,Upper> UMatrixType;
    using Base::m_flags;
    using Base::m_status;

  public:
    typedef _MatrixType MatrixType;
    typedef typename MatrixType::Index Index;

    SparseLU(int flags = NaturalOrdering)
      : Base(flags), m_numeric(0)
    {
    }

    SparseLU(const MatrixType& matrix, int flags = NaturalOrdering)
      : Base(flags), m_numeric(0)
    {
      compute(matrix);
    }

    ~SparseLU()
    {
      if (m_numeric)
        umfpack_free_numeric(&m_numeric,Scalar());
    }

    inline const LMatrixType& matrixL() const
    {
      if (m_extractedDataAreDirty) extractData();
      return m_l;
    }

    inline const UMatrixType& matrixU() const
    {
      if (m_extractedDataAreDirty) extractData();
      return m_u;
    }

    inline const IntColVectorType& permutationP() const
    {
      if (m_extractedDataAreDirty) extractData();
      return m_p;
    }

    inline const IntRowVectorType& permutationQ() const
    {
      if (m_extractedDataAreDirty) extractData();
      return m_q;
    }

    Scalar determinant() const;

    template<typename BDerived, typename XDerived>
    bool solve(const MatrixBase<BDerived> &b, MatrixBase<XDerived>* x) const;

    template<typename Rhs>
      inline const internal::solve_retval<SparseLU<MatrixType, UmfPack>, Rhs>
    solve(const MatrixBase<Rhs>& b) const
    {
      eigen_assert(true && "SparseLU is not initialized.");
      return internal::solve_retval<SparseLU<MatrixType, UmfPack>, Rhs>(*this, b.derived());
    }

    void compute(const MatrixType& matrix);

    inline Index cols() const { return m_matrixRef->cols(); }
    inline Index rows() const { return m_matrixRef->rows(); }

    inline const MatrixType& matrixLU() const
    {
      //eigen_assert(m_isInitialized && "LU is not initialized.");
      return *m_matrixRef;
    }

    const void* numeric() const
    {
      return m_numeric;
    }

  protected:

    void extractData() const;
  
  protected:
    // cached data:
    void* m_numeric;
    const MatrixType* m_matrixRef;
    mutable LMatrixType m_l;
    mutable UMatrixType m_u;
    mutable IntColVectorType m_p;
    mutable IntRowVectorType m_q;
    mutable bool m_extractedDataAreDirty;
};

namespace internal {

template<typename _MatrixType, typename Rhs>
  struct solve_retval<SparseLU<_MatrixType, UmfPack>, Rhs>
  : solve_retval_base<SparseLU<_MatrixType, UmfPack>, Rhs>
{
  typedef SparseLU<_MatrixType, UmfPack> SpLUDecType;
  EIGEN_MAKE_SOLVE_HELPERS(SpLUDecType,Rhs)

  template<typename Dest> void evalTo(Dest& dst) const
  {
    const int rhsCols = rhs().cols();

    eigen_assert((Rhs::Flags&RowMajorBit)==0 && "UmfPack backend does not support non col-major rhs yet");
    eigen_assert((Dest::Flags&RowMajorBit)==0 && "UmfPack backend does not support non col-major result yet");

    void* numeric = const_cast<void*>(dec().numeric());

    EIGEN_UNUSED int errorCode = 0;
    for (int j=0; j<rhsCols; ++j)
    {
      errorCode = umfpack_solve(UMFPACK_A,
                                dec().matrixLU()._outerIndexPtr(), dec().matrixLU()._innerIndexPtr(), dec().matrixLU()._valuePtr(),
                                &dst.col(j).coeffRef(0), &rhs().const_cast_derived().col(j).coeffRef(0), numeric, 0, 0);
      eigen_assert(!errorCode && "UmfPack could not solve the system.");
    }
  }
    
};

} // end namespace internal

template<typename MatrixType>
void SparseLU<MatrixType,UmfPack>::compute(const MatrixType& a)
{
  typedef typename MatrixType::Index Index;
  const Index rows = a.rows();
  const Index cols = a.cols();
  eigen_assert((MatrixType::Flags&RowMajorBit)==0 && "Row major matrices are not supported yet");

  m_matrixRef = &a;

  if (m_numeric)
    umfpack_free_numeric(&m_numeric,Scalar());

  void* symbolic;
  int errorCode = 0;
  errorCode = umfpack_symbolic(rows, cols, a._outerIndexPtr(), a._innerIndexPtr(), a._valuePtr(),
                                  &symbolic, 0, 0);
  if (errorCode==0)
    errorCode = umfpack_numeric(a._outerIndexPtr(), a._innerIndexPtr(), a._valuePtr(),
                                   symbolic, &m_numeric, 0, 0);

  umfpack_free_symbolic(&symbolic,Scalar());

  m_extractedDataAreDirty = true;

  Base::m_succeeded = (errorCode==0);
}

template<typename MatrixType>
void SparseLU<MatrixType,UmfPack>::extractData() const
{
  if (m_extractedDataAreDirty)
  {
    // get size of the data
    int lnz, unz, rows, cols, nz_udiag;
    umfpack_get_lunz(&lnz, &unz, &rows, &cols, &nz_udiag, m_numeric, Scalar());

    // allocate data
    m_l.resize(rows,std::min(rows,cols));
    m_l.resizeNonZeros(lnz);
    
    m_u.resize(std::min(rows,cols),cols);
    m_u.resizeNonZeros(unz);

    m_p.resize(rows);
    m_q.resize(cols);

    // extract
    umfpack_get_numeric(m_l._outerIndexPtr(), m_l._innerIndexPtr(), m_l._valuePtr(),
                        m_u._outerIndexPtr(), m_u._innerIndexPtr(), m_u._valuePtr(),
                        m_p.data(), m_q.data(), 0, 0, 0, m_numeric);
    
    m_extractedDataAreDirty = false;
  }
}

template<typename MatrixType>
typename SparseLU<MatrixType,UmfPack>::Scalar SparseLU<MatrixType,UmfPack>::determinant() const
{
  Scalar det;
  umfpack_get_determinant(&det, 0, m_numeric, 0);
  return det;
}

template<typename MatrixType>
template<typename BDerived,typename XDerived>
bool SparseLU<MatrixType,UmfPack>::solve(const MatrixBase<BDerived> &b, MatrixBase<XDerived> *x) const
{
  //const int size = m_matrix.rows();
  const int rhsCols = b.cols();
//   eigen_assert(size==b.rows());
  eigen_assert((BDerived::Flags&RowMajorBit)==0 && "UmfPack backend does not support non col-major rhs yet");
  eigen_assert((XDerived::Flags&RowMajorBit)==0 && "UmfPack backend does not support non col-major result yet");

  int errorCode;
  for (int j=0; j<rhsCols; ++j)
  {
    errorCode = umfpack_solve(UMFPACK_A,
        m_matrixRef->_outerIndexPtr(), m_matrixRef->_innerIndexPtr(), m_matrixRef->_valuePtr(),
        &x->col(j).coeffRef(0), &b.const_cast_derived().col(j).coeffRef(0), m_numeric, 0, 0);
    if (errorCode!=0)
      return false;
  }
//   errorCode = umfpack_di_solve(UMFPACK_A,
//       m_matrixRef._outerIndexPtr(), m_matrixRef._innerIndexPtr(), m_matrixRef._valuePtr(),
//       x->derived().data(), b.derived().data(), m_numeric, 0, 0);

  return true;
}

#endif // EIGEN_UMFPACKSUPPORT_H
