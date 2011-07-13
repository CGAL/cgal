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

#ifndef EIGEN_CHOLMODSUPPORT_LEGACY_H
#define EIGEN_CHOLMODSUPPORT_LEGACY_H

namespace internal {

template<typename Scalar, typename CholmodType>
void cholmod_configure_matrix_legacy(CholmodType& mat)
{
  if (internal::is_same<Scalar,float>::value)
  {
    mat.xtype = CHOLMOD_REAL;
    mat.dtype = CHOLMOD_SINGLE;
  }
  else if (internal::is_same<Scalar,double>::value)
  {
    mat.xtype = CHOLMOD_REAL;
    mat.dtype = CHOLMOD_DOUBLE;
  }
  else if (internal::is_same<Scalar,std::complex<float> >::value)
  {
    mat.xtype = CHOLMOD_COMPLEX;
    mat.dtype = CHOLMOD_SINGLE;
  }
  else if (internal::is_same<Scalar,std::complex<double> >::value)
  {
    mat.xtype = CHOLMOD_COMPLEX;
    mat.dtype = CHOLMOD_DOUBLE;
  }
  else
  {
    eigen_assert(false && "Scalar type not supported by CHOLMOD");
  }
}

template<typename _MatrixType>
cholmod_sparse cholmod_map_eigen_to_sparse(_MatrixType& mat)
{
  typedef typename _MatrixType::Scalar Scalar;
  cholmod_sparse res;
  res.nzmax   = mat.nonZeros();
  res.nrow    = mat.rows();;
  res.ncol    = mat.cols();
  res.p       = mat._outerIndexPtr();
  res.i       = mat._innerIndexPtr();
  res.x       = mat._valuePtr();
  res.xtype   = CHOLMOD_REAL;
  res.itype   = CHOLMOD_INT;
  res.sorted  = 1;
  res.packed  = 1;
  res.dtype   = 0;
  res.stype   = -1;

  internal::cholmod_configure_matrix_legacy<Scalar>(res);


  if (_MatrixType::Flags & SelfAdjoint)
  {
    if (_MatrixType::Flags & Upper)
      res.stype = 1;
    else if (_MatrixType::Flags & Lower)
      res.stype = -1;
    else
      res.stype = 0;
  }
  else
    res.stype = -1; // by default we consider the lower part

  return res;
}

template<typename Derived>
cholmod_dense cholmod_map_eigen_to_dense(MatrixBase<Derived>& mat)
{
  EIGEN_STATIC_ASSERT((internal::traits<Derived>::Flags&RowMajorBit)==0,THIS_METHOD_IS_ONLY_FOR_COLUMN_MAJOR_MATRICES);
  typedef typename Derived::Scalar Scalar;

  cholmod_dense res;
  res.nrow   = mat.rows();
  res.ncol   = mat.cols();
  res.nzmax  = res.nrow * res.ncol;
  res.d      = Derived::IsVectorAtCompileTime ? mat.derived().size() : mat.derived().outerStride();
  res.x      = mat.derived().data();
  res.z      = 0;

  internal::cholmod_configure_matrix_legacy<Scalar>(res);

  return res;
}

template<typename Scalar, int Flags, typename Index>
MappedSparseMatrix<Scalar,Flags,Index> map_cholmod_sparse_to_eigen(cholmod_sparse& cm)
{
  return MappedSparseMatrix<Scalar,Flags,Index>
         (cm.nrow, cm.ncol, reinterpret_cast<Index*>(cm.p)[cm.ncol],
          reinterpret_cast<Index*>(cm.p), reinterpret_cast<Index*>(cm.i),reinterpret_cast<Scalar*>(cm.x) );
}

} // namespace internal

template<typename _MatrixType>
class SparseLLT<_MatrixType, Cholmod> : public SparseLLT<_MatrixType>
{
  protected:
    typedef SparseLLT<_MatrixType> Base;
    typedef typename Base::Scalar Scalar;
    typedef typename Base::RealScalar RealScalar;
    typedef typename Base::CholMatrixType CholMatrixType;
    using Base::MatrixLIsDirty;
    using Base::SupernodalFactorIsDirty;
    using Base::m_flags;
    using Base::m_matrix;
    using Base::m_status;

  public:
    typedef _MatrixType MatrixType;
    typedef typename MatrixType::Index Index;

    SparseLLT(int flags = 0)
      : Base(flags), m_cholmodFactor(0)
    {
      cholmod_start(&m_cholmod);
    }

    SparseLLT(const MatrixType& matrix, int flags = 0)
      : Base(flags), m_cholmodFactor(0)
    {
      cholmod_start(&m_cholmod);
      compute(matrix);
    }

    ~SparseLLT()
    {
      if (m_cholmodFactor)
        cholmod_free_factor(&m_cholmodFactor, &m_cholmod);
      cholmod_finish(&m_cholmod);
    }

    inline const CholMatrixType& matrixL() const;

    template<typename Derived>
    bool solveInPlace(MatrixBase<Derived> &b) const;

    template<typename Rhs>
    inline const internal::solve_retval<SparseLLT<MatrixType, Cholmod>, Rhs>
    solve(const MatrixBase<Rhs>& b) const
    {
      eigen_assert(true && "SparseLLT is not initialized.");
      return internal::solve_retval<SparseLLT<MatrixType, Cholmod>, Rhs>(*this, b.derived());
    }

    void compute(const MatrixType& matrix);

    inline Index cols() const { return m_matrix.cols(); }
    inline Index rows() const { return m_matrix.rows(); }

    inline const cholmod_factor* cholmodFactor() const
    { return m_cholmodFactor; }

    inline cholmod_common* cholmodCommon() const
    { return &m_cholmod; }

    bool succeeded() const;

  protected:
    mutable cholmod_common m_cholmod;
    cholmod_factor* m_cholmodFactor;
};


namespace internal {

template<typename _MatrixType, typename Rhs>
  struct solve_retval<SparseLLT<_MatrixType, Cholmod>, Rhs>
  : solve_retval_base<SparseLLT<_MatrixType, Cholmod>, Rhs>
{
  typedef SparseLLT<_MatrixType, Cholmod> SpLLTDecType;
  EIGEN_MAKE_SOLVE_HELPERS(SpLLTDecType,Rhs)

  template<typename Dest> void evalTo(Dest& dst) const
  {
    //Index size = dec().cholmodFactor()->n;
    eigen_assert((Index)dec().cholmodFactor()->n==rhs().rows());
    
    cholmod_factor* cholmodFactor = const_cast<cholmod_factor*>(dec().cholmodFactor());
    cholmod_common* cholmodCommon = const_cast<cholmod_common*>(dec().cholmodCommon());
    // this uses Eigen's triangular sparse solver
    // if (m_status & MatrixLIsDirty)
    //   matrixL();
    // Base::solveInPlace(b);
    // as long as our own triangular sparse solver is not fully optimal,
    // let's use CHOLMOD's one:
    cholmod_dense cdb = internal::cholmod_map_eigen_to_dense(rhs().const_cast_derived());
    cholmod_dense* x = cholmod_solve(CHOLMOD_A, cholmodFactor, &cdb, cholmodCommon);

    dst = Matrix<typename Base::Scalar,Dynamic,1>::Map(reinterpret_cast<typename Base::Scalar*>(x->x), rhs().rows());  

    cholmod_free_dense(&x, cholmodCommon);

  }
    
};

} // namespace internal



template<typename _MatrixType>
void SparseLLT<_MatrixType,Cholmod>::compute(const _MatrixType& a)
{
  if (m_cholmodFactor)
  {
    cholmod_free_factor(&m_cholmodFactor, &m_cholmod);
    m_cholmodFactor = 0;
  }

  cholmod_sparse A = internal::cholmod_map_eigen_to_sparse(const_cast<_MatrixType&>(a));
//   m_cholmod.supernodal = CHOLMOD_AUTO;
  // TODO
//   if (m_flags&IncompleteFactorization)
//   {
//     m_cholmod.nmethods = 1;
//     m_cholmod.method[0].ordering = CHOLMOD_NATURAL;
//     m_cholmod.postorder = 0;
//   }
//   else
//   {
//     m_cholmod.nmethods = 1;
//     m_cholmod.method[0].ordering = CHOLMOD_NATURAL;
//     m_cholmod.postorder = 0;
//   }
//   m_cholmod.final_ll = 1;
  m_cholmodFactor = cholmod_analyze(&A, &m_cholmod);
  cholmod_factorize(&A, m_cholmodFactor, &m_cholmod);

  this->m_status = (this->m_status & ~Base::SupernodalFactorIsDirty) | Base::MatrixLIsDirty;
}


// TODO
template<typename _MatrixType>
bool SparseLLT<_MatrixType,Cholmod>::succeeded() const
{ return true; }



template<typename _MatrixType>
inline const typename SparseLLT<_MatrixType,Cholmod>::CholMatrixType&
SparseLLT<_MatrixType,Cholmod>::matrixL() const
{
  if (this->m_status & Base::MatrixLIsDirty)
  {
    eigen_assert(!(this->m_status & Base::SupernodalFactorIsDirty));

    cholmod_sparse* cmRes = cholmod_factor_to_sparse(m_cholmodFactor, &m_cholmod);
    const_cast<typename Base::CholMatrixType&>(this->m_matrix) = 
      internal::map_cholmod_sparse_to_eigen<Scalar,ColMajor,Index>(*cmRes);
    free(cmRes);

    this->m_status = (this->m_status & ~Base::MatrixLIsDirty);
  }
  return this->m_matrix;
}




template<typename _MatrixType>
template<typename Derived>
bool SparseLLT<_MatrixType,Cholmod>::solveInPlace(MatrixBase<Derived> &b) const
{
  //Index size = m_cholmodFactor->n;
  eigen_assert((Index)m_cholmodFactor->n==b.rows());

  // this uses Eigen's triangular sparse solver
  //   if (m_status & MatrixLIsDirty)
  //     matrixL();
  //   Base::solveInPlace(b);
  // as long as our own triangular sparse solver is not fully optimal,
  // let's use CHOLMOD's one:
  cholmod_dense cdb = internal::cholmod_map_eigen_to_dense(b);

  cholmod_dense* x = cholmod_solve(CHOLMOD_A, m_cholmodFactor, &cdb, &m_cholmod);
  eigen_assert(x && "Eigen: cholmod_solve failed.");

  b = Matrix<typename Base::Scalar,Dynamic,1>::Map(reinterpret_cast<typename Base::Scalar*>(x->x),b.rows());
  cholmod_free_dense(&x, &m_cholmod);
  return true;
}











template<typename _MatrixType>
class SparseLDLT<_MatrixType,Cholmod> : public SparseLDLT<_MatrixType>
{
  protected:
    typedef SparseLDLT<_MatrixType> Base;
    typedef typename Base::Scalar Scalar;
    typedef typename Base::RealScalar RealScalar;
    using Base::MatrixLIsDirty;
    using Base::SupernodalFactorIsDirty;
    using Base::m_flags;
    using Base::m_matrix;
    using Base::m_status;

  public:
    typedef _MatrixType MatrixType;
    typedef typename MatrixType::Index Index;

    SparseLDLT(int flags = 0)
      : Base(flags), m_cholmodFactor(0)
    {
      cholmod_start(&m_cholmod);
    }

    SparseLDLT(const _MatrixType& matrix, int flags = 0)
      : Base(flags), m_cholmodFactor(0)
    {
      cholmod_start(&m_cholmod);
      compute(matrix);
    }

    ~SparseLDLT()
    {
      if (m_cholmodFactor)
        cholmod_free_factor(&m_cholmodFactor, &m_cholmod);
      cholmod_finish(&m_cholmod);
    }

    inline const typename Base::CholMatrixType& matrixL(void) const;

    template<typename Derived>
    void solveInPlace(MatrixBase<Derived> &b) const;

    template<typename Rhs>
    inline const internal::solve_retval<SparseLDLT<MatrixType, Cholmod>, Rhs>
    solve(const MatrixBase<Rhs>& b) const
    {
      eigen_assert(true && "SparseLDLT is not initialized.");
      return internal::solve_retval<SparseLDLT<MatrixType, Cholmod>, Rhs>(*this, b.derived());
    }

    void compute(const _MatrixType& matrix);

    inline Index cols() const { return m_matrix.cols(); }
    inline Index rows() const { return m_matrix.rows(); }

    inline const cholmod_factor* cholmodFactor() const
    { return m_cholmodFactor; }

    inline cholmod_common* cholmodCommon() const
    { return &m_cholmod; }

    bool succeeded() const;

  protected:
    mutable cholmod_common m_cholmod;
    cholmod_factor* m_cholmodFactor;
};



namespace internal {

template<typename _MatrixType, typename Rhs>
  struct solve_retval<SparseLDLT<_MatrixType, Cholmod>, Rhs>
  : solve_retval_base<SparseLDLT<_MatrixType, Cholmod>, Rhs>
{
  typedef SparseLDLT<_MatrixType, Cholmod> SpLDLTDecType;
  EIGEN_MAKE_SOLVE_HELPERS(SpLDLTDecType,Rhs)

  template<typename Dest> void evalTo(Dest& dst) const
  {
    //Index size = dec().cholmodFactor()->n;
    eigen_assert((Index)dec().cholmodFactor()->n==rhs().rows());
    
    cholmod_factor* cholmodFactor = const_cast<cholmod_factor*>(dec().cholmodFactor());
    cholmod_common* cholmodCommon = const_cast<cholmod_common*>(dec().cholmodCommon());
    // this uses Eigen's triangular sparse solver
    // if (m_status & MatrixLIsDirty)
    //   matrixL();
    // Base::solveInPlace(b);
    // as long as our own triangular sparse solver is not fully optimal,
    // let's use CHOLMOD's one:
    cholmod_dense cdb = internal::cholmod_map_eigen_to_dense(rhs().const_cast_derived());
    cholmod_dense* x = cholmod_solve(CHOLMOD_LDLt, cholmodFactor, &cdb, cholmodCommon);

    dst = Matrix<typename Base::Scalar,Dynamic,1>::Map(reinterpret_cast<typename Base::Scalar*>(x->x), rhs().rows());  
    cholmod_free_dense(&x, cholmodCommon);

  }
    
};


} // namespace internal

template<typename _MatrixType>
void SparseLDLT<_MatrixType,Cholmod>::compute(const _MatrixType& a)
{
  if (m_cholmodFactor)
  {
    cholmod_free_factor(&m_cholmodFactor, &m_cholmod);
    m_cholmodFactor = 0;
  }

  cholmod_sparse A = internal::cholmod_map_eigen_to_sparse(const_cast<_MatrixType&>(a));
 
  //m_cholmod.supernodal = CHOLMOD_AUTO;
  m_cholmod.supernodal = CHOLMOD_SIMPLICIAL;
  //m_cholmod.supernodal = CHOLMOD_SUPERNODAL;
  // TODO
  if (this->m_flags & IncompleteFactorization)
  {
    m_cholmod.nmethods = 1;
    //m_cholmod.method[0].ordering = CHOLMOD_NATURAL;
    m_cholmod.method[0].ordering = CHOLMOD_COLAMD;
    m_cholmod.postorder = 1;
  }
  else
  {
    m_cholmod.nmethods = 1;
    m_cholmod.method[0].ordering = CHOLMOD_NATURAL;
    m_cholmod.postorder = 0;
  }
  m_cholmod.final_ll = 0;
  m_cholmodFactor = cholmod_analyze(&A, &m_cholmod);
  cholmod_factorize(&A, m_cholmodFactor, &m_cholmod);
  
  this->m_status = (this->m_status & ~Base::SupernodalFactorIsDirty) | Base::MatrixLIsDirty;
}


// TODO
template<typename _MatrixType>
bool SparseLDLT<_MatrixType,Cholmod>::succeeded() const
{ return true; }


template<typename _MatrixType>
inline const typename SparseLDLT<_MatrixType>::CholMatrixType&
SparseLDLT<_MatrixType,Cholmod>::matrixL() const
{
  if (this->m_status & Base::MatrixLIsDirty)
  {
    eigen_assert(!(this->m_status & Base::SupernodalFactorIsDirty));

    cholmod_sparse* cmRes = cholmod_factor_to_sparse(m_cholmodFactor, &m_cholmod);
    const_cast<typename Base::CholMatrixType&>(this->m_matrix) = MappedSparseMatrix<Scalar>(*cmRes);
    free(cmRes);

    this->m_status = (this->m_status & ~Base::MatrixLIsDirty);
  }
  return this->m_matrix;
}






template<typename _MatrixType>
template<typename Derived>
void SparseLDLT<_MatrixType,Cholmod>::solveInPlace(MatrixBase<Derived> &b) const
{
  //Index size = m_cholmodFactor->n;
  eigen_assert((Index)m_cholmodFactor->n == b.rows());

  // this uses Eigen's triangular sparse solver
  //   if (m_status & MatrixLIsDirty)
  //     matrixL();
  //   Base::solveInPlace(b);
  // as long as our own triangular sparse solver is not fully optimal,
  // let's use CHOLMOD's one:
  cholmod_dense cdb = internal::cholmod_map_eigen_to_dense(b);
  cholmod_dense* x = cholmod_solve(CHOLMOD_A, m_cholmodFactor, &cdb, &m_cholmod);
  b = Matrix<typename Base::Scalar,Dynamic,1>::Map(reinterpret_cast<typename Base::Scalar*>(x->x),b.rows());
  cholmod_free_dense(&x, &m_cholmod);
}






#endif // EIGEN_CHOLMODSUPPORT_LEGACY_H
