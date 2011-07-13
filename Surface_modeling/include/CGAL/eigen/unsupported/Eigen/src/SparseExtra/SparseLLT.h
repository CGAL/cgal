// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2008 Gael Guennebaud <gael.guennebaud@inria.fr>
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

#ifndef EIGEN_SPARSELLT_H
#define EIGEN_SPARSELLT_H

/** \ingroup Sparse_Module
  *
  * \class SparseLLT
  *
  * \brief LLT Cholesky decomposition of a sparse matrix and associated features
  *
  * \param MatrixType the type of the matrix of which we are computing the LLT Cholesky decomposition
  *
  * \sa class LLT, class LDLT
  */
template<typename _MatrixType, typename Backend = DefaultBackend>
class SparseLLT
{
  protected:
    typedef typename _MatrixType::Scalar Scalar;
    typedef typename NumTraits<typename _MatrixType::Scalar>::Real RealScalar;

    enum {
      SupernodalFactorIsDirty      = 0x10000,
      MatrixLIsDirty               = 0x20000
    };

  public:
    typedef SparseMatrix<Scalar> CholMatrixType;
    typedef _MatrixType MatrixType;
    typedef typename MatrixType::Index Index;

    /** Creates a dummy LLT factorization object with flags \a flags. */
    SparseLLT(int flags = 0)
      : m_flags(flags), m_status(0)
    {
      m_precision = RealScalar(0.1) * Eigen::NumTraits<RealScalar>::dummy_precision();
    }

    /** Creates a LLT object and compute the respective factorization of \a matrix using
      * flags \a flags. */
    SparseLLT(const MatrixType& matrix, int flags = 0)
      : m_matrix(matrix.rows(), matrix.cols()), m_flags(flags), m_status(0)
    {
      m_precision = RealScalar(0.1) * Eigen::NumTraits<RealScalar>::dummy_precision();
      compute(matrix);
    }

    /** Sets the relative threshold value used to prune zero coefficients during the decomposition.
      *
      * Setting a value greater than zero speeds up computation, and yields to an imcomplete
      * factorization with fewer non zero coefficients. Such approximate factors are especially
      * useful to initialize an iterative solver.
      *
      * \warning if precision is greater that zero, the LLT factorization is not guaranteed to succeed
      * even if the matrix is positive definite.
      *
      * Note that the exact meaning of this parameter might depends on the actual
      * backend. Moreover, not all backends support this feature.
      *
      * \sa precision() */
    void setPrecision(RealScalar v) { m_precision = v; }

    /** \returns the current precision.
      *
      * \sa setPrecision() */
    RealScalar precision() const { return m_precision; }

    /** Sets the flags. Possible values are:
      *  - CompleteFactorization
      *  - IncompleteFactorization
      *  - MemoryEfficient          (hint to use the memory most efficient method offered by the backend)
      *  - SupernodalMultifrontal   (implies a complete factorization if supported by the backend,
      *                              overloads the MemoryEfficient flags)
      *  - SupernodalLeftLooking    (implies a complete factorization  if supported by the backend,
      *                              overloads the MemoryEfficient flags)
      *
      * \sa flags() */
    void setFlags(int f) { m_flags = f; }
    /** \returns the current flags */
    int flags() const { return m_flags; }

    /** Computes/re-computes the LLT factorization */
    void compute(const MatrixType& matrix);

    /** \returns the lower triangular matrix L */
    inline const CholMatrixType& matrixL(void) const { return m_matrix; }

    template<typename Derived>
    bool solveInPlace(MatrixBase<Derived> &b) const;

    template<typename Rhs>
    inline const internal::solve_retval<SparseLLT<MatrixType>, Rhs>
    solve(const MatrixBase<Rhs>& b) const
    {
      eigen_assert(true && "SparseLLT is not initialized.");
      return internal::solve_retval<SparseLLT<MatrixType>, Rhs>(*this, b.derived());
    }

    inline Index cols() const { return m_matrix.cols(); }
    inline Index rows() const { return m_matrix.rows(); }

    /** \returns true if the factorization succeeded */
    inline bool succeeded(void) const { return m_succeeded; }

  protected:
    CholMatrixType m_matrix;
    RealScalar m_precision;
    int m_flags;
    mutable int m_status;
    bool m_succeeded;
};


namespace internal {

template<typename _MatrixType, typename Rhs>
struct solve_retval<SparseLLT<_MatrixType>, Rhs>
  : solve_retval_base<SparseLLT<_MatrixType>, Rhs>
{
  typedef SparseLLT<_MatrixType> SpLLTDecType;
  EIGEN_MAKE_SOLVE_HELPERS(SpLLTDecType,Rhs)

  template<typename Dest> void evalTo(Dest& dst) const
  {
    const Index size = dec().matrixL().rows();
    eigen_assert(size==rhs().rows());
    
    Rhs b(rhs().rows(), rhs().cols());
    b = rhs();
    
    dec().matrixL().template triangularView<Lower>().solveInPlace(b);
    dec().matrixL().adjoint().template triangularView<Upper>().solveInPlace(b);
    
    dst = b;
    
  }
    
};

} // end namespace internal


/** Computes / recomputes the LLT decomposition of matrix \a a
  * using the default algorithm.
  */
template<typename _MatrixType, typename Backend>
void SparseLLT<_MatrixType,Backend>::compute(const _MatrixType& a)
{
  assert(a.rows()==a.cols());
  const Index size = a.rows();
  m_matrix.resize(size, size);

  // allocate a temporary vector for accumulations
  AmbiVector<Scalar,Index> tempVector(size);
  RealScalar density = a.nonZeros()/RealScalar(size*size);

  // TODO estimate the number of non zeros
  m_matrix.setZero();
  m_matrix.reserve(a.nonZeros()*2);
  for (Index j = 0; j < size; ++j)
  {
    Scalar x = internal::real(a.coeff(j,j));

    // TODO better estimate of the density !
    tempVector.init(density>0.001? IsDense : IsSparse);
    tempVector.setBounds(j+1,size);
    tempVector.setZero();
    // init with current matrix a
    {
      typename _MatrixType::InnerIterator it(a,j);
      eigen_assert(it.index()==j &&
        "matrix must has non zero diagonal entries and only the lower triangular part must be stored");
      ++it; // skip diagonal element
      for (; it; ++it)
        tempVector.coeffRef(it.index()) = it.value();
    }
    for (Index k=0; k<j+1; ++k)
    {
      typename CholMatrixType::InnerIterator it(m_matrix, k);
      while (it && it.index()<j)
        ++it;
      if (it && it.index()==j)
      {
        Scalar y = it.value();
        x -= internal::abs2(y);
        ++it; // skip j-th element, and process remaining column coefficients
        tempVector.restart();
        for (; it; ++it)
        {
          tempVector.coeffRef(it.index()) -= it.value() * y;
        }
      }
    }
    // copy the temporary vector to the respective m_matrix.col()
    // while scaling the result by 1/real(x)
    RealScalar rx = internal::sqrt(internal::real(x));
    m_matrix.insert(j,j) = rx; // FIXME use insertBack
    Scalar y = Scalar(1)/rx;
    for (typename AmbiVector<Scalar,Index>::Iterator it(tempVector, m_precision*rx); it; ++it)
    {
      // FIXME use insertBack
      m_matrix.insert(it.index(), j) = it.value() * y;
    }
  }
  m_matrix.finalize();
}

/** Computes b = L^-T L^-1 b */
template<typename _MatrixType, typename Backend>
template<typename Derived>
bool SparseLLT<_MatrixType, Backend>::solveInPlace(MatrixBase<Derived> &b) const
{
  const Index size = m_matrix.rows();
  eigen_assert(size==b.rows());

  m_matrix.template triangularView<Lower>().solveInPlace(b);
  m_matrix.adjoint().template triangularView<Upper>().solveInPlace(b);

  return true;
}

#endif // EIGEN_SPARSELLT_H
