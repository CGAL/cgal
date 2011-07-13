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

/*

NOTE: the _symbolic, and _numeric functions has been adapted from
      the LDL library:

LDL Copyright (c) 2005 by Timothy A. Davis.  All Rights Reserved.

LDL License:

    Your use or distribution of LDL or any modified version of
    LDL implies that you agree to this License.

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; either
    version 2.1 of the License, or (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public
    License along with this library; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301
    USA

    Permission is hereby granted to use or copy this program under the
    terms of the GNU LGPL, provided that the Copyright, this License,
    and the Availability of the original version is retained on all copies.
    User documentation of any code that uses this code or any modified
    version of this code must cite the Copyright, this License, the
    Availability note, and "Used by permission." Permission to modify
    the code and to distribute modified code is granted, provided the
    Copyright, this License, and the Availability note are retained,
    and a notice that the code was modified is included.
 */

#ifndef EIGEN_SPARSELDLT_LEGACY_H
#define EIGEN_SPARSELDLT_LEGACY_H

/** \ingroup Sparse_Module
  *
  * \class SparseLDLT
  *
  * \brief LDLT Cholesky decomposition of a sparse matrix and associated features
  *
  * \param MatrixType the type of the matrix of which we are computing the LDLT Cholesky decomposition
  *
  * \warning the upper triangular part has to be specified. The rest of the matrix is not used. The input matrix must be column major.
  *
  * \sa class LDLT, class LDLT
  */
template<typename _MatrixType, typename Backend = DefaultBackend>
class SparseLDLT
{
  protected:
    typedef typename _MatrixType::Scalar Scalar;
    typedef typename NumTraits<typename _MatrixType::Scalar>::Real RealScalar;
    
    typedef Matrix<Scalar,_MatrixType::ColsAtCompileTime,1> VectorType;

    enum {
      SupernodalFactorIsDirty      = 0x10000,
      MatrixLIsDirty               = 0x20000
    };

  public:
    typedef SparseMatrix<Scalar> CholMatrixType;
    typedef _MatrixType MatrixType;
    typedef typename MatrixType::Index Index;


    /** Creates a dummy LDLT factorization object with flags \a flags. */
    SparseLDLT(int flags = 0)
      : m_flags(flags), m_status(0)
    {
      eigen_assert((MatrixType::Flags&RowMajorBit)==0);
      m_precision = RealScalar(0.1) * Eigen::NumTraits<RealScalar>::dummy_precision();
    }

    /** Creates a LDLT object and compute the respective factorization of \a matrix using
      * flags \a flags. */
    SparseLDLT(const MatrixType& matrix, int flags = 0)
      : m_matrix(matrix.rows(), matrix.cols()), m_flags(flags), m_status(0)
    {
      eigen_assert((MatrixType::Flags&RowMajorBit)==0);
      m_precision = RealScalar(0.1) * Eigen::NumTraits<RealScalar>::dummy_precision();
      compute(matrix);
    }

    /** Sets the relative threshold value used to prune zero coefficients during the decomposition.
      *
      * Setting a value greater than zero speeds up computation, and yields to an imcomplete
      * factorization with fewer non zero coefficients. Such approximate factors are especially
      * useful to initialize an iterative solver.
      *
      * \warning if precision is greater that zero, the LDLT factorization is not guaranteed to succeed
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
    void settags(int f) { m_flags = f; }
    /** \returns the current flags */
    int flags() const { return m_flags; }

    /** Computes/re-computes the LDLT factorization */
    void compute(const MatrixType& matrix);

    /** Perform a symbolic factorization */
    void _symbolic(const MatrixType& matrix);
    /** Perform the actual factorization using the previously
      * computed symbolic factorization */
    bool _numeric(const MatrixType& matrix);

    /** \returns the lower triangular matrix L */
    inline const CholMatrixType& matrixL(void) const { return m_matrix; }

    /** \returns the coefficients of the diagonal matrix D */
    inline VectorType vectorD(void) const { return m_diag; }

    template<typename Derived>
    bool solveInPlace(MatrixBase<Derived> &b) const;

    template<typename Rhs>
    inline const internal::solve_retval<SparseLDLT<MatrixType>, Rhs>
    solve(const MatrixBase<Rhs>& b) const
    {
      eigen_assert(true && "SparseLDLT is not initialized.");
      return internal::solve_retval<SparseLDLT<MatrixType>, Rhs>(*this, b.derived());
    }

    inline Index cols() const { return m_matrix.cols(); }
    inline Index rows() const { return m_matrix.rows(); }

    inline const VectorType& diag() const { return m_diag; }

    /** \returns true if the factorization succeeded */
    inline bool succeeded(void) const { return m_succeeded; }

  protected:
    CholMatrixType m_matrix;
    VectorType m_diag;
    VectorXi m_parent; // elimination tree
    VectorXi m_nonZerosPerCol;
//     VectorXi m_w; // workspace
    PermutationMatrix<Dynamic> m_P;
    PermutationMatrix<Dynamic> m_Pinv;
    RealScalar m_precision;
    int m_flags;
    mutable int m_status;
    bool m_succeeded;
};

namespace internal {

template<typename _MatrixType, typename Rhs>
struct solve_retval<SparseLDLT<_MatrixType>, Rhs>
  : solve_retval_base<SparseLDLT<_MatrixType>, Rhs>
{
  typedef SparseLDLT<_MatrixType> SpLDLTDecType;
  EIGEN_MAKE_SOLVE_HELPERS(SpLDLTDecType,Rhs)

  template<typename Dest> void evalTo(Dest& dst) const
  {
    //Index size = dec().matrixL().rows();
    eigen_assert(dec().matrixL().rows()==rhs().rows());

    Rhs b(rhs().rows(), rhs().cols());
    b = rhs();

    if (dec().matrixL().nonZeros()>0) // otherwise L==I
      dec().matrixL().template triangularView<UnitLower>().solveInPlace(b);

    b = b.cwiseQuotient(dec().diag());
    if (dec().matrixL().nonZeros()>0) // otherwise L==I
      dec().matrixL().adjoint().template triangularView<UnitUpper>().solveInPlace(b);
    
    dst = b;

  }
    
};

} // end namespace internal

/** Computes / recomputes the LDLT decomposition of matrix \a a
  * using the default algorithm.
  */
template<typename _MatrixType, typename Backend>
void SparseLDLT<_MatrixType,Backend>::compute(const _MatrixType& a)
{
  _symbolic(a);
  m_succeeded = _numeric(a);
}

template<typename _MatrixType, typename Backend>
void SparseLDLT<_MatrixType,Backend>::_symbolic(const _MatrixType& a)
{
  assert(a.rows()==a.cols());
  const Index size = a.rows();
  m_matrix.resize(size, size);
  m_parent.resize(size);
  m_nonZerosPerCol.resize(size);
  Index * tags = ei_aligned_stack_new(Index, size);

  const Index* Ap = a._outerIndexPtr();
  const Index* Ai = a._innerIndexPtr();
  Index* Lp = m_matrix._outerIndexPtr();
  
  const Index* P = 0;
  Index* Pinv = 0;

  if(P)
  {
    m_P.indices()     = VectorXi::Map(P,size);
    m_Pinv = m_P.inverse();
    Pinv = m_Pinv.indices().data();
  }
  else
  {
    m_P.resize(0);
    m_Pinv.resize(0);
  }

  for (Index k = 0; k < size; ++k)
  {
    /* L(k,:) pattern: all nodes reachable in etree from nz in A(0:k-1,k) */
    m_parent[k] = -1;             /* parent of k is not yet known */
    tags[k] = k;                  /* mark node k as visited */
    m_nonZerosPerCol[k] = 0;      /* count of nonzeros in column k of L */
    Index kk = P ? P[k] : k;      /* kth original, or permuted, column */
    Index p2 = Ap[kk+1];
    for (Index p = Ap[kk]; p < p2; ++p)
    {
      /* A (i,k) is nonzero (original or permuted A) */
      Index i = Pinv ? Pinv[Ai[p]] : Ai[p];
      if (i < k)
      {
        /* follow path from i to root of etree, stop at flagged node */
        for (; tags[i] != k; i = m_parent[i])
        {
          /* find parent of i if not yet determined */
          if (m_parent[i] == -1)
            m_parent[i] = k;
          ++m_nonZerosPerCol[i];        /* L (k,i) is nonzero */
          tags[i] = k;                  /* mark i as visited */
        }
      }
    }
  }
  /* construct Lp index array from m_nonZerosPerCol column counts */
  Lp[0] = 0;
  for (Index k = 0; k < size; ++k)
    Lp[k+1] = Lp[k] + m_nonZerosPerCol[k];

  m_matrix.resizeNonZeros(Lp[size]);
  ei_aligned_stack_delete(Index, tags, size);
}

template<typename _MatrixType, typename Backend>
bool SparseLDLT<_MatrixType,Backend>::_numeric(const _MatrixType& a)
{
  assert(a.rows()==a.cols());
  const Index size = a.rows();
  assert(m_parent.size()==size);
  assert(m_nonZerosPerCol.size()==size);

  const Index* Ap = a._outerIndexPtr();
  const Index* Ai = a._innerIndexPtr();
  const Scalar* Ax = a._valuePtr();
  const Index* Lp = m_matrix._outerIndexPtr();
  Index* Li = m_matrix._innerIndexPtr();
  Scalar* Lx = m_matrix._valuePtr();
  m_diag.resize(size);

  Scalar * y = ei_aligned_stack_new(Scalar, size);
  Index * pattern = ei_aligned_stack_new(Index, size);
  Index * tags = ei_aligned_stack_new(Index, size);
  
  Index* P = 0;
  Index* Pinv = 0;
  
  if(m_P.size()==size)
  {
    P = m_P.indices().data();
    Pinv = m_Pinv.indices().data();
  }
  
  bool ok = true;

  for (Index k = 0; k < size; ++k)
  {
    /* compute nonzero pattern of kth row of L, in topological order */
    y[k] = 0.0;                     /* Y(0:k) is now all zero */
    Index top = size;               /* stack for pattern is empty */
    tags[k] = k;                    /* mark node k as visited */
    m_nonZerosPerCol[k] = 0;        /* count of nonzeros in column k of L */
    Index kk = (P) ? (P[k]) : (k);  /* kth original, or permuted, column */
    Index p2 = Ap[kk+1];
    for (Index p = Ap[kk]; p < p2; ++p)
    {
      Index i = Pinv ? Pinv[Ai[p]] : Ai[p]; /* get A(i,k) */
      if (i <= k)
      {
        y[i] += internal::conj(Ax[p]);            /* scatter A(i,k) into Y (sum duplicates) */
        Index len;
        for (len = 0; tags[i] != k; i = m_parent[i])
        {
          pattern[len++] = i;     /* L(k,i) is nonzero */
          tags[i] = k;            /* mark i as visited */
        }
        while (len > 0)
          pattern[--top] = pattern[--len];
      }
    }

    /* compute numerical values kth row of L (a sparse triangular solve) */
    m_diag[k] = y[k];                       /* get D(k,k) and clear Y(k) */
    y[k] = 0.0;
    for (; top < size; ++top)
    {
      Index i = pattern[top];       /* pattern[top:n-1] is pattern of L(:,k) */
      Scalar yi = (y[i]);             /* get and clear Y(i) */
      y[i] = 0.0;
      Index p2 = Lp[i] + m_nonZerosPerCol[i];
      Index p;
      for (p = Lp[i]; p < p2; ++p)
        y[Li[p]] -= internal::conj(Lx[p]) * (yi);
      Scalar l_ki = yi / m_diag[i];       /* the nonzero entry L(k,i) */
      m_diag[k] -= l_ki * internal::conj(yi);
      Li[p] = k;                          /* store L(k,i) in column form of L */
      Lx[p] = (l_ki);
      ++m_nonZerosPerCol[i];              /* increment count of nonzeros in col i */
    }
    if (m_diag[k] == 0.0)
    {
      ok = false;                         /* failure, D(k,k) is zero */
      break;
    }
  }

  ei_aligned_stack_delete(Scalar, y, size);
  ei_aligned_stack_delete(Index, pattern, size);
  ei_aligned_stack_delete(Index, tags, size);

  return ok;  /* success, diagonal of D is all nonzero */
}

/** Computes b = L^-T D^-1 L^-1 b */
template<typename _MatrixType, typename Backend>
template<typename Derived>
bool SparseLDLT<_MatrixType, Backend>::solveInPlace(MatrixBase<Derived> &b) const
{
  //Index size = m_matrix.rows();
  eigen_assert(m_matrix.rows()==b.rows());
  if (!m_succeeded)
    return false;
  
  if(m_P.size()>0)
    b = m_Pinv * b;

  if (m_matrix.nonZeros()>0) // otherwise L==I
    m_matrix.template triangularView<UnitLower>().solveInPlace(b);
  b = b.cwiseQuotient(m_diag);
  if (m_matrix.nonZeros()>0) // otherwise L==I
    m_matrix.adjoint().template triangularView<UnitUpper>().solveInPlace(b);

  if(m_P.size()>0)
    b = m_P * b;
  
  return true;
}

#endif // EIGEN_SPARSELDLT_LEGACY_H
