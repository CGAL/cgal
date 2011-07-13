// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2008-2010 Gael Guennebaud <gael.guennebaud@inria.fr>
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

#ifndef EIGEN_SIMPLICIAL_CHOLESKY_H
#define EIGEN_SIMPLICIAL_CHOLESKY_H

enum SimplicialCholeskyMode {
  SimplicialCholeskyLLt,
  SimplicialCholeskyLDLt
};

/** \brief A direct sparse Cholesky factorization
  *
  * This class allows to solve for A.X = B sparse linear problems via a LL^T or LDL^T Cholesky factorization.
  * The sparse matrix A must be selfadjoint and positive definite. The vectors or matrices
  * X and B can be either dense or sparse.
  *
  * \tparam _MatrixType the type of the sparse matrix A, it must be a SparseMatrix<>
  * \tparam _UpLo the triangular part that will be used for the computations. It can be Lower
  *               or Upper. Default is Lower.
  *
  */
template<typename _MatrixType, int _UpLo = Lower>
class SimplicialCholesky
{
  public:
    typedef _MatrixType MatrixType;
    enum { UpLo = _UpLo };
    typedef typename MatrixType::Scalar Scalar;
    typedef typename MatrixType::RealScalar RealScalar;
    typedef typename MatrixType::Index Index;
    typedef SparseMatrix<Scalar,ColMajor,Index> CholMatrixType;
    typedef Matrix<Scalar,MatrixType::ColsAtCompileTime,1> VectorType;

  public:

    SimplicialCholesky()
      : m_info(Success), m_isInitialized(false), m_LDLt(true)
    {}

    SimplicialCholesky(const MatrixType& matrix)
      : m_info(Success), m_isInitialized(false), m_LDLt(true)
    {
      compute(matrix);
    }

    ~SimplicialCholesky()
    {
    }
    
    inline Index cols() const { return m_matrix.cols(); }
    inline Index rows() const { return m_matrix.rows(); }
  
    SimplicialCholesky& setMode(SimplicialCholeskyMode mode)
    {
      switch(mode)
      {
        case SimplicialCholeskyLLt:
          m_LDLt = false;
          break;
        case SimplicialCholeskyLDLt:
          m_LDLt = true;
          break;
        default:
          break;
      }
      
      return *this;
    }
    
    /** \brief Reports whether previous computation was successful.
      *
      * \returns \c Success if computation was succesful,
      *          \c NumericalIssue if the matrix.appears to be negative.
      */
    ComputationInfo info() const
    {
      eigen_assert(m_isInitialized && "Decomposition is not initialized.");
      return m_info;
    }

    /** Computes the sparse Cholesky decomposition of \a matrix */
    SimplicialCholesky& compute(const MatrixType& matrix)
    {
      analyzePattern(matrix);
      factorize(matrix);
      return *this;
    }
    
    /** \returns the solution x of \f$ A x = b \f$ using the current decomposition of A.
      *
      * \sa compute()
      */
    template<typename Rhs>
    inline const internal::solve_retval<SimplicialCholesky, Rhs>
    solve(const MatrixBase<Rhs>& b) const
    {
      eigen_assert(m_isInitialized && "SimplicialCholesky is not initialized.");
      eigen_assert(rows()==b.rows()
                && "SimplicialCholesky::solve(): invalid number of rows of the right hand side matrix b");
      return internal::solve_retval<SimplicialCholesky, Rhs>(*this, b.derived());
    }
    
    /** \returns the solution x of \f$ A x = b \f$ using the current decomposition of A.
      *
      * \sa compute()
      */
//     template<typename Rhs>
//     inline const internal::sparse_solve_retval<CholmodDecomposition, Rhs>
//     solve(const SparseMatrixBase<Rhs>& b) const
//     {
//       eigen_assert(m_isInitialized && "SimplicialCholesky is not initialized.");
//       eigen_assert(rows()==b.rows()
//                 && "SimplicialCholesky::solve(): invalid number of rows of the right hand side matrix b");
//       return internal::sparse_solve_retval<SimplicialCholesky, Rhs>(*this, b.derived());
//     }
    
    /** Performs a symbolic decomposition on the sparcity of \a matrix.
      *
      * This function is particularly useful when solving for several problems having the same structure.
      * 
      * \sa factorize()
      */
    void analyzePattern(const MatrixType& a);
    
    
    /** Performs a numeric decomposition of \a matrix
      *
      * The given matrix must has the same sparcity than the matrix on which the symbolic decomposition has been performed.
      *
      * \sa analyzePattern()
      */
    void factorize(const MatrixType& a);
    
    /** \returns the permutation P
      * \sa permutationPinv() */
    const PermutationMatrix<Dynamic>& permutationP() const
    { return m_P; }
    
    /** \returns the inverse P^-1 of the permutation P
      * \sa permutationP() */
    const PermutationMatrix<Dynamic>& permutationPinv() const
    { return m_Pinv; }
    
    #ifndef EIGEN_PARSED_BY_DOXYGEN
    /** \internal */
    template<typename Rhs,typename Dest>
    void _solve(const MatrixBase<Rhs> &b, MatrixBase<Dest> &dest) const
    {
      eigen_assert(m_factorizationIsOk && "The decomposition is not in a valid state for solving, you must first call either compute() or symbolic()/numeric()");
      eigen_assert(m_matrix.rows()==b.rows());
      
      if(m_info!=Success)
        return;
    
      if(m_P.size()>0)
        dest = m_Pinv * b;
      else
        dest = b;
      
      if(m_LDLt)
      {
        if(m_matrix.nonZeros()>0) // otherwise L==I
          m_matrix.template triangularView<UnitLower>().solveInPlace(dest);
      
        dest = m_diag.asDiagonal().inverse() * dest;
      
        if (m_matrix.nonZeros()>0) // otherwise L==I
          m_matrix.adjoint().template triangularView<UnitUpper>().solveInPlace(dest);
      }
      else
      {
        if(m_matrix.nonZeros()>0) // otherwise L==I
          m_matrix.template triangularView<Lower>().solveInPlace(dest);
      
        if (m_matrix.nonZeros()>0) // otherwise L==I
          m_matrix.adjoint().template triangularView<Upper>().solveInPlace(dest);
      }
      
      if(m_P.size()>0)
        dest = m_P * dest;
    }
    
    /** \internal */
    /*
    template<typename RhsScalar, int RhsOptions, typename RhsIndex, typename DestScalar, int DestOptions, typename DestIndex>
    void _solve(const SparseMatrix<RhsScalar,RhsOptions,RhsIndex> &b, SparseMatrix<DestScalar,DestOptions,DestIndex> &dest) const
    {
      // TODO
    }
    */
    #endif // EIGEN_PARSED_BY_DOXYGEN
    
    template<typename Stream>
    void dumpMemory(Stream& s)
    {
      int total = 0;
      s << "  L:        " << ((total+=(m_matrix.cols()+1) * sizeof(int) + m_matrix.nonZeros()*(sizeof(int)+sizeof(Scalar))) >> 20) << "Mb" << "\n";
      s << "  diag:     " << ((total+=m_diag.size() * sizeof(Scalar)) >> 20) << "Mb" << "\n";
      s << "  tree:     " << ((total+=m_parent.size() * sizeof(int)) >> 20) << "Mb" << "\n";
      s << "  nonzeros: " << ((total+=m_nonZerosPerCol.size() * sizeof(int)) >> 20) << "Mb" << "\n";
      s << "  perm:     " << ((total+=m_P.size() * sizeof(int)) >> 20) << "Mb" << "\n";
      s << "  perm^-1:  " << ((total+=m_Pinv.size() * sizeof(int)) >> 20) << "Mb" << "\n";
      s << "  TOTAL:    " << (total>> 20) << "Mb" << "\n";
    }

  protected:
    /** keeps off-diagonal entries; drops diagonal entries */
    struct keep_diag {
      inline bool operator() (const Index& row, const Index& col, const Scalar&) const
      {
        return row!=col;
      }
    };

    mutable ComputationInfo m_info;
    bool m_isInitialized;
    bool m_factorizationIsOk;
    bool m_analysisIsOk;
    bool m_LDLt;
    
    CholMatrixType m_matrix;
    VectorType m_diag;                  // the diagonal coefficients in case of a LDLt decomposition
    VectorXi m_parent;                  // elimination tree
    VectorXi m_nonZerosPerCol;
    PermutationMatrix<Dynamic> m_P;     // the permutation
    PermutationMatrix<Dynamic> m_Pinv;  // the inverse permutation
};

template<typename _MatrixType, int _UpLo>
void SimplicialCholesky<_MatrixType,_UpLo>::analyzePattern(const MatrixType& a)
{
  eigen_assert(a.rows()==a.cols());
  const Index size = a.rows();
  m_matrix.resize(size, size);
  m_parent.resize(size);
  m_nonZerosPerCol.resize(size);
  
  Index* tags = ei_aligned_stack_new(Index, size);
  
  // TODO allows to configure the permutation
  {
    CholMatrixType C;
    C = a.template selfadjointView<UpLo>();
    // remove diagonal entries:
    C.prune(keep_diag());
    internal::minimum_degree_ordering(C, m_P);
  }
  
  if(m_P.size()>0)
    m_Pinv  = m_P.inverse();
  else
    m_Pinv.resize(0);
  
  SparseMatrix<Scalar,ColMajor,Index> ap(size,size);
  ap.template selfadjointView<Upper>() = a.template selfadjointView<UpLo>().twistedBy(m_Pinv);
  
  for(Index k = 0; k < size; ++k)
  {
    /* L(k,:) pattern: all nodes reachable in etree from nz in A(0:k-1,k) */
    m_parent[k] = -1;             /* parent of k is not yet known */
    tags[k] = k;                  /* mark node k as visited */
    m_nonZerosPerCol[k] = 0;      /* count of nonzeros in column k of L */
    for(typename CholMatrixType::InnerIterator it(ap,k); it; ++it)
    {
      Index i = it.index();
      if(i < k)
      {
        /* follow path from i to root of etree, stop at flagged node */
        for(; tags[i] != k; i = m_parent[i])
        {
          /* find parent of i if not yet determined */
          if (m_parent[i] == -1)
            m_parent[i] = k;
          m_nonZerosPerCol[i]++;        /* L (k,i) is nonzero */
          tags[i] = k;                  /* mark i as visited */
        }
      }
    }
  }
  
  // release workspace
  ei_aligned_stack_delete(Index, tags, size);
  
  /* construct Lp index array from m_nonZerosPerCol column counts */
  Index* Lp = m_matrix._outerIndexPtr();
  Lp[0] = 0;
  for(Index k = 0; k < size; ++k)
    Lp[k+1] = Lp[k] + m_nonZerosPerCol[k] + (m_LDLt ? 0 : 1);

  m_matrix.resizeNonZeros(Lp[size]);
  
  m_isInitialized     = true;
  m_info              = Success;
  m_analysisIsOk      = true;
  m_factorizationIsOk = false;
}


template<typename _MatrixType, int _UpLo>
void SimplicialCholesky<_MatrixType,_UpLo>::factorize(const MatrixType& a)
{
  eigen_assert(m_analysisIsOk && "You must first call analyzePattern()");
  eigen_assert(a.rows()==a.cols());
  const Index size = a.rows();
  eigen_assert(m_parent.size()==size);
  eigen_assert(m_nonZerosPerCol.size()==size);

  const Index* Lp = m_matrix._outerIndexPtr();
  Index* Li = m_matrix._innerIndexPtr();
  Scalar* Lx = m_matrix._valuePtr();

  Scalar* y       = ei_aligned_stack_new(Scalar, size);
  Index* pattern  = ei_aligned_stack_new(Index, size);
  Index* tags     = ei_aligned_stack_new(Index, size);

  SparseMatrix<Scalar,ColMajor,Index> ap(size,size);
  ap.template selfadjointView<Upper>() = a.template selfadjointView<UpLo>().twistedBy(m_Pinv);
  
  bool ok = true;
  m_diag.resize(m_LDLt ? size : 0);
  
  for(Index k = 0; k < size; ++k)
  {
    // compute nonzero pattern of kth row of L, in topological order
    y[k] = 0.0;                     // Y(0:k) is now all zero
    Index top = size;               // stack for pattern is empty
    tags[k] = k;                    // mark node k as visited
    m_nonZerosPerCol[k] = 0;        // count of nonzeros in column k of L
    for(typename MatrixType::InnerIterator it(ap,k); it; ++it)
    {
      Index i = it.index();
      if(i <= k)
      {
        y[i] += internal::conj(it.value());            /* scatter A(i,k) into Y (sum duplicates) */
        Index len;
        for(len = 0; tags[i] != k; i = m_parent[i])
        {
          pattern[len++] = i;     /* L(k,i) is nonzero */
          tags[i] = k;            /* mark i as visited */
        }
        while(len > 0)
          pattern[--top] = pattern[--len];
      }
    }

    /* compute numerical values kth row of L (a sparse triangular solve) */
    Scalar d = y[k];                  // get D(k,k) and clear Y(k)
    y[k] = 0.0;
    for(; top < size; ++top)
    {
      Index i = pattern[top];       /* pattern[top:n-1] is pattern of L(:,k) */
      Scalar yi = y[i];             /* get and clear Y(i) */
      y[i] = 0.0;
      
      /* the nonzero entry L(k,i) */
      Scalar l_ki;
      if(m_LDLt)
        l_ki = yi / m_diag[i];       
      else
        yi = l_ki = yi / Lx[Lp[i]];
      
      Index p2 = Lp[i] + m_nonZerosPerCol[i];
      Index p;
      for(p = Lp[i] + (m_LDLt ? 0 : 1); p < p2; ++p)
        y[Li[p]] -= internal::conj(Lx[p]) * yi;
      d -= l_ki * internal::conj(yi);
      Li[p] = k;                          /* store L(k,i) in column form of L */
      Lx[p] = l_ki;
      ++m_nonZerosPerCol[i];              /* increment count of nonzeros in col i */
    }
    if(m_LDLt)
      m_diag[k] = d;
    else
    {
      Index p = Lp[k]+m_nonZerosPerCol[k]++;
      Li[p] = k ;                /* store L(k,k) = sqrt (d) in column k */
      Lx[p] = internal::sqrt(d) ;
    }
    if(d == Scalar(0))
    {
      ok = false;                         /* failure, D(k,k) is zero */
      break;
    }
  }

  // release workspace
  ei_aligned_stack_delete(Scalar, y, size);
  ei_aligned_stack_delete(Index, pattern, size);
  ei_aligned_stack_delete(Index, tags, size);
  
  m_info = ok ? Success : NumericalIssue;
  m_factorizationIsOk = true;
}

namespace internal {
  
template<typename _MatrixType, int _UpLo, typename Rhs>
struct solve_retval<SimplicialCholesky<_MatrixType,_UpLo>, Rhs>
  : solve_retval_base<SimplicialCholesky<_MatrixType,_UpLo>, Rhs>
{
  typedef SimplicialCholesky<_MatrixType,_UpLo> Dec;
  EIGEN_MAKE_SOLVE_HELPERS(Dec,Rhs)

  template<typename Dest> void evalTo(Dest& dst) const
  {
    dec()._solve(rhs(),dst);
  }
};

template<typename _MatrixType, int _UpLo, typename Rhs>
struct sparse_solve_retval<SimplicialCholesky<_MatrixType,_UpLo>, Rhs>
  : sparse_solve_retval_base<SimplicialCholesky<_MatrixType,_UpLo>, Rhs>
{
  typedef SimplicialCholesky<_MatrixType,_UpLo> Dec;
  EIGEN_MAKE_SPARSE_SOLVE_HELPERS(Dec,Rhs)

  template<typename Dest> void evalTo(Dest& dst) const
  {
    dec()._solve(rhs(),dst);
  }
};

}

#endif // EIGEN_SIMPLICIAL_CHOLESKY_H
