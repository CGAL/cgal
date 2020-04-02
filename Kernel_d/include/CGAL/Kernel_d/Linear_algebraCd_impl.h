// Copyright (c) 1999
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Herve.Bronnimann@sophia.inria.fr
//                 Michael.Seel@mpi-sb.mpg.de

#ifndef CGAL_LINEAR_ALGEBRACD_C
#define CGAL_LINEAR_ALGEBRACD_C

#include <algorithm>
#include <functional>

namespace CGAL {

template < class FT, class AL >
typename Linear_algebraCd<FT,AL>::Matrix
Linear_algebraCd<FT,AL>::
transpose(const Matrix &M)
{
  Matrix P(transpose(M.dimension()));
  int i, j;
  for (i=0; i<P.row_dimension(); ++i)
    for (j=0; j<P.column_dimension(); ++j)
      P[i][j] = M[j][i];
  return P;
}

template < class FT, class AL >
inline // in order to facilitate the optimization with unused variables
void
Linear_algebraCd<FT,AL>::
Gaussian_elimination(const Matrix &M,
                     // return parameters
                     Matrix &L, Matrix &U,
                     std::vector<int> &row_permutation,
                     std::vector<int> &column_permutation,
                     FT &det, int &rank, Vector &c)
{
  // Method: we use Gaussian elimination with division at each step
  // We do not use the variant by Bareiss (because we are on a field)
  // We obtain L, U, q, such that LM' = U, and U is upper diagonal,
  // where M' is M whose rows and columns are permuted
  // Picked up on the way:
  // det, rank, non-trivial solution c if system is degenerate (c^t M = 0)

  int i, j, k;
  int dim  = M.row_dimension(), cdim = M.column_dimension();
  // All the parameters are already initialized (as in C++)
  int sign = 1;
  // First create a copy of M into U, and set L and permutations to identity
  U = M;
  typename Matrix::Identity IDENTITY;
  L = Matrix(dim,IDENTITY);
  for (i=0; i<dim; ++i) row_permutation.push_back(i);
  for (i=0; i<cdim; ++i) column_permutation.push_back(i);
  // Main loop : invariant is that L * M[q] = U
  // M[q] stands for M with row i permuted with row q[i]
    CGAL_KD_TRACEN("START GAUSS ELIMINATION");
  det = 1;
  for (k=0; k<dim; ++k) {
    // Total pivoting, without looking for the maximum entry
    for (i=k,j=k;
         j<cdim && U[i][j] == FT(0);
         (void)((++i==dim)? ++j,i=k : 0 )) {}
      CGAL_KD_TRACEN("before swap [k="<<k<<"] :");
      CGAL_KD_TRACEN(" found i="<<i<<" and j="<<j);
      CGAL_KD_TRACEV(U);
    if (j==cdim) break;
    if (i!=k) {
      CGAL_KD_TRACEN("swap row i="<<i<<" and k="<<k);
      U.swap_rows(i,k); L.swap_rows(i,k);
      std::swap(row_permutation[k], row_permutation[i]);
      sign = -sign;
    }
    if (j!=k) {
      CGAL_KD_TRACEN("swap column j="<<j<<" and k="<<k);
      U.swap_columns(j,k);
      std::swap(column_permutation[j],column_permutation[k]);
      sign = -sign;
    }
      CGAL_KD_TRACEN("after swap: "<<U);
    FT pivot = U[k][k];
    CGAL_assertion(pivot != FT(0));
    det *= pivot;
    for (i=k+1; i<dim; ++i) {
      FT temp = U[i][k] / pivot;
      for (j=0; j<dim; ++j)
        L[i][j] -= L[k][j] * temp;
      U[i][k] = FT(0);
      for (j=k+1; j<cdim; ++j)
        U[i][j] -= U[k][j] * temp;
    }
  }

    CGAL_KD_TRACEN("END GAUSS ELIMINATION"); CGAL_KD_TRACEV(L);CGAL_KD_TRACEV(U);
    // By invariant, L * M[q] = U and det(M) = det
  rank = k;
  if (rank == dim) {
    det *= FT(sign);
  } else {
    det = FT(0);
    // A vector c such that M[q] * c == 0 is obtained by L.row(dim-1)
    c = L.row(dim-1);
  }
}

template < class FT, class AL >
bool
Linear_algebraCd<FT,AL>::
Triangular_system_solver(const Matrix &U, const Matrix& L, const Vector &b,
                         int rank, Vector &x, FT &D)
{
  // METHOD: solve system Ux=b, returning solution (x/D)
  // back substitution of x[rdim], x[rdim-1], etc.
  // depends on "free" variables x[rdim+1], etc. x[cdim]
  CGAL_kernel_assertion( U.row_dimension() == b.dimension());
    CGAL_KD_TRACEN("Triangular_system_solver");CGAL_KD_TRACEV(U);CGAL_KD_TRACEV(b);
  D = FT(1); int i;
  for (i = rank; i < U.row_dimension(); ++i)
      if ( b[i] != FT(0) ) {
          x = L.row(i);
          return false;
      }

  x = Vector(U.column_dimension());
  for (i = rank-1; i>=0; --i) {
    x[i] = b[i];
    for (int j = i+1; j<rank; ++j)
      x[i] -= U[i][j] * x[j];
    x[i] /= U[i][i];
  }
  return true;
}

template < class FT, class AL >
inline // in order to facilitate the optimization with unused variables
void
Linear_algebraCd<FT,AL>::
Triangular_left_inverse(const Matrix &U, Matrix &Uinv)
{
  int i, j, k;
  CGAL_kernel_precondition(U.dimension() == transpose(Uinv.dimension()));
    CGAL_KD_TRACEN("system : " << U);
  for (i=U.row_dimension()-1; i>=0; --i) {
    Uinv[i][i] = FT(1)/U[i][i];
    for (j=i+1; j<U.column_dimension(); ++j) {
      for (k=i; k<j; ++k)
        Uinv[i][j] -= Uinv[i][k] * U[k][j];
      Uinv[i][j] /= U[j][j];
    }
  }
    CGAL_KD_TRACEN("finally : " << Uinv);
}

template < class FT, class AL >
bool
Linear_algebraCd<FT,AL>::
inverse(const Matrix &M, Matrix &I, FT &D, Vector &c)
{
  CGAL_kernel_precondition(M.row_dimension()==M.column_dimension());
  int rank;
  Matrix L,U;
  std::vector<int> rq, cq;
  Gaussian_elimination(M, L, U, rq, cq, D, rank, c);
  if (D == FT(0)) return false; // c holds the witness
  // Otherwise, compute the inverse of U
  Matrix Uinv(M.column_dimension(),M.row_dimension());
  Triangular_left_inverse(U,Uinv);
  Uinv = Uinv * L;
  // Don't forget to permute the rows of M back

    CGAL_KD_TRACEN("inverse before permutation : "<<I);
  I = Matrix(M.column_dimension(),M.row_dimension());
  typename Matrix::row_iterator rit, wit;
  for (rank=0; rank<I.column_dimension(); ++rank)
    for ( wit = I.row_begin(rank), rit = Uinv.row_begin(cq[rank]);
          rit != Uinv.row_end(cq[rank]); ++wit, ++rit ) *wit = *rit;
  /* does not work with MS:
    std::copy(Uinv.row_begin(cq[rank]), Uinv.row_end(cq[rank]),
              I.row_begin(rank)); */
  D = FT(1);
  return true;
}

template < class FT, class AL >
inline
typename Linear_algebraCd<FT,AL>::Matrix
Linear_algebraCd<FT,AL>::
inverse(const Matrix &M, FT &D)
{
  CGAL_kernel_precondition(M.row_dimension()==M.column_dimension());
  Matrix I; Vector c;
  CGAL_kernel_assertion_code( bool result = ) inverse(M,I,D,c);
  CGAL_kernel_assertion( result );
  return I;
}

template < class FT, class AL >
typename Linear_algebraCd<FT,AL>::FT
Linear_algebraCd<FT,AL>::
determinant(const Matrix &M, Matrix &L, Matrix &U,
            std::vector<int> &q, Vector &c)
{
  FT det; int rank;
  std::vector<int> cq;
  Gaussian_elimination(M, L, U, q, cq, det, rank, c);
  return det;
}

template < class FT, class AL >
inline
typename Linear_algebraCd<FT,AL>::FT
Linear_algebraCd<FT,AL>::
determinant(const Matrix &M)
{
  Matrix L,U; Vector c;
  std::vector<int> q;
  return determinant(M,L,U,q,c);
}

template < class FT, class AL >
inline
Sign
Linear_algebraCd<FT,AL>::
sign_of_determinant(const Matrix &M)
{ return CGAL_NTS sign(determinant(M)); }

template < class FT, class AL >
bool
Linear_algebraCd<FT,AL>::
verify_determinant(const Matrix & /*M*/,
                   const FT & /*D*/,
                   const Matrix & /*L*/,
                   const Matrix & /*U*/,
                   const std::vector<int> & /*q*/,
                   const Vector & /*c*/)
{
  // TODO: verify_determinant
  return true;
}

template < class FT, class AL >
bool
Linear_algebraCd<FT,AL>::
linear_solver(const Matrix &M, const Vector &b,
              Vector &x, FT &D, Vector &c)
{
  CGAL_kernel_precondition( b.dimension() == M.row_dimension() );
  Matrix L,U;
  int rank;
  std::vector<int> dummy, var;
    CGAL_KD_TRACEN("linear_solver");CGAL_KD_TRACEV(M); CGAL_KD_TRACEV(b);
  Gaussian_elimination(M, L, U, dummy, var, D, rank, c);
  // Compute a solution by solving triangular system
  // Since LM=U, and x is a solution of Mx=b, then Ux=Lb
  // Temporary store the solution in c
  if ( !Triangular_system_solver(U, L, L*b, rank, c, D) )
    return false;
  // Don't forget to permute the rows of M back
  x = Vector(M.column_dimension());
  for (int i=0; i<U.column_dimension(); ++i)
    x[ var[i] ] = c[i];
#ifdef CGAL_LA_SELFTEST
  CGAL_assertion( (M*x) == b );
#endif

  return true;
}


template < class FT, class AL >
bool
Linear_algebraCd<FT,AL>::
linear_solver(const Matrix &M, const Vector &b,
              Vector &x, FT &D, Matrix &spanning_vectors,
              Vector &c)
{
  CGAL_kernel_precondition( b.dimension() == M.row_dimension() );
  Matrix L,U;
  int rank;
  std::vector<int> dummy, var;
    CGAL_KD_TRACEN("linear_solver");CGAL_KD_TRACEV(M); CGAL_KD_TRACEV(b);
  Gaussian_elimination(M, L, U, dummy, var, D, rank, c);
  // Compute a solution by solving triangular system
  // Since LM=U, and x is a solution of Mx=b, then Ux=Lb
  // Temporary store the solution in c
  if ( !Triangular_system_solver(U, L, L*b, rank, c, D) )
    return false;
  // Don't forget to permute the rows of M back:
  x = Vector(M.column_dimension());
  for (int i=0; i<U.column_dimension(); ++i)
    x[ var[i] ] = c[i];

#ifdef CGAL_LA_SELFTEST
  CGAL_assertion( (M*x) == b );
#endif
  int defect = M.column_dimension() - rank;
  spanning_vectors = Matrix(M.column_dimension(),defect);
  if (defect > 0) {
    for(int l=0; l < defect; ++l) {
      spanning_vectors(var[rank + l],l)=FT(1);
      for(int i = rank - 1; i >= 0 ; i--) {
        FT h = - U(i,rank+l);
        for (int j= i + 1; j<rank; ++j)
          h -= U(i,j)*spanning_vectors(var[j],l);
        spanning_vectors(var[i],l)= h / U(i,i);
      }
      CGAL_KD_TRACEV(spanning_vectors.column(l));

#ifdef CGAL_LA_SELFTEST
      CGAL_assertion( (M*spanning_vectors.column(l)).is_zero() );
#endif
    }
  }
  return true;
}




template < class FT, class AL >
bool
Linear_algebraCd<FT,AL>::
linear_solver(const Matrix &M, const Vector &b, Vector &x, FT &D)
{ Vector c;
  return linear_solver(M, b, x, D, c);
}

template < class FT, class AL >
bool
Linear_algebraCd<FT,AL>::
is_solvable(const Matrix &M,  const Vector &b)
{ Vector x; FT D;
  return linear_solver(M, b, x, D);
}

template < class FT, class AL >
int
Linear_algebraCd<FT,AL>::
homogeneous_linear_solver(const Matrix &M, Matrix &spanning_vectors)
{
  Matrix L,U; Vector c, b(M.row_dimension());
  std::vector<int> dummy,var;
  FT D; int rank,i;
  Gaussian_elimination(M, L, U, dummy, var, D, rank, c);

#ifdef CGAL_LA_SELFTEST
  Vector x;
  Triangular_system_solver(U, L, b, rank, c, D);
  CGAL_KD_TRACEV(M);CGAL_KD_TRACEV(U);CGAL_KD_TRACEV(b);CGAL_KD_TRACEV(rank);CGAL_KD_TRACEV(c);CGAL_KD_TRACEV(D);
  x = Vector(M.column_dimension());
  for (i=0; i<U.row_dimension(); ++i)
    x[ var[i] ] = c[i];
  CGAL_assertion( (M*x).is_zero() );
#endif

  int defect = M.column_dimension() - rank;
  spanning_vectors = Matrix(M.column_dimension(),defect);

  if (defect > 0) {
    /* In the $l$-th spanning vector, $0 \le l < |defect|$ we set
       variable |var[rank + l]| to $1 = |denom|/|denom|$ and then the
       dependent variables as dictated by the $|rank| + l$ - th column of
       |C|.*/

    for(int l=0; l < defect; ++l) {
      spanning_vectors(var[rank + l],l)=FT(1);
      for(i = rank - 1; i >= 0 ; i--) {
        FT h = - U(i,rank+l);
        for (int j= i + 1; j<rank; ++j)
          h -= U(i,j)*spanning_vectors(var[j],l);
        spanning_vectors(var[i],l)= h / U(i,i);
      }
      CGAL_KD_TRACEV(spanning_vectors.column(l));

#ifdef CGAL_LA_SELFTEST
      /* we check whether the $l$ - th spanning vector is a solution
         of the homogeneous system */
      if ( !(M*spanning_vectors.column(l)).is_zero() )
        std::cerr << M*spanning_vectors.column(l) << std::endl;
      CGAL_assertion( (M*spanning_vectors.column(l)).is_zero() );
#endif
    }
  }
  return defect;
}

template < class FT, class AL >
bool
Linear_algebraCd<FT,AL>::
homogeneous_linear_solver(const Matrix &M, Vector &x)
{
  x = Vector(M.row_dimension());
  Matrix spanning_vectors;
  int defect = homogeneous_linear_solver(M,spanning_vectors);
  if (defect == 0) return false;

  /* return first column of |spanning_vectors| */
  for (int i = 0; i < spanning_vectors.row_dimension(); i++)
    x[i] = spanning_vectors(i,0);
  return true;

}


template < class FT, class AL >
int
Linear_algebraCd<FT,AL>::
independent_columns(const Matrix &M, std::vector<int> &q)
{
  int rank;
  std::vector<int> dummy;
  Matrix Dummy1,Dummy2; Vector Dummy3; FT null(0);
  Gaussian_elimination(M, Dummy1, Dummy2, dummy, q, null, rank, Dummy3);
  return rank;
}

template < class FT, class AL >
int
Linear_algebraCd<FT,AL>::
rank(const Matrix &M)
{ std::vector<int> q;
  return independent_columns(M,q);
}

} //namespace CGAL
#endif // CGAL_LINEAR_ALGEBRACD_C
