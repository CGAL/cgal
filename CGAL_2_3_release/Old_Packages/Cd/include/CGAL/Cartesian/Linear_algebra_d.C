// ======================================================================
//
// Copyright (c) 2000 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       :
// release_date  :
//
// file          : include/CGAL/Cartesian/Linear_algebra_d.C
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Herve.Bronnimann@sophia.inria.fr
// coordinator   : INRIA Sophia-Antipolis (Mariette.Yvinec@sophia.inria.fr)
//
// ======================================================================

#ifndef CGAL_CARTESIAN_LINEAR_ALGEBRA_D_C
#define CGAL_CARTESIAN_LINEAR_ALGEBRA_D_C

#include <CGAL/Cartesian/Linear_algebra_d.h>
#include <CGAL/Cartesian/Linear_algebra_vector.C>
#include <CGAL/Cartesian/Linear_algebra_matrix.C>
#include <algorithm>
#include <functional>

CGAL_BEGIN_NAMESPACE

template < class FT >
std::pair<int,int>
Linear_algebraCd<FT>::
transpose(std::pair<int,int> p) const
{
  std::swap(p.first,p.second);
  return p;
}

template < class FT >
Linear_algebraCd<FT>::Matrix
Linear_algebraCd<FT>::
transpose(const Linear_algebraCd<FT>::Matrix &M) const
{
  Matrix P(transpose(M.dimension()));
  int i, j;
  for (i=0; i<P.row_dimension(); ++i)
    for (j=0; j<P.column_dimension(); ++j)
      P[i][j] = M[j][i];
  return P;
}

template < class FT >
inline // in order to facilitate the optimization with unused variables
void
Linear_algebraCd<FT>::
Gaussian_elimination(const Linear_algebraCd<FT>::Matrix &M,
                     // return parameters
                     Linear_algebraCd<FT>::Matrix &L,
          	     Linear_algebraCd<FT>::Matrix &U,
                     std::vector<int> &row_permutation,
                     std::vector<int> &column_permutation,
         	     Linear_algebraCd<FT>::FT &det,
	             int &rank,
         	     Linear_algebraCd<FT>::Vector &c) const
{
  // Method: we use Gaussian elimination with division at each step
  // We do not use the variant by Bareiss (because we are on a field)
  // We obtain L, U, q, such that LM' = U, and U is upper diagonal,
  // where M' is M whose rows and columns are permuted
  // Picked up on the way:
  // det, rank, non-trivial solution c if system is degenerate (c^t M = 0)
  // Preconditions:
  CGAL_kernel_precondition( M.row_dimension() <= M.column_dimension() );
  CGAL_kernel_precondition( U.dimension() == M.dimension() );
  CGAL_kernel_precondition( L.row_dimension() == U.row_dimension() );
  CGAL_kernel_precondition( L.column_dimension() == U.row_dimension() );
  CGAL_kernel_precondition( M.row_dimension() == c.dimension() );
  // Temporaries
  int i, j, k;
  int dim  = M.row_dimension(), cdim = M.column_dimension();
  // All the parameters are already initialized (as in C++)
  int sign = 1;
  // First create a copy of M into U, and set L and permutations to identity
  std::copy(M.begin(), M.end(), U.begin());
  std::fill(L.begin(), L.end(), FT(0)); // should be unnecessary
  for (i=0; i<dim; ++i) {
    L[i][i] = FT(1);
    row_permutation.push_back(i);
  }
  for (i=0; i<cdim; ++i)
    column_permutation.push_back(i);
  // Main loop : invariant is that L * M[q] = U
  // M[q] stands for M with row i permuted with row q[i]
  //DEBUG: std::cerr << "START GAUSS ELIMINATION" << std::endl;
  det = 1;
  for (k=0; k<dim; ++k) {
    // Total pivoting, without looking for the maximum entry
    for (i=k,j=k;
         j<cdim && U[i][j] == FT(0);
         (++i==dim)? ++j,i=k : 0 );
    //DEBUG: std::cerr << "before swap [k="<<k<<"] :";
    //std::cerr << " found i="<<i<<" and j="<<j<<std::endl;
    //DEBUG: std::cerr << U << std::endl;
    if (j==cdim) break;
    if (i!=k) {
      //DEBUG: std::cerr << "swap row i="<<i<<" and k="<<k<<std::endl;
      std::swap_ranges(U.row_begin(k), U.row_end(k), U.row_begin(i));
      std::swap_ranges(L.row_begin(k), L.row_end(k), L.row_begin(i));
      std::swap(row_permutation[k], row_permutation[i]);
      sign = -sign;
    }
    if (j!=k) {
      //DEBUG: std::cerr << "swap column j="<<j<<" and k="<<k<<std::endl;
      U.swap_columns(j,k);
      L.swap_columns(j,k);
      std::swap(column_permutation[j],column_permutation[k]);
      sign = -sign;
    }
    //DEBUG: std::cerr << "after swap: " << std::endl;
    //DEBUG: std::cerr << U << std::endl;
    // At this point, U[k][k] is not zero
    FT pivot = U[k][k];
    det *= pivot;
    for (i=k+1; i<dim; ++i) {
      FT temp = U[i][k] / pivot;
      for (j=0; j<cdim; ++j)
        L[i][j] -= L[k][j] * temp;
      U[i][k] = FT(0);
      for (j=k+1; j<cdim; ++j)
        U[i][j] -= U[k][j] * temp;
    }
  }
  //DEBUG: std::cerr << "END GAUSS ELIMINATION" << std::endl;
  //std::cerr << "finally : "<< std::endl << L << std::endl << U << std::endl;
  // By invariant, L * M[q] = U and det(M) = det
  rank = k;
  if (rank == dim) {
    //DEBUG: std::cerr << "det=" << det << " sign=" << sign << std::endl; 
    det *= sign;
  } else {
    det = FT(0);
    // A vector c such that M[q] * c == 0 is obtained by L.row(dim-1)
    std::copy(L.row_begin(dim-1),L.row_end(dim-1),c.begin());
  }
}

template < class FT >
inline // in order to facilitate the optimization with unused variables
void
Linear_algebraCd<FT>::
Triangular_system_solver(const Linear_algebraCd<FT>::Matrix &U,
         	         const Linear_algebraCd<FT>::Vector &c,
                         // return parameters
         	         Linear_algebraCd<FT>::Vector &x,
         	         Linear_algebraCd<FT>::FT &D) const
{
  // METHOD: solve system Ux=c, returning solution (x/D)
  // back substitution of x[rdim], x[rdim-1], etc.
  // depends on "free" variables x[rdim+1], etc. x[cdim]
  CGAL_kernel_assertion( U.column_dimension() == x.dimension());
  CGAL_kernel_assertion( U.row_dimension() == c.dimension());
  //DEBUG: std::cerr << "system : " << U << std::endl << c << std::endl;
  //DEBUG: std::cerr << "left hand : " << c << std::endl;
  int i, j;
  for (i=U.row_dimension()-1; i>=0; --i) {
    x[i] = c[i];
    for (j=i+1; j<U.column_dimension(); ++j) 
      x[i] -= U[i][j] * x[j];
    x[i] /= U[i][i];
  }
  D = FT(1);
  //DEBUG: std::cerr << "finally : " << x << "/" << D << std::endl;
  // std::cerr << "sanity check : " << U*x << " == " << c*D << std::endl;
}

template < class FT >
inline // in order to facilitate the optimization with unused variables
void
Linear_algebraCd<FT>::
Triangular_left_inverse(const Linear_algebraCd<FT>::Matrix &U,
                        // return parameters
         	        Linear_algebraCd<FT>::Matrix &Uinv) const
{
  int i, j, k;
  CGAL_kernel_precondition( U.dimension() == transpose(Uinv.dimension()) );
  //DEBUG: std::cerr << "system : " << U << std::endl;
  // std::fill(Uinv.begin(), Uinv.end(), FT(0));
  for (i=U.row_dimension()-1; i>=0; --i) {
    Uinv[i][i] = FT(1)/U[i][i];
    for (j=i+1; j<U.column_dimension(); ++j) {
      for (k=i; k<j; ++k)
        Uinv[i][j] -= Uinv[i][k] * U[k][j];
      Uinv[i][j] /= U[j][j];
    }
  }
  //DEBUG: std::cerr << "finally : " << Uinv << std::endl;
  // std::cerr << "sanity check : " << U*Uinv << std::endl;
}

template < class FT >
bool
Linear_algebraCd<FT>::
inverse(const Linear_algebraCd<FT>::Matrix &M,
        Linear_algebraCd<FT>::Matrix &I, FT &D,
        Linear_algebraCd<FT>::Vector &c) const
{
  Matrix L(M.dimension());
  Matrix U(M.dimension());
  Matrix Uinv(M.column_dimension(),M.row_dimension());
  int rank;
  std::vector<int> rq, cq;
  Gaussian_elimination(M, L, U, rq, cq, D, rank, c);
  if (D == FT(0)) return false; // c holds the witness
  // Otherwise, compute the inverse of U
  Triangular_left_inverse(U,Uinv);
  Uinv = Uinv * L;
  // Don't forget to permute the rows of M back
  //DEBUG: std::cerr << "inverse before permutation : " << I << std::endl;
  for (rank=0; rank<I.column_dimension(); ++rank)
    std::copy(Uinv.row_begin(cq[rank]), Uinv.row_end(cq[rank]),
              I.row_begin(rank));
  D = FT(1);
  return true;
}

template < class FT >
inline
Linear_algebraCd<FT>::Matrix
Linear_algebraCd<FT>::
inverse(const Linear_algebraCd<FT>::Matrix &M,
        Linear_algebraCd<FT>::FT &D) const
{
  CGAL_kernel_precondition( M.row_dimension() == M.column_dimension() );
  Matrix I(M.column_dimension(),M.row_dimension());
  Vector c(M.row_dimension());
  bool result = inverse(M,I,D,c);
  CGAL_kernel_precondition( result );
  return I;
}

template < class FT >
Linear_algebraCd<FT>::FT
Linear_algebraCd<FT>::
determinant(const Linear_algebraCd<FT>::Matrix &M,
            Linear_algebraCd<FT>::Matrix &L,
            Linear_algebraCd<FT>::Matrix &U,
            std::vector<int> &q,
            Linear_algebraCd<FT>::Vector &c) const
{
  FT det;
  int rank;
  std::vector<int> cq;
  Gaussian_elimination(M, L, U, q, cq, det, rank, c);
  return det;
}

template < class FT >
inline
Linear_algebraCd<FT>::FT
Linear_algebraCd<FT>::
determinant(const Linear_algebraCd<FT>::Matrix &M) const
{
  Matrix L(M.dimension());
  Matrix U(M.dimension());
  std::vector<int> q;
  Vector c(M.column_dimension());
  return determinant(M,L,U,q,c);
}

template < class FT >
inline
Sign
Linear_algebraCd<FT>::
sign_of_determinant(const Linear_algebraCd<FT>::Matrix &M) const
{
  return CGAL::sign(determinant(M));
}

template < class FT >
bool
Linear_algebraCd<FT>::
verify_determinant(const Linear_algebraCd<FT>::Matrix & /*M*/,
                   const Linear_algebraCd<FT>::Matrix & /*L*/,
		   const Linear_algebraCd<FT>::Matrix & /*U*/,
                   const FT & /*D*/,
		   const std::vector<int> & /*q*/,
		   const Linear_algebraCd<FT>::Vector & /*c*/) const
{
  // TODO: verify_determinant
  CGAL_kernel_assertion( false );
  return false;
}

template < class FT >
bool
Linear_algebraCd<FT>::
linear_solver(const Linear_algebraCd<FT>::Matrix &/*M*/,
              const Linear_algebraCd<FT>::Vector &/*b*/,
              Linear_algebraCd<FT>::Vector &/*x*/, FT &/*D*/,
	      Linear_algebraCd<FT>::Matrix &/*spanning_vectors*/,
              Linear_algebraCd<FT>::Vector &/*c*/) const
{
  // TODO: full linear_solver() with spanning vectors
  CGAL_kernel_assertion( false );
  return false;
}

template < class FT >
inline
bool
Linear_algebraCd<FT>::
linear_solver(const Linear_algebraCd<FT>::Matrix &M,
              const Linear_algebraCd<FT>::Vector &b,
              Linear_algebraCd<FT>::Vector &x, FT &D,
              Linear_algebraCd<FT>::Vector &c) const
{
  CGAL_kernel_precondition( M.row_dimension() <= M.column_dimension() );
  CGAL_kernel_precondition( b.dimension() == M.row_dimension() );
  CGAL_kernel_precondition( x.dimension() == M.column_dimension() );
  CGAL_kernel_precondition( c.dimension() == M.column_dimension() );
  Matrix L(M.dimension());
  Matrix U(M.dimension());
  FT det;
  int rank;
  std::vector<int> rq, cq;

  //DEBUG: std::cerr << "system : " << M << std::endl << b << std::endl;
  Gaussian_elimination(M, L, U, rq, cq, D, rank, c);
  //DEBUG: std::cerr << "transformed : " << L << std::endl << U << std::endl;
  if (rank != U.row_dimension()) return false;
  // Compute a solution by solving triangular system
  // Since LM=U, and x is a solution of Mx=b, then Ux=Lb
  // Temporary store the solution in c
  Triangular_system_solver(U, L*b, c, D);
  // Don't forget to permute the rows of M back
  //DEBUG: std::cerr << "solution before permutation : " << c << std::endl;
  for (rank=0; rank<U.row_dimension(); ++rank)
    x[ cq[rank] ] = c[rank];
  return true;
}

template < class FT >
inline
bool
Linear_algebraCd<FT>::
linear_solver(const Linear_algebraCd<FT>::Matrix &M,
              const Linear_algebraCd<FT>::Vector &b,
              Linear_algebraCd<FT>::Vector &x, FT &D) const
{
  Vector c(M.column_dimension());
  return linear_solver(M, b, x, D, c);
}

template < class FT >
inline
bool
Linear_algebraCd<FT>::
is_solvable(const Linear_algebraCd<FT>::Matrix &M,
            const Linear_algebraCd<FT>::Vector &b) const
{
  Vector x(M.column_dimension());
  FT D;
  return linear_solver(M, b, x, D);
}

template < class FT >
bool
Linear_algebraCd<FT>::
homogeneous_linear_solver(const Linear_algebraCd<FT>::Matrix &M,
                          Linear_algebraCd<FT>::Vector &x) const
{
  // TODO: homogeneous_linear_solver
}

template < class FT >
int
Linear_algebraCd<FT>::
homogeneous_linear_solver(const Linear_algebraCd<FT>::Matrix &M,
                          Linear_algebraCd<FT>::Matrix &spanning_vectors) const
{
  // TODO: full homogeneous_linear_solver() with spanning vectors
}

template < class FT >
int
Linear_algebraCd<FT>::
rank(const Linear_algebraCd<FT>::Matrix &M, std::vector<int> &q) const
{
  int rank;
  std::vector<int> cq;
  Gaussian_elimination(M, Matrix(M.dimension()), Matrix(M.dimension()),
                       q, cq, FT(0), rank, Vector(M.row_dimension()));
  return rank;
}

template < class FT >
inline
int
Linear_algebraCd<FT>::
rank(const Linear_algebraCd<FT>::Matrix &M)
{
  std::vector<int> q;
  return rank(M,q);
}

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_LINEAR_ALGEBRA_D_C
