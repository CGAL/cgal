// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Herve.Bronnimann@sophia.inria.fr

#ifndef CGAL_CARTESIAN_LINEAR_ALGEBRA_D_C
#define CGAL_CARTESIAN_LINEAR_ALGEBRA_D_C

#include <CGAL/Cartesian/Linear_algebra_d.h>
#include <CGAL/Cartesian/Linear_algebra_vector.C>
#include <CGAL/Cartesian/Linear_algebra_matrix.C>
#include <algorithm>
#include <functional>

CGAL_BEGIN_NAMESPACE

template < class FT >
Linear_algebraCd<FT>::Matrix
Linear_algebraCd<FT>::
transpose(const Linear_algebraCd<FT>::Matrix &M) const
{
  Matrix P(M.column_dimension(),M.row_dimension());
  int i, j;
  for (i=0; i<M.row_dimension(); ++i)
    for (j=0; j<M.column_dimension(); ++j)
      P[j][i] = M[i][j];
}

template < class FT >
inline // in order to facilitate the optimization with unused variables
void
Linear_algebraCd<FT>::
Gaussian_elimination(const Linear_algebraCd<FT>::Matrix &M,
                     // return parameters
                     Linear_algebraCd<FT>::Matrix &L,
          	     Linear_algebraCd<FT>::Matrix &U,
         	     Linear_algebraCd<FT>::FT &det,
	             int &rank,
                     std::vector<int> &q,
         	     Linear_algebraCd<FT>::Vector &c) const
{
  // Method: we use Gaussian elimination with division at each step
  // We do not use the variant by Bareiss (because we are on a field)
  // TODO treat case when M is not square matrix
  // Preconditions:
  CGAL_kernel_precondition( M.row_dimension() == M.column_dimension() );
  CGAL_kernel_precondition( L.row_dimension() == L.column_dimension() );
  CGAL_kernel_precondition( U.row_dimension() == U.column_dimension() );
  CGAL_kernel_precondition( M.dimension() == U.dimension() );
  CGAL_kernel_precondition( M.dimension() == L.dimension() );
  CGAL_kernel_precondition( M.row_dimension() == c.dimension() );
  // Temporaries
  int i, j, k;
  int dim  = M.row_dimension();
  // All the parameters are already initialized (as in C++)
  int sign = 1;
  // First create a copy of M into U, and set L to identity
  std::copy(M.begin(), M.end(), U.begin());
  std::fill(L.begin(), L.end(), FT(0)); // should be unnecessary
  for (i=0; i<dim; ++i) L[i][i] = FT(1);
  // Main loop : invariant is that L * M[q] = U
  // M[q] stands for M with row i permuted with row q[i]
  //DEBUG: std::cerr << "START GAUSS ELIMINATION" << std::endl;
  det = 1;
  for (k=0; k<dim; ++k) {
    // Partial pivoting, without looking for the maximum entry
    for (i=k; i<dim && (U[i][k] == FT(0)); ++i) ;
    //DEBUG: std::cerr << "before swap [k="<<k<<"] : found i="<<i<<std::endl;
    //DEBUG: std::cerr << U << std::endl;
    if (i==dim) break;
    if (i!=k) {
      std::swap_ranges(U.row_begin(k), U.row_end(k), U.row_begin(i));
      std::swap_ranges(L.row_begin(k), L.row_end(k), L.row_begin(i));
      q.push_back(i);
      sign = -sign;
    }
    //DEBUG: std::cerr << "after swap: " << std::endl;
    //DEBUG: std::cerr << U << std::endl;
    // At this point, U[k][k] is not zero
    FT pivot = U[k][k];
    det *= pivot;
    for (i=k+1; i<dim; ++i) {
      FT temp = U[i][k] / pivot;
      for (j=0; j<dim; ++j)
        L[i][j] -= L[k][j] * temp;
      for (j=k+1; j<dim; ++j)
        U[i][j] -= U[k][j] * temp;
      U[i][k] = FT(0);
    }
  }
  //DEBUG: std::cerr << "END GAUSS ELIMINATION" << std::endl;
  // By invariant, denom * L * M[q] = U, and det(U)/det(L) = U[dim-1][dim-1]
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
bool
Linear_algebraCd<FT>::
inverse(const Linear_algebraCd<FT>::Matrix &M,
        Linear_algebraCd<FT>::Matrix &I, FT &det,
        Linear_algebraCd<FT>::Vector &c) const
{
  Matrix L(M.dimension());
  Matrix U(M.dimension());
  int rank;
  Gaussian_elimination(M, L, U, det, rank, q, c);
  // TODO
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
  Matrix I(M.dimension());
  CGAL_kernel_precondition( inverse(M,I,D,Vector()) );
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
  Gaussian_elimination(M, L, U, det, rank, q, c);
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
  FT det = determinant(M,L,U,q,c);
  return det;
}

template < class FT >
inline
Sign
Linear_algebraCd<FT>::
sign_of_determinant(const Linear_algebraCd<FT>::Matrix &M) const
{
  return CGAL::sign(M.determinant());
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
  // TODO
  CGAL_kernel_assertion( false );
  return false;
}

template < class FT >
bool
Linear_algebraCd<FT>::
linear_solver(const Linear_algebraCd<FT>::Matrix &M,
              const Linear_algebraCd<FT>::Vector &b,
              Linear_algebraCd<FT>::Vector &x, FT &D,
	      Linear_algebraCd<FT>::Matrix &spanning_vectors,
              Linear_algebraCd<FT>::Vector &c) const
{
  // TODO
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
  return linear_solver(M, b, x, D, Matrix(M.dimension()), c);
}

template < class FT >
inline
bool
Linear_algebraCd<FT>::
linear_solver(const Linear_algebraCd<FT>::Matrix &M,
              const Linear_algebraCd<FT>::Vector &b,
              Linear_algebraCd<FT>::Vector &x, FT &D) const
{
  return linear_solver(M, b, x, D,
		       Matrix(M.dimension()),
		       Vector(M.column_dimension()));
}

template < class FT >
inline
bool
Linear_algebraCd<FT>::
is_solvable(const Linear_algebraCd<FT>::Matrix &M,
            const Linear_algebraCd<FT>::Vector &b) const
{
  return linear_solver(M, b, Vector(column_dimension()), FT(0),
		       Matrix(M.dimension()),
		       Vector(M.column_dimension()));
}

template < class FT >
bool
Linear_algebraCd<FT>::
homogeneous_linear_solver(const Linear_algebraCd<FT>::Matrix &M,
                          Linear_algebraCd<FT>::Vector &x) const
{
  // TODO
}

template < class FT >
int
Linear_algebraCd<FT>::
homogeneous_linear_solver(const Linear_algebraCd<FT>::Matrix &M,
                          Linear_algebraCd<FT>::Matrix &spanning_vectors) const
{
  // TODO
}

template < class FT >
int
Linear_algebraCd<FT>::
rank(const Linear_algebraCd<FT>::Matrix &M, std::vector<int> &q) const
{
  int rank;
  Gaussian_elimination(M, Matrix(M.dimension()), Matrix(M.dimension()),
                       FT(0), rank, q, Vector(M.row_dimension()));
  return rank;
}

template < class FT >
inline
int
Linear_algebraCd<FT>::
rank(const Linear_algebraCd<FT>::Matrix &M)
{
  return rank(M,std::vector<int>());
}

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_LINEAR_ALGEBRA_D_C
