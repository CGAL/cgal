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
#include <vector>

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
      P(j,i) = M(i,j);
}

template < class FT >
bool
Linear_algebraCd<FT>::
inverse(const Linear_algebraCd<FT>::Matrix &M,
        Linear_algebraCd<FT>::Matrix &I, FT &D,
        Linear_algebraCd<FT>::Vector &c) const
{
}

template < class FT >
Linear_algebraCd<FT>::Matrix
Linear_algebraCd<FT>::
inverse(const Linear_algebraCd<FT>::Matrix M, FT &D) const
{
}

template < class FT >
FT
Linear_algebraCd<FT>::
determinant(const Linear_algebraCd<FT>::Matrix &M,
            Linear_algebraCd<FT>::Matrix &L,
	    Linear_algebraCd<FT>::Matrix &U,
            std::vector<int> &q,
	    Linear_algebraCd<FT>::Vector &c) const
{
}

template < class FT >
bool
Linear_algebraCd<FT>::
verify_determinant(const Linear_algebraCd<FT>::Matrix &M,
                   const Linear_algebraCd<FT>::Matrix &L,
		   const Linear_algebraCd<FT>::Matrix &U,
                   const FT &D,
		   const std::vector<int> &q,
		   const Linear_algebraCd<FT>::Vector &c) const
{
}

template < class FT >
FT
Linear_algebraCd<FT>::
determinant(const Linear_algebraCd<FT>::Matrix &M) const
{
}

template < class FT >
Sign
Linear_algebraCd<FT>::
sign_of_determinant(const Linear_algebraCd<FT>::Matrix &M) const
{
}

template < class FT >
bool
Linear_algebraCd<FT>::
linear_solver(const Linear_algebraCd<FT>::Matrix &M, const Linear_algebraCd<FT>::Vector &b,
                   Linear_algebraCd<FT>::Vector &x, FT &D, Linear_algebraCd<FT>::Matrix &spanning_vectors,
                   Linear_algebraCd<FT>::Vector &c) const
{
}

template < class FT >
bool
Linear_algebraCd<FT>::
linear_solver(const Linear_algebraCd<FT>::Matrix &M,
              const Linear_algebraCd<FT>::Vector &b,
              Linear_algebraCd<FT>::Vector &x, FT &D,
              Linear_algebraCd<FT>::Vector &c) const
{
}

template < class FT >
bool
Linear_algebraCd<FT>::
linear_solver(const Linear_algebraCd<FT>::Matrix &M,
              const Linear_algebraCd<FT>::Vector &b,
              Linear_algebraCd<FT>::Vector &x, FT &D) const
{
}

template < class FT >
bool
Linear_algebraCd<FT>::
is_solvable(const Linear_algebraCd<FT>::Matrix &M,
            const Linear_algebraCd<FT>::Vector &b) const
{
}

template < class FT >
bool
Linear_algebraCd<FT>::
homogeneous_linear_solver(const Linear_algebraCd<FT>::Matrix &M,
                          Linear_algebraCd<FT>::Vector &x) const
{
}

template < class FT >
int
Linear_algebraCd<FT>::
homogeneous_linear_solver(const Linear_algebraCd<FT>::Matrix &M,
                          Linear_algebraCd<FT>::Matrix &spanning_vectors) const
{
}

template < class FT >
int
Linear_algebraCd<FT>::
rank(const Linear_algebraCd<FT>::Matrix &M)
{
}

template < class FT >
int
Linear_algebraCd<FT>::
rank(const Linear_algebraCd<FT>::Matrix &M, std::vector<int> &q) const
{
}

CGAL_END_NAMESPACE

#ifndef CGAL_CARTESIAN_CLASS_DEFINED
#include <CGAL/Cartesian/Linear_algebra_d.C>
#endif 

#endif // CGAL_CARTESIAN_VECTOR_D_H
