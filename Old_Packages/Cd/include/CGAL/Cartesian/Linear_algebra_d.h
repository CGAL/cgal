// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Herve.Bronnimann@sophia.inria.fr

#ifndef CGAL_CARTESIAN_LINEAR_ALGEBRA_D_H
#define CGAL_CARTESIAN_LINEAR_ALGEBRA_D_H

#include <memory>
#include <vector>
#include <CGAL/Cartesian/redefine_names_d.h>
#include <CGAL/Cartesian/Linear_algebra_vector.h>
#include <CGAL/Cartesian/Linear_algebra_matrix.h>

CGAL_BEGIN_NAMESPACE

template < class _FT >
class Linear_algebraCd
{
public:
  typedef _FT                              FT;
  typedef _FT                              RT;
  typedef Linear_algebraCd<FT>             Self;
  typedef LA_vectorCd<Self>                Vector;
  typedef LA_matrixCd<Self>                Matrix;
  typedef const FT*                        const_iterator;
  typedef FT*                              iterator;
  
  Matrix transpose(const Matrix &M) const;

  bool   inverse(const Matrix &M, Matrix &I, FT &D, Vector &c) const;
  Matrix inverse(const Matrix M, RT &D) const;
  
  FT     determinant(const Matrix &M, Matrix &L, Matrix &U,
	     std::vector<int> &q, Vector &c) const;
  bool   verify_determinant(const Matrix &M,
             const Matrix &L, const Matrix &U, const RT &D,
             const std::vector<int> &q, const Vector &c) const;
  FT     determinant(const Matrix &M) const;
  Sign   sign_of_determinant(const Matrix &M) const;
 
  bool   linear_solver(const Matrix &M, const Vector &b,
             Vector &x, RT &D, Matrix &spanning_vectors, Vector &c) const;
  bool   linear_solver(const Matrix &M, const Vector &b,
             Vector &x, RT &D, Vector &c) const;
  bool   linear_solver(const Matrix &M, const Vector &b,
             Vector &x, RT &D) const;
  bool   is_solvable(const Matrix &M, const Vector &b) const;

  bool   homogeneous_linear_solver(const Matrix &M, Vector &x) const;
  int    homogeneous_linear_solver(const Matrix &M,
             Matrix &spanning_vectors) const;

  int    rank(const Matrix &M);
  int    rank(const Matrix &M, std::vector<int> &q) const;
};

CGAL_END_NAMESPACE

#ifndef CGAL_CARTESIAN_CLASS_DEFINED
#include <CGAL/Cartesian/Linear_algebra_d.C>
#endif 

#endif // CGAL_CARTESIAN_VECTOR_D_H
