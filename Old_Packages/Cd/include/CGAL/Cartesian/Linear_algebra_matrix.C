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
// file          : include/CGAL/Cartesian/Linear_algebra_matrix.C
// revision      : $Revision$
// revision_date : $Date$
// author        : Herve Brönnimann
// coordinator   : INRIA Sophia-Antipolis (Mariette.Yvinec@sophia.inria.fr)
//
// ======================================================================

#ifndef CGAL_CARTESIAN_LINEAR_ALGEBRA_MATRIX_D_C
#define CGAL_CARTESIAN_LINEAR_ALGEBRA_MATRIX_D_C

#include <CGAL/Cartesian/d_utils.h>
#include <CGAL/Cartesian/Linear_algebra_matrix.h>
#include <iterator>
#include <algorithm>
#include <numeric>

#ifdef CGAL_CFG_TYPENAME_BUG
#define typename
#endif

CGAL_BEGIN_NAMESPACE

template < class LA >
inline
LA_matrixCd<LA>::LA_matrixCd(int n)
{
  new_rep(n,n);
}

template < class LA >
inline
LA_matrixCd<LA>::LA_matrixCd(int n, int m)
{
  new_rep(n,m);
}

template < class LA >
inline
LA_matrixCd<LA>::LA_matrixCd(const std::pair<int,int> &d)
{
  new_rep(d.first,d.second);
}

template < class LA >
inline
LA_matrixCd<LA>::LA_matrixCd(const LA_matrixCd<LA> &v)
  : Handle(v), _rdim(v.row_dimension()), _cdim(v.column_dimension())
{}

template < class LA >
inline
LA_matrixCd<LA>::~LA_matrixCd()
{}

template < class LA >
inline
LA_matrixCd<LA> &
LA_matrixCd<LA>::operator=(const LA_matrixCd<LA> &v)
{
  Handle::operator=((const Handle &)v);
  _rdim = v.row_dimension();
  _cdim = v.column_dimension();
  return *this;
}

template < class LA >
inline
bool
LA_matrixCd<LA>::operator==(const LA_matrixCd<LA> &M) const
{
  if (dimension() != M.dimension()) return false;
  if (ptr() == M.ptr()) return true; // identical
  return std::equal(begin(),end(),M.begin());
}

template < class LA >
inline
bool
LA_matrixCd<LA>::operator!=(const LA_matrixCd<LA> &v) const
{
  return !(*this == v);
}

template < class LA >
inline
std::pair<int,int>
LA_matrixCd<LA>::dimension() const
{
  return std::make_pair(row_dimension(),column_dimension());
}

template < class LA >
inline
typename LA_matrixCd<LA>::Vector
LA_matrixCd<LA>::row(int i) const
{
  CGAL_kernel_precondition( i>=0 && i<row_dimension());
  return Vector(row_begin(i), row_end(i));
}

template < class LA >
inline
typename LA_matrixCd<LA>::Vector
LA_matrixCd<LA>::column(int j) const
{
  CGAL_kernel_precondition( j>=0 && j<column_dimension());
  Vector v(row_dimension()); // beware: row_dimension == *number* of raws
  // We could avoid all this if we had an adapter for iterator with step
  // Then we would just use std::copy_n()... what a shame!
  const_iterator i;
  typename Vector::iterator k;
  for (i=column_begin(j), k=v.begin();
       i!=column_end(j);
       ++k, i+=column_dimension())
    *k = *i;
  return v;
}

template < class LA >
LA_matrixCd<LA>
LA_matrixCd<LA>::operator+(const LA_matrixCd<LA> &M) const
{
  CGAL_kernel_precondition( dimension() == M.dimension() );
  Self w(dimension());
  std::transform(begin(),end(),M.begin(),w.begin(),std::plus<FT>());
  return w;
}

template < class LA >
LA_matrixCd<LA>
LA_matrixCd<LA>::operator-(const LA_matrixCd<LA> &M) const
{
  CGAL_kernel_precondition( dimension() == M.dimension() );
  Self w(dimension());
  std::transform(begin(),end(),M.begin(),w.begin(),std::minus<FT>());
  return w;
}

template < class LA >
LA_matrixCd<LA>
LA_matrixCd<LA>::operator-() const
{
  Self v(dimension());
  std::transform(begin(),end(),v.begin(),std::negate<FT>());
  return v;
}

template < class LA >
inline
LA_matrixCd<LA>
LA_matrixCd<LA>::operator*(const typename LA_matrixCd<LA>::FT &c) const
{
  Self v(dimension());
  std::transform(begin(),end(),v.begin(),std::bind2nd(std::multiplies<FT>(),c));
  return v;
}

template < class LA >
inline
LA_matrixCd<LA>
LA_matrixCd<LA>::operator/(const typename LA_matrixCd<LA>::FT &c) const
{
  Self v(dimension());
  std::transform(begin(),end(),v.begin(),std::bind2nd(std::divides<FT>(),c));
  return v;
}

template < class LA >
LA_matrixCd<LA>
LA_matrixCd<LA>::operator*(const LA_matrixCd<LA> &M) const
{
  CGAL_kernel_precondition( column_dimension() == M.row_dimension() );
  Self w( row_dimension(), M.column_dimension() );
  // TODO: Iterator-based matrix-multiplication
  int i, j, k;
  for (i=0; i<row_dimension(); ++i)
    for (j=0; j<M.column_dimension(); ++j) {
      w[i][j] = FT(0);
      for (k=0; k<column_dimension(); ++k)
        w[i][j] += (*this)[i][k] * M[k][j];
    }
  return w;
}

template < class LA >
typename LA_matrixCd<LA>::Vector
LA_matrixCd<LA>::operator*(const typename LA_matrixCd<LA>::Vector &v) const
{
  CGAL_kernel_precondition( column_dimension() == v.dimension() );
  Vector w( row_dimension() ); 
  const_iterator i;
  typename Vector::iterator j;
  for (j=w.begin(), i=begin(); j!=w.end(); ++j,i+=column_dimension())
    *j = std::inner_product(i, i+column_dimension(), v.begin(), FT(0));
  return w;
}

template < class LA >
void
LA_matrixCd<LA>::swap_columns(int j, int k)
{
  iterator jit, kit;
  for (jit=column_begin(j),kit=column_begin(k);
       jit != column_end(j);
       jit+=column_dimension(),kit+=column_dimension())
    std::swap(*jit,*kit);
}

#ifndef CGAL_CARTESIAN_NO_OSTREAM_INSERT_LA_VECTORCD
template < class LA >
std::ostream &
operator<<(std::ostream &os, const LA_matrixCd<LA> &M)
{
  typedef typename LA::FT FT;
  print_d<FT> prt(&os);
  if (os.iword(IO::mode)==IO::PRETTY) os << "LA_Matrix(";
  prt(M.row_dimension());
  prt(M.column_dimension());
  if (os.iword(IO::mode)==IO::PRETTY) { os << ", ("; prt.reset(); }
  std::for_each(M.begin(),M.end(),prt);
  if (os.iword(IO::mode)==IO::PRETTY) os << "))";
  return os;
}
#endif // CGAL_CARTESIAN_NO_OSTREAM_INSERT_LA_VECTORCD

#ifndef CGAL_CARTESIAN_NO_ISTREAM_EXTRACT_LA_VECTORCD
template < class LA >
std::istream &
operator>>(std::istream &is, LA_matrixCd<LA> &M)
{
  typedef typename LA_matrixCd<LA>::FT FT;
  int cdim, rdim, dim;
  FT* w;
  switch(is.iword(IO::mode)) {
    case IO::ASCII :
    case IO::BINARY :
      is >> rdim >> cdim;
      dim = cdrim*rdim;
      M = LA_matrixCd<LA>(rdim, cdim);
      std::copy_n(std::istream_iterator<FT>(is),dim, M.begin());
      break;
    default:
      std::cerr << "" << std::endl;
      std::cerr << "Stream must be in ascii or binary mode" << std::endl;
      break;
  }
  return is;
}
#endif // CGAL_CARTESIAN_NO_ISTREAM_EXTRACT_LA_VECTORCD

CGAL_END_NAMESPACE

#ifdef CGAL_CFG_TYPENAME_BUG
#undef typename
#endif

#endif // CGAL_CARTESIAN_MATRIX_D_C
