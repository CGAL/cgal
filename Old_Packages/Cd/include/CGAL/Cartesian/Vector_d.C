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
// file          : include/CGAL/Cartesian/Vector_d.C
// revision      : $Revision$
// revision_date : $Date$
// author        : Herve Brönnimann
// coordinator   : INRIA Sophia-Antipolis (Mariette.Yvinec@sophia.inria.fr)
//
// ======================================================================

#ifndef CGAL_CARTESIAN_VECTOR_D_C
#define CGAL_CARTESIAN_VECTOR_D_C

#include <CGAL/Cartesian/Vector_d.h>
#include <CGAL/Cartesian/Direction_d.h>
#include <CGAL/Cartesian/d_utils.h>
#include <iterator>
#include <algorithm>
#include <numeric>

#ifndef CGAL_CTAG
#define CGAL_CTAG
#endif

#ifdef CGAL_CFG_TYPENAME_BUG
#define typename
#endif

CGAL_BEGIN_NAMESPACE

template < class R >
VectorCd<R CGAL_CTAG>::
VectorCd(int d)
{
  PTR = new _d_tuple<FT>(d);
}

template < class R >
VectorCd<R CGAL_CTAG>::
VectorCd(const VectorCd<R CGAL_CTAG> &v)
  : Handle(v)
{}

template < class R >
inline
VectorCd<R CGAL_CTAG>::
VectorCd(const typename VectorCd<R CGAL_CTAG>::Point_d &p)
  : Handle((const Handle&)p)
{
}

template < class R >
inline VectorCd<R CGAL_CTAG>::
VectorCd(const typename VectorCd<R CGAL_CTAG>::Direction_d &d)
  : Handle((const Handle&)d)
{
}

template < class R >
VectorCd<R CGAL_CTAG>::
VectorCd(int dim, const Null_vector &)
{
  CGAL_kernel_precondition( dim > 0);
  PTR = new _d_tuple<FT>(dim);
  std::fill(begin(), end(), FT(0));
}

template < class R >
VectorCd<R CGAL_CTAG>::~VectorCd()
{}

template < class R >
VectorCd<R CGAL_CTAG> &
VectorCd<R CGAL_CTAG>::operator=(const VectorCd<R CGAL_CTAG> &v)
{
  Handle::operator=((const Handle &)v);
  return *this;
}

template < class R >
bool
VectorCd<R CGAL_CTAG>::operator==(const VectorCd<R CGAL_CTAG> &v) const
{
  if (dimension() != v.dimension()) return false;
  if (ptr() == v.ptr()) return true; // identical
  return std::equal(begin(),end(),v.begin());
}

template < class R >
inline
bool
VectorCd<R CGAL_CTAG>::operator!=(const VectorCd<R CGAL_CTAG> &v) const
{
  return !(*this == v);
}

template < class R >
bool VectorCd<R CGAL_CTAG>::operator==(const Null_vector &) const
{
  const_iterator non_zero;
  non_zero = find_if(begin(),end(),bind1st(not_equal_to<FT>(),FT(0)));
  return non_zero == end();
}

template < class R >
inline bool VectorCd<R CGAL_CTAG>::operator!=(const Null_vector &v) const
{
  return !(*this == v);
}

template < class R >
inline
typename VectorCd<R CGAL_CTAG>::FT
VectorCd<R CGAL_CTAG>::cartesian(int i) const
{
  CGAL_kernel_precondition( (i>=0) && (i<dimension()) );
  return *(begin()+i);
}

template < class R >
inline
typename VectorCd<R CGAL_CTAG>::FT
VectorCd<R CGAL_CTAG>::operator[](int i) const
{
  return cartesian(i);
}

template < class R >
typename VectorCd<R CGAL_CTAG>::FT
VectorCd<R CGAL_CTAG>::homogeneous(int i) const
{
  CGAL_kernel_precondition( (i>=0) && (i<=dimension()) );
  return (i==dimension()) ? FT(1) : cartesian(i);
}

template < class R >
inline
VectorCd<R CGAL_CTAG>
VectorCd<R CGAL_CTAG>::operator+(const VectorCd<R CGAL_CTAG> &v) const
{
  CGAL_kernel_precondition( dimension() == v.dimension() );
  Self w(dimension());
  std::transform(begin(),end(),v.begin(),w.begin(),std::plus<FT>());
  return w;
}

template < class R >
inline
VectorCd<R CGAL_CTAG>
VectorCd<R CGAL_CTAG>::operator-(const VectorCd<R CGAL_CTAG> &v) const
{
  CGAL_kernel_precondition( dimension() == v.dimension() );
  Self w(dimension());
  std::transform(begin(),end(),v.begin(),w.begin(),std::minus<FT>());
  return w;
}

template < class R >
inline
VectorCd<R CGAL_CTAG>
VectorCd<R CGAL_CTAG>::operator-() const
{
  Self v(dimension());
  std::transform(begin(),end(),v.begin(),std::negate<FT>());
  return v;
}

template < class R >
inline
typename VectorCd<R CGAL_CTAG>::FT
VectorCd<R CGAL_CTAG>::operator*(const VectorCd<R CGAL_CTAG> &v) const
{
  CGAL_kernel_precondition( dimension() == v.dimension() );
  return std::inner_product(begin(),end(),v.begin(),FT(0));
}

template < class R >
inline
VectorCd<R CGAL_CTAG>
VectorCd<R CGAL_CTAG>::
operator*(const typename VectorCd<R CGAL_CTAG>::FT &c) const
{
  Self v(dimension());
  std::transform(begin(),end(),v.begin(),std::bind2nd(std::multiplies<FT>(),c));
  return v;
}

template < class R >
inline
VectorCd<R CGAL_CTAG>
VectorCd<R CGAL_CTAG>::
operator/(const typename VectorCd<R CGAL_CTAG>::FT &c) const
{
  Self v(dimension());
  std::transform(begin(),end(),v.begin(),std::bind2nd(std::divides<FT>(),c));
  return v;
}

template < class R >
inline
typename VectorCd<R CGAL_CTAG>::Direction_d
VectorCd<R CGAL_CTAG>::direction() const
{
  return Direction_d(*this);
}

template < class R >
VectorCd<R CGAL_CTAG>
VectorCd<R CGAL_CTAG>::
transform(const typename VectorCd<R CGAL_CTAG>::Aff_transformation_d &t) const
{
  return t.transform(*this);
}

#ifndef CGAL_CARTESIAN_NO_OSTREAM_INSERT_VECTORCD
template < class R >
std::ostream &
operator<<(std::ostream &os, const VectorCd<R CGAL_CTAG> &v)
{
  typedef typename VectorCd<R CGAL_CTAG>::FT FT;
  print_d<FT> prt(os);
  if (os.iword(IO::mode) == IO::PRETTY) os << "VectorCd(";
  prt(v.dimension());
  if (os.iword(IO::mode) == IO::PRETTY) { os << ", ("; prt.reset(); }
  std::for_each(v.begin(),v.end(),prt(os));
  if (os.iword(IO::mode) == IO::PRETTY) os << "))";
  return os;
}
#endif // CGAL_CARTESIAN_NO_OSTREAM_INSERT_VECTORCD

#ifndef CGAL_CARTESIAN_NO_ISTREAM_EXTRACT_VECTORCD
template < class R >
std::istream &
operator>>(std::istream &is, VectorCd<R CGAL_CTAG> &v)
{
  typedef typename VectorCd<R CGAL_CTAG>::FT FT;
  int dim;
  FT* w;
  switch(is.iword(IO::mode)) {
    case IO::ASCII :
    case IO::BINARY :
      is >> dim;
      w = new FT[dim];
      std::copy_n(std::istream_iterator<FT>(is),dim, w);
      break;
    default:
      std::cerr << "" << std::endl;
      std::cerr << "Stream must be in ascii or binary mode" << std::endl;
      break;
  }
  if (is)
      v = VectorCd<R CGAL_CTAG>(dim, w, w+dim);
  return is;
}
#endif // CGAL_CARTESIAN_NO_ISTREAM_EXTRACT_VECTORCD

CGAL_END_NAMESPACE

#ifdef CGAL_CFG_TYPENAME_BUG
#undef typename
#endif

#endif // CGAL_CARTESIAN_VECTOR_D_C
