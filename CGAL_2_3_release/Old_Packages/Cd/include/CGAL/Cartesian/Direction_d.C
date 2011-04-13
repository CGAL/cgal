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
// file          : include/CGAL/Cartesian/Direction_d.C
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Herve Brönnimann
// coordinator   : INRIA Sophia-Antipolis (Mariette.Yvinec@sophia.inria.fr)
//
// ======================================================================

#ifndef CGAL_CARTESIAN_DIRECTION_D_C
#define CGAL_CARTESIAN_DIRECTION_D_C

#include <CGAL/Cartesian/redefine_names_d.h>
#include <CGAL/Cartesian/Direction_d.h>
#include <CGAL/Cartesian/predicates_on_directions_d.h>
#include <CGAL/Cartesian/d_utils.h>

#ifndef CGAL_CARTESIAN_REDEFINE_NAMES_D_H
#define CGAL_CTAG
#endif

#ifdef CGAL_CFG_TYPENAME_BUG
#define typename
#endif

CGAL_BEGIN_NAMESPACE

template < class R >
DirectionCd<R CGAL_CTAG>::
DirectionCd(int d)
{
  PTR = new _d_tuple<RT>(d);
}

template < class R >
DirectionCd<R CGAL_CTAG>::
DirectionCd(const DirectionCd<R CGAL_CTAG> &d)
  : Handle((const Handle&)d)
{}

template < class R >
DirectionCd<R CGAL_CTAG>::
DirectionCd(const typename DirectionCd<R CGAL_CTAG>::Vector_d &v)
  : Handle((const Handle&)v)
{}

template < class R >
DirectionCd<R CGAL_CTAG>::~DirectionCd()
{
}

template < class R >
DirectionCd<R CGAL_CTAG> &
DirectionCd<R CGAL_CTAG>::operator=(const DirectionCd<R CGAL_CTAG> &d)
{
  Handle::operator=((const Handle &)d);
  return *this;
}

template < class R >
bool
DirectionCd<R CGAL_CTAG>::operator==(const DirectionCd<R CGAL_CTAG> &d) const
{
  if (dimension() != d.dimension()) return false;
  if (ptr() == d.ptr()) return true; // identical
  return equal_direction(*this,d);
}

template < class R >
inline
bool
DirectionCd<R CGAL_CTAG>::operator!=(const DirectionCd<R CGAL_CTAG> &d) const
{
  return !(*this == d);
}

template < class R >
inline
long
DirectionCd<R CGAL_CTAG>::id() const
{
  return (long) PTR;
}

template < class R >
inline
typename DirectionCd<R CGAL_CTAG>::Vector_d
DirectionCd<R CGAL_CTAG>::to_vector() const
{
  return Vector_d(*this);
}

template < class R >
inline
bool
DirectionCd<R CGAL_CTAG>::is_degenerate() const
{
  return to_vector() == Vector_d(dimension(), NULL_VECTOR);
}

template < class R >
inline
DirectionCd<R CGAL_CTAG>
DirectionCd<R CGAL_CTAG>::
transform
  (const typename DirectionCd<R CGAL_CTAG>::Aff_transformation_d &t) const
{
  return t.transform(*this);
}

template < class R >
inline
DirectionCd<R CGAL_CTAG> 
DirectionCd<R CGAL_CTAG>::operator-() const
{
  Self v(dimension());
  std::transform(begin(),end(),v.begin(),std::negate<RT>());
  return v;
}

template < class R >
typename DirectionCd<R CGAL_CTAG>::RT 
DirectionCd<R CGAL_CTAG>::delta(int i) const
{
  CGAL_kernel_precondition( i >= 0 && i < dimension() );
  return *(begin()+i);
}

#ifndef CGAL_NO_OSTREAM_INSERT_DIRECTIONCd
template < class R >
std::ostream &operator<<(std::ostream &os, const DirectionCd<R CGAL_CTAG> &d)
{
  typedef typename DirectionCd<R CGAL_CTAG>::RT RT;
  typename DirectionCd<R CGAL_CTAG>::Vector_d v = d.to_vector();
  switch(os.iword(IO::mode)) {
    case IO::ASCII :
    case IO::BINARY :
      std::for_each(d.begin(),d.end(),print_d<RT>(os));
      break;
    default:
      os << "DirectionCd(" << d.dimension() << ", ";
      std::for_each(v.begin(),d.end(),print_d<RT>(os));
      os << ")";
  }
  return os;
}
#endif // CGAL_NO_OSTREAM_INSERT_DIRECTIONCd

#ifndef CGAL_NO_ISTREAM_EXTRACT_DIRECTIONCd
template < class R >
std::istream &operator>>(std::istream &is, DirectionCd<R CGAL_CTAG> &d)
{
  typedef typename DirectionCd<R CGAL_CTAG>::RT RT;
  int dim;
  RT *v;
  switch(is.iword(IO::mode)) {
    case IO::ASCII :
    case IO::BINARY :
      is >> dim;
      v = new RT[dim];
      std::copy_n(istream_iterator<RT>(is),dim, v);
      break;
    default:
      std::cerr << "" << std::endl;
      std::cerr << "Stream must be in ascii or binary mode" << std::endl;
      break;
  }
  if (is)
      d = DirectionCd<R CGAL_CTAG>(dim, v, v+dim);
  return is;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_DIRECTIONCd

CGAL_END_NAMESPACE

#ifdef CGAL_CFG_TYPENAME_BUG
#undef typename
#endif

#endif // CGAL_CARTESIAN_DIRECTION_D_C
