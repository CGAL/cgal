// ======================================================================
//
// Copyright (c) 1999 The CGAL Consortium
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
// file          : Direction_2.h
// package       : _2
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ======================================================================
 

#ifndef CGAL_DIRECTION_2_H
#define CGAL_DIRECTION_2_H

#ifndef CGAL_REP_CLASS_DEFINED
#error  no representation class defined
#endif

#ifdef CGAL_HOMOGENEOUS_H
#include <CGAL/DirectionH2.h>
#endif

#if defined CGAL_CARTESIAN_H || defined CGAL_SIMPLE_CARTESIAN_H
#include <CGAL/Cartesian/Direction_2.h>
#endif

#include <CGAL/Vector_2.h>

CGAL_BEGIN_NAMESPACE

template <class R_>
class Direction_2 : public R_::Direction_2_base
{
public:
  typedef  R_   R;
  typedef typename R::RT                    RT;
  typedef typename R::FT                    FT;
  typedef typename R::Vector_2_base  RVector_2;
  typedef typename R::Direction_2_base  RDirection_2;
  // typedef typename R::Aff_transformation_2  Aff_transformation_2;

  Direction_2()
  {}

  Direction_2(const CGAL::Direction_2<R> &d)
    : RDirection_2(static_cast<const RDirection_2&>(d))
  {}


  Direction_2(const RDirection_2& d)
    : RDirection_2(d)
  {}


  Direction_2(const RVector_2& v)
    : RDirection_2(v)
  {}

  Direction_2(const RT &x, const RT &y)
    :  RDirection_2(x,y)
  {}

  bool
  operator==(const CGAL::Direction_2<R> &d) const
  { return RDirection_2::operator==(d); }

  bool
  operator!=(const CGAL::Direction_2<R> &d) const
  { return !(*this == d); }

  bool
  operator>=(const CGAL::Direction_2<R> &d) const
  { return RDirection_2::operator>=(d); }

  bool
  operator<=(const CGAL::Direction_2<R> &d) const
  { return RDirection_2::operator<=(d); }

  bool
  operator>(const CGAL::Direction_2<R> &d) const
  { return RDirection_2::operator>(d); }

  bool
  operator<(const CGAL::Direction_2<R> &d) const
  { return RDirection_2::operator<(d); }

  bool
  counterclockwise_in_between(const CGAL::Direction_2<R> &d1,
                                   const CGAL::Direction_2<R> &d2) const
  { return RDirection_2::counterclockwise_in_between(d1,d2); }

  CGAL::Vector_2<R>
  vector() const
  { return static_cast<CGAL::Vector_2<R> >(RDirection_2::to_vector()); }

  CGAL::Vector_2<R>
  to_vector() const
  { return static_cast<CGAL::Vector_2<R> >(RDirection_2::to_vector()); }

  CGAL::Direction_2<R>
  transform(const CGAL::Aff_transformation_2<R> &t) const
  { return RDirection_2::transform(t); }

  CGAL::Direction_2<R>
  operator-() const
  { return RDirection_2::operator-(); }

  RT
  delta(int i) const
  { return RDirection_2::delta(i); }

  RT
  dx() const
  { return RDirection_2::dx(); }

  RT
  dy() const
  { return RDirection_2::dy(); }
};


#ifndef CGAL_NO_OSTREAM_INSERT_DIRECTION_2
template < class R >
std::ostream &
operator<<(std::ostream &os, const Direction_2<R> &d)
{
  typedef typename  R::Direction_2_base  RDirection_2;
  return os << static_cast<const RDirection_2&>(d);
}
#endif // CGAL_NO_OSTREAM_INSERT_DIRECTION_2

#ifndef CGAL_NO_ISTREAM_EXTRACT_DIRECTION_2
template < class R >
std::istream &
operator>>(std::istream &is, Direction_2<R> &p)
{
  typedef typename  R::Direction_2_base  RDirection_2;
  return is >> static_cast<RDirection_2&>(p);
}
#endif // CGAL_NO_ISTREAM_EXTRACT_DIRECTION_2

CGAL_END_NAMESPACE

#endif // CGAL_DIRECTION_2_H
