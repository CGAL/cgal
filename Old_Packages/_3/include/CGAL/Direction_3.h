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
// file          : Direction_3.h
// package       : _3
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas Fabri
//                 Stefan Schirra
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ======================================================================
 
#ifndef CGAL_DIRECTION_3_H
#define CGAL_DIRECTION_3_H

#include <CGAL/Vector_3.h>
#include <CGAL/Line_3.h>
#include <CGAL/Ray_3.h>
#include <CGAL/Segment_3.h>

CGAL_BEGIN_NAMESPACE

template <class R_>
class Direction_3 : public R_::Direction_3_base
{
public:
  typedef          R_                       R;
  typedef typename R::RT                    RT;
  typedef typename R::FT                    FT;
  typedef typename R::Direction_3_base  RDirection_3;
  typedef typename R::Vector_3_base  RVector_3;
  typedef typename R::Line_3_base           RLine_3;
  typedef typename R::Ray_3_base            RRay_3;
  typedef typename R::Segment_3_base        RSegment_3;

  Direction_3()
  {}

  Direction_3(const CGAL::Direction_3<R>& d)
    : RDirection_3( static_cast<const RDirection_3&>(d) )
  {}

  Direction_3(const RDirection_3&  d)
    : RDirection_3(d)
  {}

  Direction_3(const RVector_3&  v)
    : RDirection_3(v)
  {}

  Direction_3(const RLine_3& l)
    : RDirection_3(l)
  {}

  Direction_3(const RRay_3& r)
    : RDirection_3(r)
  {}

  Direction_3(const RSegment_3& s)
    : RDirection_3(s)
  {}

  Direction_3(const RT& hx, const RT& hy, const RT& hz)
    : RDirection_3(hx, hy, hz)
  {}

  CGAL::Vector_3<R> vector() const
  { return static_cast<CGAL::Vector_3<R> >(RDirection_3::to_vector()); }

  CGAL::Vector_3<R> to_vector() const
  { return static_cast<CGAL::Vector_3<R> >(RDirection_3::to_vector()); }

  CGAL::Direction_3<R> transform(const CGAL::Aff_transformation_3<R> & t) const
  { return RDirection_3::transform(t); }

  CGAL::Direction_3<R> operator-() const
  { return RDirection_3::operator-(); }
};

#ifndef CGAL_NO_OSTREAM_INSERT_DIRECTION_3
template < class R >
std::ostream& operator<<(std::ostream& os, const Direction_3<R>& d)
{
  typedef typename  R::Direction_3_base  RDirection_3;
  return os << static_cast<const RDirection_3&>(d); }
#endif // CGAL_NO_OSTREAM_INSERT_DIRECTION_3

#ifndef CGAL_NO_ISTREAM_EXTRACT_DIRECTION_3
template < class R >
std::istream& operator>>(std::istream& is, Direction_3<R>& p)
{
  typedef typename  R::Direction_3_base  RDirection_3;
  return is >> static_cast<RDirection_3&>(p); }
#endif // CGAL_NO_ISTREAM_EXTRACT_DIRECTION_3

CGAL_END_NAMESPACE

#endif // CGAL_DIRECTION_3_H
