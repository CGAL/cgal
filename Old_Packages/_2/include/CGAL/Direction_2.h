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

CGAL_BEGIN_NAMESPACE

template <class R_>
class Direction_2 : public R_::Direction_2_base
{
  typedef typename R_::RT                    RT;
  typedef typename R_::Vector_2              Vector_2;
  typedef typename R_::Line_2                Line_2;
  typedef typename R_::Ray_2                 Ray_2;
  typedef typename R_::Segment_2             Segment_2;
  typedef typename R_::Direction_2_base      RDirection_2;
public:
  typedef  R_   R;

  Direction_2()
  {}

  Direction_2(const CGAL::Direction_2<R> &d)
    : RDirection_2(static_cast<const RDirection_2&>(d)) {}

  Direction_2(const RDirection_2& d)
    : RDirection_2(d) {}

  Direction_2(const Vector_2& v)
    : RDirection_2(v) {}

  Direction_2(const Line_2& l)
    : RDirection_2(l) {}

  Direction_2(const Ray_2& r)
    : RDirection_2(r) {}

  Direction_2(const Segment_2& s)
    : RDirection_2(s) {}

  Direction_2(const RT &x, const RT &y)
    :  RDirection_2(x,y) {}
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
