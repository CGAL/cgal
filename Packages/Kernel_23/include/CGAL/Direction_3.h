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
// author(s)     : Andreas Fabri, Stefan Schirra
//
// coordinator   : MPI, Saarbruecken
// ======================================================================
 
#ifndef CGAL_DIRECTION_3_H
#define CGAL_DIRECTION_3_H

CGAL_BEGIN_NAMESPACE

template <class R_>
class Direction_3 : public R_::Direction_3_base
{
  typedef typename R_::RT                    RT;
  typedef typename R_::Vector_3              Vector_3;
  typedef typename R_::Line_3                Line_3;
  typedef typename R_::Ray_3                 Ray_3;
  typedef typename R_::Segment_3             Segment_3;
  typedef typename R_::Direction_3_base  RDirection_3;
public:
  typedef          R_                       R;

  Direction_3()
  {}

  Direction_3(const CGAL::Direction_3<R>& d)
    : RDirection_3( static_cast<const RDirection_3&>(d) ) {}

  Direction_3(const RDirection_3& d)
    : RDirection_3(d) {}

  Direction_3(const Vector_3& v)
    : RDirection_3(v) {}

  Direction_3(const Line_3& l)
    : RDirection_3(l) {}

  Direction_3(const Ray_3& r)
    : RDirection_3(r) {}

  Direction_3(const Segment_3& s)
    : RDirection_3(s) {}

  Direction_3(const RT& hx, const RT& hy, const RT& hz)
    : RDirection_3(hx, hy, hz) {}
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
