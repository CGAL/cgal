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
// release_date  : 2000, July 30
// 
// source        : Triangle_2.fw
// file          : Triangle_2.h
// package       : _2 (3.6)
// maintainer    : Stefan Schirra <stschirr@mpi-sb.mpg.de>
// revision      : 3.6
// revision_date : 30 Jul 2000 
// author(s)     : Andreas Fabri
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ======================================================================
 

#ifndef CGAL_TRIANGLE_2_H
#define CGAL_TRIANGLE_2_H

#ifndef CGAL_REP_CLASS_DEFINED
#error  no representation class defined
#endif  // CGAL_REP_CLASS_DEFINED


#ifdef CGAL_HOMOGENEOUS_H
#ifndef CGAL_TRIANGLEH2_H
#include <CGAL/TriangleH2.h>
#endif // CGAL_TRIANGLEH2_H
#endif // CGAL_HOMOGENEOUS_H

#ifdef CGAL_CARTESIAN_H
#ifndef CGAL_TRIANGLEC2_H
#include <CGAL/Cartesian/Triangle_2.h>
#endif // CGAL_TRIANGLEC2_H
#endif // CGAL_CARTESIAN_H

#ifdef CGAL_SIMPLE_CARTESIAN_H
#include <CGAL/SimpleCartesian/TriangleS2.h>
#endif // CGAL_SIMPLE_CARTESIAN_H


#ifndef CGAL_POINT_2_H
#include <CGAL/Point_2.h>
#endif // CGAL_POINT_2_H

CGAL_BEGIN_NAMESPACE

template <class R_>
  class Triangle_2 : public R_::Triangle_2_base
{
public:
  typedef  R_                          R;
  typedef typename R::RT               RT;
  typedef typename R::FT               FT;
  typedef typename R::Triangle_2_base  RTriangle_2;

  Triangle_2() : RTriangle_2() {}

  Triangle_2(const CGAL::Triangle_2<R> &t) : RTriangle_2((RTriangle_2&)t) {}

  Triangle_2(const RTriangle_2& t) : RTriangle_2(t) {}

  Triangle_2(const CGAL::Point_2<R> &p,
             const CGAL::Point_2<R> &q,
             const CGAL::Point_2<R> &r) : RTriangle_2(p,q,r) {}


  bool
  operator==(const CGAL::Triangle_2<R> &t) const
  { return RTriangle_2::operator==(t); }

  bool
  operator!=(const CGAL::Triangle_2<R> &t) const
  { return !(*this == t); }

  CGAL::Point_2<R>
  vertex(int i) const
  { return RTriangle_2::vertex(i); }

  CGAL::Point_2<R>
  operator[](int i) const
  { return vertex(i); }

  CGAL::Triangle_2<R>
  transform(const CGAL::Aff_transformation_2<R> &t) const
  { return  RTriangle_2::transform(t); }

  CGAL::Triangle_2<R>  opposite() const
  { return  CGAL::Triangle_2<R>(vertex(0), vertex(2), vertex(1)); }

  Orientation
  orientation() const
  { return RTriangle_2::orientation(); }

  Oriented_side
  oriented_side(const CGAL::Point_2<R> &p) const
  { return RTriangle_2::oriented_side(p); }

  Bounded_side
  bounded_side(const CGAL::Point_2<R> &p) const
  { return RTriangle_2::bounded_side(p); }

  bool
  has_on_boundary(const CGAL::Point_2<R> &p) const
  { return RTriangle_2::has_on_boundary(p); }

  bool
  has_on_positive_side(const CGAL::Point_2<R> &p) const
  { return RTriangle_2::has_on_positive_side(p); }

  bool
  has_on_negative_side(const CGAL::Point_2<R> &p) const
  { return RTriangle_2::has_on_negative_side(p); }

  bool
  has_on_unbounded_side(const CGAL::Point_2<R> &p) const
  { return RTriangle_2::has_on_unbounded_side(p); }

  bool
  has_on_bounded_side(const CGAL::Point_2<R> &p) const
  { return RTriangle_2::has_on_bounded_side(p); }

  bool
  is_degenerate() const
  { return RTriangle_2::is_degenerate(); }

  Bbox_2
  bbox() const
  { return vertex(0).bbox() + vertex(1).bbox() + vertex(2).bbox(); }
};

#ifndef NO_OSTREAM_INSERT_TRIANGLE_2
template < class R >
std::ostream &
operator<<(std::ostream &os, const Triangle_2<R> &t)
{
  typedef typename  R::Triangle_2_base  RTriangle_2;
  return os << (const RTriangle_2&)t;
}
#endif // NO_OSTREAM_INSERT_TRIANGLE_2

#ifndef NO_ISTREAM_EXTRACT_TRIANGLE_2
template < class R >
std::istream &
operator>>(std::istream &is, Triangle_2<R> &t)
{
  typedef typename  R::Triangle_2_base  RTriangle_2;
  return is >> (RTriangle_2&)t;
}
#endif // NO_ISTREAM_EXTRACT_TRIANGLE_2

CGAL_END_NAMESPACE


#endif // CGAL_TRIANGLE_2_H
