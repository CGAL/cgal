// ============================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
// 
// release       :
// release_date  :
// 
// source        : Ray_3.fw
// file          : Ray_3.h
// revision      : 2.4
// revision_date : 24 Aug 1999 
// author(s)     : Andreas Fabri
//                 Stefan Schirra
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ============================================================================
 

#ifndef CGAL_RAY_3_H
#define CGAL_RAY_3_H

#ifndef CGAL_REP_CLASS_DEFINED
#error  no representation class defined
#endif  // CGAL_REP_CLASS_DEFINED

#ifdef CGAL_HOMOGENEOUS_H
#ifndef CGAL_RAYH3_H
#include <CGAL/RayH3.h>
#endif // CGAL_RAYH3_H
#endif // CGAL_HOMOGENEOUS_H

#ifdef CGAL_CARTESIAN_H
#ifndef CGAL_RAYC3_H
#include <CGAL/RayC3.h>
#endif // CGAL_RAYC3_H
#endif // CGAL_CARTESIAN_H

CGAL_BEGIN_NAMESPACE

template <class _R>
class Ray_3 : public _R::Ray_3_base
{
public:
  typedef          _R                       R;
  typedef typename R::RT                    RT;
  typedef typename R::FT                    FT;
  typedef typename R::Ray_3_base  RRay_3;

  Ray_3() : RRay_3()
  {}
  Ray_3(const Ray_3<R>& r) : RRay_3(r)
  {}
  Ray_3(const RRay_3&  r) : RRay_3(r)
  {}
  Ray_3(const Point_3<R>& sp,
            const Point_3<R>& secondp)
    : RRay_3(sp, secondp)
  {}
  Ray_3(const Point_3<R>& sp,
            const Direction_3<R>& d)
    : RRay_3(sp, d)
  {}

  Ray_3<R>&      operator=(const Ray_3<R>& r)
  {
      RRay_3::operator=(r);
      return *this;
  }
  bool                operator==(const Ray_3<R>& r) const
  { return RRay_3::operator==(r); }
  bool                operator!=(const Ray_3<R>& r) const
  { return !(*this == r); }

  int                 id() const  /* XXX */
  { return (int)  PTR ; }

  Point_3<R>     start() const
  { return RRay_3::start(); }
  Point_3<R>     source() const
  { return RRay_3::source(); }
  Point_3<R>     second_point() const
  { return RRay_3::second_point(); }
  Point_3<R>     point(int i) const
  { return RRay_3::point(i); }
  Direction_3<R> direction() const
  { return RRay_3::direction(); }
  Line_3<R>      supporting_line() const
  { return RRay_3::supporting_line(); }
  Ray_3<R>       opposite() const
  { return RRay_3::opposite(); }
  Ray_3<R>       transform(const Aff_transformation_3<R>& t) const
  { return RRay_3::transform(t); }
  bool                is_degenerate() const
  { return RRay_3::is_degenerate(); }
  bool                has_on(const Point_3<R>& p) const
  { return RRay_3::has_on(p); }
};

#ifndef NO_OSTREAM_INSERT_RAY_3
template < class R >
std::ostream&
operator<<(std::ostream& os, const Ray_3<R>& r)
{
  typedef typename  R::Ray_3_base  RRay_3;
  return os << (const RRay_3& )r;
}
#endif // NO_OSTREAM_INSERT_RAY_3

#ifndef NO_ISTREAM_EXTRACT_RAY_3
template < class R >
std::istream&
operator>>(std::istream& is, Ray_3<R>& r)
{
  typedef typename  R::Ray_3_base  RRay_3;
  return is >> (RRay_3& )r;
}
#endif // NO_ISTREAM_EXTRACT_RAY_3


CGAL_END_NAMESPACE


#ifndef CGAL_LINE_3_H
#include <CGAL/Line_3.h>
#endif // CGAL_LINE_3_H

#endif // CGAL_RAY_3_H
