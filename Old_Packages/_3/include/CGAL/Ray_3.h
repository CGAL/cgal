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
// file          : Ray_3.h
// package       : _3
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas Fabri
//                 Stefan Schirra
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ======================================================================
 
#ifndef CGAL_RAY_3_H
#define CGAL_RAY_3_H

#include <CGAL/Point_3.h>
#include <CGAL/Direction_3.h>
#include <CGAL/Line_3.h>

CGAL_BEGIN_NAMESPACE

template <class R_>
class Ray_3 : public R_::Ray_3_base
{
public:
  typedef          R_                       R;
  typedef typename R::RT                    RT;
  typedef typename R::FT                    FT;
  typedef typename R::Ray_3_base  RRay_3;

  Ray_3() : RRay_3()
  {}

  Ray_3(const CGAL::Ray_3<R>& r) : RRay_3(r)
  {}

  Ray_3(const RRay_3&  r) : RRay_3(r)
  {}

  Ray_3(const CGAL::Point_3<R>& sp, const CGAL::Point_3<R>& secondp)
    : RRay_3(sp, secondp)
  {}

  Ray_3(const CGAL::Point_3<R>& sp, const CGAL::Direction_3<R>& d)
    : RRay_3(sp, d)
  {}

  CGAL::Point_3<R>     start() const
  { return RRay_3::start(); }
  CGAL::Point_3<R>     source() const
  { return RRay_3::source(); }
  CGAL::Point_3<R>     second_point() const
  { return RRay_3::second_point(); }
  CGAL::Point_3<R>     point(int i) const
  { return RRay_3::point(i); }
  CGAL::Direction_3<R> direction() const
  { return RRay_3::direction(); }
  CGAL::Line_3<R>      supporting_line() const
  { return RRay_3::supporting_line(); }
  CGAL::Ray_3<R>       opposite() const
  { return RRay_3::opposite(); }
  CGAL::Ray_3<R>       transform(const CGAL::Aff_transformation_3<R>& t) const
  { return RRay_3::transform(t); }
};

#ifndef CGAL_NO_OSTREAM_INSERT_RAY_3
template < class R >
std::ostream&
operator<<(std::ostream& os, const Ray_3<R>& r)
{
  typedef typename  R::Ray_3_base  RRay_3;
  return os << static_cast<const RRay_3&>(r);
}
#endif // CGAL_NO_OSTREAM_INSERT_RAY_3

#ifndef CGAL_NO_ISTREAM_EXTRACT_RAY_3
template < class R >
std::istream&
operator>>(std::istream& is, Ray_3<R>& r)
{
  typedef typename  R::Ray_3_base  RRay_3;
  return is >> static_cast<RRay_3&>(r);
}
#endif // CGAL_NO_ISTREAM_EXTRACT_RAY_3

CGAL_END_NAMESPACE

#endif // CGAL_RAY_3_H
