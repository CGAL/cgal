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
// file          : Ray_2.h
// package       : _2
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas Fabri
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ======================================================================

#ifndef CGAL_RAY_2_H
#define CGAL_RAY_2_H

#include <CGAL/Point_2.h>

CGAL_BEGIN_NAMESPACE

template <class R_>
class Ray_2 : public R_::Ray_2_base
{
public:
  typedef  R_   R;
  typedef typename R::RT                    RT;
  typedef typename R::FT                    FT;
  typedef typename R::Ray_2_base  RRay_2;

  Ray_2()
    : RRay_2()
  {}

  Ray_2(const CGAL::Ray_2<R> &r)
    : RRay_2(static_cast<const RRay_2&>(r))
  {}

  Ray_2(const RRay_2& r)
    : RRay_2(r)
  {}

  Ray_2(const CGAL::Point_2<R> &sp, const CGAL::Point_2<R> &secondp)
    : RRay_2(sp, secondp)
  {}

  Ray_2(const CGAL::Point_2<R> &sp, const CGAL::Direction_2<R> &d)
    : RRay_2(sp, d)
  {}

  CGAL::Point_2<R> start() const
  { return RRay_2::start(); }

  CGAL::Point_2<R> source() const
  { return RRay_2::source(); }

  CGAL::Point_2<R> second_point() const
  { return RRay_2::second_point(); }

  CGAL::Point_2<R> point(int i) const
  { return RRay_2::point(i); }

  CGAL::Direction_2<R> direction() const
  { return RRay_2::direction(); }

  CGAL::Line_2<R> supporting_line() const
  { return RRay_2::supporting_line(); }

  CGAL::Ray_2<R> opposite() const
  { return RRay_2::opposite(); }

  CGAL::Ray_2<R> transform(const CGAL::Aff_transformation_2<R> &t) const
  { return RRay_2::transform(t); }
};

#ifndef CGAL_NO_OSTREAM_INSERT_RAY_2
template < class R >
std::ostream &
operator<<(std::ostream &os, const Ray_2<R> &r)
{
  typedef typename  R::Ray_2_base  RRay_2;
  return os << static_cast<const RRay_2&>(r);
}
#endif // CGAL_NO_OSTREAM_INSERT_RAY_2

#ifndef CGAL_NO_ISTREAM_EXTRACT_RAY_2
template < class R >
std::istream &
operator>>(std::istream &is, Ray_2<R> &r)
{
  typedef typename  R::Ray_2_base  RRay_2;
  return is >> static_cast<RRay_2&>(r);
}
#endif // CGAL_NO_ISTREAM_EXTRACT_RAY_2

CGAL_END_NAMESPACE

#endif  // CGAL_RAY_2_H
