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

CGAL_BEGIN_NAMESPACE

template <class R_>
class Ray_2 : public R_::Ray_2_base
{
  typedef typename R_::RT                    RT;
  typedef typename R_::Point_2               Point_2;
  typedef typename R_::Direction_2           Direction_2;
  typedef typename R_::Ray_2_base  RRay_2;
public:
  typedef  R_   R;

  Ray_2()
    : RRay_2() {}

  Ray_2(const CGAL::Ray_2<R> &r)
    : RRay_2(static_cast<const RRay_2&>(r)) {}

  Ray_2(const RRay_2& r)
    : RRay_2(r) {}

  Ray_2(const Point_2 &sp, const Point_2 &secondp)
    : RRay_2(sp, secondp) {}

  Ray_2(const Point_2 &sp, const Direction_2 &d)
    : RRay_2(sp, d) {}
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
