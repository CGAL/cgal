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
// author(s)     : Andreas Fabri, Stefan Schirra
//
// coordinator   : MPI, Saarbruecken
// ======================================================================
 
#ifndef CGAL_RAY_3_H
#define CGAL_RAY_3_H

CGAL_BEGIN_NAMESPACE

template <class R_>
class Ray_3 : public R_::Ray_3_base
{
  typedef typename R_::RT                    RT;
  typedef typename R_::Point_3               Point_3;
  typedef typename R_::Direction_3           Direction_3;
  typedef typename R_::Ray_3_base  RRay_3;
public:
  typedef          R_                       R;

  Ray_3()
      : RRay_3() {}

  Ray_3(const CGAL::Ray_3<R>& r)
      : RRay_3(r) {}

  Ray_3(const RRay_3& r)
      : RRay_3(r) {}

  Ray_3(const Point_3& sp, const Point_3& secondp)
    : RRay_3(sp, secondp) {}

  Ray_3(const Point_3& sp, const Direction_3& d)
    : RRay_3(sp, d) {}
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
