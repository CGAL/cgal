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
// file          : Plane_3.h
// package       : _3
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas Fabri
//                 Stefan Schirra
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ======================================================================
 
#ifndef CGAL_PLANE_3_H
#define CGAL_PLANE_3_H

#include <CGAL/Line_3.h>
#include <CGAL/Point_2.h>

CGAL_BEGIN_NAMESPACE

template <class R_>
class Plane_3 : public R_::Plane_3_base
{
public:
  typedef          R_                       R;
  typedef typename R::RT                    RT;
  typedef typename R::FT                    FT;
  typedef typename R::Plane_3_base  RPlane_3;

  Plane_3() : RPlane_3()
  {}

  Plane_3(const CGAL::Plane_3<R>& p) : RPlane_3(p)
  {}

  Plane_3(const RPlane_3&  p) : RPlane_3(p)
  {}

  Plane_3(const CGAL::Point_3<R>& p,
          const CGAL::Point_3<R>& q,
          const CGAL::Point_3<R>& r)
    : RPlane_3(p,q,r)
  {}

  Plane_3(const CGAL::Point_3<R>& p, const CGAL::Direction_3<R>& d)
    : RPlane_3(p,d)
  {}

  Plane_3(const CGAL::Point_3<R>& p, const CGAL::Vector_3<R>& v)
    : RPlane_3(p,v)
  {}

  Plane_3(const RT& a, const RT& b, const RT& c, const RT& d)
    : RPlane_3(a,b,c,d)
  {}

  Plane_3(const CGAL::Line_3<R>& l, const CGAL::Point_3<R>& p)
    : RPlane_3(l,p)
  {}

  Plane_3(const CGAL::Segment_3<R>& s, const CGAL::Point_3<R>& p)
    : RPlane_3(s,p)
  {}

  Plane_3(CGAL::Ray_3<R>& r, const CGAL::Point_3<R>& p)
    : RPlane_3(r,p)
  {}

  CGAL::Line_3<R>       perpendicular_line(const CGAL::Point_3<R>& p) const
  { return RPlane_3::perpendicular_line(p); }

  CGAL::Plane_3<R>      opposite() const
  { return RPlane_3::opposite(); }

  CGAL::Point_3<R>      projection(const CGAL::Point_3<R>& p) const
  { return RPlane_3::projection(p); }

  CGAL::Point_3<R>      point() const
  { return RPlane_3::point(); }

  CGAL::Vector_3<R>     orthogonal_vector() const
  { return RPlane_3::orthogonal_vector(); }

  CGAL::Direction_3<R>  orthogonal_direction() const
  { return RPlane_3::orthogonal_direction(); }

  CGAL::Vector_3<R>     base1() const
  { return RPlane_3::base1(); }

  CGAL::Vector_3<R>     base2() const
  { return RPlane_3::base2(); }

  CGAL::Point_2<R>      to_2d(const CGAL::Point_3<R>& p) const
  { return RPlane_3::to_2d(p); }

  CGAL::Point_3<R>      to_3d(const CGAL::Point_2<R>& p) const
  { return RPlane_3::to_3d(p); }

  CGAL::Plane_3<R>      transform( CGAL::Aff_transformation_3<R>& t) const
  { return CGAL::Plane_3<R>( RPlane_3::transform(t) ); }
};

#ifndef CGAL_NO_OSTREAM_INSERT_PLANE_3
template < class R >
std::ostream&
operator<<(std::ostream& os, const Plane_3<R>& p)
{
  typedef typename  R::Plane_3_base  RPlane_3;
  return os << static_cast<const RPlane_3&>(p);
}
#endif // CGAL_NO_OSTREAM_INSERT_PLANE_3

#ifndef CGAL_NO_ISTREAM_EXTRACT_PLANE_3
template < class R >
std::istream&
operator>>(std::istream& is, Plane_3<R>& t)
{
  typedef typename  R::Plane_3_base  RPlane_3;
  return is >> static_cast<RPlane_3&>(t);
}
#endif // CGAL_NO_ISTREAM_EXTRACT_PLANE_3

CGAL_END_NAMESPACE

#endif // CGAL_PLANE_3_H
