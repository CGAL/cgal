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

#ifndef CGAL_REP_CLASS_DEFINED
#error  no representation class defined
#endif  // CGAL_REP_CLASS_DEFINED

#if defined(CGAL_CFG_INCOMPLETE_TYPE_BUG_1) && \
   !defined(CGAL_NO_PLANE_TRANSFORM_IN_AT)
#define CGAL_NO_PLANE_TRANSFORM_IN_AT
#endif // CGAL_CFG_INCOMPLETE_TYPE_BUG_1

#if defined CGAL_HOMOGENEOUS_H || defined CGAL_SIMPLE_HOMOGENEOUS_H
#include <CGAL/PlaneH3.h>
#endif

#if defined CGAL_CARTESIAN_H || defined CGAL_SIMPLE_CARTESIAN_H
#include <CGAL/Cartesian/Plane_3.h>
#endif

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

  bool                  operator==(const CGAL::Plane_3<R>& p) const
  { return RPlane_3::operator==(p); }

  bool                  operator!=(const CGAL::Plane_3<R>& p) const
  { return !(*this == p); }


  RT a() const
  { return RPlane_3::a(); }

  RT                    b() const
  { return RPlane_3::b(); }

  RT                    c() const
  { return RPlane_3::c(); }

  RT                    d() const
  { return RPlane_3::d(); }

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

  Oriented_side   oriented_side(const CGAL::Point_3<R>& p) const
  { return RPlane_3::oriented_side(p); }

  bool                 has_on(const  CGAL::Point_3<R>& p) const
  { return RPlane_3::has_on_boundary(p); }

  bool                 has_on(const  CGAL::Line_3<R>& l) const
  { return RPlane_3::has_on_boundary(l); }

  bool                 has_on_boundary(const  CGAL::Point_3<R>& p) const
  { return RPlane_3::has_on_boundary(p); }

  bool                 has_on_boundary(const  CGAL::Line_3<R>& l) const
  { return RPlane_3::has_on_boundary(l); }

  bool                 has_on_positive_side(const  CGAL::Point_3<R>& p) const
  { return RPlane_3::has_on_positive_side(p); }

  bool                 has_on_negative_side(const  CGAL::Point_3<R>& p) const
  { return RPlane_3::has_on_negative_side(p); }

  bool                 is_degenerate() const
  { return RPlane_3::is_degenerate(); }
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
