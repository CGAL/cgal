// Copyright (c) 2000  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved. 
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0+
// 
//
// Author(s)     : Andreas Fabri

#ifndef CGAL_CARTESIAN_PLANE_3_H
#define CGAL_CARTESIAN_PLANE_3_H

#include <CGAL/array.h>
#include <CGAL/Handle_for.h>
#include <CGAL/Cartesian/solve_3.h>
#include <CGAL/Cartesian/plane_constructions_3.h>

namespace CGAL {

template <class R_>
class PlaneC3
{
  typedef typename R_::FT                   FT;
  typedef typename R_::Point_2              Point_2;
  typedef typename R_::Point_3              Point_3;
  typedef typename R_::Vector_3             Vector_3;
  typedef typename R_::Direction_3          Direction_3;
  typedef typename R_::Line_3               Line_3;
  typedef typename R_::Ray_3                Ray_3;
  typedef typename R_::Segment_3            Segment_3;
  typedef typename R_::Plane_3              Plane_3;
  typedef typename R_::Circle_3             Circle_3;
  typedef typename R_::Construct_point_3    Construct_point_3;
  typedef typename R_::Construct_point_2    Construct_point_2;

  typedef cpp11::array<FT, 4>               Rep;
  typedef typename R_::template Handle<Rep>::type  Base;

  Base base;

public:

  typedef R_                                     R;

  PlaneC3() {}

  PlaneC3(const Point_3 &p, const Point_3 &q, const Point_3 &r)
  { *this = plane_from_points<R>(p, q, r); }

  PlaneC3(const Point_3 &p, const Direction_3 &d)
  { *this = plane_from_point_direction<R>(p, d); }

  PlaneC3(const Point_3 &p, const Vector_3 &v)
  { *this = plane_from_point_direction<R>(p, v.direction()); }

  PlaneC3(const FT &a, const FT &b, const FT &c, const FT &d)
    : base(CGAL::make_array(a, b, c, d)) {}

  PlaneC3(const Line_3 &l, const Point_3 &p)
  { *this = plane_from_points<R>(l.point(),
	                      l.point()+l.direction().to_vector(),
			      p); }

  PlaneC3(const Segment_3 &s, const Point_3 &p)
  { *this = plane_from_points<R>(s.start(), s.end(), p); }

  PlaneC3(const Ray_3 &r, const Point_3 &p)
  { *this = plane_from_points<R>(r.start(), r.second_point(), p); }

  typename R::Boolean   operator==(const PlaneC3 &p) const;
  typename R::Boolean   operator!=(const PlaneC3 &p) const;

  const FT & a() const
  {
      return get_pointee_or_identity(base)[0];
  }
  const FT & b() const
  {
      return get_pointee_or_identity(base)[1];
  }
  const FT & c() const
  {
      return get_pointee_or_identity(base)[2];
  }
  const FT & d() const
  {
      return get_pointee_or_identity(base)[3];
  }

  Line_3       perpendicular_line(const Point_3 &p) const;
  Plane_3      opposite() const;

  Point_3      point() const;
  Point_3      projection(const Point_3 &p) const;
  Vector_3     orthogonal_vector() const;
  Direction_3  orthogonal_direction() const;
  Vector_3     base1() const;
  Vector_3     base2() const;

  Point_3      to_plane_basis(const Point_3 &p) const;

  Point_2      to_2d(const Point_3 &p) const;
  Point_3      to_3d(const Point_2 &p) const;

  typename R::Oriented_side     oriented_side(const Point_3 &p) const;
  typename R::Boolean           has_on_positive_side(const Point_3 &l) const;
  typename R::Boolean           has_on_negative_side(const Point_3 &l) const;
  typename R::Boolean           has_on(const Point_3 &p) const
  {
    return oriented_side(p) == ON_ORIENTED_BOUNDARY;
  }
  typename R::Boolean           has_on(const Line_3 &l) const
  {
    return has_on(l.point())
       &&  has_on(l.point() + l.direction().to_vector());
  }
  typename R::Boolean           has_on(const Circle_3 &circle) const
  {
    if(circle.squared_radius() != FT(0)) {
      const Plane_3& p = circle.supporting_plane();
      if(is_zero(a())) {
        if(!is_zero(p.a())) return false;
        if(is_zero(b())) {
          if(!is_zero(p.b())) return false;
          return c() * p.d() == d() * p.c();
        }
        return (p.c() * b() == c() * p.b()) &&
               (p.d() * b() == d() * p.b());
      }
      return (p.b() * a() == b() * p.a()) &&
             (p.c() * a() == c() * p.a()) &&
             (p.d() * a() == d() * p.a());
    } else return has_on(circle.center());
  }

  typename R::Boolean           is_degenerate() const;
};

template < class R >
CGAL_KERNEL_INLINE
typename R::Boolean
PlaneC3<R>::operator==(const PlaneC3<R> &p) const
{
  if (CGAL::identical(base, p.base))
      return true;
  return equal_plane(*this, p);
}

template < class R >
inline
typename R::Boolean
PlaneC3<R>::operator!=(const PlaneC3<R> &p) const
{
  return !(*this == p);
}

template < class R >
inline
typename PlaneC3<R>::Point_3
PlaneC3<R>::point() const
{
  return point_on_plane(*this);
}

template < class R >
inline
typename PlaneC3<R>::Point_3
PlaneC3<R>::
projection(const typename PlaneC3<R>::Point_3 &p) const
{
  return projection_plane(p, *this);
}

template < class R >
inline
typename PlaneC3<R>::Vector_3
PlaneC3<R>::orthogonal_vector() const
{
  return R().construct_orthogonal_vector_3_object()(*this);
}

template < class R >
inline
typename PlaneC3<R>::Direction_3
PlaneC3<R>::orthogonal_direction() const
{
  return Direction_3(a(), b(), c());
}

template < class R >
typename PlaneC3<R>::Vector_3
PlaneC3<R>::base1() const
{
  return R().construct_base_vector_3_object()(*this, 1);
}

template < class R >
typename PlaneC3<R>::Vector_3
PlaneC3<R>::base2() const
{
  return R().construct_base_vector_3_object()(*this, 2);
}

template < class R >
typename PlaneC3<R>::Point_3
PlaneC3<R>::
to_plane_basis(const typename PlaneC3<R>::Point_3 &p) const
{
  FT alpha, beta, gamma;
  Construct_point_3 construct_point_3;
  Cartesian_internal::solve(base1(), base2(), orthogonal_vector(), p - point(),
	alpha, beta, gamma);

  return construct_point_3(alpha, beta, gamma);
}

template < class R >
typename PlaneC3<R>::Point_2
PlaneC3<R>::
to_2d(const typename PlaneC3<R>::Point_3 &p) const
{
  FT alpha, beta, gamma;
  Construct_point_2 construct_point_2;

  Cartesian_internal::solve(base1(), base2(), orthogonal_vector(), p - point(),
	alpha, beta, gamma);

  return construct_point_2(alpha, beta);
}

template < class R >
inline
typename PlaneC3<R>::Point_3
PlaneC3<R>::
to_3d(const typename PlaneC3<R>::Point_2 &p) const
{
  return R().construct_lifted_point_3_object()(*this, p);
}

template < class R >
inline
typename PlaneC3<R>::Line_3
PlaneC3<R>::
perpendicular_line(const typename PlaneC3<R>::Point_3 &p) const
{
  return Line_3(p, orthogonal_direction());
}

template < class R >
inline
typename PlaneC3<R>::Plane_3
PlaneC3<R>::opposite() const
{
  return PlaneC3<R>(-a(), -b(), -c(), -d());
}

template < class R >
inline
typename R::Oriented_side
PlaneC3<R>::
oriented_side(const typename PlaneC3<R>::Point_3 &p) const
{
  return side_of_oriented_plane(*this, p);
}

template < class R >
inline
typename R::Boolean
PlaneC3<R>::
has_on_positive_side(const  typename PlaneC3<R>::Point_3 &p) const
{
  return oriented_side(p) == ON_POSITIVE_SIDE;
}

template < class R >
inline
typename R::Boolean
PlaneC3<R>::
has_on_negative_side(const  typename PlaneC3<R>::Point_3 &p) const
{
  return oriented_side(p) == ON_NEGATIVE_SIDE;
}

template < class R >
inline
typename R::Boolean
PlaneC3<R>::
is_degenerate() const
{ // FIXME : predicate
  return CGAL_NTS is_zero(a()) && CGAL_NTS is_zero(b()) &&
         CGAL_NTS is_zero(c());
}

} //namespace CGAL

#endif // CGAL_CARTESIAN_PLANE_3_H
