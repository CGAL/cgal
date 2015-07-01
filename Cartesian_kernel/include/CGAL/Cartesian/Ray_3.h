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
// 
//
// Author(s)     : Andreas Fabri

#ifndef CGAL_CARTESIAN_RAY_3_H
#define CGAL_CARTESIAN_RAY_3_H

#include <CGAL/array.h>
#include <CGAL/Handle_for.h>

namespace CGAL {

template < class R_ >
class RayC3
{
  typedef typename R_::FT                   FT;
  typedef typename R_::Point_3              Point_3;
  typedef typename R_::Direction_3          Direction_3;
  typedef typename R_::Vector_3             Vector_3;
  typedef typename R_::Line_3               Line_3;
  typedef typename R_::Ray_3                Ray_3;

  typedef cpp11::array<Point_3, 2>          Rep;
  typedef typename R_::template Handle<Rep>::type  Base;

  Base base;

public:
  typedef R_                                     R;

  RayC3() {}

  RayC3(const Point_3 &sp, const Point_3 &secondp)
    : base(CGAL::make_array(sp, secondp)) {}

  RayC3(const Point_3 &sp, const Vector_3 &v)
    : base(CGAL::make_array(sp, sp + v)) {}

  RayC3(const Point_3 &sp, const Direction_3 &d)
    : base(CGAL::make_array(sp, sp + d.to_vector())) {}

  RayC3(const Point_3 &sp, const Line_3 &l)
    : base(CGAL::make_array(sp, sp + l.to_vector())) {}

  typename R::Boolean          operator==(const RayC3 &r) const;
  typename R::Boolean          operator!=(const RayC3 &r) const;

  const Point_3 &   source() const
  {
      return get_pointee_or_identity(base)[0];
  }
  const Point_3 &   second_point() const
  {
      return get_pointee_or_identity(base)[1];
  }
  Point_3     point(int i) const;

  Direction_3 direction() const;
  Vector_3    to_vector() const;
  Line_3      supporting_line() const;
  Ray_3       opposite() const;

  typename R::Boolean          is_degenerate() const;
  typename R::Boolean          has_on(const Point_3 &p) const;
  typename R::Boolean          collinear_has_on(const Point_3 &p) const;
};

template < class R >
inline
typename R::Boolean
RayC3<R>::operator==(const RayC3<R> &r) const
{
    if (CGAL::identical(base, r.base))
	return true;
    return source() == r.source() && direction() == r.direction();
}

template < class R >
inline
typename R::Boolean
RayC3<R>::operator!=(const RayC3<R> &r) const
{
  return !(*this == r);
}

template < class R >
CGAL_KERNEL_INLINE
typename RayC3<R>::Point_3
RayC3<R>::point(int i) const
{
  CGAL_kernel_precondition( i >= 0 );
  if (i == 0) return source();
  if (i == 1) return second_point();
  return source() + FT(i) * (second_point() - source());
}

template < class R >
inline
typename RayC3<R>::Vector_3
RayC3<R>::to_vector() const
{
  return second_point() - source();
}

template < class R >
inline
typename RayC3<R>::Direction_3
RayC3<R>::direction() const
{
  return Direction_3( second_point() - source() );
}

template < class R >
inline
typename RayC3<R>::Line_3
RayC3<R>::supporting_line() const
{
  return Line_3(*this);
}

template < class R >
inline
typename RayC3<R>::Ray_3
RayC3<R>::opposite() const
{
  return RayC3<R>( source(), - direction() );
}

template < class R >
typename R::Boolean
RayC3<R>::
has_on(const typename RayC3<R>::Point_3 &p) const
{
  return (p == source()) ||
         ( collinear(source(), p, second_point())
           && ( Direction_3(p - source()) == direction() ));
}

template < class R >
inline
typename R::Boolean
RayC3<R>::is_degenerate() const
{
  return source() == second_point();
}

template < class R >
inline
typename R::Boolean
RayC3<R>::
collinear_has_on(const typename RayC3<R>::Point_3 &p) const
{
  CGAL_kernel_exactness_precondition( collinear(source(), p, second_point()) );

  typename R::Comparison_result cx = compare_x(source(), second_point());
  if (cx != EQUAL)
    return cx != compare_x(p, source());

  typename R::Comparison_result cy = compare_y(source(), second_point());
  if (cy != EQUAL)
    return cy != compare_y(p, source());

  typename R::Comparison_result cz = compare_z(source(), second_point());
  if (cz != EQUAL)
    return cz != compare_z(p, source());

  return true; // p == source()
}

} //namespace CGAL

#endif // CGAL_CARTESIAN_RAY_3_H
