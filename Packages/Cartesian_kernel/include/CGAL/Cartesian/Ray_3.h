// Copyright (c) 2000  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbruecken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Andreas Fabri

#ifndef CGAL_CARTESIAN_RAY_3_H
#define CGAL_CARTESIAN_RAY_3_H

#include <CGAL/Twotuple.h>

CGAL_BEGIN_NAMESPACE

template < class R_ >
class RayC3
  : public R_::template Handle<Twotuple<typename R_::Point_3> >::type
{
CGAL_VC7_BUG_PROTECTED
  typedef typename R_::FT                   FT;
  typedef typename R_::Point_3              Point_3;
  typedef typename R_::Direction_3          Direction_3;
  typedef typename R_::Vector_3             Vector_3;
  typedef typename R_::Line_3               Line_3;
  typedef typename R_::Ray_3                Ray_3;
  typedef typename R_::Aff_transformation_3 Aff_transformation_3;

  typedef Twotuple<Point_3>                        rep;
  typedef typename R_::template Handle<rep>::type  base;

  const base& Base() const { return *this; }
  base& Base() { return *this; }

public:
  typedef R_                                     R;

  RayC3() {}

  RayC3(const Point_3 &sp, const Point_3 &secondp)
    : base(sp, secondp) {}

  RayC3(const Point_3 &sp, const Vector_3 &v)
    : base(sp, sp + v) {}

  RayC3(const Point_3 &sp, const Direction_3 &d)
    : base(sp, sp + d.to_vector()) {}

  RayC3(const Point_3 &sp, const Line_3 &l)
    : base(sp, sp + l.to_vector()) {}

  bool        operator==(const RayC3 &r) const;
  bool        operator!=(const RayC3 &r) const;

  const Point_3 &   start() const;
  const Point_3 &   source() const
  {
      return get(Base()).e0;
  }
  const Point_3 &   second_point() const
  {
      return get(Base()).e1;
  }
  Point_3     point(int i) const;

  Direction_3 direction() const;
  Vector_3    to_vector() const;
  Line_3      supporting_line() const;
  Ray_3       opposite() const;

  Ray_3       transform(const Aff_transformation_3 &t) const
  {
    return RayC3<R>(t.transform(source()), t.transform(second_point()));
  }

  bool        is_degenerate() const;
  bool        has_on(const Point_3 &p) const;
  bool        collinear_has_on(const Point_3 &p) const;
};

template < class R >
inline
bool
RayC3<R>::operator==(const RayC3<R> &r) const
{
    if (CGAL::identical(Base(), r.Base()))
	return true;
    return source() == r.source() && direction() == r.direction();
}

template < class R >
inline
bool
RayC3<R>::operator!=(const RayC3<R> &r) const
{
  return !(*this == r);
}

template < class R >
inline
const typename RayC3<R>::Point_3 &
RayC3<R>::start() const
{
  return source();
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
bool
RayC3<R>::
has_on(const typename RayC3<R>::Point_3 &p) const
{
  return (p == source()) ||
         ( collinear(source(), p, second_point())
           && ( Direction_3(p - source()) == direction() ));
}

template < class R >
inline
bool
RayC3<R>::is_degenerate() const
{
  return source() == second_point();
}

template < class R >
inline
bool
RayC3<R>::
collinear_has_on(const typename RayC3<R>::Point_3 &p) const
{
  CGAL_kernel_exactness_precondition( collinear(source(), p, second_point()) );

  Comparison_result cx = compare_x(source(), second_point());
  if (cx != EQUAL)
    return cx != compare_x(p, source());

  Comparison_result cy = compare_y(source(), second_point());
  if (cy != EQUAL)
    return cy != compare_y(p, source());

  Comparison_result cz = compare_z(source(), second_point());
  if (cz != EQUAL)
    return cz != compare_z(p, source());

  return true; // p == source()
}

#ifndef CGAL_NO_OSTREAM_INSERT_RAYC3
template < class R >
std::ostream &
operator<<(std::ostream &os, const RayC3<R> &r)
{
    switch(os.iword(IO::mode)) {
    case IO::ASCII :
        return os << r.start() << ' ' << r.direction();
    case IO::BINARY :
        return os<< r.start() << r.direction();
    default:
        return os << "RayC3(" << r.start() <<  ", " << r.direction() << ")";
    }
}
#endif // CGAL_NO_OSTREAM_INSERT_RAYC3

#ifndef CGAL_NO_ISTREAM_EXTRACT_RAYC3
template < class R >
std::istream &
operator>>(std::istream &is, RayC3<R> &r)
{
    typename R::Point_3 p;
    typename R::Direction_3 d;

    is >> p >> d;

    if (is)
	r = RayC3<R>(p, d);
    return is;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_RAYC3

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_RAY_3_H
