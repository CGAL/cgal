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
// Author(s)     : Andreas Fabri, Herve Bronnimann

#ifndef CGAL_CARTESIAN_RAY_2_H
#define CGAL_CARTESIAN_RAY_2_H

#include <CGAL/Twotuple.h>

CGAL_BEGIN_NAMESPACE

template < class R_ >
class RayC2
{
  typedef typename R_::FT                   FT;
  typedef typename R_::Point_2              Point_2;
  typedef typename R_::Direction_2          Direction_2;
  typedef typename R_::Vector_2             Vector_2;
  typedef typename R_::Line_2               Line_2;
  typedef typename R_::Ray_2                Ray_2;
  typedef typename R_::Aff_transformation_2 Aff_transformation_2;

  typedef Twotuple<Point_2>                        Rep;
  typedef typename R_::template Handle<Rep>::type  Base;

  Base base;

public:
  typedef R_                                     R;
  typename R::Construct_translated_point_2 Construct_translated_point;

  RayC2() {}

  RayC2(const Point_2 &sp, const Point_2 &secondp)
    : base(sp, secondp) {}

  RayC2(const Point_2 &sp, const Direction_2 &d)
    : base(sp, Construct_translated_point( sp, d.to_vector())){}

  RayC2(const Point_2 &sp, const Vector_2 &v)
    : base(sp, Construct_translated_point(sp,  v)){}

  RayC2(const Point_2 &sp, const Line_2 &l)
    : base(sp, Construct_translated_point(sp, l.to_vector())){}

  bool        operator==(const RayC2 &r) const;
  bool        operator!=(const RayC2 &r) const;

  const Point_2 &     start() const;
  const Point_2 &     source() const
  {
      return get(base).e0;
  }
  Point_2     point(int i) const;
  const Point_2 &     second_point() const
  {
      return get(base).e1;
  }

  Direction_2 direction() const;
  Vector_2    to_vector() const;
  Line_2      supporting_line() const;
  Ray_2       opposite() const;

  Ray_2       transform(const Aff_transformation_2 &t) const
  {
    return RayC2<R>(t.transform(source()), t.transform(second_point()));
  }

  bool        is_horizontal() const;
  bool        is_vertical() const;
  bool        is_degenerate() const;
  bool        has_on(const Point_2 &p) const;
  bool        collinear_has_on(const Point_2 &p) const;
};

template < class R >
CGAL_KERNEL_INLINE
bool
RayC2<R>::operator==(const RayC2<R> &r) const
{
  if (CGAL::identical(base, r.base))
      return true;
  return source() == r.source() && direction() == r.direction();
}

template < class R >
inline
bool
RayC2<R>::operator!=(const RayC2<R> &r) const
{
  return !(*this == r);
}

template < class R >
inline
const typename RayC2<R>::Point_2 &
RayC2<R>::start() const
{
  return source();
}

template < class R >
CGAL_KERNEL_INLINE
typename RayC2<R>::Point_2
RayC2<R>::point(int i) const
{
  CGAL_kernel_precondition( i >= 0 );

  typename R::Construct_vector_2 construct_vector;
  typename R::Construct_scaled_vector_2 construct_scaled_vector;
  typename R::Construct_translated_point_2 construct_translated_point;
  if (i == 0) return source();
  if (i == 1) return second_point();
  return construct_translated_point(source(),
				    construct_scaled_vector(construct_vector(source(), 
									     second_point()),
							    FT(i)));
}

template < class R >
inline
typename RayC2<R>::Vector_2
RayC2<R>::to_vector() const
{
  typename R::Construct_vector_2 construct_vector;
  return construct_vector(source(), second_point());
}

template < class R >
inline
typename RayC2<R>::Direction_2
RayC2<R>::direction() const
{
  typename R::Construct_vector_2 construct_vector;
  return Direction_2( construct_vector(source(), second_point()) );
}

template < class R >
inline
typename RayC2<R>::Line_2
RayC2<R>::supporting_line() const
{
  return Line_2(*this);
}

template < class R >
inline
typename RayC2<R>::Ray_2
RayC2<R>::opposite() const
{
  return RayC2<R>( source(), - direction() );
}

template < class R >
CGAL_KERNEL_INLINE
bool RayC2<R>::is_horizontal() const
{
  return R().equal_y_2_object()(source(), second_point());
}

template < class R >
CGAL_KERNEL_INLINE
bool RayC2<R>::is_vertical() const
{
  return R().equal_x_2_object()(source(), second_point());
}

template < class R >
CGAL_KERNEL_INLINE
bool RayC2<R>::is_degenerate() const
{
  return source() == second_point();
}

template < class R >
CGAL_KERNEL_INLINE
bool
RayC2<R>::has_on(const typename RayC2<R>::Point_2 &p) const
{
  typename R::Construct_vector_2  construct_vector;
  return p == source()
      || R().collinear_2_object()(source(), p, second_point())
      && Direction_2(construct_vector( source(), p)) == direction();
}

template < class R >
inline
bool
RayC2<R>::
collinear_has_on(const typename RayC2<R>::Point_2 &p) const
{
  return R().collinear_has_on_2_object()
               (static_cast<const typename R::Ray_2>(*this), p);
}

#ifndef CGAL_NO_OSTREAM_INSERT_RAYC2
template < class R >
std::ostream &
operator<<(std::ostream &os, const RayC2<R> &r)
{
    switch(os.iword(IO::mode)) {
    case IO::ASCII :
        return os << r.source() << ' ' << r.direction();
    case IO::BINARY :
        return os << r.source() << r.direction();
    default:
        return os << "RayC2(" << r.source() <<  ", " << r.direction() << ")";
    }
}
#endif // CGAL_NO_OSTREAM_INSERT_RAYC2

#ifndef CGAL_NO_ISTREAM_EXTRACT_RAYC2
template < class R >
std::istream &
operator>>(std::istream &is, RayC2<R> &r)
{
    typename R::Point_2 p;
    typename R::Direction_2 d;

    is >> p >> d;

    if (is)
	r = RayC2<R>(p, d);
    return is;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_RAYC2

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_RAY_2_H
