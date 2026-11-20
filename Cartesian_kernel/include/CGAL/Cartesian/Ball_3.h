// Copyright (c) 2000
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Herve Bronnimann

#ifndef CGAL_CARTESIAN_BALL_3_H
#define CGAL_CARTESIAN_BALL_3_H

#include <CGAL/Handle_for.h>
#include <CGAL/Interval_nt.h>
#include <boost/tuple/tuple.hpp>
#include <CGAL/Kernel/global_functions_3.h>

namespace CGAL {

template <class R_>
class BallC3
{
  typedef typename R_::Boolean              Boolean;
  typedef typename R_::Bounded_side         Bounded_side;
  typedef typename R_::FT                   FT;
// https://doc.cgal.org/latest/Manual/devman_code_format.html#secprogramming_conventions
  typedef typename R_::Point_3              Point_3_;
  typedef typename R_::Vector_3             Vector_3;
  typedef typename R_::Ball_3             Ball_3;
  typedef typename R_::Circle_3             Circle_3;

  typedef boost::tuple<Point_3_, FT, Orientation>   Rep;
  typedef typename R_::template Handle<Rep>::type  Base;

  Base base;

public:
  typedef R_                                     R;

  BallC3() {}

  BallC3(const Point_3_ &center, const FT &squared_radius,
         const Orientation &o = COUNTERCLOCKWISE)
  {
    CGAL_kernel_precondition( (squared_radius >= FT(0)) &
                              (o != COLLINEAR) );

    base = Rep(center, squared_radius, o);
  }

  // Ball passing through and oriented by p,q,r,s
  BallC3(const Point_3_ &p, const Point_3_ &q,
         const Point_3_ &r, const Point_3_ &s)
  {
    Orientation orient = make_certain(CGAL::orientation(p, q, r, s));
    Point_3_ center = CGAL::circumcenter(p, q, r, s);
    FT      squared_radius = squared_distance(p, center);

    base = Rep(center, squared_radius, orient);
  }

  // Ball with great circle passing through p,q,r, oriented by o
  BallC3(const Point_3_ &p, const Point_3_ &q, const Point_3_ &r,
         const Orientation &o = COUNTERCLOCKWISE)
  {
    CGAL_kernel_precondition(o != COLLINEAR);

    Point_3_ center = CGAL::circumcenter(p, q, r);
    FT      squared_radius = squared_distance(p, center);

    base = Rep(center, squared_radius, o);
  }

  // Ball with diameter pq and orientation o
  BallC3(const Point_3_ &p, const Point_3_ &q,
         const Orientation &o = COUNTERCLOCKWISE)
  {
    CGAL_kernel_precondition(o != COLLINEAR);

    Point_3_ center = midpoint(p, q);
    FT      squared_radius = squared_distance(p, center);

    base = Rep(center, squared_radius, o);
  }

  explicit BallC3(const Point_3_ &center,
           const Orientation& o = COUNTERCLOCKWISE)
  {
    CGAL_kernel_precondition(o != COLLINEAR);

    base = Rep(center, FT(0), o);
  }

  typename R::Boolean   operator==(const BallC3 &) const;
  typename R::Boolean   operator!=(const BallC3 &) const;

  const Point_3_ & center() const
  {
      return get_pointee_or_identity(base).template get<0>();
  }
  const FT & squared_radius() const
  {
      // Returns the square of the radius (instead of the radius itself,
      // which would require square roots)
      return get_pointee_or_identity(base).template get<1>();
  }
  Orientation orientation() const
  {
      return get_pointee_or_identity(base).template get<2>();
  }

  // A circle is degenerate if its (squared) radius is null or negative
  typename R::Boolean   is_degenerate() const;

  // Returns a circle with opposite orientation
  Ball_3 opposite() const;

  typename R_::Oriented_side  oriented_side(const Point_3_ &p) const;
  //! precond: ! x.is_degenerate() (when available)
  // Returns R::ON_POSITIVE_SIDE, R::ON_ORIENTED_BOUNDARY or
  // R::ON_NEGATIVE_SIDE
  Boolean has_on(const Circle_3 &p) const;
  Boolean has_on(const Point_3_ &p) const;
  Boolean has_on_boundary(const Point_3_ &p) const;
  Boolean has_on_positive_side(const Point_3_ &p) const;
  Boolean has_on_negative_side(const Point_3_ &p) const;

  Bounded_side bounded_side(const Point_3_ &p) const;
  //! precond: ! x.is_degenerate() (when available)
  // Returns R::ON_BOUNDED_SIDE, R::ON_BOUNDARY or R::ON_UNBOUNDED_SIDE
  Boolean has_on_bounded_side(const Point_3_ &p) const;
  Boolean has_on_unbounded_side(const Point_3_ &p) const;
};

template < class R >
CGAL_KERNEL_INLINE
typename R::Boolean
BallC3<R>::operator==(const BallC3<R> &t) const
{
  if (CGAL::identical(base, t.base))
      return true;
  return center() == t.center() &&
         squared_radius() == t.squared_radius() &&
         orientation() == t.orientation();
}

template < class R >
inline
typename R::Boolean
BallC3<R>::operator!=(const BallC3<R> &t) const
{
  return !(*this == t);
}

template < class R >
CGAL_KERNEL_MEDIUM_INLINE
typename R::Oriented_side
BallC3<R>::
oriented_side(const typename BallC3<R>::Point_3_ &p) const
{
  typedef typename R::Oriented_side Oriented_side;
  return enum_cast<Oriented_side>(bounded_side(p)) * orientation();
}

template < class R >
CGAL_KERNEL_INLINE
typename R::Bounded_side
BallC3<R>::
bounded_side(const typename BallC3<R>::Point_3_ &p) const
{
  typedef typename R::Bounded_side Bounded_side;
  return enum_cast<Bounded_side>(compare(squared_radius(),
                                         squared_distance(center(), p)));
}

template < class R >
inline
typename R::Boolean
BallC3<R>::
has_on(const typename BallC3<R>::Circle_3 &c) const
{
  typedef typename BallC3<R>::Point_3_ Point_3_;
  typedef typename BallC3<R>::FT      FT;
  Point_3_ proj = c.supporting_plane().projection(center());
  if(!(proj == c.center())) return false;
  const FT d2 = squared_distance(center(),c.center());
  return ((squared_radius() - d2) == c.squared_radius());
}

template < class R >
inline
typename R::Boolean
BallC3<R>::
has_on(const typename BallC3<R>::Point_3_ &p) const
{
  return has_on_boundary(p);
}

template < class R >
inline
typename R::Boolean
BallC3<R>::
has_on_boundary(const typename BallC3<R>::Point_3_ &p) const
{
    // FIXME: it's a predicate...
  return squared_distance(center(),p) == squared_radius();
  // NB: J'ai aussi trouve ailleurs :
  // return oriented_side(p)==ON_ORIENTED_BOUNDARY;
  // a voir...
}

template < class R >
CGAL_KERNEL_INLINE
typename R::Boolean
BallC3<R>::
has_on_negative_side(const typename BallC3<R>::Point_3_ &p) const
{
  if (orientation() == COUNTERCLOCKWISE)
    return has_on_unbounded_side(p);
  return has_on_bounded_side(p);
  // NB: J'ai aussi trouve ailleurs :
  // return oriented_side(p)==ON_NEGATIVE_SIDE;
}

template < class R >
CGAL_KERNEL_INLINE
typename R::Boolean
BallC3<R>::
has_on_positive_side(const typename BallC3<R>::Point_3_ &p) const
{
  if (orientation() == COUNTERCLOCKWISE)
    return has_on_bounded_side(p);
  return has_on_unbounded_side(p);
  // NB: J'ai aussi trouve ailleurs :
  // return oriented_side(p)==ON_POSITIVE_SIDE;
}

template < class R >
inline
typename R::Boolean
BallC3<R>::
has_on_bounded_side(const typename BallC3<R>::Point_3_ &p) const
{
    // FIXME: it's a predicate...
  return squared_distance(center(),p) < squared_radius();
  // NB: J'ai aussi trouve ailleurs :
  // return bounded_side(p)==ON_BOUNDED_SIDE;
}

template < class R >
inline
typename R::Boolean
BallC3<R>::
has_on_unbounded_side(const typename BallC3<R>::Point_3_ &p) const
{
    // FIXME: it's a predicate...
  return squared_distance(center(),p) > squared_radius();
  // NB: J'ai aussi trouve ailleurs :
  // return bounded_side(p)==ON_UNBOUNDED_SIDE;
}

template < class R >
inline
typename R::Boolean
BallC3<R>::
is_degenerate() const
{
    // FIXME: it's a predicate (?)
  return CGAL_NTS is_zero(squared_radius());
}

template < class R >
inline
typename BallC3<R>::Ball_3
BallC3<R>::opposite() const
{
  return BallC3<R>(center(), squared_radius(),
                               CGAL::opposite(orientation()) );
}

} //namespace CGAL

#endif // CGAL_CARTESIAN_BALL_3_H
