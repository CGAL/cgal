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
// Author(s)     : Herve Bronnimann

#ifndef CGAL_CARTESIAN_SPHERE_3_H
#define CGAL_CARTESIAN_SPHERE_3_H

#include <CGAL/Handle_for.h>
#include <CGAL/Interval_nt.h>
#include <boost/tuple/tuple.hpp>
#include <CGAL/Kernel/global_functions_3.h>

namespace CGAL {

template <class R_>
class SphereC3
{
  typedef typename R_::FT                   FT;
// http://www.cgal.org/Members/Manual_test/LAST/Developers_internal_manual/Developers_manual/Chapter_code_format.html#sec:programming_conventions
  typedef typename R_::Point_3              Point_3_;
  typedef typename R_::Vector_3             Vector_3;
  typedef typename R_::Sphere_3             Sphere_3;
  typedef typename R_::Circle_3             Circle_3;

  typedef boost::tuple<Point_3_, FT, Orientation>   Rep;
  typedef typename R_::template Handle<Rep>::type  Base;

  Base base;

public:
  typedef R_                                     R;

  SphereC3() {}

  SphereC3(const Point_3_ &center, const FT &squared_radius,
           const Orientation &o = COUNTERCLOCKWISE)
  {
    CGAL_kernel_precondition( (squared_radius >= FT(0)) &
                              (o != COLLINEAR) );

    base = Rep(center, squared_radius, o);
  }

  // Sphere passing through and oriented by p,q,r,s
  SphereC3(const Point_3_ &p, const Point_3_ &q,
           const Point_3_ &r, const Point_3_ &s)
  {
    Orientation orient = make_certain(CGAL::orientation(p, q, r, s));
    Point_3_ center = circumcenter(p, q, r, s);
    FT      squared_radius = squared_distance(p, center);

    base = Rep(center, squared_radius, orient);
  }

  // Sphere with great circle passing through p,q,r, oriented by o
  SphereC3(const Point_3_ &p, const Point_3_ &q, const Point_3_ &r,
	   const Orientation &o = COUNTERCLOCKWISE)
  {
    CGAL_kernel_precondition(o != COLLINEAR);

    Point_3_ center = circumcenter(p, q, r);
    FT      squared_radius = squared_distance(p, center);

    base = Rep(center, squared_radius, o);
  }

  // Sphere with diameter pq and orientation o
  SphereC3(const Point_3_ &p, const Point_3_ &q,
           const Orientation &o = COUNTERCLOCKWISE)
  {
    CGAL_kernel_precondition(o != COLLINEAR);

    Point_3_ center = midpoint(p, q);
    FT      squared_radius = squared_distance(p, center);

    base = Rep(center, squared_radius, o);
  }

  explicit SphereC3(const Point_3_ &center,
           const Orientation& o = COUNTERCLOCKWISE)
  {
    CGAL_kernel_precondition(o != COLLINEAR);

    base = Rep(center, FT(0), o);
  }

  typename R::Boolean   operator==(const SphereC3 &) const;
  typename R::Boolean   operator!=(const SphereC3 &) const;

  const Point_3_ & center() const
  {
      return get(base).template get<0>();
  }
  const FT & squared_radius() const
  {
      // Returns the square of the radius (instead of the radius itself,
      // which would require square roots)
      return get(base).template get<1>();
  }
  Orientation orientation() const
  {
      return get(base).template get<2>();
  }

  // A circle is degenerate if its (squared) radius is null or negative
  typename R::Boolean   is_degenerate() const;

  // Returns a circle with opposite orientation
  Sphere_3 opposite() const;

  typename R_::Oriented_side  oriented_side(const Point_3_ &p) const;
  //! precond: ! x.is_degenerate() (when available)
  // Returns R::ON_POSITIVE_SIDE, R::ON_ORIENTED_BOUNDARY or
  // R::ON_NEGATIVE_SIDE
  typename R::Boolean   has_on(const Circle_3 &p) const;
  typename R::Boolean   has_on(const Point_3_ &p) const;
  typename R::Boolean   has_on_boundary(const Point_3_ &p) const;
  typename R::Boolean   has_on_positive_side(const Point_3_ &p) const;
  typename R::Boolean   has_on_negative_side(const Point_3_ &p) const;

  typename R_::Bounded_side bounded_side(const Point_3_ &p) const;
  //! precond: ! x.is_degenerate() (when available)
  // Returns R::ON_BOUNDED_SIDE, R::ON_BOUNDARY or R::ON_UNBOUNDED_SIDE
  typename R::Boolean   has_on_bounded_side(const Point_3_ &p) const;
  typename R::Boolean   has_on_unbounded_side(const Point_3_ &p) const;
};

template < class R >
CGAL_KERNEL_INLINE
typename R::Boolean
SphereC3<R>::operator==(const SphereC3<R> &t) const
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
SphereC3<R>::operator!=(const SphereC3<R> &t) const
{
  return !(*this == t);
}

template < class R >
CGAL_KERNEL_MEDIUM_INLINE
typename R::Oriented_side
SphereC3<R>::
oriented_side(const typename SphereC3<R>::Point_3_ &p) const
{
  return enum_cast<Oriented_side>(bounded_side(p)) * orientation();
}

template < class R >
CGAL_KERNEL_INLINE
typename R::Bounded_side
SphereC3<R>::
bounded_side(const typename SphereC3<R>::Point_3_ &p) const
{
  return enum_cast<Bounded_side>(compare(squared_radius(),
                                         squared_distance(center(), p)));
}

template < class R >
inline
typename R::Boolean
SphereC3<R>::
has_on(const typename SphereC3<R>::Circle_3 &c) const
{
  typedef typename SphereC3<R>::Point_3_ Point_3_;
  typedef typename SphereC3<R>::FT      FT;
  Point_3_ proj = c.supporting_plane().projection(center());
  if(!(proj == c.center())) return false;
  const FT d2 = squared_distance(center(),c.center());
  return ((squared_radius() - d2) == c.squared_radius());
}

template < class R >
inline
typename R::Boolean
SphereC3<R>::
has_on(const typename SphereC3<R>::Point_3_ &p) const
{
  return has_on_boundary(p);
}

template < class R >
inline
typename R::Boolean
SphereC3<R>::
has_on_boundary(const typename SphereC3<R>::Point_3_ &p) const
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
SphereC3<R>::
has_on_negative_side(const typename SphereC3<R>::Point_3_ &p) const
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
SphereC3<R>::
has_on_positive_side(const typename SphereC3<R>::Point_3_ &p) const
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
SphereC3<R>::
has_on_bounded_side(const typename SphereC3<R>::Point_3_ &p) const
{
    // FIXME: it's a predicate...
  return squared_distance(center(),p) < squared_radius();
  // NB: J'ai aussi trouve ailleurs :
  // return bounded_side(p)==ON_BOUNDED_SIDE;
}

template < class R >
inline
typename R::Boolean
SphereC3<R>::
has_on_unbounded_side(const typename SphereC3<R>::Point_3_ &p) const
{
    // FIXME: it's a predicate...
  return squared_distance(center(),p) > squared_radius();
  // NB: J'ai aussi trouve ailleurs :
  // return bounded_side(p)==ON_UNBOUNDED_SIDE;
}

template < class R >
inline
typename R::Boolean
SphereC3<R>::
is_degenerate() const
{
    // FIXME: it's a predicate (?)
  return CGAL_NTS is_zero(squared_radius());
}

template < class R >
inline
typename SphereC3<R>::Sphere_3
SphereC3<R>::opposite() const
{
  return SphereC3<R>(center(), squared_radius(),
                               CGAL::opposite(orientation()) );
}

} //namespace CGAL

#endif // CGAL_CARTESIAN_SPHERE_3_H
