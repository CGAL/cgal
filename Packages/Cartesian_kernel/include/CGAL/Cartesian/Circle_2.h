// Copyright (c) 1997-2004  Utrecht University (The Netherlands),
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

#ifndef CGAL_CARTESIAN_CIRCLE_2_H
#define CGAL_CARTESIAN_CIRCLE_2_H

#include <CGAL/utility.h>
#include <CGAL/Interval_arithmetic.h>
#include <CGAL/Cartesian/predicates_on_points_2.h>

CGAL_BEGIN_NAMESPACE

template <class R_ >
class CircleC2
{
  typedef typename R_::FT                   FT;
  typedef typename R_::Circle_2             Circle_2;
  typedef typename R_::Point_2              Point_2;
  typedef typename R_::Aff_transformation_2 Aff_transformation_2;

  typedef Triple<Point_2, FT, Orientation>         Rep;
  typedef typename R_::template Handle<Rep>::type  Base;

  Base base;

public:
  typedef R_                                     R;

  CircleC2() {}

  CircleC2(const Point_2 &center, const FT &squared_radius = FT(0),
           const Orientation &orient = COUNTERCLOCKWISE) // Is this new?
  {
    CGAL_kernel_precondition( ( squared_radius >= FT(0) ) &&
                              ( orient    != COLLINEAR) );

    base = Rep(center, squared_radius, orient);
  }

  CircleC2(const Point_2 &center, const Orientation &orient) // Is this new?
  {
    CGAL_kernel_precondition( orient != COLLINEAR );

    base = Rep(center, FT(0), orient);
  }

  CircleC2(const Point_2 &p, const Point_2 &q,
           const Orientation &orient = COUNTERCLOCKWISE) // And this too?
  { // FIXME : construction
    CGAL_kernel_precondition( orient != COLLINEAR);

    typename R::Compute_squared_distance_2 squared_distance;
    typename R::Construct_midpoint_2 midpoint;
    if (p != q) {
      Point_2 center = midpoint(p, q);
      base = Rep(center, squared_distance(p, center), orient);
    } else
      base = Rep(p, FT(0), orient);
  }

  CircleC2(const Point_2 &p, const Point_2 &q, const Point_2 &r)
  { // FIXME : construction
    typename R::Orientation_2 orientation;
    typename R::Compute_squared_distance_2 squared_distance;
    typename R::Construct_circumcenter_2 circumcenter;
    Orientation orient = orientation(p, q, r);
    CGAL_kernel_precondition( orient != COLLINEAR);

    Point_2 center = circumcenter(p, q, r);
    base = Rep(center, squared_distance(p, center), orient);
  }

  bool           operator==(const CircleC2 &s) const;
  bool           operator!=(const CircleC2 &s) const;

  const Point_2 & center() const
  {
    return get(base).first;
  }

  const FT & squared_radius() const
  {
    return get(base).second;
  }

  Orientation orientation() const
  {
    return get(base).third;
  }

  Circle_2           opposite() const;

  Circle_2           orthogonal_transform(const Aff_transformation_2 &t) const;

  Oriented_side  oriented_side(const Point_2 &p) const;
  Bounded_side   bounded_side(const Point_2 &p) const;

  bool           has_on_boundary(const Point_2 &p) const;
  bool           has_on_negative_side(const Point_2 &p) const;
  bool           has_on_positive_side(const Point_2 &p) const;

  bool           has_on_bounded_side(const Point_2 &p) const;
  bool           has_on_unbounded_side(const Point_2 &p) const;

  bool           is_degenerate() const;

  Bbox_2         bbox() const;
};

template < class R >
CGAL_KERNEL_INLINE
bool
CircleC2<R>::operator==(const CircleC2<R> &c) const
{ // FIXME : predicate
  if (CGAL::identical(base, c.base))
      return true;
  return center() == c.center() &&
         squared_radius() == c.squared_radius() &&
         orientation() == c.orientation();
}

template < class R >
inline
bool
CircleC2<R>::operator!=(const CircleC2<R> &c) const
{
  return !(*this == c);
}

template < class R >
CGAL_KERNEL_MEDIUM_INLINE
Oriented_side
CircleC2<R>::
oriented_side(const typename CircleC2<R>::Point_2 &p) const
{
  return Oriented_side(bounded_side(p) * orientation());
}

template < class R >
CGAL_KERNEL_INLINE
Bounded_side
CircleC2<R>::
bounded_side(const typename CircleC2<R>::Point_2 &p) const
{
  typename R::Compute_squared_distance_2 squared_distance;
  return Bounded_side(CGAL_NTS compare(squared_radius(),
                                       squared_distance(center(),p)));
}

template < class R >
inline
bool
CircleC2<R>::
has_on_boundary(const typename CircleC2<R>::Point_2 &p) const
{
  return bounded_side(p) == ON_BOUNDARY;
}

template < class R >
inline
bool
CircleC2<R>::
has_on_bounded_side(const typename CircleC2<R>::Point_2 &p) const
{
  return bounded_side(p) == ON_BOUNDED_SIDE;
}

template < class R >
inline
bool
CircleC2<R>::
has_on_unbounded_side(const typename CircleC2<R>::Point_2 &p) const
{
  return bounded_side(p) == ON_UNBOUNDED_SIDE;
}

template < class R >
CGAL_KERNEL_INLINE
bool
CircleC2<R>::
has_on_negative_side(const typename CircleC2<R>::Point_2 &p) const
{
  if (orientation() == COUNTERCLOCKWISE)
    return has_on_unbounded_side(p);
  return has_on_bounded_side(p);
}

template < class R >
CGAL_KERNEL_INLINE
bool
CircleC2<R>::
has_on_positive_side(const typename CircleC2<R>::Point_2 &p) const
{
  if (orientation() == COUNTERCLOCKWISE)
    return has_on_bounded_side(p);
  return has_on_unbounded_side(p);
}

template < class R >
inline
bool
CircleC2<R>::is_degenerate() const
{
  return CGAL_NTS is_zero(squared_radius());
}

template < class R >
inline
typename CircleC2<R>::Circle_2
CircleC2<R>::opposite() const
{
  return CircleC2<R>(center(),
                               squared_radius(),
                               CGAL::opposite(orientation()) );
}

template < class R >
CGAL_KERNEL_INLINE
Bbox_2
CircleC2<R>::bbox() const
{
  typename R::Construct_bbox_2 construct_bbox_2;
  Bbox_2 b = construct_bbox_2(center());

  Interval_nt<> x (b.xmin(), b.xmax());
  Interval_nt<> y (b.ymin(), b.ymax());

  Interval_nt<> sqr = CGAL_NTS to_interval(squared_radius());
  Interval_nt<> r = CGAL::sqrt(sqr);
  Interval_nt<> minx = x-r;
  Interval_nt<> maxx = x+r;
  Interval_nt<> miny = y-r;
  Interval_nt<> maxy = y+r;

  return Bbox_2(minx.inf(), miny.inf(), maxx.sup(), maxy.sup());
}

template < class R >
CGAL_KERNEL_INLINE
typename CircleC2<R>::Circle_2
CircleC2<R>::orthogonal_transform
  (const typename CircleC2<R>::Aff_transformation_2 &t) const
{
  typename R::Vector_2 vec(FT(1), FT(0) );   // unit vector
  vec = vec.transform(t);             // transformed
  FT sq_scale = vec.squared_length();       // squared scaling factor

  return CircleC2<R>(t.transform(center()),
                               sq_scale * squared_radius(),
                               t.is_even() ? orientation()
                                           : CGAL::opposite(orientation()));
}

#ifndef CGAL_NO_OSTREAM_INSERT_CIRCLEC2
template < class R >
CGAL_KERNEL_INLINE
std::ostream &
operator<<(std::ostream &os, const CircleC2<R> &c)
{
    switch(os.iword(IO::mode)) {
    case IO::ASCII :
        os << c.center() << ' ' << c.squared_radius() << ' '
           << static_cast<int>(c.orientation());
        break;
    case IO::BINARY :
        os << c.center();
        write(os, c.squared_radius());
        write(os, static_cast<int>(c.orientation()));
        break;
    default:
        os << "CircleC2(" << c.center() <<  ", " << c.squared_radius() ;
        switch (c.orientation()) {
        case CLOCKWISE:
            os << ", clockwise)";
            break;
        case COUNTERCLOCKWISE:
            os << ", counterclockwise)";
            break;
        default:
            os << ", collinear)";
            break;
        }
        break;
    }
    return os;
}
#endif // CGAL_NO_OSTREAM_INSERT_CIRCLEC2

#ifndef CGAL_NO_ISTREAM_EXTRACT_CIRCLEC2
template < class R >
CGAL_KERNEL_INLINE
std::istream&
operator>>(std::istream &is, CircleC2<R> &c)
{
    typename R::Point_2 center;
    typename R::FT squared_radius;
    int o;
    switch(is.iword(IO::mode)) {
    case IO::ASCII :
        is >> center >> squared_radius >> o;
        break;
    case IO::BINARY :
        is >> center;
        read(is, squared_radius);
        is >> o;
        break;
    default:
        std::cerr << "" << std::endl;
        std::cerr << "Stream must be in ascii or binary mode" << std::endl;
        break;
    }
    if (is)
	c = CircleC2<R>(center, squared_radius,
		                  static_cast<Orientation>(o));
    return is;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_CIRCLEC2

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_CIRCLE_2_H
