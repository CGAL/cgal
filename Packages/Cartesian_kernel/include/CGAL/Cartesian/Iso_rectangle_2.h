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

#ifndef CGAL_CARTESIAN_ISO_RECTANGLE_2_H
#define CGAL_CARTESIAN_ISO_RECTANGLE_2_H

#include <CGAL/Twotuple.h>

CGAL_BEGIN_NAMESPACE

template <class R_>
class Iso_rectangleC2
{
  typedef typename R_::FT                   FT;
  typedef typename R_::Point_2              Point_2;
  typedef typename R_::Iso_rectangle_2      Iso_rectangle_2;
  typedef typename R_::Aff_transformation_2 Aff_transformation_2;
  typedef typename R_::Construct_point_2    Construct_point_2;

  typedef Twotuple<Point_2>                        Rep;
  typedef typename R_::template Handle<Rep>::type  Base;

  Base base;

public:
  typedef R_                                     R;

  Iso_rectangleC2() {}

  Iso_rectangleC2(const Point_2 &p, const Point_2 &q)
  {
    FT minx, maxx, miny, maxy;
    if (p.x() < q.x()) { minx = p.x(); maxx = q.x(); }
    else               { minx = q.x(); maxx = p.x(); }
    if (p.y() < q.y()) { miny = p.y(); maxy = q.y(); }
    else               { miny = q.y(); maxy = p.y(); }
    Construct_point_2 construct_point_2;
    base = Rep(construct_point_2(minx, miny),
	       construct_point_2(maxx, maxy));
  }

  Iso_rectangleC2(const Point_2 &left, const Point_2 &right,
                  const Point_2 &bottom, const Point_2 &top)
    : base(Construct_point_2()(left.x(), bottom.y()),
           Construct_point_2()(right.x(), top.y()))
  {
    CGAL_kernel_assertion_code(typename R::Less_x_2 less_x;)
    CGAL_kernel_assertion_code(typename R::Less_y_2 less_y;)
    CGAL_kernel_assertion(!less_x(right, left));
    CGAL_kernel_assertion(!less_y(top, bottom));
  }

  Iso_rectangleC2(const FT& min_x, const FT& min_y, 
                  const FT& max_x, const FT& max_y)
    : base(Construct_point_2()(min_x, min_y),
           Construct_point_2()(max_x, max_y))
  {
    CGAL_kernel_precondition(min_x <= max_x);
    CGAL_kernel_precondition(min_y <= max_y);
  }

  Iso_rectangleC2(const FT& min_hx, const FT& min_hy, 
                  const FT& max_hx, const FT& max_hy, const FT& hw)
  {
    Construct_point_2 construct_point_2;
    if (hw == FT(1))
       base = Rep(construct_point_2(min_hx, min_hy),
	          construct_point_2(max_hx, max_hy));
    else
       base = Rep(construct_point_2(min_hx/hw, min_hy/hw),
	          construct_point_2(max_hx/hw, max_hy/hw));
  }

  bool            operator==(const Iso_rectangleC2 &s) const;
  bool            operator!=(const Iso_rectangleC2 &s) const;

  const Point_2 & min() const
  {
      return get(base).e0;
  }
  const Point_2 & max() const
  {
      return get(base).e1;
  }
  Point_2 vertex(int i) const;
  Point_2 operator[](int i) const;

  Iso_rectangle_2 transform(const Aff_transformation_2 &t) const
  {
    // FIXME : We need a precondition like this!!!
    // CGAL_kernel_precondition(t.is_axis_preserving());
    return Iso_rectangleC2<R>(t.transform(vertex(0)), t.transform(vertex(2)));
  }

  Bounded_side    bounded_side(const Point_2 &p) const;
  bool            has_on_boundary(const Point_2 &p) const;
  bool            has_on_bounded_side(const Point_2 &p) const;
  bool            has_on_unbounded_side(const Point_2 &p) const;

  bool            is_degenerate() const;

  Bbox_2          bbox() const;

  const FT &      xmin() const;
  const FT &      ymin() const;
  const FT &      xmax() const;
  const FT &      ymax() const;
  const FT &      min_coord(int i) const;
  const FT &      max_coord(int i) const;

  FT              area() const;
};

template < class R >
inline
bool
Iso_rectangleC2<R>::
operator==(const Iso_rectangleC2<R> &r) const
{
  if (CGAL::identical(base, r.base))
      return true;
  return vertex(0) == r.vertex(0) && vertex(2) == r.vertex(2);
}

template < class R >
inline
bool
Iso_rectangleC2<R>::
operator!=(const Iso_rectangleC2<R> &r) const
{
  return !(*this == r);
}

template < class R >
inline
const typename Iso_rectangleC2<R>::FT &
Iso_rectangleC2<R>::xmin() const
{
  return min().x();
}

template < class R >
inline
const typename Iso_rectangleC2<R>::FT &
Iso_rectangleC2<R>::ymin() const
{
  return min().y();
}

template < class R >
inline
const typename Iso_rectangleC2<R>::FT &
Iso_rectangleC2<R>::xmax() const
{
  return max().x();
}

template < class R >
inline
const typename Iso_rectangleC2<R>::FT &
Iso_rectangleC2<R>::ymax() const
{
  return max().y();
}

template < class R >
inline
const typename Iso_rectangleC2<R>::FT &
Iso_rectangleC2<R>::min_coord(int i) const
{
  CGAL_kernel_precondition( i == 0 || i == 1 );
  if (i == 0)
     return xmin();
  else
     return ymin();
}

template < class R >
inline
const typename Iso_rectangleC2<R>::FT &
Iso_rectangleC2<R>::max_coord(int i) const
{
  CGAL_kernel_precondition( i == 0 || i == 1 );
  if (i == 0)
     return xmax();
  else
     return ymax();
}

template < class R >
typename Iso_rectangleC2<R>::Point_2
Iso_rectangleC2<R>::vertex(int i) const
{
  Construct_point_2 construct_point_2;
  switch (i%4) {
  case 0: return min();
  case 1: return construct_point_2(xmax(), ymin());
  case 2: return max();
  default: return construct_point_2(xmin(), ymax());
  }
}

template < class R >
inline
typename Iso_rectangleC2<R>::Point_2
Iso_rectangleC2<R>::operator[](int i) const
{
  return vertex(i);
}

template < class R >
inline
typename Iso_rectangleC2<R>::FT
Iso_rectangleC2<R>::area() const
{
  return (xmax()-xmin()) * (ymax()-ymin());
}

template < class R >
CGAL_KERNEL_MEDIUM_INLINE
Bounded_side
Iso_rectangleC2<R>::
bounded_side(const typename Iso_rectangleC2<R>::Point_2 &p) const
{ // FIXME : predicate
  bool x_incr = (xmin() < p.x()) && (p.x() < xmax()),
       y_incr = (ymin() < p.y()) && (p.y() < ymax());
  if (x_incr)
    {
      if (y_incr)
          return ON_BOUNDED_SIDE;
      if ( (p.y() == ymin()) || (ymax() == p.y()) )
          return ON_BOUNDARY;
    }
  if ( (p.x() == xmin()) || (xmax() == p.x()) )
      if ( y_incr || (p.y() == ymin()) || (ymax() == p.y()) )
          return ON_BOUNDARY;

  return ON_UNBOUNDED_SIDE;
}

template < class R >
inline
bool
Iso_rectangleC2<R>::
has_on_boundary(const typename Iso_rectangleC2<R>::Point_2 &p) const
{
  return bounded_side(p) == ON_BOUNDARY;
}

template < class R >
inline
bool
Iso_rectangleC2<R>::
has_on_bounded_side(const typename Iso_rectangleC2<R>::Point_2 &p)
    const
{
  return bounded_side(p) == ON_BOUNDED_SIDE;
}

template < class R >
inline
bool
Iso_rectangleC2<R>::
has_on_unbounded_side(const typename Iso_rectangleC2<R>::Point_2 &p)
    const
{
  return bounded_side(p) == ON_UNBOUNDED_SIDE;
}

template < class R >
inline
bool
Iso_rectangleC2<R>::is_degenerate() const
{
  return (xmin() == xmax()) || (ymin() == ymax());
}

template < class R >
inline
Bbox_2
Iso_rectangleC2<R>::bbox() const
{ 
  typename R::Construct_bbox_2 construct_bbox_2;
  return construct_bbox_2(min()) + construct_bbox_2(max());
}

#ifndef CGAL_NO_OSTREAM_INSERT_ISO_RECTANGLEC2
template < class R >
std::ostream &
operator<<(std::ostream &os, const Iso_rectangleC2<R> &r)
{
    switch(os.iword(IO::mode)) {
    case IO::ASCII :
        return os << r[0] << ' ' << r[2];
    case IO::BINARY :
        return os << r[0] << r[2];
    default:
        return os << "Iso_rectangleC2(" << r[0] << ", " << r[2] << ")";
    }
}
#endif // CGAL_NO_OSTREAM_INSERT_ISO_RECTANGLEC2

#ifndef CGAL_NO_ISTREAM_EXTRACT_ISO_RECTANGLEC2
template < class R >
CGAL_KERNEL_MEDIUM_INLINE
std::istream &
operator>>(std::istream &is, Iso_rectangleC2<R> &r)
{
    typename R::Point_2 p, q;

    is >> p >> q;

    if (is)
	r = Iso_rectangleC2<R>(p, q);
    return is;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_ISO_RECTANGLEC2

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_ISO_RECTANGLE_2_H
