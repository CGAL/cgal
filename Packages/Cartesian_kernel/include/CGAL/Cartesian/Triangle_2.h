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

#ifndef CGAL_CARTESIAN_TRIANGLE_2_H
#define CGAL_CARTESIAN_TRIANGLE_2_H

#include <CGAL/Threetuple.h>
#include <CGAL/Cartesian/predicates_on_points_2.h>

CGAL_BEGIN_NAMESPACE

template <class R_>
class TriangleC2
{
  typedef typename R_::FT                   FT;
  typedef typename R_::Point_2              Point_2;
  typedef typename R_::Vector_2             Vector_2;
  typedef typename R_::Triangle_2           Triangle_2;
  typedef typename R_::Aff_transformation_2 Aff_transformation_2;

  typedef Threetuple<Point_2>	                   Rep;
  typedef typename R_::template Handle<Rep>::type  Base;

  Base base;

public:
  typedef R_                                    R;

  TriangleC2() {}

  TriangleC2(const Point_2 &p, const Point_2 &q, const Point_2 &r)
    : base(p, q, r) {}

  bool           operator==(const TriangleC2 &s) const;
  bool           operator!=(const TriangleC2 &s) const;

  const Point_2 & vertex(int i) const;
  const Point_2 & operator[](int i) const;

  Triangle_2           opposite() const;
  Triangle_2           transform(const Aff_transformation_2 &t) const
  {
    return TriangleC2<R>(t.transform(vertex(0)),
                t.transform(vertex(1)),
                t.transform(vertex(2)));
  }

  Orientation    orientation() const;
  Oriented_side  oriented_side(const Point_2 &p) const;
  Bounded_side   bounded_side(const Point_2 &p) const;

  bool           has_on_boundary(const Point_2 &p) const;

  bool           has_on_bounded_side(const Point_2 &p) const;
  bool           has_on_unbounded_side(const Point_2 &p) const;

  bool           has_on_positive_side(const Point_2 &p) const;
  bool           has_on_negative_side(const Point_2 &p) const;

  bool           is_degenerate() const;

  Bbox_2         bbox() const;

  FT             area() const;
};

template < class R >
CGAL_KERNEL_MEDIUM_INLINE
bool
TriangleC2<R>::operator==(const TriangleC2<R> &t) const
{
  if (CGAL::identical(base, t.base))
      return true;

  int i;
  for(i=0; i<3; i++)
    if ( vertex(0) == t.vertex(i) )
      break;

  return (i<3) && vertex(1) == t.vertex(i+1) && vertex(2) == t.vertex(i+2);
}

template < class R >
inline
bool
TriangleC2<R>::operator!=(const TriangleC2<R> &t) const
{
  return !(*this == t);
}

template < class R >
CGAL_KERNEL_MEDIUM_INLINE
const typename TriangleC2<R>::Point_2 &
TriangleC2<R>::vertex(int i) const
{
  if (i>2) i = i%3;
  else if (i<0) i = (i%3) + 3;
  return (i==0) ? get(base).e0 :
         (i==1) ? get(base).e1 :
                  get(base).e2;
}

template < class R >
inline
const typename TriangleC2<R>::Point_2 &
TriangleC2<R>::operator[](int i) const
{
  return vertex(i);
}

template < class R >
CGAL_KERNEL_MEDIUM_INLINE
typename TriangleC2<R>::FT
TriangleC2<R>::area() const
{
  typename R::Compute_area_2 compute_area;
  return compute_area(vertex(0), vertex(1), vertex(2));
}

template < class R >
inline
Orientation
TriangleC2<R>::orientation() const
{
  typename R::Orientation_2 orientation;
  return orientation(vertex(0), vertex(1), vertex(2));
}

template < class R >
CGAL_KERNEL_LARGE_INLINE
Bounded_side
TriangleC2<R>::
bounded_side(const typename TriangleC2<R>::Point_2 &p) const
{
  typename R::Collinear_are_ordered_along_line_2 
    collinear_are_ordered_along_line;
  typename R::Orientation_2 orientation;
  Orientation o1 = orientation(vertex(0), vertex(1), p),
              o2 = orientation(vertex(1), vertex(2), p),
              o3 = orientation(vertex(2), vertex(3), p);

  if (o2 == o1 && o3 == o1)
    return ON_BOUNDED_SIDE;
  return
     (o1 == COLLINEAR
       && collinear_are_ordered_along_line(vertex(0), p, vertex(1))) ||
     (o2 == COLLINEAR
       && collinear_are_ordered_along_line(vertex(1), p, vertex(2))) ||
     (o3 == COLLINEAR
       && collinear_are_ordered_along_line(vertex(2), p, vertex(3)))
     ? ON_BOUNDARY
     : ON_UNBOUNDED_SIDE;
}

template < class R >
CGAL_KERNEL_LARGE_INLINE
Oriented_side
TriangleC2<R>::
oriented_side(const typename TriangleC2<R>::Point_2 &p) const
{
  typename R::Collinear_are_ordered_along_line_2 
    collinear_are_ordered_along_line;
  typename R::Orientation_2 orientation;
  // depends on the orientation of the vertices
  Orientation o1 = orientation(vertex(0), vertex(1), p),
              o2 = orientation(vertex(1), vertex(2), p),
              o3 = orientation(vertex(2), vertex(3), p),
              ot = orientation(vertex(0), vertex(1), vertex(2));

  if (o1 == ot && o2 == ot && o3 == ot) // ot cannot be COLLINEAR
    return Oriented_side(ot);
  return
     (o1 == COLLINEAR
       && collinear_are_ordered_along_line(vertex(0), p, vertex(1))) ||
     (o2 == COLLINEAR
       && collinear_are_ordered_along_line(vertex(1), p, vertex(2))) ||
     (o3 == COLLINEAR
       && collinear_are_ordered_along_line(vertex(2), p, vertex(3)))
     ? ON_ORIENTED_BOUNDARY
     : Oriented_side(-ot);
}

template < class R >
CGAL_KERNEL_LARGE_INLINE
bool
TriangleC2<R>::
has_on_bounded_side(const typename TriangleC2<R>::Point_2 &p) const
{
  return bounded_side(p) == ON_BOUNDED_SIDE;
}

template < class R >
CGAL_KERNEL_LARGE_INLINE
bool
TriangleC2<R>::
has_on_unbounded_side(const typename TriangleC2<R>::Point_2 &p) const
{
  return bounded_side(p) == ON_UNBOUNDED_SIDE;
}

template < class R >
inline
bool
TriangleC2<R>::
has_on_boundary(const typename TriangleC2<R>::Point_2 &p) const
{
  return bounded_side(p) == ON_BOUNDARY;
}

template < class R >
inline
bool
TriangleC2<R>::
has_on_negative_side(const typename TriangleC2<R>::Point_2 &p) const
{
  return oriented_side(p) == ON_NEGATIVE_SIDE;
}

template < class R >
inline
bool
TriangleC2<R>::
has_on_positive_side(const typename TriangleC2<R>::Point_2 &p) const
{
  return oriented_side(p) == ON_POSITIVE_SIDE;
}

template < class R >
inline
bool
TriangleC2<R>::is_degenerate() const
{
  typename R::Collinear_2 collinear;
  return collinear(vertex(0), vertex(1), vertex(2));
}

template < class R >
inline
Bbox_2
TriangleC2<R>::bbox() const
{
  typename R::Construct_bbox_2 construct_bbox_2;
  return construct_bbox_2(vertex(0)) 
       + construct_bbox_2(vertex(1)) 
       + construct_bbox_2(vertex(2));
}

template < class R >
inline
typename TriangleC2<R>::Triangle_2
TriangleC2<R>::opposite() const
{
  return TriangleC2<R>(vertex(0), vertex(2), vertex(1));
}

#ifndef CGAL_NO_OSTREAM_INSERT_TRIANGLEC2
template < class R >
std::ostream &
operator<<(std::ostream &os, const TriangleC2<R> &t)
{
    switch(os.iword(IO::mode)) {
    case IO::ASCII :
        return os << t[0] << ' ' << t[1] << ' ' << t[2];
    case IO::BINARY :
        return os << t[0] << t[1]  << t[2];
    default:
        return os<< "TriangleC2(" << t[0] << ", " 
		 << t[1] << ", " << t[2] <<")";
    }
}
#endif // CGAL_NO_OSTREAM_INSERT_TRIANGLEC2

#ifndef CGAL_NO_ISTREAM_EXTRACT_TRIANGLEC2
template < class R >
std::istream &
operator>>(std::istream &is, TriangleC2<R> &t)
{
    typename R::Point_2 p, q, r;

    is >> p >> q >> r;

    if (is)
	t = TriangleC2<R>(p, q, r);
    return is;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_TRIANGLEC2

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_TRIANGLE_2_H
