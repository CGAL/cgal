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

#ifndef CGAL_CARTESIAN_SEGMENT_2_H
#define CGAL_CARTESIAN_SEGMENT_2_H

#include <CGAL/Twotuple.h>
#include <CGAL/Cartesian/predicates_on_points_2.h>

CGAL_BEGIN_NAMESPACE

template < class R_ >
class SegmentC2
{
  typedef typename R_::FT                   FT;
  typedef typename R_::Point_2              Point_2;
  typedef typename R_::Vector_2             Vector_2;
  typedef typename R_::Direction_2          Direction_2;
  typedef typename R_::Line_2               Line_2;
  typedef typename R_::Segment_2            Segment_2;
  typedef typename R_::Aff_transformation_2 Aff_transformation_2;

  typedef Twotuple<Point_2>                        Rep;
  typedef typename R_::template Handle<Rep>::type  Base;

  Base base;

public:
  typedef R_                                     R;

  SegmentC2() {}

  SegmentC2(const Point_2 &sp, const Point_2 &ep)
    : base(sp, ep) {}

  bool        is_horizontal() const;
  bool        is_vertical() const;
  bool        has_on(const Point_2 &p) const;
  bool        collinear_has_on(const Point_2 &p) const;

  bool        operator==(const SegmentC2 &s) const;
  bool        operator!=(const SegmentC2 &s) const;

  const Point_2 &   source() const
  {
      return get(base).e0;
  }
  const Point_2 &   target() const
  {
      return get(base).e1;
  }
  
  const Point_2 &    start() const;
  const Point_2 &    end() const;

  const Point_2 &   min() const;
  const Point_2 &   max() const;
  const Point_2 &   vertex(int i) const;
  const Point_2 &   point(int i) const;
  const Point_2 &   operator[](int i) const;

  FT          squared_length() const;

  Direction_2 direction() const;
  Vector_2    to_vector() const;
  Line_2      supporting_line() const;
  Segment_2        opposite() const;
  Segment_2        transform(const Aff_transformation_2 &t) const
  {
    return SegmentC2<R>(t.transform(source()), t.transform(target()));
  }

  bool        is_degenerate() const;
  Bbox_2      bbox() const;
};

template < class R >
inline
bool
SegmentC2<R>::operator==(const SegmentC2<R> &s) const
{
  if (CGAL::identical(base, s.base))
      return true;
  return source() == s.source() && target() == s.target();
}

template < class R >
inline
bool
SegmentC2<R>::operator!=(const SegmentC2<R> &s) const
{
  return !(*this == s);
}

template < class R >
inline
const typename SegmentC2<R>::Point_2 &
SegmentC2<R>::start() const
{
  return source();
}

template < class R >
inline
const typename SegmentC2<R>::Point_2 &
SegmentC2<R>::end() const
{
  return target();
}

template < class R >
CGAL_KERNEL_INLINE
const typename SegmentC2<R>::Point_2 &
SegmentC2<R>::min() const
{
  typename R::Less_xy_2 less_xy; 
  return less_xy(source(),target()) ? source() : target();
}

template < class R >
CGAL_KERNEL_INLINE
const typename SegmentC2<R>::Point_2 &
SegmentC2<R>::max() const
{
  typename R::Less_xy_2 less_xy; 
  return less_xy(source(),target()) ? target() : source();
}

template < class R >
CGAL_KERNEL_INLINE
const typename SegmentC2<R>::Point_2 &
SegmentC2<R>::vertex(int i) const
{
  return (i%2 == 0) ? source() : target();
}

template < class R >
inline
const typename SegmentC2<R>::Point_2 &
SegmentC2<R>::point(int i) const
{
  return vertex(i);
}

template < class R >
inline
const typename SegmentC2<R>::Point_2 &
SegmentC2<R>::operator[](int i) const
{
  return vertex(i);
}

template < class R >
CGAL_KERNEL_INLINE
typename SegmentC2<R>::FT
SegmentC2<R>::squared_length() const
{
  return squared_distance(source(), target());
}

template < class R >
CGAL_KERNEL_INLINE
typename SegmentC2<R>::Direction_2
SegmentC2<R>::direction() const
{
  typename R::Construct_vector_2 construct_vector;
  return Direction_2( construct_vector( source(), target()));
}

template < class R >
CGAL_KERNEL_INLINE
typename SegmentC2<R>::Vector_2
SegmentC2<R>::to_vector() const
{
  typename R::Construct_vector_2 construct_vector;
  return construct_vector( source(), target());
}

template < class R >
inline
typename SegmentC2<R>::Line_2
SegmentC2<R>::supporting_line() const
{
  typename R::Construct_line_2 construct_line;

  return construct_line(*this);
}

template < class R >
inline
typename SegmentC2<R>::Segment_2
SegmentC2<R>::opposite() const
{
  return Segment_2(target(), source());
}

template < class R >
CGAL_KERNEL_INLINE
Bbox_2
SegmentC2<R>::bbox() const
{
  typename R::Construct_bbox_2 construct_bbox_2;
  return construct_bbox_2(source()) + construct_bbox_2(target());
}

template < class R >
inline
bool
SegmentC2<R>::is_degenerate() const
{
  return R().equal_2_object()(source(), target());
}

template < class R >
CGAL_KERNEL_INLINE
bool
SegmentC2<R>::is_horizontal() const
{
  return R().equal_y_2_object()(source(), target());
}

template < class R >
CGAL_KERNEL_INLINE
bool
SegmentC2<R>::is_vertical() const
{
  return R().equal_x_2_object()(source(), target());
}

template < class R >
CGAL_KERNEL_INLINE
bool
SegmentC2<R>::
has_on(const typename SegmentC2<R>::Point_2 &p) const
{
  return R().are_ordered_along_line_2_object()(source(), 
					       p, 
					       target());
}

template < class R >
inline
bool
SegmentC2<R>::
collinear_has_on(const typename SegmentC2<R>::Point_2 &p) const
{
  return R().collinear_has_on_2_object()
               (static_cast<const typename R::Segment_2>(*this), p);
}

#ifndef CGAL_NO_OSTREAM_INSERT_SEGMENTC2
template < class R >
std::ostream &
operator<<(std::ostream &os, const SegmentC2<R> &s)
{
    switch(os.iword(IO::mode)) {
    case IO::ASCII :
        return os << s.source() << ' ' << s.target();
    case IO::BINARY :
        return os << s.source() << s.target();
    default:
        return os << "SegmentC2(" << s.source() <<  ", " << s.target() << ")";
    }
}
#endif // CGAL_NO_OSTREAM_INSERT_SEGMENTC2

#ifndef CGAL_NO_ISTREAM_EXTRACT_SEGMENTC2
template < class R >
std::istream &
operator>>(std::istream &is, SegmentC2<R> &s)
{
    typename R::Point_2 p, q;

    is >> p >> q;

    if (is)
	s = SegmentC2<R>(p, q);
    return is;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_SEGMENTC2

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_SEGMENT_2_H
