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
// Author(s)     : Andreas Fabri

#ifndef CGAL_CARTESIAN_SEGMENT_3_H
#define CGAL_CARTESIAN_SEGMENT_3_H

#include <CGAL/array.h>
#include <CGAL/Handle_for.h>

namespace CGAL {

template < class R_ >
class SegmentC3
{
  typedef typename R_::Point_3              Point_3;
  typedef typename R_::Direction_3          Direction_3;
  typedef typename R_::Vector_3             Vector_3;
  typedef typename R_::Line_3               Line_3;
  typedef typename R_::Segment_3            Segment_3;

  typedef std::array<Point_3, 2>          Rep;
  typedef typename R_::template Handle<Rep>::type  Base;

  Base base;

public:
  typedef R_                                     R;

  SegmentC3() {}

  SegmentC3(const Point_3 &sp, const Point_3 &ep)
    : base(CGAL::make_array(sp, ep)) {}

  bool        has_on(const Point_3 &p) const;
  bool        collinear_has_on(const Point_3 &p) const;

  bool        operator==(const SegmentC3 &s) const;
  bool        operator!=(const SegmentC3 &s) const;

  const Point_3 &   source() const
  {
      return get_pointee_or_identity(base)[0];
  }
  const Point_3 &   target() const
  {
      return get_pointee_or_identity(base)[1];
  }

  const Point_3 &   start() const;
  const Point_3 &   end() const;

  const Point_3 &   min BOOST_PREVENT_MACRO_SUBSTITUTION () const;
  const Point_3 &   max BOOST_PREVENT_MACRO_SUBSTITUTION () const;
  const Point_3 &   vertex(int i) const;
  const Point_3 &   point(int i) const;
  const Point_3 &   operator[](int i) const;

  Direction_3 direction() const;
  Vector_3    to_vector() const;
  Line_3      supporting_line() const;
  Segment_3   opposite() const;

  bool        is_degenerate() const;
};

template < class R >
inline
bool
SegmentC3<R>::operator==(const SegmentC3<R> &s) const
{
  if (CGAL::identical(base, s.base))
      return true;
  return source() == s.source() && target() == s.target();
}

template < class R >
inline
bool
SegmentC3<R>::operator!=(const SegmentC3<R> &s) const
{
  return !(*this == s);
}

template < class R >
const typename SegmentC3<R>::Point_3 &
SegmentC3<R>::start() const
{
  return source();
}

template < class R >
const typename SegmentC3<R>::Point_3 &
SegmentC3<R>::end() const
{
  return target();
}

template < class R >
inline
const typename SegmentC3<R>::Point_3 &
SegmentC3<R>::min BOOST_PREVENT_MACRO_SUBSTITUTION () const
{
  return lexicographically_xyz_smaller(source(),target()) ? source()
                                                          : target();
}

template < class R >
inline
const typename SegmentC3<R>::Point_3 &
SegmentC3<R>::max BOOST_PREVENT_MACRO_SUBSTITUTION () const
{
  return lexicographically_xyz_smaller(source(),target()) ? target()
                                                          : source();
}

template < class R >
inline
const typename SegmentC3<R>::Point_3 &
SegmentC3<R>::vertex(int i) const
{
  return (i%2 == 0) ? source() : target();
}

template < class R >
inline
const typename SegmentC3<R>::Point_3 &
SegmentC3<R>::point(int i) const
{
  return vertex(i);
}

template < class R >
inline
const typename SegmentC3<R>::Point_3 &
SegmentC3<R>::operator[](int i) const
{
  return vertex(i);
}

template < class R >
inline
typename SegmentC3<R>::Vector_3
SegmentC3<R>::to_vector() const
{
  return target() - source();
}

template < class R >
inline
typename SegmentC3<R>::Direction_3
SegmentC3<R>::direction() const
{
  return Direction_3( target() - source() );
}

template < class R >
inline
typename SegmentC3<R>::Line_3
SegmentC3<R>::supporting_line() const
{
  return Line_3(*this);
}

template < class R >
inline
typename SegmentC3<R>::Segment_3
SegmentC3<R>::opposite() const
{
  return SegmentC3<R>(target(), source());
}

template < class R >
inline
bool
SegmentC3<R>::is_degenerate() const
{
  return source() == target();
}

template < class R >
inline
bool
SegmentC3<R>::
has_on(const typename SegmentC3<R>::Point_3 &p) const
{
  return are_ordered_along_line(source(), p, target());
}

template < class R >
inline
bool
SegmentC3<R>::
collinear_has_on(const typename SegmentC3<R>::Point_3 &p) const
{
  return collinear_are_ordered_along_line(source(), p, target());
}

} //namespace CGAL

#endif // CGAL_CARTESIAN_SEGMENT_3_H
