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

#ifndef CGAL_CARTESIAN_SEGMENT_3_H
#define CGAL_CARTESIAN_SEGMENT_3_H

#include <CGAL/Twotuple.h>

CGAL_BEGIN_NAMESPACE

template < class R_ >
class SegmentC3
  : public R_::template Handle<Twotuple<typename R_::Point_3> >::type
{
  typedef typename R_::FT                   FT;
  typedef typename R_::Point_3              Point_3;
  typedef typename R_::Direction_3          Direction_3;
  typedef typename R_::Vector_3             Vector_3;
  typedef typename R_::Line_3               Line_3;
  typedef typename R_::Segment_3            Segment_3;
  typedef typename R_::Aff_transformation_3 Aff_transformation_3;

  typedef Twotuple<Point_3>                        rep;
  typedef typename R_::template Handle<rep>::type  base;

  const base& Base() const { return *this; }
  base& Base() { return *this; }

public:
  typedef R_                                     R;

  SegmentC3() {}

  SegmentC3(const Point_3 &sp, const Point_3 &ep)
    : base(sp, ep) {}

  bool        has_on(const Point_3 &p) const;
  bool        collinear_has_on(const Point_3 &p) const;

  bool        operator==(const SegmentC3 &s) const;
  bool        operator!=(const SegmentC3 &s) const;

  const Point_3 &   source() const
  {
      return get(Base()).e0;
  }
  const Point_3 &   target() const
  {
      return get(Base()).e1;
  }

  const Point_3 &   start() const;
  const Point_3 &   end() const;

  const Point_3 &   min() const;
  const Point_3 &   max() const;
  const Point_3 &   vertex(int i) const;
  const Point_3 &   point(int i) const;
  const Point_3 &   operator[](int i) const;

  FT          squared_length() const;

  Direction_3 direction() const;
  Vector_3    to_vector() const;
  Line_3      supporting_line() const;
  Segment_3   opposite() const;
  Segment_3   transform(const Aff_transformation_3 &t) const
  {
    return SegmentC3<R>(t.transform(source()), t.transform(target()));
  }

  bool        is_degenerate() const;
  Bbox_3      bbox() const;
};

template < class R >
inline
bool
SegmentC3<R>::operator==(const SegmentC3<R> &s) const
{
  if (CGAL::identical(Base(), s.Base()))
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
SegmentC3<R>::min() const
{
  return lexicographically_xyz_smaller(source(),target()) ? source()
                                                          : target();
}

template < class R >
inline
const typename SegmentC3<R>::Point_3 &
SegmentC3<R>::max() const
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
  return (i%2 == 0) ? source() : target();
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
typename SegmentC3<R>::FT
SegmentC3<R>::squared_length() const
{
  return squared_distance(target(), source());
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
Bbox_3
SegmentC3<R>::bbox() const
{
  typename R::Construct_bbox_3 construct_bbox_3;
  return construct_bbox_3(source()) + construct_bbox_3(target());
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

#ifndef CGAL_NO_OSTREAM_INSERT_SEGMENTC3
template < class R >
std::ostream &
operator<<(std::ostream &os, const SegmentC3<R> &s)
{
    switch(os.iword(IO::mode)) {
    case IO::ASCII :
        return os << s.source() << ' ' << s.target();
    case IO::BINARY :
        return os << s.source() << s.target();
    default:
        return os << "SegmentC3(" << s.source() <<  ", " << s.target() << ")";
    }
}
#endif // CGAL_NO_OSTREAM_INSERT_SEGMENTC3

#ifndef CGAL_NO_ISTREAM_EXTRACT_SEGMENTC3
template < class R >
std::istream &
operator>>(std::istream &is, SegmentC3<R> &s)
{
    typename R::Point_3 p, q;

    is >> p >> q;

    if (is)
	s = SegmentC3<R>(p, q);
    return is;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_SEGMENTC3

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_SEGMENT_3_H
