// ======================================================================
//
// Copyright (c) 2000 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       :
// release_date  :
//
// file          : include/CGAL/Cartesian/Segment_2.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas Fabri, Herve Bronnimann
// coordinator   : INRIA Sophia-Antipolis (Mariette.Yvinec@sophia.inria.fr)
//
// ======================================================================

#ifndef CGAL_CARTESIAN_SEGMENT_2_H
#define CGAL_CARTESIAN_SEGMENT_2_H

#include <CGAL/Twotuple.h>
#include <CGAL/Cartesian/predicates_on_points_2.h>

CGAL_BEGIN_NAMESPACE

template < class R_ >
class SegmentC2
  : public R_::template Handle<Twotuple<typename R_::Point_2> >::type
{
CGAL_VC7_BUG_PROTECTED
  typedef typename R_::FT                   FT;
  typedef typename R_::Point_2              Point_2;
  typedef typename R_::Direction_2          Direction_2;
  typedef typename R_::Line_2               Line_2;
  typedef typename R_::Segment_2            Segment_2;
  typedef typename R_::Aff_transformation_2 Aff_transformation_2;

  typedef Twotuple<Point_2>                        rep;
  typedef typename R_::template Handle<rep>::type  base;

public:
  typedef R_                                     R;

  SegmentC2()
    : base(rep()) {}

  SegmentC2(const Point_2 &sp, const Point_2 &ep)
    : base(rep(sp, ep)) {}

  bool        is_horizontal() const;
  bool        is_vertical() const;
  bool        has_on(const Point_2 &p) const;
  bool        collinear_has_on(const Point_2 &p) const;

  bool        operator==(const SegmentC2 &s) const;
  bool        operator!=(const SegmentC2 &s) const;

  const Point_2 &   source() const
  {
      return Ptr()->e0;
  }
  const Point_2 &   target() const
  {
      return Ptr()->e1;
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
  if (identical(s))
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
  return lexicographically_xy_smaller(source(),target()) ? source() : target();
}

template < class R >
CGAL_KERNEL_INLINE
const typename SegmentC2<R>::Point_2 &
SegmentC2<R>::max() const
{
  return lexicographically_xy_smaller(source(),target()) ? target() : source();
}

template < class R >
CGAL_KERNEL_INLINE
const typename SegmentC2<R>::Point_2 &
SegmentC2<R>::vertex(int i) const
{
  return (i%2 == 0) ? source() : target();
}

template < class R >
CGAL_KERNEL_INLINE
const typename SegmentC2<R>::Point_2 &
SegmentC2<R>::point(int i) const
{
  return (i%2 == 0) ? source() : target();
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
  return Direction_2( target() - source() );
}

template < class R >
inline
typename SegmentC2<R>::Line_2
SegmentC2<R>::supporting_line() const
{
  return Line_2(*this);
}

template < class R >
inline
typename SegmentC2<R>::Segment_2
SegmentC2<R>::opposite() const
{
  return SegmentC2<R>(target(), source());
}

template < class R >
CGAL_KERNEL_INLINE
Bbox_2
SegmentC2<R>::bbox() const
{
  return source().bbox() + target().bbox();
}

template < class R >
inline
bool
SegmentC2<R>::is_degenerate() const
{
  return source() == target();
}

template < class R >
CGAL_KERNEL_INLINE
bool
SegmentC2<R>::is_horizontal() const
{
  return y_equal(source(), target());
}

template < class R >
CGAL_KERNEL_INLINE
bool
SegmentC2<R>::is_vertical() const
{
  return x_equal(source(), target());
}

template < class R >
CGAL_KERNEL_INLINE
bool
SegmentC2<R>::
has_on(const typename SegmentC2<R>::Point_2 &p) const
{
  return are_ordered_along_line(source(), p, target());
}

template < class R >
CGAL_KERNEL_MEDIUM_INLINE
bool
SegmentC2<R>::
collinear_has_on(const typename SegmentC2<R>::Point_2 &p) const
{
    return collinear_are_ordered_along_line(source(), p, target());
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
