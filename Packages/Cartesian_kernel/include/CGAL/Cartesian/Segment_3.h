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
// file          : include/CGAL/Cartesian/Segment_3.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas Fabri
// coordinator   : INRIA Sophia-Antipolis (Mariette.Yvinec@sophia.inria.fr)
//
// ======================================================================

#ifndef CGAL_CARTESIAN_SEGMENT_3_H
#define CGAL_CARTESIAN_SEGMENT_3_H

#include <CGAL/Twotuple.h>

CGAL_BEGIN_NAMESPACE

template < class R_ >
class SegmentC3
  : public R_::template Handle<Twotuple<typename R_::Point_3> >::type
{
CGAL_VC7_BUG_PROTECTED
  typedef typename R_::FT                   FT;
  typedef typename R_::Point_3              Point_3;
  typedef typename R_::Direction_3          Direction_3;
  typedef typename R_::Line_3               Line_3;
  typedef typename R_::Segment_3            Segment_3;
  typedef typename R_::Aff_transformation_3 Aff_transformation_3;

  typedef Twotuple<Point_3>                        rep;
  typedef typename R_::template Handle<rep>::type  base;

public:
  typedef R_                                     R;

  SegmentC3()
    : base(rep()) {}

  SegmentC3(const Point_3 &sp, const Point_3 &ep)
    : base(rep(sp, ep)) {}

  bool        has_on(const Point_3 &p) const;
  bool        collinear_has_on(const Point_3 &p) const;

  bool        operator==(const SegmentC3 &s) const;
  bool        operator!=(const SegmentC3 &s) const;

  const Point_3 &   source() const
  {
      return Ptr()->e0;
  }
  const Point_3 &   target() const
  {
      return Ptr()->e1;
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
  Line_3      supporting_line() const;
  Segment_3   opposite() const;
  Segment_3   transform(const Aff_transformation_3 &t) const
  {
    return SegmentC3<R>(t.transform(source()), t.transform(target()));
  }

  bool        is_degenerate() const;
  Bbox_3      bbox() const;
};

#ifdef CGAL_CFG_TYPENAME_BUG
#define typename
#endif

template < class R >
inline
bool
SegmentC3<R>::operator==(const SegmentC3<R> &s) const
{
  if (identical(s))
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
  return source().bbox() + target().bbox();
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

#ifdef CGAL_CFG_TYPENAME_BUG
#undef typename
#endif

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_SEGMENT_3_H
