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

#include <CGAL/Cartesian/redefine_names_3.h>

CGAL_BEGIN_NAMESPACE

template < class R_ >
class SegmentC3 CGAL_ADVANCED_KERNEL_PARTIAL_SPEC
  : public R_::Segment_handle_3
{
public:
  typedef R_                                    R;
  typedef typename R::FT                        FT;
  typedef typename R::RT                        RT;

  typedef typename R::Segment_handle_3          Segment_handle_3_;
  typedef typename Segment_handle_3_::element_type Segment_ref_3;

#ifndef CGAL_CFG_NO_ADVANCED_KERNEL
  typedef SegmentC3<R CGAL_CTAG>                Self;
  typedef typename R::Point_3                   Point_3;
  typedef typename R::Direction_3               Direction_3;
  typedef typename R::Line_3                    Line_3;
  typedef typename R::Aff_transformation_3      Aff_transformation_3;
#else
  typedef SegmentC3<R>                          Self;
  typedef typename R::Point_3_base              Point_3;
  typedef typename R::Direction_3_base          Direction_3;
  typedef typename R::Line_3_base               Line_3;
  typedef typename R::Aff_transformation_3_base Aff_transformation_3;
#endif

  SegmentC3()
    : Segment_handle_3_(Segment_ref_3()) {}

  SegmentC3(const Point_3 &sp, const Point_3 &ep)
    : Segment_handle_3_(Segment_ref_3(sp, ep)) {}

  bool        has_on(const Point_3 &p) const;
  bool        collinear_has_on(const Point_3 &p) const;

  bool        operator==(const Self &s) const;
  bool        operator!=(const Self &s) const;

  Point_3     source() const
  {
      return Ptr()->e0;
  }
  Point_3     target() const
  {
      return Ptr()->e1;
  }

  Point_3     start() const;
  Point_3     end() const;

  Point_3     min() const;
  Point_3     max() const;
  Point_3     vertex(int i) const;
  Point_3     point(int i) const;
  Point_3     operator[](int i) const;

  FT          squared_length() const;

  Direction_3 direction() const;
  Line_3      supporting_line() const;
  Self        opposite() const;
  Self        transform(const Aff_transformation_3 &t) const
  {
    return Self(t.transform(source()), t.transform(target()));
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
SegmentC3<R CGAL_CTAG>::operator==(const SegmentC3<R CGAL_CTAG> &s) const
{
  if (identical(s))
      return true;
  return source() == s.source() && target() == s.target();
}

template < class R >
inline
bool
SegmentC3<R CGAL_CTAG>::operator!=(const SegmentC3<R CGAL_CTAG> &s) const
{
  return !(*this == s);
}

template < class R >
typename SegmentC3<R CGAL_CTAG>::Point_3
SegmentC3<R CGAL_CTAG>::start() const
{
  return source();
}

template < class R >
typename SegmentC3<R CGAL_CTAG>::Point_3
SegmentC3<R CGAL_CTAG>::end() const
{
  return target();
}

template < class R >
inline
typename SegmentC3<R CGAL_CTAG>::Point_3
SegmentC3<R CGAL_CTAG>::min() const
{
  return lexicographically_xyz_smaller(source(),target()) ? source()
                                                          : target();
}

template < class R >
inline
typename SegmentC3<R CGAL_CTAG>::Point_3
SegmentC3<R CGAL_CTAG>::max() const
{
  return lexicographically_xyz_smaller(source(),target()) ? target()
                                                          : source();
}

template < class R >
inline
typename SegmentC3<R CGAL_CTAG>::Point_3
SegmentC3<R CGAL_CTAG>::vertex(int i) const
{
  return (i%2 == 0) ? source() : target();
}

template < class R >
inline
typename SegmentC3<R CGAL_CTAG>::Point_3
SegmentC3<R CGAL_CTAG>::point(int i) const
{
  return (i%2 == 0) ? source() : target();
}

template < class R >
inline
typename SegmentC3<R CGAL_CTAG>::Point_3
SegmentC3<R CGAL_CTAG>::operator[](int i) const
{
  return vertex(i);
}

template < class R >
inline
typename SegmentC3<R CGAL_CTAG>::FT
SegmentC3<R CGAL_CTAG>::squared_length() const
{
  return squared_distance(target(), source());
}

template < class R >
inline
typename SegmentC3<R CGAL_CTAG>::Direction_3
SegmentC3<R CGAL_CTAG>::direction() const
{
  return Direction_3( target() - source() );
}

template < class R >
inline
typename SegmentC3<R CGAL_CTAG>::Line_3
SegmentC3<R CGAL_CTAG>::supporting_line() const
{
  return Line_3(*this);
}

template < class R >
inline
SegmentC3<R CGAL_CTAG>
SegmentC3<R CGAL_CTAG>::opposite() const
{
  return SegmentC3<R CGAL_CTAG>(target(), source());
}

template < class R >
inline
bool
SegmentC3<R CGAL_CTAG>::is_degenerate() const
{
  return source() == target();
}

template < class R >
inline
Bbox_3
SegmentC3<R CGAL_CTAG>::bbox() const
{
  return source().bbox() + target().bbox();
}

template < class R >
inline
bool
SegmentC3<R CGAL_CTAG>::
has_on(const typename SegmentC3<R CGAL_CTAG>::Point_3 &p) const
{
  return are_ordered_along_line(source(), p, target());
}

template < class R >
inline
bool
SegmentC3<R CGAL_CTAG>::
collinear_has_on(const typename SegmentC3<R CGAL_CTAG>::Point_3 &p) const
{
  return collinear_are_ordered_along_line(source(), p, target());
}

#ifndef CGAL_NO_OSTREAM_INSERT_SEGMENTC3
template < class R >
std::ostream &
operator<<(std::ostream &os, const SegmentC3<R CGAL_CTAG> &s)
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
operator>>(std::istream &is, SegmentC3<R CGAL_CTAG> &s)
{
    typename SegmentC3<R CGAL_CTAG>::Point_3 p, q;

    is >> p >> q;

    if (is)
	s = SegmentC3<R CGAL_CTAG>(p, q);
    return is;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_SEGMENTC3

#ifdef CGAL_CFG_TYPENAME_BUG
#undef typename
#endif

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_SEGMENT_3_H
