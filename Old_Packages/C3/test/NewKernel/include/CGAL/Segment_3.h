// ============================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
// 
// release       :
// release_date  :
// 
// source        : Segment_3.fw
// file          : Segment_3.h
// revision      : 2.4
// revision_date : 24 Aug 1999 
// author(s)     : Andreas Fabri
//                 Stefan Schirra
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ============================================================================
 

#ifndef CGAL_SEGMENT_3_H
#define CGAL_SEGMENT_3_H

#ifndef CGAL_REP_CLASS_DEFINED
#error  no representation class defined
#endif  // CGAL_REP_CLASS_DEFINED

#ifdef CGAL_HOMOGENEOUS_H
#ifndef CGAL_SEGMENTH3_H
#include <CGAL/SegmentH3.h>
#endif // CGAL_SEGMENTH3_H
#endif // CGAL_HOMOGENEOUS_H

#ifdef CGAL_CARTESIAN_H
#ifndef CGAL_SEGMENTC3_H
#include <CGAL/SegmentC3.h>
#endif // CGAL_SEGMENTC3_H
#endif // CGAL_CARTESIAN_H

#ifndef CGAL_LINE_3_H
#include <CGAL/Line_3.h>
#endif // CGAL_LINE_3_H

CGAL_BEGIN_NAMESPACE

template <class _R>
class Segment_3 : public _R::Segment_3_base
{
public:
  typedef          _R                       R;
  typedef typename R::RT                    RT;
  typedef typename R::FT                    FT;
  typedef typename R::Segment_3_base  RSegment_3;

  Segment_3() : RSegment_3()
  {}
  Segment_3(const Segment_3<_R>& s) : RSegment_3(s)
  {}
  Segment_3(const Point_3<_R>& sp, const Point_3<_R>& ep)
    : RSegment_3(sp,ep)
  {}
  Segment_3(const RSegment_3&  s) : RSegment_3(s)
  {}

  Segment_3<_R>&   operator=(const Segment_3<_R>& s)
  {
    RSegment_3::operator=(s);
    return *this;
  }
  bool                 has_on(const Point_3<_R>& p) const
  { return RSegment_3::has_on(p); }
  bool                 operator==(const Segment_3<_R>& s) const
  { return RSegment_3::operator==(s); }
  bool                 operator!=(const Segment_3<_R>& s) const
  { return !(*this == s); }
  int                  id() const   /* XXX */
  { return (int) PTR ; }
  Point_3<_R>     start() const
  { return RSegment_3::start(); }
  Point_3<_R>     end() const
  { return RSegment_3::end(); }
  Point_3<_R>     source() const
  { return RSegment_3::source(); }
  Point_3<_R>     target() const
  { return RSegment_3::target(); }
  Point_3<_R>     min() const
  { return RSegment_3::min(); }
  Point_3<_R>     max() const
  { return RSegment_3::max(); }
  Point_3<_R>     vertex(int i) const
  { return RSegment_3::vertex(i); }
  Point_3<_R>     operator[](int i) const
  { return vertex(i); }
  FT                   squared_length() const
  { return RSegment_3::squared_length(); }
  Direction_3<_R> direction() const
  { return RSegment_3::direction(); }
  Segment_3<_R>  opposite() const
  { return Segment_3<_R>(target(),source()); }
  Segment_3<_R>  transform(const Aff_transformation_3<_R>& t) const
  { return RSegment_3::transform(t); }
  Line_3<_R>     supporting_line() const
  { return RSegment_3::supporting_line(); }
  bool                is_degenerate() const
  { return RSegment_3::is_degenerate(); }
  Bbox_3         bbox() const
  { return source().bbox() + target().bbox(); }
};


#ifndef NO_OSTREAM_INSERT_SEGMENT_3
template < class R>
std::ostream&
operator<<(std::ostream& os, const Segment_3<R>& s)
{
  typedef typename  R::Segment_3_base  RSegment_3;
  return os << (const RSegment_3& )s;
}
#endif // NO_OSTREAM_INSERT_SEGMENT_3

#ifndef NO_ISTREAM_EXTRACT_SEGMENT_3
template < class R>
std::istream&
operator>>(std::istream& is, Segment_3<R>& s)
{
  typedef typename  R::Segment_3_base  RSegment_3;
  return is >> (RSegment_3& )s;
}
#endif // NO_ISTREAM_EXTRACT_SEGMENT_3


CGAL_END_NAMESPACE


#endif // CGAL_SEGMENT_3_H
