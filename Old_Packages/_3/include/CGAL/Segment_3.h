// ======================================================================
//
// Copyright (c) 1999 The CGAL Consortium
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
// file          : Segment_3.h
// package       : _3
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas Fabri
//                 Stefan Schirra
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ======================================================================
 

#ifndef CGAL_SEGMENT_3_H
#define CGAL_SEGMENT_3_H

#ifndef CGAL_REP_CLASS_DEFINED
#error  no representation class defined
#endif  // CGAL_REP_CLASS_DEFINED

#ifdef CGAL_HOMOGENEOUS_H
#include <CGAL/SegmentH3.h>
#endif

#if defined CGAL_CARTESIAN_H || defined CGAL_SIMPLE_CARTESIAN_H
#include <CGAL/Cartesian/Segment_3.h>
#endif

#include <CGAL/Line_3.h>

CGAL_BEGIN_NAMESPACE

template <class R_>
class Segment_3 : public R_::Segment_3_base
{
public:
  typedef          R_                       R;
  typedef typename R::RT                    RT;
  typedef typename R::FT                    FT;
  typedef typename R::Segment_3_base  RSegment_3;

  Segment_3() : RSegment_3()
  {}
  Segment_3(const CGAL::Segment_3<R>& s) : RSegment_3(s)
  {}
  Segment_3(const CGAL::Point_3<R>& sp, const CGAL::Point_3<R>& ep)
    : RSegment_3(sp,ep)
  {}
  Segment_3(const RSegment_3&  s) : RSegment_3(s)
  {}

  bool                 has_on(const CGAL::Point_3<R>& p) const
  { return RSegment_3::has_on(p); }
  bool                 operator==(const CGAL::Segment_3<R>& s) const
  { return RSegment_3::operator==(s); }
  bool                 operator!=(const CGAL::Segment_3<R>& s) const
  { return !(*this == s); }
  CGAL::Point_3<R>     start() const
  { return RSegment_3::start(); }
  CGAL::Point_3<R>     end() const
  { return RSegment_3::end(); }
  CGAL::Point_3<R>     source() const
  { return RSegment_3::source(); }
  CGAL::Point_3<R>     target() const
  { return RSegment_3::target(); }
  CGAL::Point_3<R>     min() const
  { return RSegment_3::min(); }
  CGAL::Point_3<R>     max() const
  { return RSegment_3::max(); }
  CGAL::Point_3<R>     vertex(int i) const
  { return RSegment_3::vertex(i); }
  CGAL::Point_3<R>     operator[](int i) const
  { return vertex(i); }
  FT                   squared_length() const
  { return RSegment_3::squared_length(); }
  CGAL::Direction_3<R> direction() const
  { return RSegment_3::direction(); }
  CGAL::Segment_3<R>  opposite() const
  { return CGAL::Segment_3<R>(target(),source()); }
  CGAL::Segment_3<R>  transform(const CGAL::Aff_transformation_3<R>& t) const
  { return RSegment_3::transform(t); }
  CGAL::Line_3<R>     supporting_line() const
  { return RSegment_3::supporting_line(); }
  bool                is_degenerate() const
  { return RSegment_3::is_degenerate(); }
  Bbox_3         bbox() const
  { return source().bbox() + target().bbox(); }
};


#ifndef CGAL_NO_OSTREAM_INSERT_SEGMENT_3
template < class R>
std::ostream&
operator<<(std::ostream& os, const Segment_3<R>& s)
{
  typedef typename  R::Segment_3_base  RSegment_3;
  return os << (const RSegment_3& )s;
}
#endif // CGAL_NO_OSTREAM_INSERT_SEGMENT_3

#ifndef CGAL_NO_ISTREAM_EXTRACT_SEGMENT_3
template < class R>
std::istream&
operator>>(std::istream& is, Segment_3<R>& s)
{
  typedef typename  R::Segment_3_base  RSegment_3;
  return is >> (RSegment_3& )s;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_SEGMENT_3

CGAL_END_NAMESPACE

#endif // CGAL_SEGMENT_3_H
