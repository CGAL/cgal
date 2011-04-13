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
// file          : Segment_2.h
// package       : _2
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas Fabri
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ======================================================================

#ifndef CGAL_SEGMENT_2_H
#define CGAL_SEGMENT_2_H

#include <CGAL/Point_2.h>
#include <CGAL/Line_2.h>
#include <CGAL/Direction_2.h>

CGAL_BEGIN_NAMESPACE

template <class R_>
class Segment_2 : public R_::Segment_2_base
{
public:
  typedef  R_                               R;
  typedef typename R::RT                    RT;
  typedef typename R::FT                    FT;
  typedef typename R::Segment_2_base  RSegment_2;

  Segment_2()     // doesn't the default constructor do the same ???
    : RSegment_2()  // does the handle stuff
  {}

  Segment_2(const CGAL::Segment_2<R>& s)
    : RSegment_2(static_cast<const RSegment_2&>(s))  // does the handle stuff
  {}

  Segment_2(const CGAL::Point_2<R> &sp, const CGAL::Point_2<R> &ep)
    :  RSegment_2(sp,ep)
  {}

  // conversion from implementation class object to interface class object
  Segment_2(const RSegment_2& s)
    : RSegment_2(s)  // does the handle stuff
  {}

  CGAL::Point_2<R>     start() const
  { return RSegment_2::start(); }

  CGAL::Point_2<R>     end() const
  { return RSegment_2::end(); }

  CGAL::Point_2<R>     source() const
  { return RSegment_2::source(); }

  CGAL::Point_2<R>     target() const
  { return RSegment_2::target(); }

  CGAL::Point_2<R>     min() const
  { return RSegment_2::min(); }

  CGAL::Point_2<R>     max() const
  { return RSegment_2::max(); }

  CGAL::Point_2<R>     vertex(int i) const
  { return RSegment_2::vertex(i); }

  CGAL::Point_2<R>     point(int i) const
  { return RSegment_2::vertex(i); }

  CGAL::Point_2<R>     operator[](int i) const
  { return vertex(i); }

  CGAL::Direction_2<R> direction() const
  { return RSegment_2::direction(); }

  CGAL::Segment_2<R>  opposite() const
  { return CGAL::Segment_2<R>(target(),source()); }

  // this makes use of the constructor of the interface class
  // taking an object of the implemetation class as argument.

  CGAL::Segment_2<R>   transform(const CGAL::Aff_transformation_2<R> &t) const
  { return  RSegment_2::transform(t); }

  CGAL::Line_2<R>      supporting_line() const
  { return RSegment_2::supporting_line(); }
};

#ifndef CGAL_NO_OSTREAM_INSERT_SEGMENT_2
template < class R>
std::ostream &
operator<<(std::ostream &os, const Segment_2<R> &s)
{
  typedef typename  R::Segment_2_base  RSegment_2;
  return os << static_cast<const RSegment_2&>(s);
}
#endif // CGAL_NO_OSTREAM_INSERT_SEGMENT_2

#ifndef CGAL_NO_ISTREAM_EXTRACT_SEGMENT_2
template < class R>
std::istream &
operator>>(std::istream &is, Segment_2<R> &s)
{
  typedef typename  R::Segment_2_base  RSegment_2;
  return is >> static_cast<RSegment_2&>(s);
}
#endif // CGAL_NO_ISTREAM_EXTRACT_SEGMENT_2

CGAL_END_NAMESPACE

#endif //  CGAL_SEGMENT_2_H
