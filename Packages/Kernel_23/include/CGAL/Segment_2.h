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

CGAL_BEGIN_NAMESPACE

template <class R_>
class Segment_2 : public R_::Segment_2_base
{
  typedef typename R_::RT                    RT;
  typedef typename R_::Point_2               Point_2;
  typedef typename R_::Segment_2_base  RSegment_2;
public:
  typedef  R_                               R;

  Segment_2()
    : RSegment_2() {}

  Segment_2(const CGAL::Segment_2<R>& s)
    : RSegment_2(static_cast<const RSegment_2&>(s)) {}

  Segment_2(const Point_2 &sp, const Point_2 &ep)
    :  RSegment_2(sp,ep) {}

  // conversion from implementation class object to interface class object
  Segment_2(const RSegment_2& s)
    : RSegment_2(s) {}
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
