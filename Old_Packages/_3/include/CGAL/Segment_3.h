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
// author(s)     : Andreas Fabri, Stefan Schirra
//
// coordinator   : MPI, Saarbruecken
// ======================================================================

#ifndef CGAL_SEGMENT_3_H
#define CGAL_SEGMENT_3_H

CGAL_BEGIN_NAMESPACE

template <class R_>
class Segment_3 : public R_::Kernel_base::Segment_3
{
  typedef typename R_::RT                    RT;
  typedef typename R_::Point_3               Point_3;
  typedef typename R_::Kernel_base::Segment_3  RSegment_3;
public:
  typedef          R_                       R;

  Segment_3()
      : RSegment_3() {}

  Segment_3(const CGAL::Segment_3<R>& s)
      : RSegment_3(s) {}

  Segment_3(const Point_3& sp, const Point_3& ep)
    : RSegment_3(sp,ep) {}

  Segment_3(const RSegment_3& s)
      : RSegment_3(s) {}
};

#ifndef CGAL_NO_OSTREAM_INSERT_SEGMENT_3
template < class R>
std::ostream&
operator<<(std::ostream& os, const Segment_3<R>& s)
{
  typedef typename  R::Kernel_base::Segment_3  RSegment_3;
  return os << static_cast<const RSegment_3&>(s);
}
#endif // CGAL_NO_OSTREAM_INSERT_SEGMENT_3

#ifndef CGAL_NO_ISTREAM_EXTRACT_SEGMENT_3
template < class R>
std::istream&
operator>>(std::istream& is, Segment_3<R>& s)
{
  typedef typename  R::Kernel_base::Segment_3  RSegment_3;
  return is >> static_cast<RSegment_3&>(s);
}
#endif // CGAL_NO_ISTREAM_EXTRACT_SEGMENT_3

CGAL_END_NAMESPACE

#endif // CGAL_SEGMENT_3_H
