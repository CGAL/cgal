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
// file          : include/CGAL/Old_style_kernel/Segment_3.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ============================================================================
 

#ifndef CGAL_OLD_STYLE_KERNEL_SEGMENT_3_H
#define CGAL_OLD_STYLE_KERNEL_SEGMENT_3_H

namespace CGAL {

template <class R_>
class Segment_3< R_, Old_style_tag> : public R_::Segment_3_base
{
 public:
  typedef  R_                        R;
  typedef typename R::RT             RT;
  typedef typename R::FT             FT;
  typedef typename R::Point_3_base   RPoint_3;
  typedef typename R::Segment_3_base RSegment_3;

  Segment_3() : RSegment_3() {}
  Segment_3(const RSegment_3& s) : RSegment_3(s) {}
  Segment_3(const RPoint_3& sp, const RPoint_3& ep) : RSegment_3(sp,ep) {}

  Line_3< R_, Old_style_tag>
  supporting_line() const
  { return RSegment_3::supporting_line(); }
};

} // namespace CGAL
#endif // CGAL_OLD_STYLE_KERNEL_SEGMENT_3_H
