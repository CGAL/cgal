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
// file          : include/CGAL/Old_style_kernel/Segment_2.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ============================================================================
 

#ifndef CGAL_OLD_STYLE_KERNEL_SEGMENT_2_H
#define CGAL_OLD_STYLE_KERNEL_SEGMENT_2_H

namespace CGAL {

template <class R_>
class Segment_2< R_, Old_style_tag> : public R_::Segment_2_base
{
 public:
  typedef  R_                        R;
  typedef typename R::RT             RT;
  typedef typename R::FT             FT;
  typedef typename R::Point_2_base   RPoint_2;
  typedef typename R::Segment_2_base RSegment_2;

  Segment_2() : RSegment_2() {}
  Segment_2(const RSegment_2& s) : RSegment_2(s) {}
  Segment_2(const RPoint_2& sp, const RPoint_2& ep) : RSegment_2(sp,ep) {}
};

} // namespace CGAL
#endif // CGAL_OLD_STYLE_KERNEL_SEGMENT_2_H
