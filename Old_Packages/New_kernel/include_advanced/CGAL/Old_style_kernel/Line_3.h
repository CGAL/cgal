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
// file          : include/CGAL/Old_style_kernel/Line_3.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ============================================================================
 

#ifndef CGAL_OLD_STYLE_KERNEL_LINE_3_H
#define CGAL_OLD_STYLE_KERNEL_LINE_3_H

namespace CGAL {

template <class R_>
class Line_3< R_, Old_style_tag> : public R_::Line_3_base
{
 public:
  typedef  R_                          R;
  typedef typename R::RT               RT;
  typedef typename R::FT               FT;
  typedef typename R::Point_3_base     RPoint_3;
  typedef typename R::Direction_3_base RDirection_3;
  typedef typename R::Line_3_base      RLine_3;
  typedef typename R::Segment_3_base   RSegment_3;
  typedef typename R::Ray_3_base       RRay_3;

  Line_3() : RLine_3() {}
  Line_3(const RLine_3& l) : RLine_3(l) {}
  Line_3(const RPoint_3 &p, const RPoint_3 &q) : RLine_3(p,q) {}
  Line_3(const RSegment_3& s) : RLine_3(s) {}
  Line_3(const RRay_3& r) : RLine_3(r) {}
  Line_3(const RPoint_3& p, const RDirection_3& d) : RLine_3(p,d) {}
};

} // namespace CGAL
#endif // CGAL_OLD_STYLE_KERNEL_LINE_3_H
