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
// file          : include/CGAL/Old_style_kernel/Line_2.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ============================================================================
 

#ifndef CGAL_OLD_STYLE_KERNEL_LINE_2_H
#define CGAL_OLD_STYLE_KERNEL_LINE_2_H

namespace CGAL {

template <class R_>
class Line_2< R_, Old_style_tag> : public R_::Line_2_base
{
 public:
  typedef  R_                          R;
  typedef typename R::RT               RT;
  typedef typename R::FT               FT;
  typedef typename R::Point_2_base     RPoint_2;
  typedef typename R::Direction_2_base RDirection_2;
  typedef typename R::Line_2_base      RLine_2;
  typedef typename R::Segment_2_base   RSegment_2;
  typedef typename R::Ray_2_base       RRay_2;

  Line_2() : RLine_2() {}
  Line_2(const RLine_2& l) : RLine_2(l) {}
  Line_2(const RPoint_2 &p, const RPoint_2 &q) : RLine_2(p,q) {}
  Line_2(const RT& a, const RT& b, const RT& c) : RLine_2(a,b,c) {}
  Line_2(const RSegment_2& s) : RLine_2(s) {}
  Line_2(const RRay_2& r) : RLine_2(r) {}
  Line_2(const RPoint_2& p, const RDirection_2& d) : RLine_2(p,d) {}
};

} // namespace CGAL
#endif // CGAL_OLD_STYLE_KERNEL_LINE_2_H
