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
// file          : include/CGAL/Old_style_kernel/Plane_3.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ============================================================================
 

#ifndef CGAL_OLD_STYLE_KERNEL_PLANE_3_H
#define CGAL_OLD_STYLE_KERNEL_PLANE_3_H

namespace CGAL {

template <class R_>
class Plane_3< R_, Old_style_tag> : public R_::Plane_3_base
{
 public:
  typedef  R_                          R;
  typedef typename R::RT               RT;
  typedef typename R::FT               FT;
  typedef typename R::Point_3_base     RPoint_3;
  typedef typename R::Vector_3_base    RVector_3;
  typedef typename R::Direction_3_base RDirection_3;
  typedef typename R::Line_3_base      RLine_3;
  typedef typename R::Plane_3_base     RPlane_3;
  typedef typename R::Segment_3_base   RSegment_3;
  typedef typename R::Ray_3_base       RRay_3;

  Plane_3() : RPlane_3() {}
  Plane_3(const RPlane_3&  p) : RPlane_3(p) {}
  Plane_3(const RPoint_3& p, const RPoint_3& q, const RPoint_3& r)
    : RPlane_3(p,q,r) {}
  Plane_3(const RPoint_3& p, const RDirection_3& d) : RPlane_3(p,d) {}
  Plane_3(const RPoint_3& p, const RVector_3& v) : RPlane_3(p,v) {}
  Plane_3(const RT& a, const RT& b, const RT& c, const RT& d)
    : RPlane_3(a,b,c,d) {}
  Plane_3(const RLine_3& l, const RPoint_3& p) : RPlane_3(l,p) {}
  Plane_3(const RSegment_3& s, const RPoint_3& p) : RPlane_3(s,p) {}
  Plane_3(RRay_3& r, const RPoint_3& p) : RPlane_3(r,p) {}
};

} // namespace CGAL
#endif // CGAL_OLD_STYLE_KERNEL_PLANE_3_H
