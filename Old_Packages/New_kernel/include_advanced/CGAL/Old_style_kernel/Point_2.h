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
// file          : include/CGAL/Old_style_kernel/Point_2.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ============================================================================
 

#ifndef CGAL_OLD_STYLE_KERNEL_POINT_2_H
#define CGAL_OLD_STYLE_KERNEL_POINT_2_H

namespace CGAL {

template <class R_>
class Point_2< R_, Old_style_tag> : public R_::Point_2_base
{
 public:
  typedef  R_                        R;
  typedef typename R::RT             RT;
  typedef typename R::Point_2_base   RPoint_2;
  typedef typename R::Vector_2_base  RVector_2;

  Point_2() {}
  Point_2(const Origin& o) : RPoint_2(o) {}
  Point_2(const RPoint_2& p) : RPoint_2(p) {}
  Point_2(const RT& hx, const RT& hy) : RPoint_2(hx, hy) {}
  Point_2(const RT& hx, const RT& hy, const RT& hw) : RPoint_2(hx, hy, hw) {}

 private:
  Point_2(const RVector_2& v) : RPoint_2(v) {}
};

} // namespace CGAL
#endif // CGAL_OLD_STYLE_KERNEL_POINT_2_H
