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
// file          : include/CGAL/Old_style_kernel/Point_3.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ============================================================================
 

#ifndef CGAL_OLD_STYLE_KERNEL_POINT_3_H
#define CGAL_OLD_STYLE_KERNEL_POINT_3_H

namespace CGAL {

template <class R_>
class Point_3< R_, Old_style_tag> : public R_::Point_3_base
{
 public:
  typedef  R_                        R;
  typedef typename R::RT             RT;
  typedef typename R::Point_3_base   RPoint_3;
  typedef typename R::Vector_3_base  RVector_3;

  Point_3() {}
  Point_3(const Origin& o) : RPoint_3(o) {}
  Point_3(const RPoint_3& p) : RPoint_3(p) {}
  Point_3(const RT& hx, const RT& hy, const RT& hz)
    : RPoint_3(hx, hy, hz) {}
  Point_3(const RT& hx, const RT& hy, const RT& hz, const RT& hw )
    : RPoint_3(hx, hy, hz, hw) {}

 private:
  Point_3(const RVector_3& v) : RPoint_3(v) {}
};

} // namespace CGAL
#endif // CGAL_OLD_STYLE_KERNEL_POINT_3_H
