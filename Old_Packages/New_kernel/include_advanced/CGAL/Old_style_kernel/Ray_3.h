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
// file          : include/CGAL/Old_style_kernel/Ray_3.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ============================================================================
 

#ifndef CGAL_OLD_STYLE_KERNEL_RAY_3_H
#define CGAL_OLD_STYLE_KERNEL_RAY_3_H

namespace CGAL {

template <class R_>
class Ray_3< R_, Old_style_tag> : public R_::Ray_3_base
{
 public:
  typedef  R_                          R;
  typedef typename R::RT               RT;
  typedef typename R::FT               FT;
  typedef typename R::Point_3_base     RPoint_3;
  typedef typename R::Ray_3_base       RRay_3;
  typedef typename R::Direction_3_base RDirection_3;

  Ray_3() : RRay_3() {}
  Ray_3(const RRay_3& r) : RRay_3(r) {}
  Ray_3(const RPoint_3& sp,
        const RPoint_3& secondp)
    : RRay_3(sp, secondp) {}
  Ray_3(const RPoint_3& sp, const RDirection_3& d)
    : RRay_3(sp, d) {}

  Line_3< R_, Old_style_tag>
  supporting_line() const
  { return RRay_3::supporting_line(); }
};

} // namespace CGAL
#endif //CGAL_OLD_STYLE_KERNEL_RAY_3_H
