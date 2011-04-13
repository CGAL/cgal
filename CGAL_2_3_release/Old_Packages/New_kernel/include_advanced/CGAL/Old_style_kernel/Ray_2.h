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
// file          : include/CGAL/Old_style_kernel/Ray_2.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ============================================================================
 

#ifndef CGAL_OLD_STYLE_KERNEL_RAY_2_H
#define CGAL_OLD_STYLE_KERNEL_RAY_2_H

namespace CGAL {

template <class R_>
class Ray_2< R_, Old_style_tag> : public R_::Ray_2_base
{
 public:
  typedef  R_                          R;
  typedef typename R::RT               RT;
  typedef typename R::FT               FT;
  typedef typename R::Point_2_base     RPoint_2;
  typedef typename R::Ray_2_base       RRay_2;
  typedef typename R::Direction_2_base RDirection_2;

  Ray_2() : RRay_2() {}
  Ray_2(const RRay_2& r) : RRay_2(r) {}
  Ray_2(const RPoint_2& sp,
        const RPoint_2& secondp)
    : RRay_2(sp, secondp) {}
  Ray_2(const RPoint_2& sp, const RDirection_2& d)
    : RRay_2(sp, d) {}
};

} // namespace CGAL
#endif //CGAL_OLD_STYLE_KERNEL_RAY_2_H
