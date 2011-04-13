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
// file          : include/CGAL/Old_style_kernel/Circle_2.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ============================================================================
 

#ifndef CGAL_OLD_STYLE_KERNEL_CIRCLE_2_H
#define CGAL_OLD_STYLE_KERNEL_CIRCLE_2_H

namespace CGAL {

template <class R_>
class Circle_2< R_, Old_style_tag> : public R_::Circle_2_base
{
 public:
  typedef  R_                        R;
  typedef typename R::RT             RT;
  typedef typename R::FT             FT;
  typedef typename R::Point_2_base   RPoint_2;
  typedef typename R::Circle_2_base  RCircle_2;

    Circle_2() : RCircle_2() {}
    Circle_2(const RCircle_2& t) : RCircle_2(t) {}
    Circle_2(const RPoint_2& center,
             const FT &squared_radius,
             const Orientation &orientation)
      : RCircle_2(center, squared_radius, orientation) {}
    Circle_2(const RPoint_2 &center,
             const FT &squared_radius)
      : RCircle_2(center, squared_radius, COUNTERCLOCKWISE) {}
    Circle_2(const RPoint_2& p,
             const RPoint_2& q,
             const RPoint_2& r)
      : RCircle_2(p,q,r) {}
    Circle_2(const RPoint_2& p,
             const RPoint_2& q,
             const Orientation &orientation)
      : RCircle_2(p,q,orientation) {}
    Circle_2(const RPoint_2& p,
             const RPoint_2& q)
      : RCircle_2(p,q,COUNTERCLOCKWISE) {}
    Circle_2(const RPoint_2& center,
             const Orientation& orientation)
      : RCircle_2(center,FT(0),orientation) {}
    Circle_2(const RPoint_2& center)
      : RCircle_2(center,FT(0),COUNTERCLOCKWISE) {}
};

} // namespace CGAL
#endif // CGAL_OLD_STYLE_KERNEL_CIRCLE_2_H
