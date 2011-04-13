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
// file          : include/CGAL/Old_style_kernel/Direction_2.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ============================================================================
 

#ifndef CGAL_OLD_STYLE_KERNEL_DIRECTION_2_H
#define CGAL_OLD_STYLE_KERNEL_DIRECTION_2_H

namespace CGAL {

template <class R_>
class Direction_2< R_, Old_style_tag> : public R_::Direction_2_base
{
 public:
  typedef  R_                           R;
  typedef typename R::RT                RT;
  typedef typename R::Vector_2_base     RVector_2;
  typedef typename R::Direction_2_base  RDirection_2;

  Direction_2() {}
  Direction_2(const RDirection_2& d) : RDirection_2(d) {}
  Direction_2(const RVector_2& v) : RDirection_2(v) {}
  Direction_2(const RT &x, const RT &y) :  RDirection_2(x,y) {}
};

} // namespace CGAL
#endif // CGAL_OLD_STYLE_KERNEL_DIRECTION_2_H
