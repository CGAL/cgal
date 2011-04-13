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
// file          : include/CGAL/Old_style_kernel/Direction_3.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ============================================================================
 

#ifndef CGAL_OLD_STYLE_KERNEL_DIRECTION_3_H
#define CGAL_OLD_STYLE_KERNEL_DIRECTION_3_H

namespace CGAL {

template <class R_>
class Direction_3< R_, Old_style_tag> : public R_::Direction_3_base
{
 public:
  typedef  R_                           R;
  typedef typename R::RT                RT;
  typedef typename R::Vector_3_base     RVector_3;
  typedef typename R::Direction_3_base  RDirection_3;

  Direction_3() {}
  Direction_3(const RDirection_3& d) : RDirection_3(d) {}
  Direction_3(const RVector_3& v) : RDirection_3(v) {}
  Direction_3(const RT& x, const RT& y, const RT& z)
    : RDirection_3(x,y,z) {}
};

} // namespace CGAL
#endif // CGAL_OLD_STYLE_KERNEL_DIRECTION_3_H
