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
// file          : include/CGAL/Old_style_kernel/Vector_2.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ============================================================================
 

#ifndef CGAL_OLD_STYLE_KERNEL_VECTOR_2_H
#define CGAL_OLD_STYLE_KERNEL_VECTOR_2_H

namespace CGAL {

template <class R_>
class Vector_2< R_, Old_style_tag> : public R_::Vector_2_base
{
 public:
  typedef  R_                        R;
  typedef typename R::RT             RT;
  typedef typename R::Point_2_base   RPoint_2;
  typedef typename R::Vector_2_base  RVector_2;

  Vector_2() {}
  Vector_2(const RVector_2& v) : RVector_2(v) {}
  Vector_2(const Null_vector &v) : RVector_2(v) {}
  Vector_2(const RT &x, const RT &y) : RVector_2(x,y) {}
  Vector_2(const RT &x, const RT &y, const RT &w) : RVector_2(x,y,w) {}

 private:
  Vector_2(const RPoint_2& p) : RVector_2(p) {}
};

} // namespace CGAL
#endif // CGAL_OLD_STYLE_KERNEL_VECTOR_2_H
