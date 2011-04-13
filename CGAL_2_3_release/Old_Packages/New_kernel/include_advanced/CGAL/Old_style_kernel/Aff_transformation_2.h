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
// file          : include/CGAL/Old_style_kernel/Aff_transformation_2.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ============================================================================
 

#ifndef CGAL_OLD_STYLE_KERNEL_AFF_TRANSFORMATION_2_H
#define CGAL_OLD_STYLE_KERNEL_AFF_TRANSFORMATION_2_H

namespace CGAL {

template <class R_>
class Aff_transformation_2< R_, Old_style_tag>
  : public R_::Aff_transformation_2_base
{
 public:
  typedef  R_                           R;
  typedef typename R::RT                RT;
  typedef typename R::Vector_2_base     RVector_2;
  typedef typename R::Direction_2_base  RDirection_2;
  typedef typename R::Aff_transformation_2_base  RAff_transformation_2;

  Aff_transformation_2()
    : RAff_transformation_2() {}
  Aff_transformation_2(const RAff_transformation_2& t)
    : RAff_transformation_2(t) {}
  Aff_transformation_2(const Identity_transformation tag)
    : RAff_transformation_2(tag) {}
  Aff_transformation_2(const Translation tag,
                       const RVector_2& v)
    : RAff_transformation_2(tag, v) {}
  Aff_transformation_2(const Rotation tag,
                       const RDirection_2& d,
                       const RT &num,
                       const RT &den = RT(1))
    : RAff_transformation_2(tag, d, num, den) {}
  Aff_transformation_2(const Rotation tag,
                       const RT &sin,
                       const RT &cos,
                       const RT &den = RT(1))
    : RAff_transformation_2(tag, sin, cos, den) {}
  Aff_transformation_2(const Scaling tag,
                       const RT &s,
                       const RT &w= RT(1))
    : RAff_transformation_2(tag, s, w) {}
  Aff_transformation_2(const RT & m11, const RT & m12, const RT & m13,
                       const RT & m21, const RT & m22, const RT & m23,
                                                       const RT &w= RT(1))
    : RAff_transformation_2(m11, m12, m13,
                            m21, m22, m23,
                                      w) {}
  Aff_transformation_2(const RT & m11, const RT & m12,
                       const RT & m21, const RT & m22,
                                       const RT &w = RT(1))
    : RAff_transformation_2(m11, m12,
                            m21, m22,
                                      w) {}
};

} // namespace CGAL
#endif // CGAL_OLD_STYLE_KERNEL_AFF_TRANSFORMATION_2_H
