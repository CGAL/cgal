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
// file          : include/CGAL/Old_style_kernel/Aff_transformation_3.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ============================================================================
 

#ifndef CGAL_OLD_STYLE_KERNEL_AFF_TRANSFORMATION_3_H
#define CGAL_OLD_STYLE_KERNEL_AFF_TRANSFORMATION_3_H

namespace CGAL {

template <class R_>
class Aff_transformation_3< R_, Old_style_tag>
  : public R_::Aff_transformation_3_base
{
 public:
  typedef  R_                           R;
  typedef typename R::RT                RT;
  typedef typename R::Vector_3_base     RVector_3;
  typedef typename R::Direction_3_base  RDirection_3;
  typedef typename R::Aff_transformation_3_base  RAff_transformation_3;

  Aff_transformation_3()
    : RAff_transformation_3() {}
  Aff_transformation_3(const RAff_transformation_3& t)
    : RAff_transformation_3(t) {}
  Aff_transformation_3(const Identity_transformation tag)
    : RAff_transformation_3(tag) {}
  Aff_transformation_3(const Translation tag,
                       const RVector_3& v)
    : RAff_transformation_3(tag, v) {}
  Aff_transformation_3(const Scaling tag,
                       const RT &s,
                       const RT &w= RT(1))
    : RAff_transformation_3(tag, s, w) {}
  Aff_transformation_3(
      const RT& m11, const RT& m12, const RT& m13, const RT& m14,
      const RT& m21, const RT& m22, const RT& m23, const RT& m24,
      const RT& m31, const RT& m32, const RT& m33, const RT& m34,
                                                   const RT& w= RT(1) )
    : RAff_transformation_3(m11, m12, m13, m14,
                            m21, m22, m23, m24,
                            m31, m32, m33, m34,
                                           w) {}
  Aff_transformation_3(
      const RT& m11, const RT& m12, const RT& m13,
      const RT& m21, const RT& m22, const RT& m23,
      const RT& m31, const RT& m32, const RT& m33,
                                                   const RT& w = RT(1) )
    : RAff_transformation_3(m11, m12, m13,
                           m21, m22, m23,
                           m31, m32, m33,
                                          w) {}

};

} // namespace CGAL
#endif // CGAL_OLD_STYLE_KERNEL_AFF_TRANSFORMATION_3_H
