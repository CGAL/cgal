// ======================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
// 
// release       : 
// release_date  : 
// 
// file          : Aff_transformation_3.h
// package       : _3
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas Fabri, Stefan Schirra
//
// coordinator   : MPI, Saarbruecken
// ======================================================================

#ifndef CGAL_AFF_TRANSFORMATION_3_H
#define CGAL_AFF_TRANSFORMATION_3_H

CGAL_BEGIN_NAMESPACE

template <class R_>
class Aff_transformation_3 : public R_::Kernel_base::Aff_transformation_3
{
  typedef typename R_::RT            RT;
  typedef typename R_::Vector_3      Vector_3;
  typedef typename R_::Kernel_base::Aff_transformation_3 RAff_transformation_3;
public:
  typedef R_                        R;

  Aff_transformation_3()
      : RAff_transformation_3() {}

  Aff_transformation_3(const CGAL::Aff_transformation_3<R>& t)
    : RAff_transformation_3(t) {}

  Aff_transformation_3(const RAff_transformation_3&  t)
    : RAff_transformation_3(t) {}

  Aff_transformation_3(const Identity_transformation& tag)
    : RAff_transformation_3(tag) {}

  Aff_transformation_3(const Translation tag,
                       const Vector_3& v)
    : RAff_transformation_3(tag, v) {}

  Aff_transformation_3(const Scaling tag,
                       const RT& s,
                       const RT& w= RT(1) )
    : RAff_transformation_3(tag, s, w) {}

  // the general case:
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

#ifndef CGAL_NO_OSTREAM_INSERT_AFF_TRANSFORMATION_3
template < class R >
std::ostream&
operator<<(std::ostream& os, const CGAL::Aff_transformation_3<R>& t)
{
  typedef typename R::Kernel_base::Aff_transformation_3 RAff_transformation_3;
  return os << static_cast<const RAff_transformation_3&>(t);
}
#endif // CGAL_NO_OSTREAM_INSERT_AFF_TRANSFORMATION_3

#ifndef CGAL_NO_ISTREAM_EXTRACT_AFF_TRANSFORMATION_3
template < class R >
std::istream&
operator>>(std::istream& is, CGAL::Aff_transformation_3<R>& t)
{
  typedef typename R::Kernel_base::Aff_transformation_3 RAff_transformation_3;
  return is >> static_cast<const RAff_transformation_3&>(t);
}
#endif // CGAL_NO_ISTREAM_EXTRACT_AFF_TRANSFORMATION_3

CGAL_END_NAMESPACE

#endif // CGAL_AFF_TRANSFORMATION_3_H
