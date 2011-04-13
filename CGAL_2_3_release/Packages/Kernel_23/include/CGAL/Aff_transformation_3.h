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
// author(s)     : Andreas Fabri
//                 Stefan Schirra
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ======================================================================

#ifndef CGAL_AFF_TRANSFORMATION_3_H
#define CGAL_AFF_TRANSFORMATION_3_H

#include <CGAL/Point_3.h>
#include <CGAL/Vector_3.h>
#include <CGAL/Direction_3.h>
#include <CGAL/Plane_3.h>

CGAL_BEGIN_NAMESPACE

template <class R_>
class Aff_transformation_3 : public R_::Aff_transformation_3_base
{
public:
  typedef R_                                R;
  typedef typename R::RT                    RT;
  typedef typename R::FT                    FT;
  typedef typename R::Plane_3_base  RPlane_3;
  typedef typename R::Aff_transformation_3_base  RAff_transformation_3;

  Aff_transformation_3() : RAff_transformation_3()
  {}

  Aff_transformation_3(const CGAL::Aff_transformation_3<R>& t)
    : RAff_transformation_3(t)
  {}

  Aff_transformation_3(const RAff_transformation_3&  t)
    : RAff_transformation_3(t)
  {}

  Aff_transformation_3(const Identity_transformation& tag)
    : RAff_transformation_3(tag)
  {}

  Aff_transformation_3(const Translation tag,
                       const CGAL::Vector_3<R>& v)
    : RAff_transformation_3(tag, v)
  {}

  Aff_transformation_3(const Scaling tag,
                       const RT& s,
                       const RT& w= RT(1) )
    : RAff_transformation_3(tag, s, w)
  {}

  // the general case:
  Aff_transformation_3(
      const RT& m11, const RT& m12, const RT& m13, const RT& m14,
      const RT& m21, const RT& m22, const RT& m23, const RT& m24,
      const RT& m31, const RT& m32, const RT& m33, const RT& m34,
                                                   const RT& w= RT(1) )
    : RAff_transformation_3(m11, m12, m13, m14,
                            m21, m22, m23, m24,
                            m31, m32, m33, m34,
                                           w)
  {}

  Aff_transformation_3(
      const RT& m11, const RT& m12, const RT& m13,
      const RT& m21, const RT& m22, const RT& m23,
      const RT& m31, const RT& m32, const RT& m33,
                                                   const RT& w = RT(1) )
    : RAff_transformation_3(m11, m12, m13,
                           m21, m22, m23,
                           m31, m32, m33,
                                          w)
  {}

  // transformations
  CGAL::Point_3<R>       transform(const CGAL::Point_3<R>& p) const
                        { return RAff_transformation_3::transform(p); }
  CGAL::Point_3<R>       operator()(const CGAL::Point_3<R>& p) const
                        { return RAff_transformation_3::transform(p); }
  CGAL::Vector_3<R>      transform(const CGAL::Vector_3<R>& v) const
                        { return RAff_transformation_3::transform(v); }
  CGAL::Vector_3<R>      operator()(const CGAL::Vector_3<R>& v) const
                        { return RAff_transformation_3::transform(v); }
  CGAL::Direction_3<R>   transform(const CGAL::Direction_3<R>& d) const
                        { return RAff_transformation_3::transform(d); }
  CGAL::Direction_3<R>   operator()(const CGAL::Direction_3<R>& d) const
                        { return RAff_transformation_3::transform(d); }
  CGAL::Plane_3<R>       transform(const CGAL::Plane_3<R>& pl) const
                        { return RAff_transformation_3::transform(pl); }
  CGAL::Plane_3<R>       operator()(const CGAL::Plane_3<R>& pl) const
                        { return transform(pl); }
  // further members
  CGAL::Aff_transformation_3<R>
                        inverse() const
                        { return RAff_transformation_3::inverse(); }
  // composition
  CGAL::Aff_transformation_3<R>
                        operator*(const CGAL::Aff_transformation_3<R>& t) const
                        {
                          return
                          static_cast<const RAff_transformation_3&>(*this) *
                          static_cast<const RAff_transformation_3&>(t) ;
                        }
};

#ifndef CGAL_NO_OSTREAM_INSERT_AFF_TRANSFORMATION_3
template < class R >
std::ostream&
operator<<(std::ostream& os, const CGAL::Aff_transformation_3<R>& t)
{
  typedef typename   R::Aff_transformation_3_base  RAff_transformation_3;
  return os << static_cast<const RAff_transformation_3&>(t);
}
#endif // CGAL_NO_OSTREAM_INSERT_AFF_TRANSFORMATION_3

#ifndef CGAL_NO_ISTREAM_EXTRACT_AFF_TRANSFORMATION_3
template < class R >
std::istream&
operator>>(std::istream& is, CGAL::Aff_transformation_3<R>& t)
{
  typedef typename   R::Aff_transformation_3_base  RAff_transformation_3;
  return is >> static_cast<const RAff_transformation_3&>(t);
}
#endif // CGAL_NO_ISTREAM_EXTRACT_AFF_TRANSFORMATION_3

CGAL_END_NAMESPACE

#endif // CGAL_AFF_TRANSFORMATION_3_H
