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
// source        : Aff_transformation_3.fw
// file          : Aff_transformation_3.h
// revision      : 2.4
// revision_date : 24 Aug 1999 
// author(s)     : Andreas Fabri
//                 Stefan Schirra
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ============================================================================
 

#ifndef CGAL_AFF_TRANSFORMATION_3_H
#define CGAL_AFF_TRANSFORMATION_3_H

#ifndef CGAL_REP_CLASS_DEFINED
#error  no representation class defined
#endif  // CGAL_REP_CLASS_DEFINED

#if defined(CGAL_CFG_INCOMPLETE_TYPE_BUG_1) && \
   !defined(CGAL_NO_PLANE_TRANSFORM_IN_AT)
#define CGAL_NO_PLANE_TRANSFORM_IN_AT
#endif // CGAL_CFG_INCOMPLETE_TYPE_BUG_1

#ifdef CGAL_HOMOGENEOUS_H
#include <CGAL/Aff_transformationH3.h>
#endif // CGAL_HOMOGENEOUS_H

#ifdef CGAL_CARTESIAN_H
#include <CGAL/Aff_transformationC3.h>
#endif // CGAL_CARTESIAN_H

#include <CGAL/Point_3.h>
#include <CGAL/Vector_3.h>
#include <CGAL/Direction_3.h>
#include <CGAL/Plane_3.h>

CGAL_BEGIN_NAMESPACE

template <class _R>
class Aff_transformation_3 : public _R::Aff_transformation_3_base
{
public:
  typedef _R                                R;
  typedef typename R::RT                    RT;
  typedef typename R::FT                    FT;
  typedef typename R::Plane_3_base  RPlane_3;
  typedef typename R::Aff_transformation_3_base  RAff_transformation_3;

  // default constructor
  Aff_transformation_3() : RAff_transformation_3()
  {}
  // copy constructor
  Aff_transformation_3(const Aff_transformation_3<R>& t)
    : RAff_transformation_3(t)
  {}
  // up cast constructor
  Aff_transformation_3(const RAff_transformation_3&  t)
    : RAff_transformation_3(t)
  {}
  // translation:
  Aff_transformation_3(const Translation tag,
                            const Vector_3<R>& v)
    : RAff_transformation_3(tag, v)
  {}
  // scaling:
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
  // dtor
  ~Aff_transformation_3()
  {}
  // transformations
  Point_3<R>       transform(const Point_3<R>& p) const
                        { return RAff_transformation_3::transform(p); }
  Point_3<R>       operator()(const Point_3<R>& p) const
                        { return RAff_transformation_3::transform(p); }
  Vector_3<R>      transform(const Vector_3<R>& v) const
                        { return RAff_transformation_3::transform(v); }
  Vector_3<R>      operator()(const Vector_3<R>& v) const
                        { return RAff_transformation_3::transform(v); }
  Direction_3<R>   transform(const Direction_3<R>& d) const
                        { return RAff_transformation_3::transform(d); }
  Direction_3<R>   operator()(const Direction_3<R>& d) const
                        { return RAff_transformation_3::transform(d); }
#ifndef CGAL_NO_PLANE_TRANSFORM_IN_AT
  Plane_3<R>       transform(const Plane_3<R>& pl) const
                        { return RAff_transformation_3::transform(pl); }
#else
  Plane_3<R>
    transform(const Plane_3<R>& pl) const
    {
      return
      (( const RPlane_3&  )pl).transform( (const RAff_transformation_3& )(*this) );
    }
#endif // CGAL_NO_PLANE_TRANSFORM_IN_AT
  Plane_3<R>       operator()(const Plane_3<R>& pl) const
                        { return transform(pl); }
  // further members
  Aff_transformation_3<R>
                        inverse() const
                        { return RAff_transformation_3::inverse(); }
  bool                  is_even() const
                        { return RAff_transformation_3::is_even(); }
  bool                  is_odd() const
                        { return !is_even(); }
  // access
  FT                    cartesian(int i, int j) const
                        { return RAff_transformation_3::cartesian(i,j); }
  RT                    homogeneous(int i, int j) const
                        { return RAff_transformation_3::homogeneous(i,j); }
  FT                    m(int i, int j) const
                        { return RAff_transformation_3::m(i,j); }
  RT                    hm(int i, int j) const
                        { return RAff_transformation_3::hm(i,j); }
  // composition
  Aff_transformation_3<R>
                        operator*(const Aff_transformation_3<R>& t) const
                        {
                          return
                          static_cast<const RAff_transformation_3&>(*this) *
                          static_cast<const RAff_transformation_3&>(t) ;
                        }
};

// I/O operators
#ifndef NO_OSTREAM_INSERT_AFF_TRANSFORMATION_3
template < class R >
std::ostream&
operator<<(std::ostream& os, const Aff_transformation_3<R>& t)
{
  typedef typename   R::Aff_transformation_3_base  RAff_transformation_3;
  return os << static_cast<const RAff_transformation_3&>(t);
}
#endif // NO_OSTREAM_INSERT_AFF_TRANSFORMATION_3

#ifndef NO_ISTREAM_EXTRACT_AFF_TRANSFORMATION_3
template < class R >
std::istream&
operator>>(std::istream& is, Aff_transformation_3<R>& t)
{
  typedef typename   R::Aff_transformation_3_base  RAff_transformation_3;
  return is >> static_cast<const RAff_transformation_3&>(t);
}
#endif // NO_ISTREAM_EXTRACT_AFF_TRANSFORMATION_3

CGAL_END_NAMESPACE


#endif // CGAL_AFF_TRANSFORMATION_3_H
