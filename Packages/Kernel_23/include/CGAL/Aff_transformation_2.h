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
// file          : Aff_transformation_2.h
// package       : _2
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas Fabri
//                 Stefan Schirra
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ======================================================================
 

#ifndef CGAL_AFF_TRANSFORMATION_2_H
#define CGAL_AFF_TRANSFORMATION_2_H

#ifndef CGAL_REP_CLASS_DEFINED
#error  no representation class defined
#endif  // CGAL_REP_CLASS_DEFINED

#if defined(CGAL_CFG_INCOMPLETE_TYPE_BUG_1) && \
   !defined(CGAL_NO_LINE_TRANSFORM_IN_AT)
#define CGAL_NO_LINE_TRANSFORM_IN_AT
#endif // CGAL_CFG_INCOMPLETE_TYPE_BUG_1

#include <CGAL/Line_2.h>

#if defined CGAL_HOMOGENEOUS_H || defined CGAL_SIMPLE_HOMOGENEOUS_H
#include <CGAL/Aff_transformationH2.h>
#endif

#if defined CGAL_CARTESIAN_H || defined CGAL_SIMPLE_CARTESIAN_H
#include <CGAL/Cartesian/Aff_transformation_2.h>
#endif


CGAL_BEGIN_NAMESPACE

template <class R_>
class Aff_transformation_2 : public R_::Aff_transformation_2_base
{
public:
  typedef  R_                               R;
  typedef typename R::RT                    RT;
  typedef typename R::FT                    FT;
  typedef typename R::Line_2_base  RLine_2;
  typedef typename R::Direction_2_base  RDirection_2;
  typedef typename R::Aff_transformation_2_base  RAff_transformation_2;

  // default constructor
  Aff_transformation_2()
    : RAff_transformation_2()
  {}

  // copy constructor
  Aff_transformation_2(const CGAL::Aff_transformation_2<R> &t)
    : RAff_transformation_2(t)
  {}

  // up cast constructor
  Aff_transformation_2(const RAff_transformation_2& t)
    : RAff_transformation_2(t)
  {}

  // Identity:
  Aff_transformation_2(const Identity_transformation tag)
    : RAff_transformation_2(tag)
  {}

  // Translation:
  Aff_transformation_2(const Translation tag,
                       const CGAL::Vector_2<R> &v)
    : RAff_transformation_2(tag, v)
  {}

  // Rational Rotation:
  Aff_transformation_2(const Rotation tag,
                       const CGAL::Direction_2<R> &d,
                       const RT &num,
                       const RT &den = RT(1))
    : RAff_transformation_2(tag, RDirection_2(d), num, den)
  {}

  Aff_transformation_2(const Rotation tag,
                       const RT &sin,
                       const RT &cos,
                       const RT &den = RT(1))
    : RAff_transformation_2(tag, sin, cos, den)
  {}

  // Reflection:
  Aff_transformation_2(const Reflection tag,
                       const CGAL::Line_2<R>& l )
    : RAff_transformation_2(tag, l)
  {}


  // Scaling:
  Aff_transformation_2(const Scaling tag,
                       const RT &s,
                       const RT &w= RT(1))
    : RAff_transformation_2(tag, s, w)
  {}

  // The general case:
  Aff_transformation_2(const RT & m11,
                       const RT & m12,
                       const RT & m13,

                       const RT & m21,
                       const RT & m22,
                       const RT & m23,

                       const RT &w= RT(1))
    : RAff_transformation_2(m11, m12, m13,
                            m21, m22, m23,
                                      w)
  {}

  Aff_transformation_2(const RT & m11, const RT & m12,
                       const RT & m21, const RT & m22,
                       const RT &w = RT(1))
    : RAff_transformation_2(m11, m12,
                            m21, m22,
                                      w)
  {}

  ~Aff_transformation_2()
  {}

  CGAL::Point_2<R>     transform(const CGAL::Point_2<R> &p) const
                      { return RAff_transformation_2::transform(p); }

  CGAL::Point_2<R>     operator()(const CGAL::Point_2<R> &p) const
                      { return transform(p); }

  CGAL::Vector_2<R>    transform(const CGAL::Vector_2<R> &v) const
                      { return RAff_transformation_2::transform(v); }

  CGAL::Vector_2<R>    operator()(const CGAL::Vector_2<R> &v) const
                      { return transform(v); }

  CGAL::Direction_2<R> transform(const CGAL::Direction_2<R> &d) const
                      { return RAff_transformation_2::transform(d); }

  CGAL::Direction_2<R> operator()(const CGAL::Direction_2<R> &d) const
                      { return transform(d); }

  CGAL::Line_2<R>      transform(const CGAL::Line_2<R> &l) const
  {
#ifndef CGAL_NO_LINE_TRANSFORM_IN_AT
                        return RAff_transformation_2::transform(l);
#else
                        return
      ((const RLine_2&)l).transform((const RAff_transformation_2&)(*this));
#endif // CGAL_NO_LINE_TRANSFORM_IN_AT
                      }

  CGAL::Line_2<R>      operator()(const CGAL::Line_2<R> &l) const
                      { return transform(l); }


  CGAL::Aff_transformation_2<R>
                      inverse() const
                      { return RAff_transformation_2::inverse(); }

  bool                is_even() const
                      { return RAff_transformation_2::is_even(); }

  bool                is_odd() const
                      { return !is_even(); }

  CGAL::Aff_transformation_2<R>
                      operator*(const CGAL::Aff_transformation_2<R> &t) const
                      { return RAff_transformation_2::operator*(t); }

  FT          cartesian(int i, int j) const
              { return RAff_transformation_2::cartesian(i,j); }
  RT          homogeneous(int i, int j) const
              { return RAff_transformation_2::homogeneous(i,j); }
  FT          m(int i, int j) const
              { return RAff_transformation_2::m(i,j); }
  RT          hm(int i, int j) const
              { return RAff_transformation_2::hm(i,j); }
  
};

#ifndef CGAL_NO_OSTREAM_INSERT_AFF_TRANSFORMATION_2
template < class R >
std::ostream &
operator<<(std::ostream &os, const CGAL::Aff_transformation_2<R> &t)
{
  typedef typename  R::Aff_transformation_2_base  RAff_transformation_2;
  return os << static_cast<const RAff_transformation_2&>(t);
}
#endif // CGAL_NO_OSTREAM_INSERT_AFF_TRANSFORMATION_2

#ifndef CGAL_NO_ISTREAM_EXTRACT_AFF_TRANSFORMATION_2
template < class R >
std::istream &
operator>>(std::istream &is, CGAL::Aff_transformation_2<R> &t)
{
  typedef typename  R::Aff_transformation_2_base  RAff_transformation_2;
  return is >> static_cast<RAff_transformation_2&>(t);
}
#endif // CGAL_NO_ISTREAM_EXTRACT_AFF_TRANSFORMATION_2

CGAL_END_NAMESPACE

#endif // CGAL_AFF_TRANSFORMATION_2_H
