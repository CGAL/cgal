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
// release_date  : 2000, October 15
// 
// source        : Point_3.fw
// file          : Point_3.h
// package       : _3 (3.9)
// maintainer    : Stefan Schirra <stschirr@mpi-sb.mpg.de>
// revision      : 3.9
// revision_date : 15 Oct 2000 
// author(s)     : Andreas Fabri
//                 Stefan Schirra
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ======================================================================
 

#ifndef CGAL_POINT_3_H
#define CGAL_POINT_3_H

#ifndef CGAL_REP_CLASS_DEFINED
#error  no representation class defined
#endif  // CGAL_REP_CLASS_DEFINED

#ifdef CGAL_HOMOGENEOUS_H
#ifndef CGAL_POINTH3_H
#include <CGAL/PointH3.h>
#endif // CGAL_POINTH3_H
#endif // CGAL_HOMOGENEOUS_H

#ifdef CGAL_CARTESIAN_H
#ifndef CGAL_POINTC3_H
#include <CGAL/Cartesian/Point_3.h>
#endif // CGAL_POINTC3_H
#endif // CGAL_CARTESIAN_H

#ifdef CGAL_SIMPLE_CARTESIAN_H
#include <CGAL/SimpleCartesian/PointS3.h>
#endif // CGAL_SIMPLE_CARTESIAN_H


#include <CGAL/point_vector_declarations_3.h>

CGAL_BEGIN_NAMESPACE

template <class R_>
class Point_3 : public R_::Point_3_base
{
public:
  typedef          R_                       R;
  typedef typename R::RT                    RT;
  typedef typename R::FT                    FT;
  typedef typename R::Point_3_base  RPoint_3;
  typedef typename R::Vector_3_base  RVector_3;

friend  CGAL_FRIEND_INLINE
        CGAL::Point_3<R>
        vector_to_point_conversion CGAL_NULL_TMPL_ARGS
                                            (const CGAL::Vector_3<R>& v);

public:
  Point_3()
  {}
  Point_3(const Origin& o) : RPoint_3(o)
  {}
  Point_3(const CGAL::Point_3<R>& p) : RPoint_3( (const RPoint_3& )p )
  {}
  Point_3(const RPoint_3&  p) : RPoint_3(p)
  {}
  Point_3(const RT& hx, const RT& hy, const RT& hz)
    : RPoint_3(hx, hy, hz)
  {}
  Point_3(const RT& hx, const RT& hy, const RT& hz, const RT& hw)
    : RPoint_3(hx, hy, hz, hw)
  {}

  bool operator==(const CGAL::Point_3<R>& p) const
  { return RPoint_3::operator==(p); }

  bool operator!=(const CGAL::Point_3<R>& p) const
  { return !(*this == p); }


  RT hx() const
  { return RPoint_3::hx(); }

  RT hy() const
  { return RPoint_3::hy(); }

  RT hz() const
  { return RPoint_3::hz(); }

  RT hw() const
  { return RPoint_3::hw(); }

  FT x() const
  { return RPoint_3::x(); }

  FT y() const
  { return RPoint_3::y(); }

  FT z() const
  { return RPoint_3::z(); }

  RT homogeneous(int i) const
  { return RPoint_3::homogeneous(i); }

  FT cartesian(int i) const
  { return RPoint_3::cartesian(i); }

  FT operator[](int i) const
  { return cartesian(i); }

  int dimension() const
  { return 3; }

  Bbox_3       bbox() const
  { return RPoint_3::bbox(); }

  CGAL::Point_3<R> transform(const CGAL::Aff_transformation_3<R>& t) const
  { return RPoint_3::transform(t); }

private:
  Point_3(const RVector_3&  v) : RPoint_3(v)
  {}
};

CGAL_END_NAMESPACE


#ifndef CGAL_VECTOR_3_H
#include <CGAL/Vector_3.h>
#endif // CGAL_VECTOR_3_H

#include <CGAL/point_vector_definitions_3.C>

#ifndef CGAL_AFF_TRANSFORMATION_3_H
#include <CGAL/Aff_transformation_3.h>
#endif // CGAL_AFF_TRANSFORMATION_3_H

CGAL_BEGIN_NAMESPACE

template <class R>
inline
bool
operator==(const Origin& o, const Point_3<R>& p)
{ return p == o; }

template <class R>
inline
bool
operator!=(const Origin& o, const Point_3<R>& p)
{ return p != o; }


#ifndef NO_OSTREAM_INSERT_POINT_3

template < class R >
std::ostream&
operator<<(std::ostream& os, const Point_3<R>& p)
{
  typedef typename  R::Point_3_base  RPoint_3;
  return os << (const RPoint_3& )p;
}
#endif // NO_OSTREAM_INSERT_POINT_3

#ifndef NO_ISTREAM_EXTRACT_POINT_3
template < class R >
std::istream& operator>>(std::istream& is, Point_3<R>& p)
{
  typedef typename  R::Point_3_base  RPoint_3;
  return is >> (RPoint_3& )p;
}
#endif // NO_ISTREAM_EXTRACT_POINT_3

CGAL_END_NAMESPACE


#endif // CGAL_POINT_3_H
