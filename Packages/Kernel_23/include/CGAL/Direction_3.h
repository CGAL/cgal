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
// source        : Direction_3.fw
// file          : Direction_3.h
// package       : _3 (3.9)
// maintainer    : Stefan Schirra <stschirr@mpi-sb.mpg.de>
// revision      : 3.9
// revision_date : 15 Oct 2000 
// author(s)     : Andreas Fabri
//                 Stefan Schirra
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ======================================================================
 

#ifndef CGAL_DIRECTION_3_H
#define CGAL_DIRECTION_3_H

#ifndef CGAL_REP_CLASS_DEFINED
#error  no representation class defined
#endif  // CGAL_REP_CLASS_DEFINED

#ifdef CGAL_HOMOGENEOUS_H
#ifndef CGAL_DIRECTIONH3_H
#include <CGAL/DirectionH3.h>
#endif // CGAL_DIRECTIONH3_H
#endif // CGAL_HOMOGENEOUS_H

#ifdef CGAL_CARTESIAN_H
#ifndef CGAL_DIRECTIONC3_H
#include <CGAL/Cartesian/Direction_3.h>
#endif // CGAL_DIRECTIONC3_H
#endif // CGAL_CARTESIAN_H

#ifdef CGAL_SIMPLE_CARTESIAN_H
#include <CGAL/SimpleCartesian/DirectionS3.h>
#endif // CGAL_SIMPLE_CARTESIAN_H


#ifndef CGAL_VECTOR_3_H
#include <CGAL/Vector_3.h>
#endif // CGAL_VECTOR_3_H

CGAL_BEGIN_NAMESPACE

template <class R_>
class Direction_3 : public R_::Direction_3_base
{
public:
  typedef          R_                       R;
  typedef typename R::RT                    RT;
  typedef typename R::FT                    FT;
  typedef typename R::Direction_3_base  RDirection_3;
  typedef typename R::Vector_3_base  RVector_3;

  Direction_3()
  {}
  Direction_3(const CGAL::Direction_3<R>& d)
    : RDirection_3( (const RDirection_3& )d )
  {}
  Direction_3(const RDirection_3&  d)
    : RDirection_3(d)
  {}
  Direction_3(const RVector_3&  v)
    : RDirection_3(v)
  {}
  Direction_3(const RT& hx, const RT& hy, const RT& hz)
    : RDirection_3(hx, hy, hz)
  {}

  bool operator==(const CGAL::Direction_3<R> & d) const
  { return RDirection_3::operator==(d); }

  bool operator!=(const CGAL::Direction_3<R> & d) const
  { return !(*this == d); }

  CGAL::Vector_3<R> vector() const
  { return (CGAL::Vector_3<R>)RDirection_3::to_vector(); }

  CGAL::Vector_3<R> to_vector() const
  { return (CGAL::Vector_3<R>)RDirection_3::to_vector(); }

  CGAL::Direction_3<R> transform(const CGAL::Aff_transformation_3<R> & t) const
  { return RDirection_3::transform(t); }

  CGAL::Direction_3<R> operator-() const
  { return RDirection_3::operator-(); }

  RT delta(int i) const
  { return RDirection_3::delta(i); }

  RT dx() const
  { return RDirection_3::dx(); }

  RT dy() const
  { return RDirection_3::dy(); }

  RT dz() const
  { return RDirection_3::dz(); }
};


#ifndef NO_OSTREAM_INSERT_DIRECTION_3
template < class R >
std::ostream& operator<<(std::ostream& os, const Direction_3<R>& d)
{
  typedef typename  R::Direction_3_base  RDirection_3;
  return os << (const RDirection_3& )d; }
#endif // NO_OSTREAM_INSERT_DIRECTION_3


#ifndef NO_ISTREAM_EXTRACT_DIRECTION_3
template < class R >
std::istream& operator>>(std::istream& is, Direction_3<R>& p)
{
  typedef typename  R::Direction_3_base  RDirection_3;
  return is >> (RDirection_3& )p; }
#endif // NO_ISTREAM_EXTRACT_DIRECTION_3

CGAL_END_NAMESPACE


#endif // CGAL_DIRECTION_3_H
