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
// file          : Direction_3.h
// package       : _3
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas Fabri
//                 Stefan Schirra
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ======================================================================
 

#ifndef CGAL_DIRECTION_3_H
#define CGAL_DIRECTION_3_H

#ifndef CGAL_REP_CLASS_DEFINED
#error  no representation class defined
#endif

#ifdef CGAL_HOMOGENEOUS_H
#include <CGAL/DirectionH3.h>
#endif

#if defined CGAL_CARTESIAN_H || defined CGAL_SIMPLE_CARTESIAN_H
#include <CGAL/Cartesian/Direction_3.h>
#endif

#include <CGAL/Vector_3.h>

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


#ifndef CGAL_NO_OSTREAM_INSERT_DIRECTION_3
template < class R >
std::ostream& operator<<(std::ostream& os, const Direction_3<R>& d)
{
  typedef typename  R::Direction_3_base  RDirection_3;
  return os << (const RDirection_3& )d; }
#endif // CGAL_NO_OSTREAM_INSERT_DIRECTION_3


#ifndef CGAL_NO_ISTREAM_EXTRACT_DIRECTION_3
template < class R >
std::istream& operator>>(std::istream& is, Direction_3<R>& p)
{
  typedef typename  R::Direction_3_base  RDirection_3;
  return is >> (RDirection_3& )p; }
#endif // CGAL_NO_ISTREAM_EXTRACT_DIRECTION_3

CGAL_END_NAMESPACE

#endif // CGAL_DIRECTION_3_H
