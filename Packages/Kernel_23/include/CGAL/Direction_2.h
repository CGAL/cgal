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
// file          : Direction_2.h
// package       : _2
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ======================================================================
 
#ifndef CGAL_DIRECTION_2_H
#define CGAL_DIRECTION_2_H

#include <CGAL/Vector_2.h>

CGAL_BEGIN_NAMESPACE

template <class R_>
class Direction_2 : public R_::Direction_2_base
{
public:
  typedef  R_   R;
  typedef typename R::RT                    RT;
  typedef typename R::FT                    FT;
  typedef typename R::Vector_2_base  RVector_2;
  typedef typename R::Direction_2_base  RDirection_2;

  Direction_2()
  {}

  Direction_2(const CGAL::Direction_2<R> &d)
    : RDirection_2(static_cast<const RDirection_2&>(d))
  {}

  Direction_2(const RDirection_2& d)
    : RDirection_2(d)
  {}

  Direction_2(const RVector_2& v)
    : RDirection_2(v)
  {}

  Direction_2(const RT &x, const RT &y)
    :  RDirection_2(x,y)
  {}

  CGAL::Vector_2<R>
  vector() const
  { return static_cast<CGAL::Vector_2<R> >(RDirection_2::to_vector()); }

  CGAL::Vector_2<R>
  to_vector() const
  { return static_cast<CGAL::Vector_2<R> >(RDirection_2::to_vector()); }

  CGAL::Direction_2<R>
  transform(const CGAL::Aff_transformation_2<R> &t) const
  { return RDirection_2::transform(t); }

  CGAL::Direction_2<R>
  operator-() const
  { return RDirection_2::operator-(); }
};

CGAL_END_NAMESPACE

#endif // CGAL_DIRECTION_2_H
