// ============================================================================
//
// Copyright (c) 1997-2001 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision$
// release_date  : $CGAL_Date$
//
// file          : include/CGAL/constructions/Polygon_offset_cons_ftC2.h
// package       : Straight_skeleton_2 (1.1.0)
//
// author(s)     : Fernando Cacciola
// maintainer    : Fernando Cacciola <fernando_cacciola@hotmail>
// coordinator   : Fernando Cacciola <fernando_cacciola@hotmail>
//
// ============================================================================
#ifndef CGAL_POLYGON_OFFSET_CONS_FTC2_H
#define CGAL_POLYGON_OFFSET_CONS_FTC2_H 1

#include <CGAL/constructions/Straight_skeleton_cons_ftC2.h>

CGAL_BEGIN_NAMESPACE


// Givenan offset distance 't' and 2 oriented lines l0:(a0,b0,c0), l1:(a1,b1,c1) returns the coordinates (x,y)
// of the intersection of their offsets at the given distance.
//
// PRECONDITIONS:
// The line coefficients must be normalized: a²+b²==1 and (a,b) being the leftward normal vector
// The offsets at the given distance do intersect in a single point.

template<class FT>
tuple<FT,FT> construct_offset_pointC2 ( FT t, tuple<FT,FT,FT> const& l0, tuple<FT,FT,FT> const& l1)
{
  FT a0,b0,c0,a1,b1,c1 ;

  tie(a0,b0,c0) = l0 ;
  tie(a1,b1,c1) = l1 ;

  FT den = a1 * b0 - a0 * b1 ;

  CGAL_assertion(! CGAL_NTS certified_is_zero(den) ) ;

  FT numX = t * b1 - t * b0 + b0 * c1 - b1 * c0 ;
  FT numY = t * a1 - t * a0 + a0 * c1 - a1 * c0 ;

  FT x = -numX / den ;
  FT y =  numY / den ;

  return make_tuple(x,y) ;
}

CGAL_END_NAMESPACE

#endif // CGAL_POLYGON_OFFSET_CONS_FTC2_H //
// EOF //

