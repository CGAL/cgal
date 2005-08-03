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
// file          : include/CGAL/constructions/Straight_skeleton_ftC2.h
// package       : Straight_skeleton_2 (1.1.0)
//
// author(s)     : Fernando Cacciola
// maintainer    : Fernando Cacciola <fernando_cacciola@hotmail>
// coordinator   : Fernando Cacciola <fernando_cacciola@hotmail>
//
// ============================================================================
#ifndef CGAL_STRAIGHT_SKELETON_CONSTRUCTIONS_FTC2_H
#define CGAL_STRAIGHT_SKELETON_CONSTRUCTIONS_FTC2_H 1

#include <CGAL/Quotient.h>

CGAL_BEGIN_NAMESPACE

// Given 3 oriented lines (a0,b0,c0),(a1,b1,c1) and (a2,b2,c2),
// returns the POSITIVE OFFSET DISTANCE at which the leftward-offsetted lines
// intersect in a single point
//
// PRECONDITIONS: 
// The line coefficients must be normalized: a²+b²==1 and (a,b) being the inward normal vector 
// There exist such a distance for which the 3 leftward-offsetted lines
// intersect in a single point.
//
template<class FT>
Quotient<FT> compute_offset_lines_isec_time ( FT const& a0, FT const& b0, FT const& c0
                                            , FT const& a1, FT const& b1, FT const& c1
                                            , FT const& a2, FT const& b2, FT const& c2
                                            )
{
  // A positive offset line is given by: a*x(t) + b*y(t) + c + t = 0
  // were 't' is the 'positive offset distance' (for t >= 0 ).
  // If 3 such offset lines intersect at the same offset distance, the intersection 't',
  // or 'time', can be computed solving for 't' in the linear system formed by 3 such equations.
  // The following rational expression is such a solution.
  
  FT num = (a2*b0*c1)-(a2*b1*c0)-(b2*a0*c1)+(b2*a1*c0)+(b1*a0*c2)-(b0*a1*c2);
  FT den = (-a2*b1)+(a2*b0)+(b2*a1)-(b2*a0)+(b1*a0)-(b0*a1);
  
  CGAL_assertion( ! CGAL_NTS certified_is_zero(den) ) ;

  return Quotient<FT>(num,den);
}

// Given 3 oriented lines (a0,b0,c0),(a1,b1,c1) and (a2,b2,c2),
// returns the POSITIVE OFFSET DISTANCE at which the leftward-offsetted lines 
// intersect in a single point and the coordinates (x,y) of such point.
//
// PRECONDITIONS: 
// The line coefficients must be normalized: a²+b²==1 and (a,b) being the inward normal vector 
// There exist such a distance for which the 3 leftward-offsetted lines
// intersect in a single point.
//
template<class FT>
Quotient<FT> construct_offset_lines_isec ( FT const& a0, FT const& b0, FT const& c0
                                         , FT const& a1, FT const& b1, FT const& c1
                                         , FT const& a2, FT const& b2, FT const& c2
                                         , Quotient<FT>& x, Quotient<FT>& y
                                         )
{
  typedef Quotient<FT>  QFT ;
  
  QFT t = compute_offset_lines_isec_time(a0,b0,c0,a1,b1,c1,a2,b2,c2);
  
  // Given 3 offset lines expressed in the form: a*x(t) + b*y(t) + c + t
  // The point of intersection for a given time (or offset distance) 't'
  // can be obtained by solving for x(t),y(t) in a system formed by 2 such offset lines.
  
  FT delta = a0 * b1 - b0 * a1 ;
  
  CGAL_assertion(! CGAL_NTS certified_is_zero(delta) ) ;
  
  QFT x0 = QFT(b0 * c1 - b1 * c0,delta);
  QFT x1 = t * QFT(b0 - b1      ,delta);
  x = x0 + x1 ;
  
  QFT y0 = QFT(a1 * c0 - a0 * c1,delta);
  QFT y1 = t * QFT(a1 - a0      ,delta);
  y = y0 + y1 ;
  
  return t;
}

// Given a point (px,py) and 3 oriented lines (a0,b0,c0),(a1,b1,c1) and (a2,b2,c2),
// such that their leftward-offsetted lines intersect in a single point (ix,iy),
// returns the squared distance between (px,py) and (ix,iy)
// The 'positive' offset of a line is to its left.
//
// PRECONDITIONS: 
// The line coefficients must be normalized: a²+b²==1 and (a,b) being the inward normal vector 
// There exist such a distance for which the 3 leftward-offsetted lines
// intersect in a single point.
//
template<class FT>
Quotient<FT> compute_offset_lines_isec_sdist_to_point ( FT const& px, FT const& py
                                                      , FT const& a0, FT const& b0, FT const& c0
                                                      , FT const& a1, FT const& b1, FT const& c1
                                                      , FT const& a2, FT const& b2, FT const& c2
                                                      )
{
  typedef Quotient<FT>  QFT ;

  QFT ix, iy ;
  
  construct_offset_lines_isec(a0,b0,c0,a1,b1,c1,a2,b2,c2,ix,iy);
   
  QFT dx  = ix - px ;
  QFT dy  = iy - py ;
  QFT dx2 = dx * dx ;
  QFT dy2 = dy * dy ;   
  
  QFT sdist = dx2 + dy2 ;
  
  return sdist;
}

CGAL_END_NAMESPACE

#endif // CGAL_STRAIGHT_SKELETON_CONSTRUCTIONS_FTC2_H //
// EOF //
 
