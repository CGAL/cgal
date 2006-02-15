// Copyright (c) 2005, 2006 Fernando Luis Cacciola Carballal. All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
// Author(s)     : Fernando Cacciola <fernando_cacciola@ciudad.com.ar>
//
#ifndef CGAL_STRAIGHT_SKELETON_CONS_FTC2_H
#define CGAL_STRAIGHT_SKELETON_CONS_FTC2_H 1


CGAL_BEGIN_NAMESPACE

inline MP_Float sqrt( MP_Float const& n )
{
  return CGAL_NTS approximate_sqrt(n);
}

namespace CGAL_SLS_i
{

// Given an oriented 2D stright line edge (px,py)->(qx,qy), computes the normalized coefficients (a,b,c) of the
// supporting line.
// POSTCONDITION: [a,b] is the leftward normal _unit_ (a²+b²=1) vector.
//
template<class FT>
Line<FT> compute_normalized_line_ceoffC2( Edge<FT> const& e )
{
  FT a,b,c ;

  if(e.s().y() == e.t().y())
  {
    a = 0 ;
    if(e.t().x() > e.s().x())
    {
      b = 1;
      c = -e.s().y();
    }
    else if(e.t().x() == e.s().x())
    {
      b = 0;
      c = 0;
    }
    else
    {
      b = -1;
      c = e.s().y();
    }

    CGAL_SSTRAITS_TRACE("Line coefficients for HORIZONTAL line:\npx=" << e.s().x() << "\npy=" << e.s().y()
                        << "\nqx=" << e.t().x() << "\nqy=" << e.t().y()
                        << "\na="<< a << "\nb=" << b << "\nc=" << c
                       ) ;
  }
  else if(e.t().x() == e.s().x())
  {
    b = 0;
    if(e.t().y() > e.s().y())
    {
      a = -1;
      c = e.s().x();
    }
    else if (e.t().y() == e.s().y())
    {
      a = 0;
      c = 0;
    }
    else
    {
      a = 1;
      c = -e.s().x();
    }

    CGAL_SSTRAITS_TRACE("Line coefficients for VERTICAL line:\npx=" << e.s().x() << "\npy=" << e.s().y()
                        << "\nqx=" << e.t().x() << "\nqy=" << e.t().y()
                        << "\na="<< a << "\nb=" << b << "\nc=" << c
                       ) ;

  }
  else
  {
    FT sa = e.s().y() - e.t().y();
    FT sb = e.t().x() - e.s().x();
    FT l  = CGAL_NTS sqrt( (sa*sa) + (sb*sb) );

    a = sa / l ;
    b = sb / l ;

    c = -e.s().x()*a - e.s().y()*b;

    CGAL_SSTRAITS_TRACE("Line coefficients for line:\npx=" << e.s().x() << "\npy=" << e.s().y() << "\nqx="
                        << e.t().x() << "\nqy=" << e.t().y()
                        << "\na="<< a << "\nb=" << b << "\nc=" << c << "\nl:" << l
                       ) ;
  }

  return Line<FT>(a,b,c);
}

// Given 3 oriented straight line segments: e0, e1, e2 [each segment is passed as (sx,sy,tx,ty)]
// returns the OFFSET DISTANCE (n/d) at which the offsetted lines
// intersect at a single point, IFF such intersection exist.
// If the lines intersect to the left, the returned distance is positive.
// If the lines intersect to the right, the returned distance is negative.
// If the lines do not intersect, for example, for collinear edges, or parallel edges but with the same orientation,
// returns 0 (the actual distance is undefined in this case, but 0 is a usefull return)
// NOTE: The result is a explicit rational number returned as a tuple (num,den); the caller must check that den!=0 manually
// (a predicate for instance should return indeterminate in this case)
template<class FT>
Rational<FT> compute_normal_offset_lines_isec_timeC2 ( SortedTriedge<FT> const& triedge )
{
  // DETAILS:
  //
  // An offset line is given by:
  //
  //   a*x(t) + b*y(t) + c - t = 0
  //
  // were 't>0' being to the left of the line.
  // If 3 such offset lines intersect at the same offset distance, the intersection 't',
  // or 'time', can be computed solving for 't' in the linear system formed by 3 such equations.
  // The following rational expression the solution.

  Line<FT> l0 = compute_normalized_line_ceoffC2(triedge.e0()) ;
  Line<FT> l1 = compute_normalized_line_ceoffC2(triedge.e1()) ;
  Line<FT> l2 = compute_normalized_line_ceoffC2(triedge.e2()) ;

  FT num = (l2.a()*l0.b()*l1.c())
          -(l2.a()*l1.b()*l0.c())
          -(l2.b()*l0.a()*l1.c())
          +(l2.b()*l1.a()*l0.c())
          +(l1.b()*l0.a()*l2.c())
          -(l0.b()*l1.a()*l2.c());

  FT den = (-l2.a()*l1.b())
           +(l2.a()*l0.b())
           +(l2.b()*l1.a())
           -(l2.b()*l0.a())
           +(l1.b()*l0.a())
           -(l0.b()*l1.a());

  CGAL_SSTRAITS_TRACE("Normal Event:\nn=" << num << "\nd=" << den  )

  return Rational<FT>(num,den) ;
}


// Given 3 oriented straight line segments: e0, e1, e2 [each segment is passed as (sx,sy,tx,ty)]
// such that e0 and e1 are collinear, not neccesarily consecutive but with the same orientaton, and e2 is NOT
// collinear with e0 and e1; returns the OFFSET DISTANCE (n/d) at which a line perpendicular to e0 (and e1) passing through
// the midpoint of e0.t and e1.s, intersects the offset line of e2
// If the lines intersect to the left of e0, the returned distance is positive.
// If the lines intersect to the right of e0, the returned distance is negative.
// If the lines do not intersect, for example, the three edges are collinear edges, or e0,e1 are not,
// returns 0 (the actual distance is undefined in this case, but 0 is a usefull return)
// NOTE: The result is a explicit rational number returned as a tuple (num,den); the caller must check that den!=0 manually
// (a predicate for instance should return indeterminate in this case)
template<class FT>
Rational<FT> compute_degenerate_offset_lines_isec_timeC2 ( SortedTriedge<FT> const& triedge )
{
  // DETAILS:
  //
  //   (1)
  //   The bisecting line of e0 and e1 (which are required to be collinear) is a line perpendicular to e0 (and e1)
  //   which passes through the midpoint of e0.target and e1.source (called q)
  //   This "degenerate" bisecting line is given by:
  //
  //     B0(t) = q + t*[l0.a,l0.b]
  //
  //   where l0.a and l0.b are the _normalized_ line coefficients for e0.
  //   Since [a,b] is a _unit_ vector pointing perpendicularly to the left of e0 (and e1);
  //   any point B0(k) is at a distance k from the line supporting e0 and e1.
  //
  //   (2)
  //   The bisecting line of e0 and e2 (which are required to be oblique) is given by the following SEL
  //
  //    l0.a*x(t) + l0.b*y(t) + l0.c + t = 0
  //    l2.a*x(t) + l2.b*y(t) + l2.c + t = 0
  //
  //   where (l0.a,l0.b,l0.c) and (l2.a,l2.b,l0.c) are the normalized line coefficientes of e0 and e2 resp.
  //
  //     B1(t)=[x(t),y(t)]
  //
  //   (3)
  //   These two bisecting lines B0(t) and B1(t) intersect (if they do) in a single point 'p' whose distance
  //   to the lines supporting the 3 edges is exactly 't' (since those expressions are precisely parametrized in a distance)
  //   Solving the following vectorial equation:
  //
  //     [x(y),y(t)] = q + t*[l0.a,l0.b]
  //
  //   for t gives the result we want.
  //
  //


  Line<FT> l0 = compute_normalized_line_ceoffC2(triedge.e0()) ;
  Line<FT> l2 = compute_normalized_line_ceoffC2(triedge.e2()) ;

  FT qx = ( triedge.e0().t().x() + triedge.e1().s().x() ) / static_cast<FT>(2.0);

  FT num = (l2.a() * l0.b() - l0.a() * l2.b() ) * qx + l0.b() * l2.c() - l2.b() * l0.c() ;
  FT den = (l0.a() * l0.a() - 1) * l2.b() + ( 1 - l2.a() * l0.a() ) * l0.b() ;

  CGAL_SSTRAITS_TRACE("Degenerate Event:\nn=" << num << "\nd=" << den  )

  return Rational<FT>(num,den) ;
}

template<class FT>
Rational<FT> compute_offset_lines_isec_timeC2 ( SortedTriedge<FT> const& triedge )
{
  return triedge.is_degenerate() ? compute_degenerate_offset_lines_isec_timeC2(triedge)
                                 : compute_normal_offset_lines_isec_timeC2    (triedge) ;
}

// Given 3 oriented lines l0:(l0.a,l0.b,l0.c), l1:(l1.a,l1.b,l1.c) and l2:(l2.a,l2.b,l2.c)
// such that their offsets at a certian distance intersect in a single point, returns the coordinates (x,y) of such a point.
//
// PRECONDITIONS:
// The line coefficients must be normalized: a²+b²==1 and (a,b) being the leftward normal vector
// The offsets at a certain distance do intersect in a single point.
//
template<class FT>
Vertex<FT> construct_normal_offset_lines_isecC2 ( SortedTriedge<FT> const& triedge )
{
  Line<FT> l0 = compute_normalized_line_ceoffC2(triedge.e0()) ;
  Line<FT> l1 = compute_normalized_line_ceoffC2(triedge.e1()) ;
  Line<FT> l2 = compute_normalized_line_ceoffC2(triedge.e2()) ;

  FT den = l0.a()*l2.b() - l0.a()*l1.b() - l1.a()*l2.b() + l2.a()*l1.b() + l0.b()*l1.a() - l0.b()*l2.a();

  CGAL_SSTRAITS_TRACE("Event Point:\n  d=" << den  )

  CGAL_assertion ( ! CGAL_NTS certified_is_zero(den) ) ;

  FT numX = l0.b()*l2.c() - l0.b()*l1.c() - l1.b()*l2.c() + l2.b()*l1.c() + l1.b()*l0.c() - l2.b()*l0.c();
  FT numY = l0.a()*l2.c() - l0.a()*l1.c() - l1.a()*l2.c() + l2.a()*l1.c() + l1.a()*l0.c() - l2.a()*l0.c();

  FT x =  numX / den ;
  FT y = -numY / den ;

  CGAL_SSTRAITS_TRACE("\n  x=" << x << "\n  y=" << y )

  return Vertex<FT>(x,y) ;
}

// Given 3 oriented lines l0:(l0.a,l0.b,l0.c), l1:(l1.a,l1.b,l1.c) and l2:(l2.a,l2.b,l2.c)
// such that their offsets at a certian distance intersect in a single point, returns the coordinates (x,y) of such a point.
//
// PRECONDITIONS:
// The line coefficients must be normalized: a²+b²==1 and (a,b) being the leftward normal vector
// The offsets at a certain distance do intersect in a single point.
//
template<class FT>
Vertex<FT> construct_degenerate_offset_lines_isecC2 ( SortedTriedge<FT> const& triedge )
{
  Line<FT> l0 = compute_normalized_line_ceoffC2(triedge.e0()) ;
  Line<FT> l1 = compute_normalized_line_ceoffC2(triedge.e1()) ;
  Line<FT> l2 = compute_normalized_line_ceoffC2(triedge.e2()) ;

  FT qx = ( triedge.e0().t().x() + triedge.e1().s().x() ) / static_cast<FT>(2.0);
  FT qy = ( triedge.e0().t().y() + triedge.e1().s().y() ) / static_cast<FT>(2.0);

  FT num = (l2.a() * l0.b() - l0.a() * l2.b() ) * qx + l0.b() * l2.c() - l2.b() * l0.c() ;
  FT den = (l0.a() * l0.a() - 1) * l2.b() + ( 1 - l2.a() * l0.a() ) * l0.b() ;

  CGAL_precondition( den != static_cast<FT>(0.0) ) ;

  FT x = qx + l0.a() * num / den  ;
  FT y = qy + l0.b() * num / den  ;

  CGAL_SSTRAITS_TRACE("\n  x=" << x << "\n  y=" << y )

  return Vertex<FT>(x,y) ;
}


template<class FT>
Vertex<FT> construct_offset_lines_isecC2 ( SortedTriedge<FT> const& triedge )
{
  return triedge.is_degenerate() ? construct_degenerate_offset_lines_isecC2(triedge)
                                 : construct_normal_offset_lines_isecC2    (triedge) ;
}

// Give a point (px,py) and 3 oriented lines l0:(l0.a,l0.b,l0.c), l1:(l1.a,l1.b,l1.c) and l2:(l2.a,l2.b,l2.c),
// such that their offsets at a certian distance intersect in a single point (ix,iy),
// returns the squared distance between (px,py) and (ix,iy)
//
// PRECONDITIONS:
// The line coefficients must be normalized: a²+b²==1 and (a,b) being the leftward normal vector
// The offsets at a certain distance do intersect in a single point.
//
template<class FT>
FT compute_offset_lines_isec_sdist_to_pointC2 ( Vertex<FT> const& p, SortedTriedge<FT> const& triedge )
{

  Vertex<FT> i = construct_offset_lines_isecC2(triedge);

  FT dx  = i.x() - p.x() ;
  FT dy  = i.y() - p.y() ;
  FT dx2 = dx * dx ;
  FT dy2 = dy * dy ;

  FT sdist = dx2 + dy2 ;

  return sdist;
}

} // namnepsace CGAIL_SLS_i

CGAL_END_NAMESPACE

#endif // CGAL_STRAIGHT_SKELETON_CONS_FTC2_H //
// EOF //

