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
// file          : include/CGAL/constructions/Straight_skeleton_cons_ftC2.h
// package       : Straight_skeleton_2 (1.1.0)
//
// author(s)     : Fernando Cacciola
// maintainer    : Fernando Cacciola <fernando_cacciola@hotmail>
// coordinator   : Fernando Cacciola <fernando_cacciola@hotmail>
//
// ============================================================================
#ifndef CGAL_STRAIGHT_SKELETON_CONS_FTC2_H
#define CGAL_STRAIGHT_SKELETON_CONS_FTC2_H 1


CGAL_BEGIN_NAMESPACE

// Given an oriented 2D stright line edge (px,py)->(qx,qy), computes the normalized coefficients (a,b,c) of the
// supporting line.
// POSTCONDITION: [a,b] is the leftward normal _unit_ (a²+b²=1) vector.
//
template<class FT>
tuple<FT,FT,FT>
compute_normalized_line_ceoffC2( tuple<FT,FT,FT,FT> const& edge )
{
  FT px,py,qx,qy ;
  FT a,b,c ;

  tie(px,py,qx,qy) = edge ;

  if(py == qy)
  {
    a = 0 ;
    if(qx > px)
    {
      b = 1;
      c = -py;
    }
    else if(qx == px)
    {
      b = 0;
      c = 0;
    }
    else
    {
      b = -1;
      c = py;
    }

    CGAL_SSTRAITS_TRACE("Line coefficients for HORIZONTAL line:\npx=" << px << "\npy=" << py << "\nqx=" << qx << "\nqy=" << qy
                        << "\na="<< a << "\nb=" << b << "\nc=" << c
                       ) ;
  }
  else if(qx == px)
  {
    b = 0;
    if(qy > py)
    {
      a = -1;
      c = px;
    }
    else if (qy == py)
    {
      a = 0;
      c = 0;
    }
    else
    {
      a = 1;
      c = -px;
    }

    CGAL_SSTRAITS_TRACE("Line coefficients for VERTICAL line:\npx=" << px << "\npy=" << py << "\nqx=" << qx << "\nqy=" << qy
                        << "\na="<< a << "\nb=" << b << "\nc=" << c
                       ) ;

  }
  else
  {
    FT sa = py - qy;
    FT sb = qx - px;
    FT l  = CGAL_NTS sqrt( (sa*sa) + (sb*sb) );

    a = sa / l ;
    b = sb / l ;

    c = -px*a - py*b;

    CGAL_SSTRAITS_TRACE("Line coefficients for line:\npx=" << px << "\npy=" << py << "\nqx=" << qx << "\nqy=" << qy
                        << "\na="<< a << "\nb=" << b << "\nc=" << c << "\nl:" << l
                       ) ;
  }


  return make_tuple(a,b,c);
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
tuple<FT,FT> compute_normal_offset_lines_isec_timeC2 ( tuple<FT,FT,FT,FT> const& e0
                                                     , tuple<FT,FT,FT,FT> const& e1
                                                     , tuple<FT,FT,FT,FT> const& e2
                                                     )
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

  FT a0,b0,c0,a1,b1,c1,a2,b2,c2 ;

  tie(a0,b0,c0) = compute_normalized_line_ceoffC2(e0) ;
  tie(a1,b1,c1) = compute_normalized_line_ceoffC2(e1) ;
  tie(a2,b2,c2) = compute_normalized_line_ceoffC2(e2) ;


  FT num = (a2*b0*c1)-(a2*b1*c0)-(b2*a0*c1)+(b2*a1*c0)+(b1*a0*c2)-(b0*a1*c2);
  FT den = (-a2*b1)+(a2*b0)+(b2*a1)-(b2*a0)+(b1*a0)-(b0*a1);

  CGAL_SSTRAITS_TRACE("Normal Event:\nn=" << num << "\nd=" << den  )

  return make_tuple(num,den) ;
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
tuple<FT,FT> compute_degenerate_offset_lines_isec_timeC2 ( tuple<FT,FT,FT,FT> const& e0
                                                         , tuple<FT,FT,FT,FT> const& e1
                                                         , tuple<FT,FT,FT,FT> const& e2
                                                         )
{
  // DETAILS:
  //
  //   (1)
  //   The bisecting line of e0 and e1 (which are required to be collinear) is a line perpendicular to e0 (and e1)
  //   which passes through the midpoint of e0.target and e1.source (called q)
  //   This "degenerate" bisecting line is given by:
  //
  //     B0(t) = q + t*[a0,b0]
  //
  //   where a0 and b0 are the _normalized_ line coefficients for e0.
  //   Since [a,b] is a _unit_ vector pointing perpendicularly to the left of e0 (and e1);
  //   any point B0(k) is at a distance k from the line supporting e0 and e1.
  //
  //   (2)
  //   The bisecting line of e0 and e2 (which are required to be oblique) is given by the following SEL
  //
  //    a0*x(t) + b0*y(t) + c0 + t = 0
  //    a2*x(t) + b2*y(t) + c2 + t = 0
  //
  //   where (a0,b0,c0) and (a2,b2,c0) are the normalized line coefficientes of e0 and e2 resp.
  //
  //     B1(t)=[x(t),y(t)]
  //
  //   (3)
  //   These two bisecting lines B0(t) and B1(t) intersect (if they do) in a single point 'p' whose distance
  //   to the lines supporting the 3 edges is exactly 't' (since those expressions are precisely parametrized in a distance)
  //   Solving the following vectorial equation:
  //
  //     [x(y),y(t)] = q + t*[a0,b0]
  //
  //   for t gives the result we want.
  //
  //
  FT a0,b0,c0,a2,b2,c2 ;

  tie(a0,b0,c0) = compute_normalized_line_ceoffC2(e0) ;
  tie(a2,b2,c2) = compute_normalized_line_ceoffC2(e2) ;

  FT e0tx = boost::get<2>(e0);
  FT e1sx = boost::get<0>(e1);
  FT qx   = ( e0tx + e1sx ) / static_cast<FT>(2.0);

  FT num = (a2 * b0 - a0 * b2 ) * qx + b0 * c2 - b2 * c0 ;
  FT den = (a0 * a0 - 1) * b2 + ( 1 - a2 * a0 ) * b0 ;

  CGAL_SSTRAITS_TRACE("Degenerate Event:\nn=" << num << "\nd=" << den  )

  return make_tuple(num,den) ;
}

template<class FT>
tuple<FT,FT> compute_offset_lines_isec_timeC2 ( tuple<FT,FT,FT,FT> const& l0
                                              , tuple<FT,FT,FT,FT> const& l1
                                              , tuple<FT,FT,FT,FT> const& l2
                                              , bool                      is_degenerate
                                              )
{
  return is_degenerate ? compute_degenerate_offset_lines_isec_timeC2(l0,l1,l2)
                       : compute_normal_offset_lines_isec_timeC2    (l0,l1,l2) ;
}

// Given 3 oriented lines l0:(a0,b0,c0), l1:(a1,b1,c1) and l2:(a2,b2,c2)
// such that their offsets at a certian distance intersect in a single point, returns the coordinates (x,y) of such a point.
//
// PRECONDITIONS:
// The line coefficients must be normalized: a²+b²==1 and (a,b) being the leftward normal vector
// The offsets at a certain distance do intersect in a single point.
//
template<class FT>
tuple<FT,FT> construct_offset_lines_isecC2 ( tuple<FT,FT,FT,FT> const& l0
                                           , tuple<FT,FT,FT,FT> const& l1
                                           , tuple<FT,FT,FT,FT> const& l2
                                           )
{
  FT a0,b0,c0,a1,b1,c1,a2,b2,c2 ;

  tie(a0,b0,c0) = compute_normalized_line_ceoffC2(l0) ;
  tie(a1,b1,c1) = compute_normalized_line_ceoffC2(l1) ;
  tie(a2,b2,c2) = compute_normalized_line_ceoffC2(l2) ;

  FT den = a0*b2 - a0*b1 - a1*b2 + a2*b1 + b0*a1 - b0*a2;

  CGAL_SSTRAITS_TRACE("Event Point:\n  d=" << den  )

  CGAL_assertion(! CGAL_NTS certified_is_zero(den) ) ;

  FT numX = b0*c2 - b0*c1 - b1*c2 + b2*c1 + b1*c0 - b2*c0;
  FT numY = a0*c2 - a0*c1 - a1*c2 + a2*c1 + a1*c0 - a2*c0;

  FT x =  numX / den ;
  FT y = -numY / den ;

  CGAL_SSTRAITS_TRACE("\n  x=" << x << "\n  y=" << y )

  return make_tuple(x,y) ;
}


// Give a point (px,py) and 3 oriented lines l0:(a0,b0,c0), l1:(a1,b1,c1) and l2:(a2,b2,c2),
// such that their offsets at a certian distance intersect in a single point (ix,iy),
// returns the squared distance between (px,py) and (ix,iy)
//
// PRECONDITIONS:
// The line coefficients must be normalized: a²+b²==1 and (a,b) being the leftward normal vector
// The offsets at a certain distance do intersect in a single point.
//
template<class FT>
FT compute_offset_lines_isec_sdist_to_pointC2 ( tuple<FT,FT>       const& p
                                              , tuple<FT,FT,FT,FT> const& l0
                                              , tuple<FT,FT,FT,FT> const& l1
                                              , tuple<FT,FT,FT,FT> const& l2
                                              )
{
  FT px,py,ix,iy ;

  tie(px,py) = p ;
  tie(ix,iy) = construct_offset_lines_isecC2(l0,l1,l2);

  FT dx  = ix - px ;
  FT dy  = iy - py ;
  FT dx2 = dx * dx ;
  FT dy2 = dy * dy ;

  FT sdist = dx2 + dy2 ;

  return sdist;
}

CGAL_END_NAMESPACE

#endif // CGAL_STRAIGHT_SKELETON_CONS_FTC2_H //
// EOF //

