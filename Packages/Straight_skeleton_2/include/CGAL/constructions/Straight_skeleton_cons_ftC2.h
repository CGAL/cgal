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

// Given 2 points (px,py) and (qx,qy) computes the normalized coefficients (a,b,c) of the line passing through both of them.
// POSTCONDITION: a²+b²=1 and (a,b) is the leftward normal vector
template<class FT>
tuple<FT,FT,FT>
compute_normalized_line_ceoffC2( tuple<FT,FT> const& p, tuple<FT,FT> const& q )
{
  FT px,py,qx,qy,a,b,c ;

  tie(px,py) = p ;
  tie(qx,qy) = q ;

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

// Given 3 oriented lines l0:(a0,b0,c0), l1:(a1,b1,c1) and l2:(a2,b2,c2),
// returns the OFFSET DISTANCE (n/d) at which the offsetted lines
// intersect at a single point, IFF such intersection exist.
// If the lines intersect to the left, the returned distance is positive.
// If the lines intersect to the right, the returned distance is negative.
// If the lines do not intersect, returns 0 (the actual distance is undefined in this case, but 0 is a usefull return)
// PRECONDITIONS:
// The line coefficients must be normalized: a²+b²==1 and (a,b) being the leftward normal vector
// NOTE: The result is a explicit rational number returned as a tuple (num,den); the caller must check that den!=0 manually
// (a predicate for instance should return indeterminate in this case)
template<class FT>
tuple<FT,FT> compute_offset_lines_isec_timeC2 ( tuple<FT,FT,FT> const& l0
                                              , tuple<FT,FT,FT> const& l1
                                              , tuple<FT,FT,FT> const& l2
                                              )
{
  FT a0,b0,c0,a1,b1,c1,a2,b2,c2 ;

  tie(a0,b0,c0) = l0 ;
  tie(a1,b1,c1) = l1 ;
  tie(a2,b2,c2) = l2 ;

  // An offset line is given by: a*x(t) + b*y(t) + c - t = 0
  // were 't>0' being to the left of the line.
  // If 3 such offset lines intersect at the same offset distance, the intersection 't',
  // or 'time', can be computed solving for 't' in the linear system formed by 3 such equations.
  // The following rational expression the solution.

  FT num = (a2*b0*c1)-(a2*b1*c0)-(b2*a0*c1)+(b2*a1*c0)+(b1*a0*c2)-(b0*a1*c2);
  FT den = (-a2*b1)+(a2*b0)+(b2*a1)-(b2*a0)+(b1*a0)-(b0*a1);

  CGAL_SSTRAITS_TRACE("Event:\nn=" << num << "\nd=" << den  )


  return make_tuple(num,den) ;
}

// Given 3 oriented lines l0:(a0,b0,c0), l1:(a1,b1,c1) and l2:(a2,b2,c2)
// such that their offsets at a certian distance intersect in a single point, returns the coordinates (x,y) of such a point.
//
// PRECONDITIONS:
// The line coefficients must be normalized: a²+b²==1 and (a,b) being the leftward normal vector
// The offsets at a certain distance do intersect in a single point.
//
template<class FT>
tuple<FT,FT> construct_offset_lines_isecC2 ( tuple<FT,FT,FT> const& l0
                                           , tuple<FT,FT,FT> const& l1
                                           , tuple<FT,FT,FT> const& l2
                                           )
{
  FT a0,b0,c0,a1,b1,c1,a2,b2,c2 ;

  tie(a0,b0,c0) = l0 ;
  tie(a1,b1,c1) = l1 ;
  tie(a2,b2,c2) = l2 ;

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
FT compute_offset_lines_isec_sdist_to_pointC2 ( tuple<FT,FT>    const& p
                                              , tuple<FT,FT,FT> const& l0
                                              , tuple<FT,FT,FT> const& l1
                                              , tuple<FT,FT,FT> const& l2
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

