// Copyright (c) 2006 Fernando Luis Cacciola Carballal. All rights reserved.
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

namespace CGAL_SS_i
{

template<class NT>
inline NT inexact_sqrt( NT const& n )
{
  return CGAL_NTS sqrt(n);
}

inline MP_Float inexact_sqrt( MP_Float const& n )
{
  double d = CGAL::to_double(n);
  
  if ( !CGAL_NTS is_finite(d) )
    d = CGAL_NTS sign(n) == NEGATIVE ? - std::numeric_limits<double>::max() 
                                     :   std::numeric_limits<double>::max() ;
       
  return MP_Float( CGAL_NTS sqrt(d) ) ;
}

inline Quotient<MP_Float> inexact_sqrt( Quotient<MP_Float> const& q )
{
  CGAL_precondition(q > 0);
  return Quotient<MP_Float>(CGAL_SS_i::inexact_sqrt(q.numerator()*q.denominator())
                           ,q.denominator()
                           );
}

// Given an oriented 2D straight line segment 'e', computes the normalized coefficients (a,b,c) of the
// supporting line.
// POSTCONDITION: [a,b] is the leftward normal _unit_ (a²+b²=1) vector.
// POSTCONDITION: In case of overflow, an empty optional<> is returned.
template<class K>
optional< Line_2<K> > compute_normalized_line_ceoffC2( Segment_2<K> const& e )
{
  bool finite = true ;
  
  typedef typename K::FT FT ;
  
  FT a (0.0),b (0.0) ,c(0.0)  ;

  if(e.source().y() == e.target().y())
  {
    a = 0 ;
    if(e.target().x() > e.source().x())
    {
      b = 1;
      c = -e.source().y();
    }
    else if(e.target().x() == e.source().x())
    {
      b = 0;
      c = 0;
    }
    else
    {
      b = -1;
      c = e.source().y();
    }

    CGAL_SSTRAITS_TRACE("Line coefficients for HORIZONTAL line:\npx=" << e.source().x() << "\npy=" << e.source().y()
                        << "\nqx=" << e.target().x() << "\nqy=" << e.target().y()
                        << "\na="<< a << "\nb=" << b << "\nc=" << c
                       ) ;
  }
  else if(e.target().x() == e.source().x())
  {
    b = 0;
    if(e.target().y() > e.source().y())
    {
      a = -1;
      c = e.source().x();
    }
    else if (e.target().y() == e.source().y())
    {
      a = 0;
      c = 0;
    }
    else
    {
      a = 1;
      c = -e.source().x();
    }

    CGAL_SSTRAITS_TRACE("Line coefficients for VERTICAL line:\npx=" << e.source().x() << "\npy=" << e.source().y()
                        << "\nqx=" << e.target().x() << "\nqy=" << e.target().y()
                        << "\na="<< a << "\nb=" << b << "\nc=" << c
                       ) ;
  }
  else
  {
    FT sa = e.source().y() - e.target().y();
    FT sb = e.target().x() - e.source().x();
    FT l2 = (sa*sa) + (sb*sb) ;

    if ( CGAL_NTS is_finite(l2) )
    {
      FT l = CGAL_SS_i :: inexact_sqrt(l2);
  
      a = sa / l ;
      b = sb / l ;
  
      c = -e.source().x()*a - e.source().y()*b;
    }
    else finite = false ;
    
    CGAL_SSTRAITS_TRACE("Line coefficients for line:\npx=" << e.source().x() << "\npy=" << e.source().y() << "\nqx="
                        << e.target().x() << "\nqy=" << e.target().y()
                        << "\na="<< a << "\nb=" << b << "\nc=" << c << "\nl:" << l
                       ) ;
  }
  
  if ( !finite )
    if ( !CGAL_NTS is_finite(a) || !CGAL_NTS is_finite(b) || !CGAL_NTS is_finite(c) ) 
      finite = false ;

  return cgal_make_optional( finite, K().construct_line_2_object()(a,b,c) ) ;
}

// Given 3 oriented straight line segments: e0, e1, e2 (passed in a SortedTriedge record)
// returns the OFFSET DISTANCE (n/d) at which the offsetted lines
// intersect at a single point, IFF such intersection exist.
// If the lines intersect to the left, the returned distance is positive.
// If the lines intersect to the right, the returned distance is negative.
// If the lines do not intersect, for example, for collinear edges, or parallel edges but with the same orientation,
// returns 0 (the actual distance is undefined in this case, but 0 is a usefull return)
//
// NOTE: The result is a explicit rational number returned as a tuple (num,den); the caller must check that den!=0 manually
// (a predicate for instance should return indeterminate in this case)
//
// PRECONDITION: None of e0, e1 and e2 are collinear (but two of them can be parallel)
//
// POSTCONDITION: In case of overflow an empty optional is returned.
//
template<class K>
optional< Rational< typename K::FT> > compute_normal_offset_lines_isec_timeC2 ( Sorted_triedge_2<K> const& triedge )
{
  typedef typename K::FT  FT ;
  
  typedef Line_2<K> Line_2 ;
  
  typedef optional<Line_2> Optional_line_2 ;
  
  
  FT num(0.0), den(0.0) ;
  
  // DETAILS:
  //
  // An offset line is given by:
  //
  //   a*x(t) + b*y(t) + c - t = 0
  //
  // were 't > 0' being to the left of the line.
  // If 3 such offset lines intersect at the same offset distance, the intersection 't',
  // or 'time', can be computed solving for 't' in the linear system formed by 3 such equations.
  // The result is :
  //
  //  t = a2*b0*c1 - a2*b1*c0 - b2*a0*c1 + b2*a1*c0 + b1*a0*c2 - b0*a1*c2
  //      ---------------------------------------------------------------
  //             -a2*b1 + a2*b0 + b2*a1 - b2*a0 + b1*a0 - b0*a1 ;

  bool ok = false ;
  
  Optional_line_2 l0 = compute_normalized_line_ceoffC2(triedge.e0()) ;
  Optional_line_2 l1 = compute_normalized_line_ceoffC2(triedge.e1()) ;
  Optional_line_2 l2 = compute_normalized_line_ceoffC2(triedge.e2()) ;

  if ( l0 && l1 && l2 )
  {
    num = (l2->a()*l0->b()*l1->c())
         -(l2->a()*l1->b()*l0->c())
         -(l2->b()*l0->a()*l1->c())
         +(l2->b()*l1->a()*l0->c())
         +(l1->b()*l0->a()*l2->c())
         -(l0->b()*l1->a()*l2->c());
      
    den = (-l2->a()*l1->b())
         +( l2->a()*l0->b())
         +( l2->b()*l1->a())
         -( l2->b()*l0->a())
         +( l1->b()*l0->a())
         -( l0->b()*l1->a());
         
    ok = CGAL_NTS is_finite(num) && CGAL_NTS is_finite(den);     
  }
  
  CGAL_SSTRAITS_TRACE("Normal Event:\nn=" << num << "\nd=" << den  )

  return cgal_make_optional(ok,Rational<FT>(num,den)) ;
}

// Given two oriented straight line segments e0 and e1 such that e-next follows e-prev, returns
// the coordinates of the midpoint of the segment between e-prev and e-next.
// NOTE: the edges can be oriented e0->e1 or e1->e0
//
// POSTCONDITION: In case of overflow an empty optional is returned.
//
template<class K>
optional< Point_2<K> > compute_oriented_midpoint ( Segment_2<K> const& e0, Segment_2<K> const& e1 )
{
  bool ok = false ;
  
  typedef typename K::FT FT ;
  
  FT delta01 = CGAL::squared_distance(e0.target(),e1.source());
  FT delta10 = CGAL::squared_distance(e1.target(),e0.source());
  
  Point_2<K> mp ;
   
  if ( CGAL_NTS is_finite(delta01) &&  CGAL_NTS is_finite(delta10) )
  {
    if ( delta01 <= delta10 )
         mp = CGAL::midpoint(e0.target(),e1.source());
    else mp = CGAL::midpoint(e1.target(),e0.source());
    
    ok = CGAL_NTS is_finite(mp.x()) && CGAL_NTS is_finite(mp.y());
  }
  
  return cgal_make_optional(ok,mp);
}

// Given 3 oriented straight line segments: e0, e1, e2 (passed in a SortedTriedge record)
// such that e0 and e1 are collinear, not neccesarily consecutive but with the same orientaton, and e2 is NOT
// collinear with e0 and e1; returns the OFFSET DISTANCE (n/d) at which a line perpendicular to e0 (and e1)
// passing through the oriented midpoint between e0 and e1 (see the function above) intersects the offset line of e2
//
// If the lines intersect to the left of e0, the returned distance is positive.
// If the lines intersect to the right of e0, the returned distance is negative.
// If the lines do not intersect, for example, the three edges are collinear edges, or e0,e1 are not,
// returns 0/0 (the actual distance is undefined in this case, but 0 is a usefull return)
//
// NOTE: The result is a explicit rational number returned as a tuple (num,den); the caller must check that den!=0 manually
// (a predicate for instance should return indeterminate in this case)
//
// POSTCONDITION: In case of overflow an empty optional is returned.
//
template<class K>
optional< Rational< typename K::FT> > compute_degenerate_offset_lines_isec_timeC2 ( Sorted_triedge_2<K> const& triedge )
{
  typedef typename K::FT FT ;
  
  typedef Point_2<K> Point_2 ;
  typedef Line_2 <K> Line_2 ;
  
  typedef optional<Point_2> Optional_point_2 ;
  typedef optional<Line_2>  Optional_line_2 ;
  
  // DETAILS:
  //
  //   (1)
  //   The bisecting line of e0 and e1 (which are required to be collinear) is a line perpendicular to e0 (and e1)
  //   which passes through the oriented midpoint, 'q', between e0->e1 (or e1->e0, depending on their orientatinon)
  //   This "degenerate" bisecting line is given by:
  //
  //     B0(t) = q + t*[l0.a,l0.b]
  //
  //   where l0.a and l0.b are the _normalized_ line coefficients for e0 (or e1 which is the same)
  //   Since [a,b] is a _unit_ vector pointing perpendicularly to the left of e0 (and e1);
  //   any point B0(k) is at a distance k from the line supporting e0 and e1.
  //
  //   (2)
  //   The bisecting line of e0 and e2 (which are required to be non-parallel) is given by the following SEL
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
  bool ok = false ;

  Optional_line_2 l0 = compute_normalized_line_ceoffC2(triedge.e0()) ;
  Optional_line_2 l2 = compute_normalized_line_ceoffC2(triedge.e2()) ;

  Optional_point_2 q = compute_oriented_midpoint(triedge.e0(),triedge.e1());
  
  FT num(0.0), den(0.0) ;

  if ( l0 && l2 && q )
  {
    if ( ! CGAL_NTS is_zero(l0->b()) ) // Non-vertical
    {
      num = (l2->a() * l0->b() - l0->a() * l2->b() ) * q->x() + l0->b() * l2->c() - l2->b() * l0->c() ;
      den = (l0->a() * l0->a() - 1) * l2->b() + ( 1 - l2->a() * l0->a() ) * l0->b() ;
      
      CGAL_SSTRAITS_TRACE("Non-vertical Degenerate Event:\nn=" << num << "\nd=" << den  )
    }
    else
    {
      num = (l2->a() * l0->b() - l0->a() * l2->b() ) * q->y() - l0->a() * l2->c() + l2->a() * l0->c() ;
      den = l0->a() * l0->b() * l2->b() - l0->b() * l0->b() * l2->a() + l2->a() - l0->a() ;
      
      CGAL_SSTRAITS_TRACE("Vertical Degenerate Event:\nn=" << num << "\nd=" << den  )
    }
    
    ok = CGAL_NTS is_finite(num) && CGAL_NTS is_finite(den);     
  }
  

  return cgal_make_optional(ok,Rational<FT>(num,den)) ;
}

//
// Calls the appropiate function depending on the collinearity of the edges.
//
template<class K>
optional< Rational< typename K::FT > > compute_offset_lines_isec_timeC2 ( Sorted_triedge_2<K> const& triedge )
{
  CGAL_precondition ( triedge.collinear_count() < 3 ) ;
  
  return triedge.collinear_count() == 0 ? compute_normal_offset_lines_isec_timeC2    (triedge)
                                        : compute_degenerate_offset_lines_isec_timeC2(triedge);
}

// Given 3 oriented line segments e0, e1 and e2 (passed in a SortedTriedge record)
// such that their offsets at a certian distance intersect in a single point, 
// returns the coordinates (x,y) of such a point.
//
// PRECONDITIONS:
// None of e0, e1 and e2 are collinear (but two of them can be parallel)
// The line coefficients must be normalized: a²+b²==1 and (a,b) being the leftward normal vector
// The offsets at a certain distance do intersect in a single point.
//
// POSTCONDITION: In case of overflow an empty optional is returned.
//
template<class K>
optional< Point_2<K> > construct_normal_offset_lines_isecC2 ( Sorted_triedge_2<K> const& triedge )
{
  typedef typename K::FT  FT ;
  
  typedef Point_2<K> Point_2 ;
  typedef Line_2<K>  Line_2 ;
  
  typedef optional<Line_2>  Optional_line_2 ;
  
  FT x(0.0),y(0.0) ;
  
  Optional_line_2 l0 = compute_normalized_line_ceoffC2(triedge.e0()) ;
  Optional_line_2 l1 = compute_normalized_line_ceoffC2(triedge.e1()) ;
  Optional_line_2 l2 = compute_normalized_line_ceoffC2(triedge.e2()) ;

  bool ok = false ;
  
  if ( l0 && l1 && l2 )
  {
    FT den = l0->a()*l2->b() - l0->a()*l1->b() - l1->a()*l2->b() + l2->a()*l1->b() + l0->b()*l1->a() - l0->b()*l2->a();
  
    CGAL_SSTRAITS_TRACE("Event Point:\n  d=" << den  )
  
    CGAL_assertion ( ! CGAL_NTS certified_is_zero(den) ) ;
  
    FT numX = l0->b()*l2->c() - l0->b()*l1->c() - l1->b()*l2->c() + l2->b()*l1->c() + l1->b()*l0->c() - l2->b()*l0->c();
    FT numY = l0->a()*l2->c() - l0->a()*l1->c() - l1->a()*l2->c() + l2->a()*l1->c() + l1->a()*l0->c() - l2->a()*l0->c();
  
    if ( CGAL_NTS is_finite(den) && CGAL_NTS is_finite(numX) && CGAL_NTS is_finite(numY)  )
    {
      ok = true ;
      
      x =  numX / den ;
      y = -numY / den ;
    }
    
  }
    
  CGAL_SSTRAITS_TRACE("\n  x=" << x << "\n  y=" << y )
    
  return cgal_make_optional(ok,K().construct_point_2_object()(x,y)) ;
}

// Given 3 oriented line segments e0, e1 and e2 (passed in a SortedTriedge record)
// such that their offsets at a certian distance intersect in a single point, 
// returns the coordinates (x,y) of such a point.
// e0 and e1 are collinear, not neccesarily consecutive but with the same orientaton,
// and e2 is NOT collinear with e0 and e1. 
//
// PRECONDITIONS:
// The line coefficients must be normalized: a²+b²==1 and (a,b) being the leftward normal vector
// The offsets at a certain distance do intersect in a single point.
//
// POSTCONDITION: In case of overflow an empty optional is returned.
//
template<class K>
optional< Point_2<K> > construct_degenerate_offset_lines_isecC2 ( Sorted_triedge_2<K> const& triedge )
{
  typedef typename K::FT FT ;
  
  typedef Point_2<K> Point_2 ;
  typedef Line_2<K>  Line_2 ;
  
  typedef optional<Point_2> Optional_point_2 ;
  typedef optional<Line_2>  Optional_line_2 ;
  
  FT x(0.0),y(0.0) ;
  
  Optional_line_2 l0 = compute_normalized_line_ceoffC2(triedge.e0()) ;
  Optional_line_2 l1 = compute_normalized_line_ceoffC2(triedge.e1()) ;
  Optional_line_2 l2 = compute_normalized_line_ceoffC2(triedge.e2()) ;
  
  Optional_point_2 q = compute_oriented_midpoint(triedge.e0(),triedge.e1());

  bool ok = false ;
  
  if ( l0 && l1 && l2 && q )
  {
    FT num, den ;
  
    if ( ! CGAL_NTS is_zero(l0->b()) ) // Non-vertical
    {
      num = (l2->a() * l0->b() - l0->a() * l2->b() ) * q->x() + l0->b() * l2->c() - l2->b() * l0->c() ;
      den = (l0->a() * l0->a() - 1) * l2->b() + ( 1 - l2->a() * l0->a() ) * l0->b() ;
    }
    else
    {
      num = (l2->a() * l0->b() - l0->a() * l2->b() ) * q->y() - l0->a() * l2->c() + l2->a() * l0->c() ;
      den = l0->a() * l0->b() * l2->b() - l0->b() * l0->b() * l2->a() + l2->a() - l0->a() ;
    }
  
    CGAL_precondition( ! CGAL_NTS certified_is_zero(den) ) ;
  
    if ( CGAL_NTS is_finite(den) && CGAL_NTS is_finite(num) )
    {
      x = q->x() + l0->a() * num / den  ;
      y = q->y() + l0->b() * num / den  ;
      
      ok = CGAL_NTS is_finite(x) && CGAL_NTS is_finite(y) ;
    }
  }
  

  CGAL_SSTRAITS_TRACE("\n  x=" << x << "\n  y=" << y )

  return cgal_make_optional(ok,K().construct_point_2_object()(x,y)) ;
}

//
// Calls the appropiate function depending on the collinearity of the edges.
//
template<class K>
optional< Point_2<K> > construct_offset_lines_isecC2 ( Sorted_triedge_2<K> const& triedge )
{
  CGAL_precondition ( triedge.collinear_count() < 3 ) ;
  
  return triedge.collinear_count() == 0 ? construct_normal_offset_lines_isecC2    (triedge)
                                        : construct_degenerate_offset_lines_isecC2(triedge) ;
}

// Give a point (px,py) and 3 oriented straight line segments e0,e1 and e2.
// such that their offsets at a certian distance intersect in a single point (ix,iy),
// returns the squared distance between (px,py) and (ix,iy)
//
// NOTE: e0,e1 _can_ be collinear, but e2 must not be collinear with e0 nor e1.
//
// PRECONDITIONS:
// The line coefficients must be normalized: a²+b²==1 and (a,b) being the leftward normal vector
// The offsets at a certain distance do intersect in a single point.
//
template<class K>
optional< typename K::FT> compute_offset_lines_isec_dist_to_pointC2 ( optional< Point_2<K> > const& p
                                                                    , Sorted_triedge_2<K>    const& triedge 
                                                                    )
{
  typedef typename K::FT FT ;
  
  typedef Point_2<K> Point_2 ;
  
  typedef optional<Point_2>  Optional_point_2 ;
  
  bool ok = false ;
  
  FT sdist(0.0) ;
  
  if ( p )
  {
    Optional_point_2 i = construct_offset_lines_isecC2(triedge);
    
    if ( i )
    {
      FT dx  = i->x() - p->x() ;
      FT dy  = i->y() - p->y() ;
      FT dx2 = dx * dx ;
      FT dy2 = dy * dy ;
    
      sdist = dx2 + dy2 ;
      
      ok = CGAL_NTS is_finite(sdist);
    }
  }

  return cgal_make_optional(ok,sdist);
}

} // namnepsace CGAIL_SS_i

CGAL_END_NAMESPACE

#endif // CGAL_STRAIGHT_SKELETON_CONS_FTC2_H //
// EOF //

