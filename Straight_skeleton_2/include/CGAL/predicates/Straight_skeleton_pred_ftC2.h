// Copyright (c) 2006 Fernando Luis Cacciola Carballal. All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s)     : Fernando Cacciola <fernando_cacciola@ciudad.com.ar>
//
#ifndef CGAL_STRAIGHT_SKELETON_PREDICATES_FTC2_H 
#define CGAL_STRAIGHT_SKELETON_PREDICATES_FTC2_H 1

#include <CGAL/license/Straight_skeleton_2.h>


#include <CGAL/constructions/Straight_skeleton_cons_ftC2.h>
#include <CGAL/Uncertain.h>
#include <CGAL/certified_quotient_predicates.h>

namespace CGAL {

namespace CGAL_SS_i
{

// Just like the uncertified collinear() returns true IFF r lies in the line p->q
// NOTE: r might be in the ray from p or q containing q or p, that is, there is no ordering implied, just that
// the three points are along the same line, in any order.
template<class K>
inline
Uncertain<bool> certified_collinearC2( Point_2<K> const& p
                                     , Point_2<K> const& q
                                     , Point_2<K> const& r 
                                     )
{
  return CGAL_NTS certified_is_equal( ( q.x() - p.x() ) * ( r.y() - p.y() )
                                    , ( r.x() - p.x() ) * ( q.y() - p.y() )
                                    );
}


// Just like the uncertified collinear_are_ordered_along_lineC2() returns true IFF, given p,q,r along the same line,
// q is in the closed segment [p,r].
template<class K>
inline
Uncertain<bool> certified_collinear_are_ordered_along_lineC2( Point_2<K> const& p
                                                            , Point_2<K> const& q
                                                            , Point_2<K> const& r 
                                                            )
{
  if ( CGAL_NTS certainly(p.x() < q.x()) ) return !(r.x() < q.x());
  if ( CGAL_NTS certainly(q.x() < p.x()) ) return !(q.x() < r.x());
  if ( CGAL_NTS certainly(p.y() < q.y()) ) return !(r.y() < q.y());
  if ( CGAL_NTS certainly(q.y() < p.y()) ) return !(q.y() < r.y());
  
  if ( CGAL_NTS certainly(p.x() == q.x()) && CGAL_NTS certainly(p.y() == q.y())  ) return true;
  
  return Uncertain<bool>::indeterminate();
}

// Returns true IFF segments e0,e1 share the same supporting line
template<class K>
Uncertain<bool> are_edges_collinearC2( Segment_2<K> const& e0, Segment_2<K> const& e1 )
{
  return   certified_collinearC2(e0.source(),e0.target(),e1.source())
         & certified_collinearC2(e0.source(),e0.target(),e1.target()) ;
}

// Returns true IFF the supporting lines for segments e0,e1 are parallel (or the same)
template<class K>
inline
Uncertain<bool> are_edges_parallelC2( Segment_2<K> const& e0, Segment_2<K> const& e1 )
{
  Uncertain<Sign> s = certified_sign_of_determinant2x2(e0.target().x() - e0.source().x()
                                                      ,e0.target().y() - e0.source().y()
                                                      ,e1.target().x() - e1.source().x()
                                                      ,e1.target().y() - e1.source().y()
                                                     ) ;
                                                     
  return ( s == Uncertain<Sign>(ZERO) ) ;
}


// Returns true IFF the parallel segments are equally oriented.
// Precondition: the segments must be parallel.
// the three points are along the same line, in any order.
template<class K>
inline
Uncertain<bool> are_parallel_edges_equally_orientedC2( Segment_2<K> const& e0, Segment_2<K> const& e1 )
{
  return CGAL_NTS certified_sign( (e0.target() - e0.source()) * (e1.target() - e1.source()) ) == POSITIVE;
}


// Returns true IFF segments e0,e1 share the same supporting line but do not overlap except at the vetices, and have the same orientation.
// NOTE: If e1 goes back over e0 (a degenerate antenna or alley) this returns false.
template<class K>
Uncertain<bool> are_edges_orderly_collinearC2( Segment_2<K> const& e0, Segment_2<K> const& e1 )
{
  return are_edges_collinearC2(e0,e1) & are_parallel_edges_equally_orientedC2(e0,e1);
}

template <class FT>
inline
Uncertain<Sign> certified_side_of_oriented_lineC2(const FT &a, const FT &b, const FT &c, const FT &x, const FT &y)
{
  return CGAL_NTS certified_sign(a*x+b*y+c);
}


//
// Constructs a Trisegment_2 which stores 3 edges (segments) such that
// if two of them are collinear, they are put first, as e0, e1.
// Stores also the number of collinear edges. which should be 0 or 2.
//
// If the collinearity test is indeterminate for any pair of edges the
// resulting sorted trisegment is itself indeterminate 
// (encoded as a collinear count of -1)
//
template<class K>
Uncertain<Trisegment_collinearity> certified_trisegment_collinearity ( Segment_2<K> const& e0
                                                                     , Segment_2<K> const& e1
                                                                     , Segment_2<K> const& e2
                                                                     )
{
  Uncertain<bool> is_01 = are_edges_orderly_collinearC2(e0,e1);
  if ( is_certain(is_01) )
  {
    Uncertain<bool> is_02 = are_edges_orderly_collinearC2(e0,e2);
    if ( is_certain(is_02) )
    {
      Uncertain<bool> is_12 = are_edges_orderly_collinearC2(e1,e2);
      if ( is_certain(is_12) )
      {
        if ( CGAL_NTS logical_and(is_01 , !is_02 , !is_12 ) )
          return TRISEGMENT_COLLINEARITY_01;
        else if ( CGAL_NTS logical_and(is_02 , !is_01 , !is_12 ) )
          return TRISEGMENT_COLLINEARITY_02;
        else if ( CGAL_NTS logical_and(is_12 , !is_01 , !is_02 ) )
          return TRISEGMENT_COLLINEARITY_12;
        else if ( CGAL_NTS logical_and(!is_01 , !is_02, !is_12  ) )
          return TRISEGMENT_COLLINEARITY_NONE;
        else 
          return TRISEGMENT_COLLINEARITY_ALL;
      }
    }
  }

  return Uncertain<Trisegment_collinearity>::indeterminate();
}



// Given 3 oriented straight line segments: e0, e1, e2 
// returns true if there exist some positive offset distance 't' for which the
// leftward-offsets of their supporting lines intersect at a single point.
//
// NOTE: This function can handle the case of collinear and/or parallel segments.
//
// If two segments are collinear but equally oriented (that is, they share a degenerate vertex) the event exists and
// is well defined, In that case, the degenerate vertex can be even a contour vertex or a skeleton node. If it is a skeleton
// node, it is properly defined by the trisegment tree that corresponds to the node.
// A trisegment tree stores not only the "current event" trisegment but also the trisegments for the left/right seed vertices,
// recursivey in case the seed vertices are skeleton nodes as well.
// Those seeds are used to determine the actual position of the degenerate vertex in case of collinear edges (since that point is
// not given by the collinear edges alone)
//
template<class K, class FT>
Uncertain<bool> exist_offset_lines_isec2 ( intrusive_ptr< Trisegment_2<K> > const& tri, optional<FT> const& aMaxTime )
{
  
  typedef Rational<FT>       Rational ;
  typedef optional<Rational> Optional_rational ;
  typedef Quotient<FT>       Quotient ;
  
  Uncertain<bool> rResult = Uncertain<bool>::indeterminate();

  if ( tri->collinearity() != TRISEGMENT_COLLINEARITY_ALL ) 
  {
    CGAL_STSKEL_TRAITS_TRACE( ( tri->collinearity() == TRISEGMENT_COLLINEARITY_NONE ? " normal edges" : " collinear edges" ) ) ;

    Optional_rational t = compute_offset_lines_isec_timeC2(tri) ;
    if ( t )
    {
      Uncertain<bool> d_is_zero = CGAL_NTS certified_is_zero(t->d()) ;
      if ( is_certain(d_is_zero) )
      {
        if ( !d_is_zero )
        {
          Quotient tq = t->to_quotient() ;
          
          rResult = CGAL_NTS certified_is_positive(tq) ;
          
          if ( aMaxTime && CGAL_NTS certainly(rResult) )
            rResult = CGAL_NTS certified_is_smaller_or_equal(tq,Quotient(*aMaxTime));
          
          CGAL_STSKEL_TRAITS_TRACE("\nEvent time: " << *t << ". Event " << ( rResult ? "exist." : "doesn't exist." ) ) ;
        }
        else
        {
          CGAL_STSKEL_TRAITS_TRACE("\nDenominator exactly zero, Event doesn't exist." ) ;
          rResult = false;
        }
      }
      else{
        CGAL_STSKEL_TRAITS_TRACE("\nDenominator is probably zero (but not exactly), event existance is indeterminate." ) ;
      }
    }
    else{
      CGAL_STSKEL_TRAITS_TRACE("\nEvent time overflowed, event existance is indeterminate." ) ;
    }
  }
  else
  {
    CGAL_STSKEL_TRAITS_TRACE("\nAll the edges are collinear. Event doesn't exist." ) ;
    rResult = false;
  }

  return rResult ;
}

// Given 2 triples of oriented straight line segments: (m0,m1,m2) and (n0,n1,n2), such that
// for each triple there exists distances 'mt' and 'nt' for which the offsets lines (at mt and nt resp.),
// (m0',m1',m2') and (n0',n1',n2') intersect each in a single point; returns the relative order of mt w.r.t nt.
// That is, indicates which offset triple intersects first (closer to the source lines)
// PRECONDITION: There exist distances mt and nt for which each offset triple intersect at a single point.
template<class K>
Uncertain<Comparison_result> compare_offset_lines_isec_timesC2 ( intrusive_ptr< Trisegment_2<K> > const& m
                                                               , intrusive_ptr< Trisegment_2<K> > const& n 
                                                               )
{
  typedef typename K::FT FT ;

  typedef Rational<FT>       Rational ;
  typedef Quotient<FT>       Quotient ;
  typedef optional<Rational> Optional_rational ;
  
  Uncertain<Comparison_result> rResult = Uncertain<Comparison_result>::indeterminate();

  Optional_rational mt_ = compute_offset_lines_isec_timeC2(m);
  Optional_rational nt_ = compute_offset_lines_isec_timeC2(n);
  
  if ( mt_ && nt_ )
  {
    Quotient mt = mt_->to_quotient();
    Quotient nt = nt_->to_quotient();
   
    if ( CGAL_NTS certified_is_positive(mt) && CGAL_NTS certified_is_positive(nt) ) 
      rResult = CGAL_NTS certified_compare(mt,nt);
  }
  
  return rResult ;

}


// Returns true if the point aP is on the positive side of the line supporting the edge
//
template<class K>
Uncertain<bool> is_edge_facing_pointC2 ( optional< Point_2<K> > const& aP, Segment_2<K> const& aEdge )
{
  typedef typename K::FT FT ;
  
  Uncertain<bool> rResult = Uncertain<bool>::indeterminate();
  if ( aP )
  {
    FT a,b,c ;
    line_from_pointsC2(aEdge.source().x(),aEdge.source().y(),aEdge.target().x(),aEdge.target().y(),a,b,c); 
    rResult = certified_side_of_oriented_lineC2(a,b,c,aP->x(),aP->y()) == make_uncertain(POSITIVE) ;
  }
  return rResult ;
}

// Given a triple of oriented straight line segments: (e0,e1,e2) such that their offsets
// at some distance intersects in a point (x,y), returns true if (x,y) is on the positive side of the line supporting aEdge
//
template<class K>
inline Uncertain<bool> is_edge_facing_offset_lines_isecC2 ( intrusive_ptr< Trisegment_2<K> > const& tri, Segment_2<K> const& aEdge )
{
  return is_edge_facing_pointC2(construct_offset_lines_isecC2(tri),aEdge);
}

// Given an event trisegment and two oriented straight line segments e0 and e1, returns the oriented side of the event point 
// w.r.t the (positive) bisector [e0,e1].
//
// The (positive) bisector [e0,e1] is a ray starting at the vertex (e0,e1) (called "v01")
// 
// If e0,e1 are consecutive in the input polygon, v01 = e0.target() = e1.source().
// If they are not consecutive in the input polygon they must had become consecutive at some known previous event. In this
// case, the point of the "v01_event" trisegment intersection determines v01 which is to the position of the very first
// vertex shared between e0,e1 (at the time they first become consecutive).
//
// A point *exactly on* the bisector is an offset vertex (e0*,e1*) (that is, belongs to both e0* and e1*). 
// A point to the positive side of the bisector belongs to e0* but not e1*
// A point to the negative side of the bisector belongs to e1* but not e0*
//
// If e0,e1 intersect in a single point the bisector is an angular bisector.
//
// One of e0 or e1 is considered the primary edge.
//
// If e0,e1 are parallel a line perpendicular to the primary edge passing through "v01" is used "as bisector". 
//
// If e0,e1 are collinear then this perpendicular line is a perpendicular bisector of the two segments.
// 
// If e0,e1 are parallel but not collinear then the actual bisector is an equidistant line parallel to e0 and e1.
// e0* and e1* overlap and are known to be connected sharing a known vertex v01, which is somewhere along the parallel
// line which is the bisector of e0 and e1. 
// Given a line perpendicular to e0 through v01, a point to its positive side belongs to e0* while a point to its negative side does not.
// Given a line perpendicular to e1 through v01, a point to its negative side belongs to e1* while a point to its positive side does not.
//
// This predicate is used to determine the validity of a split or edge event.
// 
// A split event is the coallision of a reflex wavefront and some opposite offset egde. Unless the three segments
// don't actually collide (there is no event), the split point is along the supporting line of the offset edge.
// Testing its validity amounts to determining if the split point is inside the closed offset segment instead of 
// the two open rays before and after the offset segment endpoints.
// The offset edge is bounded by its previous and next adjacent edges at the time of the event. Thus, the bisectors
// of this edge and its previous/next adjacent edges (at the time of the event) detemine the offset vertices that
// bound the opposite edge.
// If the opposite edge is 'e' and its previous/next edges are "preve"/"nexte" then the split point is inside the offset
// egde if it is NOT to the positive side of [preve,e] *and* NOT to the negative side o [e,nexte].
// (so this predicate answer half the question, at one and other side independenty).
// If the split point is exacty over any of this bisectors then the split point ocurres exactly and one (or both) endpoints
// of the opposite edge (so it is a pseudo-split event since the opposite edge is not itself split in two halfeves)
// When this predicate is called to test (prev,e), e is the primary edge but since it is  pass as e1, primary_is_0=false. 
// This causes the case of parallel but not collinear edges to return positive when the split point is before the source point of e*
// (a positive result means invalid split).
// Likewise, primary_is_0 must be true when testing (e,nexte) to return negative if the split point is past the target endpoint of e*.
// (in the other cases there is no need to discrminate which is 'e' in the call since the edjes do not overlap).
//
// An edge event is a coallision of three *consecutive* edges, say, e1,e2 and e3.
// The coallision causes e2 (the edge in the middle) to collapse and e1,e3 to become consecutive and form a new vertex.
// In all cases there is an edge before e1, say e0, and after e3, say e4.
// Testing for the validity of an edge event amounts to determine that (e1,e3) (the new vertex) is not before (e0,e1) nor
// past (e3,e4). 
// Thus, and edge event is valid if the new vertex NOT to the positive side of [e0,e1] *and* NOT to the negative side o [e3,e4].
//
// PRECONDITIONS:
//   There exist a single point 'p' corresponding to the event as given by the trisegment
//   e0 and e1 are known to be consectuve at the time of the event (even if they are not consecutive in the input polygon)
//   If e0 and e1 are not consecutive in the input, v01_event is the event that defined they very first offset vertex.
//   If e0 and e1 are consecutive, v01_event is null.
//
template<class K>
Uncertain<Oriented_side>
oriented_side_of_event_point_wrt_bisectorC2 ( intrusive_ptr< Trisegment_2<K> > const& event
                                            , Segment_2<K>                     const& e0
                                            , Segment_2<K>                     const& e1
                                            , intrusive_ptr< Trisegment_2<K> > const& v01_event // can be null
                                            , bool                                    primary_is_0
                                            )
{
  typedef typename K::FT FT ;
  
  typedef Point_2     <K> Point_2 ;
  typedef Line_2      <K> Line_2 ;
  
  Uncertain<Oriented_side> rResult = Uncertain<Oriented_side>::indeterminate();
  
  try
  {
    Point_2 p = validate(construct_offset_lines_isecC2(event));
    
    Line_2 l0 = validate(compute_normalized_line_ceoffC2(e0)) ;
    Line_2 l1 = validate(compute_normalized_line_ceoffC2(e1)) ;
  
    CGAL_STSKEL_TRAITS_TRACE("Getting oriented side of point " << p2str(p)
                            << " w.r.t bisector ["
                            << s2str(e0) << ( primary_is_0 ? "*" : "" )
                            << "," 
                            << s2str(e1) << ( primary_is_0 ? "" : "*" )
                            << "]"
                            ) ; 
      
    // Degenerate bisector?   
    if ( certainly( are_edges_parallelC2(e0,e1) ) )
    { 
      CGAL_STSKEL_TRAITS_TRACE("Bisector is not angular." ) ;
      
      // b01 is degenerate so we don't have an *angular bisector* but a *perpendicular* bisector.
      // We need to compute the actual bisector line.
      CGAL_assertion( v01_event || ( !v01_event && e0.target() == e1.source() ) ) ;
      
      Point_2 v01 = v01_event ? validate( construct_offset_lines_isecC2(v01_event) ) : e1.source() ;
      
      CGAL_STSKEL_TRAITS_TRACE("v01=" << p2str(v01) << ( v01_event ? " (from skelton node)" : "" ) ) ;      
      
      // (a,b,c) is a line perpedincular to the primary edge through v01.
      // If e0 and e1 are collinear this line is the actual perpendicular bisector.
      //
      // If e0 and e1 are parallel but not collinear (then neccesarrily facing each other) this line
      // is NOT the bisector, but the serves to determine the side of the point (projected along the primary ege) w.r.t vertex v01.
      
      FT a, b, c ;
      perpendicular_through_pointC2( primary_is_0 ? l0.a() : l1.a()
                                   , primary_is_0 ? l0.b() : l1.b()
                                   , v01.x(), v01.y()
                                   , a, b, c
                                   );
      
      rResult = certified_side_of_oriented_lineC2(a,b,c,p.x(),p.y());
      
      CGAL_STSKEL_TRAITS_TRACE("Point is at " << rResult << " side of degenerate bisector through v01 " << p2str(v01)) ;
    }
    else // Valid (non-degenerate) angular bisector
    {
      // Scale distance from to the lines.
      FT sd_p_l0 = validate(l0.a() * p.x() + l0.b() * p.y() + l0.c()) ;
      FT sd_p_l1 = validate(l1.a() * p.x() + l1.b() * p.y() + l1.c()) ;
      
      CGAL_STSKEL_TRAITS_TRACE("sd_p_l1=" << n2str(sd_p_l1) ) ;
      CGAL_STSKEL_TRAITS_TRACE("sd_p_l0=" << n2str(sd_p_l0) ) ;
        
      Uncertain<bool> equal = CGAL_NTS certified_is_equal(sd_p_l0,sd_p_l1) ;   
      if ( is_certain(equal) )
      {
        if ( equal )
        {
          CGAL_STSKEL_TRAITS_TRACE("Point is exactly at bisector"); 
          
          rResult = ON_ORIENTED_BOUNDARY;
        }
        else
        {
          Uncertain<bool> smaller = CGAL_NTS certified_is_smaller( validate(l0.a()*l1.b()), validate(l1.a()*l0.b()) ) ;
          if ( is_certain(smaller) )
          {
            // Reflex bisector?
            if ( smaller )
            {
              rResult = CGAL_NTS certified_is_smaller(sd_p_l0,sd_p_l1) ? ON_NEGATIVE_SIDE 
                                                                       : ON_POSITIVE_SIDE ;
                                                                      
              CGAL_STSKEL_TRAITS_TRACE("\nEvent point is at " << rResult << " side of reflex bisector" ) ;
            }
            else
            {
              rResult = CGAL_NTS certified_is_larger (sd_p_l0,sd_p_l1) ? ON_NEGATIVE_SIDE
                                                                       : ON_POSITIVE_SIDE ; 
                              
              CGAL_STSKEL_TRAITS_TRACE("\nEvent point is at " << rResult << " side of convex bisector" ) ;
            }
          }
        }
      }
    }
  }
  catch( std::overflow_error const& )
  {
    CGAL_STSKEL_TRAITS_TRACE("Unable to compute value due to overflow.");
  }
  catch( CGAL::Uncertain_conversion_exception const& )
  {
    CGAL_STSKEL_TRAITS_TRACE("Indeterminate boolean expression.");
  }
  
  return rResult ;

}


// Given 2 triples of oriented straight line segments (l0,l1,l2) and (r0,r1,r2), such that 
// the offsets at time 'tl' for triple 'l' intersects in a point (lx,ly) and 
// the offsets at time 'tr' for triple 'r' intersects in a point (rx,ry) 
// returns true if "tl==tr" and "(lx,ly)==(rx,ry)" 
// PRECONDITIONS:
//   There exist single points at which the offset lines for 'l' and 'r' at 'tl', 'tr' intersect.
//
template<class K>
Uncertain<bool> are_events_simultaneousC2 ( intrusive_ptr< Trisegment_2<K> > const& l, intrusive_ptr< Trisegment_2<K> > const& r )
{
  typedef typename K::FT FT ;
  
  typedef Point_2<K> Point_2 ;
  
  typedef Rational<FT> Rational ;
  typedef Quotient<FT> Quotient ;
  
  typedef optional<Rational> Optional_rational ;
  typedef optional<Point_2>  Optional_point_2 ;
  
  Uncertain<bool> rResult = Uncertain<bool>::indeterminate();

  Optional_rational lt_ = compute_offset_lines_isec_timeC2(l);
  Optional_rational rt_ = compute_offset_lines_isec_timeC2(r);
  
  if ( lt_ && rt_ )
  {
    Quotient lt = lt_->to_quotient();
    Quotient rt = rt_->to_quotient();

    if ( CGAL_NTS certified_is_positive(lt) && CGAL_NTS certified_is_positive(rt) ) 
    {
      Uncertain<bool> equal_times = CGAL_NTS certified_is_equal(lt,rt);
      
      if ( is_certain(equal_times) )
      {
        if ( equal_times )
        {
          Optional_point_2 li = construct_offset_lines_isecC2(l);
          Optional_point_2 ri = construct_offset_lines_isecC2(r);
  
          if ( li && ri )
            rResult = CGAL_NTS logical_and( CGAL_NTS certified_is_equal(li->x(),ri->x())
                                          , CGAL_NTS certified_is_equal(li->y(),ri->y())
                                          ) ;
        }
        else rResult = false;
      }
    }

  }
  return rResult;
}

} // namespace CGAL_SS_i

} // end namespace CGAL

#endif // CGAL_STRAIGHT_SKELETON_PREDICATES_FTC2_H //
