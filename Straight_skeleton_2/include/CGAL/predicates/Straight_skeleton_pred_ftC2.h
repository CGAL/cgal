// Copyright (c) 2006 Fernando Luis Cacciola Carballal. All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Fernando Cacciola <fernando_cacciola@ciudad.com.ar>
//
#ifndef CGAL_STRAIGHT_SKELETON_PREDICATES_FTC2_H
#define CGAL_STRAIGHT_SKELETON_PREDICATES_FTC2_H 1

#include <CGAL/license/Straight_skeleton_2.h>

#include <CGAL/certified_numeric_predicates.h>
#include <CGAL/certified_quotient_predicates.h>
#include <CGAL/Straight_skeleton_2/Straight_skeleton_aux.h>
#include <CGAL/Straight_skeleton_2/Straight_skeleton_builder_traits_2_aux.h>
#include <CGAL/constructions/Straight_skeleton_cons_ftC2.h>

#include <CGAL/Point_2.h>
#include <CGAL/Quotient.h>
#include <CGAL/Uncertain.h>

#include <optional>

#include <stdexcept>

namespace CGAL {

namespace CGAL_SS_i {

// Just like the uncertified collinear() returns true IFF r lies in the line p->q
// NOTE: r might be in the ray from p or q containing q or p, that is, there is no ordering implied, just that
// the three points are along the same line, in any order.
template<class K>
inline
Uncertain<bool> certified_collinearC2( typename K::Point_2 const& p,
                                       typename K::Point_2 const& q,
                                       typename K::Point_2 const& r )
{
  return CGAL_NTS certified_is_equal( ( q.x() - p.x() ) * ( r.y() - p.y() )
                                    , ( r.x() - p.x() ) * ( q.y() - p.y() )
                                    );
}


// Just like the uncertified collinear_are_ordered_along_lineC2() returns true IFF, given p,q,r along the same line,
// q is in the closed segment [p,r].
template<class K>
inline
Uncertain<bool> certified_collinear_are_ordered_along_lineC2( typename K::Point_2 const& p,
                                                              typename K::Point_2 const& q,
                                                              typename K::Point_2 const& r )
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
Uncertain<bool> are_edges_collinearC2( Segment_2_with_ID<K> const& e0, Segment_2_with_ID<K> const& e1 )
{
  return   certified_collinearC2(e0.source(),e0.target(),e1.source())
         & certified_collinearC2(e0.source(),e0.target(),e1.target()) ;
}

// Returns true IFF the supporting lines for segments e0,e1 are parallel (or the same)
template<class K>
inline
Uncertain<bool> are_edges_parallelC2( Segment_2_with_ID<K> const& e0, Segment_2_with_ID<K> const& e1 )
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
Uncertain<bool> are_parallel_edges_equally_orientedC2( Segment_2_with_ID<K> const& e0, Segment_2_with_ID<K> const& e1 )
{
  return CGAL_NTS certified_sign( (e0.target() - e0.source()) * (e1.target() - e1.source()) ) == POSITIVE;
}


// Returns true IFF segments e0,e1 share the same supporting line but do not overlap except at the vertices, and have the same orientation.
// NOTE: If e1 goes back over e0 (a degenerate antenna or alley) this returns false.
template<class K>
Uncertain<bool> are_edges_orderly_collinearC2( Segment_2_with_ID<K> const& e0, Segment_2_with_ID<K> const& e1 )
{
  return are_edges_collinearC2(e0,e1) & are_parallel_edges_equally_orientedC2(e0,e1);
}

template <class FT>
inline
Uncertain<Sign> certified_side_of_oriented_lineC2(const FT &a, const FT &b, const FT &c, const FT &x, const FT &y)
{
  return CGAL_NTS certified_sign(a*x+b*y+c);
}

// Given 3 oriented straight line segments: e0, e1, e2
// returns true if there exists some positive offset distance 't' for which the
// leftward-offsets of their supporting lines intersect at a single point.
//
// NOTE: This function can handle the case of collinear and/or parallel segments.
//
// If two segments are collinear but equally oriented (that is, they share a degenerate vertex) the event exists and
// is well defined, In that case, the degenerate vertex can be even a contour vertex or a skeleton node. If it is a skeleton
// node, it is properly defined by the trisegment tree that corresponds to the node.
// A trisegment tree stores not only the "current event" trisegment but also the trisegments for the left/right seed vertices,
// recursively in case the seed vertices are skeleton nodes as well.
// Those seeds are used to determine the actual position of the degenerate vertex in case of collinear edges (since that point is
// not given by the collinear edges alone)
//
template<class K, class FT, class Caches>
Uncertain<bool> exist_offset_lines_isec2 ( Trisegment_2_ptr< Trisegment_2<K, Segment_2_with_ID<K> > > const& tri,
                                           std::optional<FT> const& aMaxTime,
                                           Caches& aCaches )
{

  typedef Rational<FT>              Rational ;
  typedef std::optional<Rational> Optional_rational ;
  typedef Quotient<FT>              Quotient ;

  Uncertain<bool> rResult = Uncertain<bool>::indeterminate();

  CGAL_STSKEL_TRAITS_TRACE( "\n~~ Checking existence of an event [" << typeid(FT).name() << "]");
  CGAL_STSKEL_TRAITS_TRACE("Event:\n" << tri);

  if ( tri->collinearity() != TRISEGMENT_COLLINEARITY_ALL )
  {
    CGAL_STSKEL_TRAITS_TRACE( ( tri->collinearity() == TRISEGMENT_COLLINEARITY_NONE ? " normal edges" : " collinear edges" ) ) ;

    Optional_rational t = compute_offset_lines_isec_timeC2(tri, aCaches) ;
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

          CGAL_STSKEL_TRAITS_TRACE("Event time: " << *t << ". Event " << ( rResult ? "exist." : "doesn't exist." ) ) ;
        }
        else
        {
          CGAL_STSKEL_TRAITS_TRACE("Denominator exactly zero, Event doesn't exist." ) ;
          rResult = false;
        }
      }
      else
      {
        CGAL_STSKEL_TRAITS_TRACE("Denominator is probably zero (but not exactly), event existence is indeterminate." ) ;
      }
    }
    else
    {
      CGAL_STSKEL_TRAITS_TRACE("Event time overflowed, event existence is indeterminate." ) ;
    }
  }
  else
  {
    CGAL_STSKEL_TRAITS_TRACE("All the edges are collinear. Event doesn't exist." ) ;
    rResult = false;
  }

  return rResult ;
}

// Given 2 triples of oriented straight line segments: (m0,m1,m2) and (n0,n1,n2), such that
// for each triple there exists distances 'mt' and 'nt' for which the offsets lines (at mt and nt resp.),
// (m0',m1',m2') and (n0',n1',n2') intersect each in a single point; returns the relative order of mt w.r.t. nt.
// That is, indicates which offset triple intersects first (closer to the source lines)
// PRECONDITION: There exists distances mt and nt for which each offset triple intersect at a single point.
template<class K, class Caches>
Uncertain<Comparison_result>
compare_offset_lines_isec_timesC2 ( Trisegment_2_ptr< Trisegment_2<K, Segment_2_with_ID<K> > > const& m,
                                    Trisegment_2_ptr< Trisegment_2<K, Segment_2_with_ID<K> > > const& n,
                                    Caches& aCaches )
{
  typedef typename K::FT FT ;

  typedef Rational<FT>       Rational ;
  typedef Quotient<FT>       Quotient ;
  typedef std::optional<Rational> Optional_rational ;

  CGAL_STSKEL_TRAITS_TRACE("compare_offset_lines_isec_timesC2(\n" << m << "\n" << n << "\n) [" << typeid(FT).name() << "]" );

  Uncertain<Comparison_result> rResult = Uncertain<Comparison_result>::indeterminate();

  Optional_rational mt_ = compute_offset_lines_isec_timeC2(m, aCaches);
  Optional_rational nt_ = compute_offset_lines_isec_timeC2(n, aCaches);

  if ( mt_ && nt_ )
  {
    Quotient mt = mt_->to_quotient();
    Quotient nt = nt_->to_quotient();

    if ( CGAL_NTS certified_is_positive(mt) && CGAL_NTS certified_is_positive(nt) )
      rResult = CGAL_NTS certified_compare(mt,nt);
  }

  return rResult ;
}

template<class K>
Uncertain<Comparison_result> compare_isec_anglesC2 ( Vector_2<K> const& aBV1
                                                   , Vector_2<K> const& aBV2
                                                   , Vector_2<K> aLV
                                                   , Vector_2<K> aRV
                                                   )
{
  typedef typename K::FT FT ;
  typedef typename K::Vector_2 Vector_2 ;

  Uncertain<Comparison_result> rResult = Uncertain<Comparison_result>::indeterminate();

  const Vector_2 lBisectorDirection = aBV2 - aBV1 ;
  const FT lLNorm = CGAL_SS_i::inexact_sqrt ( K().compute_scalar_product_2_object()( aLV, aLV ) ) ;
  const FT lRNorm = CGAL_SS_i::inexact_sqrt ( K().compute_scalar_product_2_object()( aRV, aRV ) ) ;

  if (! CGAL_NTS certified_is_positive( lLNorm ) ||
      ! CGAL_NTS certified_is_positive( lRNorm ) )
    return rResult ;

  aLV = aLV / lLNorm ;
  aRV = aRV / lRNorm ;

  const FT lLSp = K().compute_scalar_product_2_object()( lBisectorDirection, aLV ) ;
  const FT lRSp = K().compute_scalar_product_2_object()( lBisectorDirection, aRV ) ;

  // Smaller if the scalar product is larger, so swapping
  rResult = CGAL_NTS certified_compare(lRSp, lLSp) ;

  return rResult;
}

// Returns true if the point aP is on the positive side of the line supporting the edge
//
template<class K>
Uncertain<bool> is_edge_facing_pointC2 ( std::optional< typename K::Point_2 > const& aP,
                                         Segment_2_with_ID<K> const& aEdge )
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
template<class K, class Caches>
inline Uncertain<bool>
is_edge_facing_offset_lines_isecC2 ( Trisegment_2_ptr< Trisegment_2<K, Segment_2_with_ID<K> > > const& tri,
                                     Segment_2_with_ID<K> const& aEdge,
                                     Caches& aCaches )
{
  return is_edge_facing_pointC2(construct_offset_lines_isecC2(tri, aCaches), aEdge);
}

// Given an event trisegment and two oriented straight line segments e0 and e1, returns the oriented side of the event point
// w.r.t. the (positive) bisector [e0,e1].
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
// If e0,e1 are parallel but in opposite directions then the bisector is an equidistant line parallel to e0 and e1.
// e0* and e1* overlap and are known to be connected sharing a known vertex v01, which is somewhere along the parallel
// line which is the bisector of e0 and e1.
// Given a line perpendicular to e0 through v01, a point to its positive side belongs to e0* while a point to its negative side does not.
// Given a line perpendicular to e1 through v01, a point to its negative side belongs to e1* while a point to its positive side does not.
//
// This predicate is used to determine the validity of a split or edge event.
//
// A split event is the collision of a reflex wavefront and some opposite offset edge. Unless the three segments
// don't actually collide (there is no event), the split point is along the supporting line of the offset edge.
// Testing its validity amounts to determining if the split point is inside the closed offset segment instead of
// the two open rays before and after the offset segment endpoints.
// The offset edge is bounded by its previous and next adjacent edges at the time of the event. Thus, the bisectors
// of this edge and its previous/next adjacent edges (at the time of the event) determine the offset vertices that
// bound the opposite edge.
// If the opposite edge is 'e' and its previous/next edges are "preve"/"nexte" then the split point is inside the offset
// edge if it is NOT to the positive side of [preve,e] *and* NOT to the negative side o [e,nexte].
// (so this predicate answers half the question, at one side).
// If the split point is exactly over any of these bisectors then the split point occurs exactly and one (or both) endpoints
// of the opposite edge (so it is a pseudo-split event since the opposite edge is not itself split in two halfedges).
// When this predicate is called to test (prev,e), e is the primary edge but since it is passed as e1, primary_is_0=false.
// This causes the case of parallel but not collinear edges to return positive when the split point is before the source point of e*
// (a positive result means invalid split).
// Likewise, primary_is_0 must be true when testing (e,nexte) to return negative if the split point is past the target endpoint of e*.
// (in the other cases there is no need to discriminate which is 'e' in the call since the edges do not overlap).
//
// An edge event is a collision of three *consecutive* edges, say, e1,e2 and e3.
// The collision causes e2 (the edge in the middle) to collapse and e1,e3 to become consecutive and form a new vertex.
// In all cases there is an edge before e1, say e0, and after e3, say e4.
// Testing for the validity of an edge event amounts to determine that (e1,e3) (the new vertex) is not before (e0,e1) nor
// past (e3,e4).
// Thus, and edge event is valid if the new vertex is NOT to the positive side of [e0,e1] *and* NOT to the negative side of [e3,e4].
//
// PRECONDITIONS:
//   There exists a single point 'p' corresponding to the event as given by the trisegment
//   e0 and e1 are known to be consecutive at the time of the event (even if they are not consecutive in the input polygon)
//   If e0 and e1 are not consecutive in the input, v01_event is the event that defined their first common offset vertex.
//   If e0 and e1 are consecutive in the input, v01_event is null.
//
template<class K, class Caches>
Uncertain<Oriented_side>
oriented_side_of_event_point_wrt_bisectorC2 ( Trisegment_2_ptr< Trisegment_2<K, Segment_2_with_ID<K> > > const& event,
                                              Segment_2_with_ID<K> const& e0,
                                              typename K::FT const& w0,
                                              Segment_2_with_ID<K> const& e1,
                                              typename K::FT const& w1,
                                              Trisegment_2_ptr< Trisegment_2<K, Segment_2_with_ID<K> > > const& v01_event, // can be null
                                              bool primary_is_0,
                                              Caches& aCaches )
{
  typedef typename K::FT FT ;

  typedef typename K::Point_2 Point_2 ;
  typedef typename K::Line_2  Line_2 ;

  Uncertain<Oriented_side> rResult = Uncertain<Oriented_side>::indeterminate();

  try
  {
    Point_2 p = validate(construct_offset_lines_isecC2(event, aCaches));

    Line_2 l0 = validate(compute_weighted_line_coeffC2(e0, w0, aCaches));
    Line_2 l1 = validate(compute_weighted_line_coeffC2(e1, w1, aCaches));

    CGAL_STSKEL_TRAITS_TRACE("\n~~ Oriented side of point [" << typeid(FT).name() << "]" );
    CGAL_STSKEL_TRAITS_TRACE("p = " << p2str(p)
                            << " w.r.t. bisector of [E" << e0.id() << " "
                            << s2str(e0) << ( primary_is_0 ? "*" : "" )
                            << ", E" << e1.id() << " "
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

      Point_2 v01 = v01_event ? validate(construct_offset_lines_isecC2(v01_event, aCaches))
                              : e1.source() ;

      CGAL_STSKEL_TRAITS_TRACE("v01=" << p2str(v01) << ( v01_event ? " (from skeleton node)" : "" ) ) ;

      // (a,b,c) is a line perpendicular to the primary edge through v01.
      // If e0 and e1 are collinear this line is the actual perpendicular bisector.
      //
      // If e0 and e1 are parallel but not collinear (then necessarily facing each other) this line
      // is NOT the bisector, but it serves to determine the side of the point (projected along
      // the primary edge) w.r.t. vertex v01.

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

      CGAL_STSKEL_TRAITS_TRACE("sd_p_l0 = " << n2str(sd_p_l0) ) ;
      CGAL_STSKEL_TRAITS_TRACE("sd_p_l1 = " << n2str(sd_p_l1) ) ;

      Uncertain<Comparison_result> lCmpResult = CGAL_NTS certified_compare(sd_p_l0,sd_p_l1) ;
      if ( is_certain(lCmpResult) )
      {
        CGAL_STSKEL_TRAITS_TRACE("compare(sd_p_l0, sd_p_l1) = " << lCmpResult ) ;
        if ( lCmpResult == EQUAL )
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
              rResult = (lCmpResult == SMALLER) ? ON_NEGATIVE_SIDE : ON_POSITIVE_SIDE ;
              CGAL_STSKEL_TRAITS_TRACE("Event point is on " << ((rResult > 0) ? "POSITIVE" : "NEGATIVE") << " side of reflex bisector" ) ;
            }
            else
            {
              rResult = (lCmpResult == LARGER) ? ON_NEGATIVE_SIDE : ON_POSITIVE_SIDE ;
              CGAL_STSKEL_TRAITS_TRACE("Event point is on " << ((rResult > 0) ? "POSITIVE" : "NEGATIVE") << " side of convex bisector" ) ;
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
//   There exists single points at which the offset lines for 'l' and 'r' at 'tl', 'tr' intersect.
//
template<class K, class Caches>
Uncertain<bool> are_events_simultaneousC2 ( Trisegment_2_ptr< Trisegment_2<K, Segment_2_with_ID<K> > > const& l,
                                            Trisegment_2_ptr< Trisegment_2<K, Segment_2_with_ID<K> > > const& r,
                                            Caches& aCaches )
{
  typedef typename K::FT FT ;

  typedef typename K::Point_2 Point_2 ;

  typedef Rational<FT> Rational ;
  typedef Quotient<FT> Quotient ;

  typedef std::optional<Rational> Optional_rational ;
  typedef std::optional<Point_2>  Optional_point_2 ;

  Uncertain<bool> rResult = Uncertain<bool>::indeterminate();

  Optional_rational lt_ = compute_offset_lines_isec_timeC2(l, aCaches);
  Optional_rational rt_ = compute_offset_lines_isec_timeC2(r, aCaches);

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
          Optional_point_2 li = construct_offset_lines_isecC2(l, aCaches);
          Optional_point_2 ri = construct_offset_lines_isecC2(r, aCaches);

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
