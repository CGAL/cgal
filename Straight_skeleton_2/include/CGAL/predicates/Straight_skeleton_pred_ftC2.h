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
#ifndef CGAL_STRAIGHT_SKELETON_PREDICATES_FTC2_H 
#define CGAL_STRAIGHT_SKELETON_PREDICATES_FTC2_H 1

#include <CGAL/constructions/Straight_skeleton_cons_ftC2.h>
#include <CGAL/Uncertain.h>
#include <CGAL/certified_quotient_predicates.h>

CGAL_BEGIN_NAMESPACE


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
// q in the closed segment [p,r].
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
  
  if ( CGAL_NTS certainly(p.x() == q.x()) && CGAL_NTS certainly(p.y() == q.y())  ) return make_uncertain(true);
  
  return Uncertain<bool>::indeterminate();
}

// Returns true IFF segments e0,e1 share the same supporting line, do not overlap except at the vetices, and have the same orientation.
// NOTE: If e1 goes back over e0 (a degenerate antenna or alley) this returns false.
template<class K>
Uncertain<bool> are_edges_orderly_collinearC2( Segment_2<K> const& e0, Segment_2<K> const& e1 )
{
  return    certified_collinearC2(e0.source(),e0.target(),e1.source())
          &
            (  certified_collinear_are_ordered_along_lineC2(e0.source(),e0.target(),e1.source())
             | certified_collinear_are_ordered_along_lineC2(e1.source(),e0.source(),e0.target())
            )
          & certified_collinearC2(e0.source(),e0.target(),e1.target()) 
          & 
            (  certified_collinear_are_ordered_along_lineC2(e0.source(),e0.target(),e1.target()) 
             | certified_collinear_are_ordered_along_lineC2(e1.target(),e0.source(),e0.target()) 
            );
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
  Uncertain<Trisegment_collinearity> rR = Uncertain<Trisegment_collinearity>::indeterminate();
  
  Uncertain<bool> is_01 = are_edges_orderly_collinearC2(e0,e1);
  if ( !CGAL_NTS is_indeterminate(is_01) )
  {
    Uncertain<bool> is_02 = are_edges_orderly_collinearC2(e0,e2);
    if ( !CGAL_NTS is_indeterminate(is_02) )
    {
      Uncertain<bool> is_12 = are_edges_orderly_collinearC2(e1,e2);
      if ( !CGAL_NTS is_indeterminate(is_12) )
      {
        if ( CGAL_NTS logical_and(is_01 , !is_02 , !is_12 ) )
          rR = make_uncertain(TRISEGMENT_COLLINEARITY_01);
        else if ( CGAL_NTS logical_and(is_02 , !is_01 , !is_12 ) )
          rR = make_uncertain(TRISEGMENT_COLLINEARITY_02);
        else if ( CGAL_NTS logical_and(is_12 , !is_01 , !is_02 ) )
          rR = make_uncertain(TRISEGMENT_COLLINEARITY_12);
        else if ( CGAL_NTS logical_and(!is_01 , !is_02, !is_12  ) )
          rR = make_uncertain(TRISEGMENT_COLLINEARITY_NONE);
        else 
          rR = make_uncertain(TRISEGMENT_COLLINEARITY_ALL);
      }
    }
  }

  return rR ;
}



// Given 3 oriented straight line segments: e0, e1, e2 
// returns true if there exist some positive offset distance 't' for which the
// leftward-offsets of their supporting lines intersect at a single point.
//
// NOTE: This function can handle the case of collinear and/or parallel segments.
//
// If two segments are collinear but equally oriented (that is, they share a degenerate vertex) the event exists and
// is well defined, In that case, the degenerate vertex can be even a contour vertex or a skeleton node. If it is a skeleton
// node, it is properly defined by the event trisegment that corresponds to the node.
// A Seeded_trisegment stores not only the "current event" trisegment but also the trisegments for the left/right seed vertices.
// Those seeds are used to determine the actual position of the degenerate vertex in case of collinear edges (since that point is
// not given by the collinear edges alone)
//
template<class K>
Uncertain<bool> exist_offset_lines_isec2 ( Seeded_trisegment_2<K> const& st )
{
  typedef typename K::FT FT ;
  
  typedef Rational<FT>       Rational ;
  typedef optional<Rational> Optional_rational ;
  
  Uncertain<bool> rResult = Uncertain<bool>::indeterminate();

  if ( st.event().collinearity() != TRISEGMENT_COLLINEARITY_ALL ) 
  {
    CGAL_STSKEL_TRAITS_TRACE( ( st.event().collinearity() == TRISEGMENT_COLLINEARITY_NONE ? " normal edges" : " collinear edges" ) ) ;

    Optional_rational t = compute_offset_lines_isec_timeC2(st) ;
    if ( t )
    {
      Uncertain<bool> d_is_zero = CGAL_NTS certified_is_zero(t->d()) ;
      if ( ! CGAL_NTS is_indeterminate(d_is_zero) )
      {
        if ( !d_is_zero )
        {
          rResult = CGAL_NTS certified_is_positive(t->to_quotient()) ;
          CGAL_STSKEL_TRAITS_TRACE("\nEvent time: " << (t->n()/t->d()) << ". Event " << ( rResult ? "exist." : "doesn't exist." ) ) ;
        }
        else
        {
          CGAL_STSKEL_TRAITS_TRACE("\nDenominator exactly zero, Event doesn't exist." ) ;
          rResult = make_uncertain(false);
        }
      }
      else
        CGAL_STSKEL_TRAITS_TRACE("\nDenominator is probably zero (but not exactly), event existance is indeterminate." ) ;
    }
    else
      CGAL_STSKEL_TRAITS_TRACE("\nEvent time overflowed, event existance is indeterminate." ) ;
  }
  else
  {
    CGAL_STSKEL_TRAITS_TRACE("\nAll the edges are collinear. Event doesn't exist." ) ;
    rResult = make_uncertain(false);
  }

  return rResult ;
}

// Given 2 triples of oriented straight line segments: (m0,m1,m2) and (n0,n1,n2), such that
// for each triple there exists distances 'mt' and 'nt' for which the offsets lines (at mt and nt resp.),
// (m0',m1',m2') and (n0',n1',n2') intersect each in a single point; returns the relative order of mt w.r.t nt.
// That is, indicates which offset triple intersects first (closer to the source lines)
// PRECONDITION: There exist distances mt and nt for which each offset triple intersect at a single point.
template<class K>
Uncertain<Comparison_result> compare_offset_lines_isec_timesC2 ( Seeded_trisegment_2<K> const& m
                                                               , Seeded_trisegment_2<K> const& n 
                                                               )
{
  typedef typename K::FT FT ;
  
  typedef Trisegment_2<K> Trisegment_2 ;
  
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
   
    CGAL_assertion ( CGAL_NTS certified_is_positive(mt) ) ;
    CGAL_assertion ( CGAL_NTS certified_is_positive(nt) ) ;

    rResult = CGAL_NTS certified_compare(mt,nt);
  }
  
  return rResult ;

}

// Given a point (px,py) and 2 triples of oriented straight line segments: (m0,m1,m2) and (n0,n1,n2),
// such that their offsets at distances 'mt' and 'nt' intersects in points (mx,my) and (nx,ny),
// returns the relative order of the distances from (px,py) to (mx,my) and (nx,ny).
// That is, indicates which offset triple intersects closer to (px,py)
// PRECONDITION: There exist single points at which the offset line triples 'm' and 'n' at 'mt' and 'nt' intersect.
template<class K>
Uncertain<Comparison_result>
compare_offset_lines_isec_dist_to_pointC2 ( optional< Point_2<K> > const& p
                                          , Seeded_trisegment_2<K> const& m
                                          , Seeded_trisegment_2<K> const& n 
                                          )
{
  typedef typename K::FT FT ;
  
  typedef Trisegment_2<K> Trisegment_2 ;
  
  typedef Rational<FT>       Rational ;
  typedef Quotient<FT>       Quotient ;
  typedef optional<Rational> Optional_rational ;
  
  Uncertain<Comparison_result> rResult = Uncertain<Comparison_result>::indeterminate();
  
  if ( p )
  {
    optional<FT> dm = compute_offset_lines_isec_dist_to_pointC2(p,m);
    optional<FT> dn = compute_offset_lines_isec_dist_to_pointC2(p,n);

    if ( dm && dn )
      rResult = CGAL_NTS certified_compare(*dm,*dn);
  }

  return rResult ;
}

// Given 3 triples of oriented straight line segments: (s0,s1,s2), (m0,m1,m2) and (n0,n1,n2),
// such that their offsets at distances 'st', 'mt' and 'nt' intersects in points (sx,sy), (mx,my) and (nx,ny),
// returns the relative order order of the distances from (sx,sy) to (mx,my) and (nx,ny).
// That is, indicates which offset triple intersects closer to (sx,sy)
// PRECONDITION: There exist single points at which the offsets at 'st', 'mt' and 'nt' intersect.
template<class K>
Uncertain<Comparison_result>
compare_offset_lines_isec_dist_to_pointC2 ( Seeded_trisegment_2<K> const& s
                                          , Seeded_trisegment_2<K> const& m
                                          , Seeded_trisegment_2<K> const& n
                                          )
{
  return compare_offset_lines_isec_dist_to_pointC2(construct_offset_lines_isecC2(s),m,n);
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
    rResult = certified_side_of_oriented_lineC2(a,b,c,aP->x(),aP->y()) == make_uncertain(POSITIVE);    
  }
  return rResult ;
}

// Given a triple of oriented straight line segments: (e0,e1,e2) such that their offsets
// at some distance intersects in a point (x,y), returns true if (x,y) is on the positive side of the line supporting aEdge
//
template<class K>
inline Uncertain<bool> is_edge_facing_offset_lines_isecC2 ( Seeded_trisegment_2<K> const& st, Segment_2<K> const& aEdge )
{
  return is_edge_facing_pointC2(construct_offset_lines_isecC2(st),aEdge);
}


// Given a triple of oriented straight line segments: (e0,e1,e2) such that their offsets
// at some distance intersects in a point (x,y), returns true if (x,y) is inside the passed offset zone.
// 
// An offset zone [z0,z1,z2], where zi is a line, is a region where the corresponding offset edges (oriented segment)
// Z0'->Z1'->Z2' remain connected in the absence of topological changes.
// It is the area (bounded or not) to the left of z1, to the right of the angular bisector (z0,z1) 
// and to the left of the angular bisector (z1,z2).
// This area represents all the possible offset edges Z1'.
//
// This predicate tells whether a split event, at (x,y), against z1, is effectively splitting the segment Z1'
// instead of hitting the supporting offseted line z1' but outside the segment.
//   
// Events are defined in terms of intersecting offset _lines_, not segments, thus if the event involves z1' 
// (that is, z1 is one of (e0,e1,e2)) then a neccesary condition for the event to actually exist is
// that the point of coallision hits a segment of z1' as bounded by the vertices shared with z0' and z2',
// and not just the line z1'. 
// This condition is equivalent to the condition that (x,y) be in the offset zone [z0,z1,z2] as described above
// (which refers to the input _lines_ and not the offset edges).
//
// This condition is neccesary but not sufficient becasue z1' might end up connected with lines other than z0' and z2' 
// after some further event, so this predicate is valid only in the context of an event at a known time 't' such
// that it is known that at time t, z1' is indeed connected to z0' and z2'.
//
// During the shrinking process, edges can anihiliate one another; that is, collide not in a single point
// but along a line segment (reach a common supporting line simultaneously).
// This is possible if and only if the edges are parallel.
// Exactly at the time when such an anhiliation event ocurrs, the two initially parallel edges become connected
// and collinear in the offset polygon (form a degenerate alley or anntenna).
// The offset zone can be "degenerate" in the sense that Z0',Z1' or Z2' can be parallel and even collinear.
// 
// Right _after_ an anhiliation event the degenerate edges collapse and dissapears from the offset polygon
// so no _subsequent_ event can involve such edges, but at the exact time of the anhilition various events
// can still involve the collapsed edges. Thus, a degnerate offset zone doesn't imply a split event cannot ocurr there.
//
// PRECONDITIONS:
//   There exist a single point at which the offset lines for e0,e1,e2 at 't' intersect.
//   'z1' must be one of (e0,e1,e2); that is, (x,y) must be exactly over the offseted z1' at time 't'.
//
template<class K>
Uncertain<bool> is_offset_lines_isec_inside_offset_zoneC2 ( Seeded_trisegment_2<K> const& st
                                                          , Seeded_trisegment_2<K> const& zone 
                                                          )
{
  typedef typename K::FT FT ;
  
  typedef Point_2<K> Point_2 ;
  typedef Line_2<K>  Line_2 ;
  
  typedef Trisegment_2<K> Trisegment_2 ;
  
  typedef optional<Point_2> Optional_point_2 ;
  typedef optional<Line_2>  Optional_line_2 ;
  
  Uncertain<bool> r = Uncertain<bool>::indeterminate();

  CGAL_assertion ( st.event().collinearity() != TRISEGMENT_COLLINEARITY_ALL ) ;
  
  Optional_line_2 zl = compute_normalized_line_ceoffC2(zone.event().e0()) ;
  Optional_line_2 zc = compute_normalized_line_ceoffC2(zone.event().e1()) ;
  Optional_line_2 zr = compute_normalized_line_ceoffC2(zone.event().e2()) ;

  // Construct intersection point (x,y)
  Optional_point_2 i = construct_offset_lines_isecC2(st);
  
  if ( zl && zc && zr && i ) // all properly computed
  {
    // Let L and R be the vertices formed by zl',zc' and zc',zr' resp.
    // Let ZC1 be the oriented segment of zc' bounded first by L then by R.
    // Let ZC0 be the ray of zc' before  L, and ZC2 the ray of zc' after R.
    // Let ZL0 be tha ray of zl' before L and ZL1 the ray of zl' after R
    // Let ZR0 be tha ray of zl' before L and ZR1 the ray of zl' after R
    
    // At the time of the event, each offseted edge zl', zc' and zr' are made of :
    //
    // zl': ZL0->L->ZL1
    // zc': ZC0->L->ZC1->R->ZC2
    // zr':         ZR0->R->ZR1
    
    // Here we need to determine if "i", known to be along zc', is inside ZC1 instead of ZC0 (and instead of ZC2)
            
    // Since all edges move at the same speed, there are 3 cases to consider:
    
    // (1) If L is convex, ZL1 is "behind" ZC1 while ZL0 is "before" ZC0.
    // (2) If L is reflex, ZL1 is "before" ZC1 while ZL0 is "behind" ZC0.
    // (3) If L is degenerate, that is, zl and zc are collinear, none if behind or before any other
    
    // If L is case (1) then "i" is inside ZC1 and not ZC0 if the signed distance to zc is smaller than to zl
    // If L is case (2) then "i" is inside ZC1 and not ZC0 if the signed distance to zc is larger  than to zl
    // If L is case (3) then "i" is inside ZC1 if its to the right of a line perpendicular to zc passing through L*
    // where L* is an offset-vertex between collinear edges
     
    // (likewise for R)
    
    // sdc : scaled (signed) distance from (x,y) to 'zc'
    FT sdc = zc->a() * i->x() + zc->b() * i->y() + zc->c() ;

    CGAL_STSKEL_TRAITS_TRACE("\nsdc=" << sdc ) ;

    // NOTE:
    //   if "i" is not on the positive side of 'zc' it isn't on it's offset zone.
    //   Also, if "i" is _over_ 'zc' (its signed distance to ec is not certainly positive) then by definition is not on its _offset_
    //   zone either.
    Uncertain<bool> cok = CGAL_NTS is_finite(sdc) ? CGAL_NTS certified_is_positive(sdc) : Uncertain<bool>::indeterminate() ;
    if ( ! CGAL_NTS is_indeterminate(cok) )
    {
      if ( cok == true )
      {
        CGAL_STSKEL_TRAITS_TRACE("\ncorrect side of zc." ) ;

        Uncertain<bool> lc_degenerate = are_edges_parallelC2(zone.event().e0(),zone.event().e1()); 
        Uncertain<bool> cr_degenerate = are_edges_parallelC2(zone.event().e1(),zone.event().e2()); 

        if ( ! CGAL_NTS is_indeterminate(lc_degenerate) && ! CGAL_NTS is_indeterminate(cr_degenerate) )
        {
          Uncertain<bool> lok = Uncertain<bool>::indeterminate() ;
          Uncertain<bool> rok = Uncertain<bool>::indeterminate() ;
          
          if ( !lc_degenerate )
          {
            // sld: scaled (signed) distances from "i" to 'zl'
            FT sdl = zl->a() * i->x() + zl->b() * i->y() + zl->c() ;
            
            if ( CGAL_NTS is_finite(sdl) )
            {
              CGAL_STSKEL_TRAITS_TRACE("\nsdl=" << sdl ) ;
              
              Uncertain<bool> lc_reflex = CGAL_NTS certified_is_smaller(zl->a()*zc->b(),zc->a()*zl->b());
              
              if ( ! CGAL_NTS is_indeterminate(lc_reflex) )
              {
                CGAL_STSKEL_TRAITS_TRACE("\nl:(zl,zc) is " << ( lc_reflex == true ? "reflex" : "non-reflex") ) ;

                lok = ( lc_reflex ? CGAL_NTS certified_is_smaller_or_equal(sdl,sdc)
                                  : CGAL_NTS certified_is_larger_or_equal (sdl,sdc) 
                      ) ;
              }
            }
            else
            {
              CGAL_STSKEL_TRAITS_TRACE("\nOverflow detected." ) ;
            }
          }
          else
          {
// std::cout << "Left zone vertex (zl,zc) is DEGENERATE." << std::endl ;
// std::cout << "  zl=" << zone.e0() << std::endl ;
// std::cout << "  zc=" << zone.e1() << std::endl ;
          CGAL_STSKEL_TRAITS_TRACE("\nl:(zl,zc) is DEGENERATE") ;
            Optional_point_2 l = compute_seed_pointC2(zone, Trisegment_2::LEFT);
            if ( l )
            {
//std::cout << "  l=(" << l->x() << "," << l->y() << ")" << std::endl ;
//std::cout << "  i=(" << i->x() << "," << i->y() << ")" << std::endl ;
              FT na, nb, nc ;
              perpendicular_through_pointC2(zc->a(),zc->b(),l->x(),l->y(),na, nb, nc);
              Uncertain<Sign> side = certified_side_of_oriented_lineC2(na,nb,nc,i->x(),i->y()) ;
              if ( !is_indeterminate(side) )
              {
//std::cout << "  Side of i=" << ((Sign)(side))<< std::endl ;              
                switch ( side )
                {
                  case POSITIVE : lok = make_uncertain(false) ; break ;
                  case NEGATIVE : lok = make_uncertain(true)  ; break ;
                  case ZERO     :     
                    lok = certified_side_of_oriented_lineC2(na,nb,nc
                                                           ,zone.event().e0().source().x()
                                                           ,zone.event().e0().source().y()
                                                           ) == make_uncertain(POSITIVE) ;
                    break ;
                }
//std::cout << "  Left vertex ok=" << (is_indeterminate(lok) ? "<don't know>" : ( (bool)lok ? "yes": "false") ) << std::endl ;
              }
//              else std::cout << " Unable to determine side of i" << std::endl ;              
            }
            else
            {
//std::cout << "  Unable to construct offset vertex" << std::endl ;

              CGAL_STSKEL_TRAITS_TRACE("\nOverflow detected." ) ;
            }
          }
          
          if ( !cr_degenerate )
          {
            // slr: scaled (signed) distances from "i" to 'zr'
            FT sdr = zr->a() * i->x() + zr->b() * i->y() + zr->c() ;
            
            if ( CGAL_NTS is_finite(sdr) )
            {
              CGAL_STSKEL_TRAITS_TRACE("\nsdr=" << sdr ) ;
              
              Uncertain<bool> cr_reflex = CGAL_NTS certified_is_smaller(zc->a()*zr->b(),zr->a()*zc->b());
              
              if ( ! CGAL_NTS is_indeterminate(cr_reflex) )
              {
                CGAL_STSKEL_TRAITS_TRACE("\nr:(zc,zr) is " << ( cr_reflex == true ? "reflex" : "non-reflex") ) ;
            
                rok = ( cr_reflex ? CGAL_NTS certified_is_smaller_or_equal(sdr,sdc)
                                  : CGAL_NTS certified_is_larger_or_equal (sdr,sdc) 
                      ) ;
              }
            }
            else
            {
              CGAL_STSKEL_TRAITS_TRACE("\nOverflow detected." ) ;
            }
          }
          else
          {
 //std::cout << "Right zone vertex (zc,zr) is DEGENERATE." << std::endl ;
 //std::cout << "  zc=" << zone.e1() << std::endl ;
 //std::cout << "  zr=" << zone.e2() << std::endl ;
            CGAL_STSKEL_TRAITS_TRACE("\nr:(zc,zr) is DEGENERATE") ;
            Optional_point_2 r = compute_seed_pointC2(zone,Trisegment_2::RIGHT);
            if ( r )
            {
//std::cout << "  r=(" << r->x() << "," << r->y() << ")" << std::endl ;
//std::cout << "  i=(" << i->x() << "," << i->y() << ")" << std::endl ;
              FT na, nb, nc ;
              perpendicular_through_pointC2(zc->a(),zc->b(),r->x(),r->y(),na, nb, nc);
              Uncertain<Sign> side = certified_side_of_oriented_lineC2(na,nb,nc,i->x(),i->y()) ;
              if ( !is_indeterminate(side) )
              {
//std::cout << "  Side of i=" << ((Sign)(side))<< std::endl ;              
                switch ( side )
                {
                  case NEGATIVE : rok = make_uncertain(false); break ;
                  case POSITIVE : rok = make_uncertain(true) ; break ;
                  case ZERO     : 
                    rok = certified_side_of_oriented_lineC2(na,nb,nc
                                                           ,zone.event().e2().target().x()
                                                           ,zone.event().e2().target().y()
                                                           )== make_uncertain(NEGATIVE);
                    break ;
                }
//std::cout << "  Right vertex ok=" << (is_indeterminate(rok) ? "<don't know>" : ( (bool)rok ? "yes": "false") ) << std::endl ;
              }
//              else std::cout << " Unable to determine side of i" << std::endl ;              
            }
            else
            {
//std::cout << "  Unable to construct offset vertex" << std::endl ;
              CGAL_STSKEL_TRAITS_TRACE("\nOverflow detected." ) ;
            }
          }
          
          CGAL_STSKEL_TRAITS_TRACE("\nlok:" << lok) ;
          CGAL_STSKEL_TRAITS_TRACE("\nrok:" << rok) ;
  
          r = CGAL_NTS logical_and(lok , rok) ;
        }
        else
        {
          CGAL_STSKEL_TRAITS_TRACE("\nUnable to reliably determine reflexivity of zone vertices." ) ;
        }
      }
      else
      {
        CGAL_STSKEL_TRAITS_TRACE("\nWRONG side of zc." ) ;
        r = make_uncertain(false);
      }
    }
    else
    {
      CGAL_STSKEL_TRAITS_TRACE("\nUnable to reliably determine side-of-line." ) ;
    }
  }
  else
  {
    CGAL_STSKEL_TRAITS_TRACE("\nOverflow detected." ) ;
  }

  return r ;
}

// Given 2 triples of oriented straight line segments (l0,l1,l2) and (r0,r1,r2), such that 
// the offsets at time 'tl' for triple 'l' intersects in a point (lx,ly) and 
// the offsets at time 'tr' for triple 'r' intersects in a point (rx,ry) 
// returns true if "tl==tr" and "(lx,ly)==(rx,ry)" 
// PRECONDITIONS:
//   There exist single points at which the offset lines for 'l' and 'r' at 'tl', 'tr' intersect.
//
template<class K>
Uncertain<bool> are_events_simultaneousC2 ( Seeded_trisegment_2<K> const& l, Seeded_trisegment_2<K> const& r )
{
  typedef typename K::FT FT ;
  
  typedef Point_2<K> Point_2 ;
  typedef Line_2<K>  Line_2 ;
  
  typedef Trisegment_2<K> Trisegment_2 ;
  
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

    CGAL_assertion ( CGAL_NTS certified_is_positive(lt) ) ;
    CGAL_assertion ( CGAL_NTS certified_is_positive(rt) ) ;

    Uncertain<bool> equal_times = CGAL_NTS certified_is_equal(lt,rt);
    
    if ( ! CGAL_NTS is_indeterminate(equal_times) )
    {
      if ( equal_times == true )
      {
        Optional_point_2 li = construct_offset_lines_isecC2(l);
        Optional_point_2 ri = construct_offset_lines_isecC2(r);

        if ( li && ri )
          rResult = CGAL_NTS logical_and( CGAL_NTS certified_is_equal(li->x(),ri->x())
                                        , CGAL_NTS certified_is_equal(li->y(),ri->y())
                                        ) ;
      }
      else rResult = make_uncertain(false);
    }
  }
  return rResult;
}


} // namespace CGAL_SS_i

CGAL_END_NAMESPACE

#endif // CGAL_STRAIGHT_SKELETON_PREDICATES_FTC2_H //
// EOF //

