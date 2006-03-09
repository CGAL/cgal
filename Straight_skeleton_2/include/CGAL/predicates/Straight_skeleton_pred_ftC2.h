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
#ifndef CGAL_STRAIGHT_SKELETON_PREDICATES_FTC2_H 
#define CGAL_STRAIGHT_SKELETON_PREDICATES_FTC2_H 1

#include <CGAL/constructions/Straight_skeleton_cons_ftC2.h>
#include <CGAL/Uncertain.h>
#include <CGAL/certified_quotient_predicates.h>

CGAL_BEGIN_NAMESPACE


namespace CGAL_SLS_i
{

template<class FT>
Uncertain<bool> certified_collinearC2( Vertex<FT> const& p, Vertex<FT> const& q, Vertex<FT> const& r )
{
  return CGAL_NTS certified_is_equal( ( q.x() - p.x() ) * ( r.y() - p.y() )
                                    , ( r.x() - p.x() ) * ( q.y() - p.y() )
                                    );
}

template<class FT>
Uncertain<bool> are_edges_collinear( Edge<FT> const& e0, Edge<FT> const& e1 )
{
  return CGAL_NTS logical_and(  certified_collinearC2(e0.s(),e0.t(),e1.s())
                             ,  certified_collinearC2(e0.s(),e0.t(),e1.t())
                             ) ;
}

template<class FT>
Uncertain<bool> are_edges_parallel( Edge<FT> const& e0, Edge<FT> const& e1 )
{
  Uncertain<Sign> s = certified_sign_of_determinant2x2(e0.t().x() - e0.s().x()
                                                      ,e0.t().y() - e0.s().y()
                                                      ,e1.t().x() - e1.s().x()
                                                      ,e1.t().y() - e1.s().y()
                                                     ) ;
  return s == Uncertain<Sign>(ZERO);  
}


template<class FT>
SortedTriedge<FT> collinear_sort ( Triedge<FT> const& triedge )
{
  int lCollinearCount = -1 ;
  
  int  idx0=0, idx1=1, idx2=2 ;

  Uncertain<bool> is_01 = are_edges_collinear(triedge.e0(),triedge.e1());
  if ( !CGAL_NTS is_indeterminate(is_01) )
  {
    Uncertain<bool> is_02 = are_edges_collinear(triedge.e0(),triedge.e2());
    if ( !CGAL_NTS is_indeterminate(is_02) )
    {
      Uncertain<bool> is_12 = are_edges_collinear(triedge.e1(),triedge.e2());
      if ( !CGAL_NTS is_indeterminate(is_12) )
      {
        if ( CGAL_NTS logical_and(is_01 , !is_02 , !is_12 ) )
        {
          idx0 = 0 ;
          idx1 = 1 ;
          idx2 = 2 ;
          lCollinearCount = 2 ;
        }
        else if ( CGAL_NTS logical_and(is_02 , !is_01 , !is_12 ) )
        {
          idx0 = 0 ;
          idx1 = 2 ;
          idx2 = 1 ;
          lCollinearCount = 2 ;
        }
        else if ( CGAL_NTS logical_and(is_12 , !is_01 , !is_02 ) )
        {
          idx0 = 1 ;
          idx1 = 2 ;
          idx2 = 0 ;
          lCollinearCount = 2 ;
        }
        else if ( CGAL_NTS logical_and(!is_01 , !is_02, !is_12  ) )
        {
          idx0 = 0 ;
          idx1 = 1 ;
          idx2 = 2 ;
          lCollinearCount = 0 ;
        }
        else
          lCollinearCount = 3 ;
      }
    }
  }

  return SortedTriedge<FT>(triedge.e(idx0),triedge.e(idx1),triedge.e(idx2),lCollinearCount);
}


// Given some input edge E1, there is a sequence of offset polygons, in some time range (possibly including time 0),
// such that E1 has E0 and E2 as adjacent edges.
// At time 0, E0 and E2 are the edges adjacent to E1 in the input polygon, but later, E0 and E2 can be any edges.
// There is some time range were E0->E1->E2 are consecutive in each offset polygon. For each specific t
// in the range, the vertices (E0,E1) and (E1,E2) move along the bisectors (E0,E1) and (E1,E2). Thus, the actual length
// and position of each E1 offset varies along the time period.
// The set of all the E1 offsets which have E0 and E2 as inmediate neighbors is called the "Offset Zone" for E1 w.r.t to E0 and E2.
// 
// Such offset zone is the part of the halfplane to the left of E1 (were _all_ offsets to E1 are) bounded between the bisectors
// (E0,E1) and (E1,E2)
//
// If E0 and/or E2 are parallel to E1, then there is just 1 E1 offset, located at the medial bisector of the E2 and its parallel(s).
//
// But E0 and E2 are the inmediate neighbors of E1 in the offset polygon (at some t), so if E0 and/or E2 are parallel to E1,
// it means that such an offset polygon has degenerated to an innermost segment. Thus, E1 can't be split further by any
// other reflex wavefront. 
// Therefore, is E0 or E2 are parallel to E1, the offset zone for E1 w.r.t to E0,E2 is said to be "degenerate".
template<class FT>
Uncertain<bool> is_offset_zone_degenerate ( Triedge<FT> const& zone )
{
  Uncertain<bool> is_01 = are_edges_parallel(zone.e0(),zone.e1());
  Uncertain<bool> is_12 = are_edges_parallel(zone.e1(),zone.e2());
  
  return CGAL_NTS logical_or(is_01,is_12);                
}

// Given 3 oriented straight line segments: l0, l1, l2 [each segment is passed as (sx,sy,tx,ty)]
// returns true if there exist some positive offset distance 't' for which the
// leftward-offsets of their supporting lines intersect at a single point.
// NOTE: This function allows l0 and l1 to be collinear if they are equally oriented.
// This allows the algorithm to handle degenerate vertices (formed by 3 collinear consecutive points)
template<class FT>
Uncertain<bool> exist_offset_lines_isec2 ( Triedge<FT> const& triedge )
{
  Uncertain<bool> rResult = Uncertain<bool>::indeterminate();

  SortedTriedge<FT> sorted = collinear_sort(triedge);

  if ( !sorted.is_indeterminate() ) // Couldn't determine collinearity
  {
    if ( sorted.collinear_count() < 3 ) // If the 3 edges are collinear thre is no event.
    {
      CGAL_SSTRAITS_TRACE( ( sorted.collinear_count() == 0 ? " normal edges" : " collinear edges" ) ) ;
  
      Rational<FT> t = compute_offset_lines_isec_timeC2(sorted) ;
  
      Uncertain<bool> d_is_zero = CGAL_NTS certified_is_zero(t.d()) ;
      if ( ! CGAL_NTS is_indeterminate(d_is_zero) )
      {
        if ( !d_is_zero )
        {
          rResult = CGAL_NTS certified_is_positive(t.to_quotient()) ;
          CGAL_SSTRAITS_TRACE("\nEvent time: " << (t.n()/t.d()) << ". Event " << ( rResult ? "exist." : "doesn't exist." ) ) ;
        }
        else
        {
          CGAL_SSTRAITS_TRACE("\nDenominator exactly zero, Event doesn't exist." ) ;
          rResult = make_uncertain(false);
        }
      }
      else
        CGAL_SSTRAITS_TRACE("\nDenominator is probably zero (but not exactly), event existance is indeterminate." ) ;
    }
    else
    {
      CGAL_SSTRAITS_TRACE("\nAll the edges are collinear. Event doesn't exist." ) ;
      rResult = make_uncertain(false);
    }
  }
  else
    CGAL_SSTRAITS_TRACE("\nEdges uncertainly collinear, event existance is indeterminate." ) ;

  return rResult ;
}

// Given 2 triples of oriented lines in _normalized_ implicit form: (m0,m1,m2) and (n0,n1,n2), such that
// for each triple there exists distances 'mt' and 'nt' for which the offsets lines (at mt and nt resp.),
// (m0',m1',m2') and (n0',n1',n2') intersect each in a single point; returns the relative order of mt w.r.t nt.
// That is, indicates which offset triple intersects first (closer to the source lines)
// PRECONDITION: There exist distances mt and nt for which each offset triple intersect at a single point.
template<class FT>
Uncertain<Comparison_result> compare_offset_lines_isec_timesC2 ( Triedge<FT> const& m, Triedge<FT> const& n )
{
  Uncertain<Comparison_result> rResult = Uncertain<Comparison_result>::indeterminate();

  SortedTriedge<FT> m_sorted = collinear_sort(m);
  SortedTriedge<FT> n_sorted = collinear_sort(n);

  if ( ! ( m_sorted.is_indeterminate() || n_sorted.is_indeterminate() ) )
  {
    CGAL_assertion ( m_sorted.collinear_count() < 3 ) ;
    CGAL_assertion ( n_sorted.collinear_count() < 3 ) ;
    
    Quotient<FT> mt = compute_offset_lines_isec_timeC2(m_sorted).to_quotient();
    Quotient<FT> nt = compute_offset_lines_isec_timeC2(n_sorted).to_quotient();
   
    CGAL_assertion ( CGAL_NTS certified_is_positive(mt) ) ;
    CGAL_assertion ( CGAL_NTS certified_is_positive(nt) ) ;

    rResult = CGAL_NTS certified_compare(mt,nt);
  }

  return rResult ;

}

// Given a point (px,py) and 2 triples of oriented lines in _normalized_ implicit form: (m0,m1,m2) and (n0,n1,n2),
// such that their offsets at distances 'mt' and 'nt' intersects in points (mx,my) and (nx,ny),
// returns the relative order order of the distances from (px,py) to (mx,my) and (nx,ny).
// That is, indicates which offset triple intersects closer to (px,py)
// PRECONDITION: There exist single points at which the offset line triples 'm' and 'n' at 'mt' and 'nt' intersect.
template<class FT>
Uncertain<Comparison_result>
compare_offset_lines_isec_dist_to_pointC2 ( Vertex<FT> const& p, Triedge<FT> const& m, Triedge<FT> const& n )
{
  Uncertain<Comparison_result> rResult = Uncertain<Comparison_result>::indeterminate();

  SortedTriedge<FT> m_sorted = collinear_sort(m);
  SortedTriedge<FT> n_sorted = collinear_sort(n);

  if ( ! ( m_sorted.is_indeterminate() || n_sorted.is_indeterminate() ) )
  {
    CGAL_assertion ( m_sorted.collinear_count() < 3 ) ;
    CGAL_assertion ( n_sorted.collinear_count() < 3 ) ;
    
    FT dm = compute_offset_lines_isec_dist_to_pointC2(p,m_sorted);
    FT dn = compute_offset_lines_isec_dist_to_pointC2(p,n_sorted);

    rResult = CGAL_NTS certified_compare(dm,dn);
  }

  return rResult ;
}

// Given 3 triples of oriented lines in _normalized_ implicit form: (s0,s1,s2), (m0,m1,m2) and (n0,n1,n2),
// such that their offsets at distances 'st', 'mt' and 'nt' intersects in points (sx,sy), (mx,my) and (nx,ny),
// returns the relative order order of the distances from (sx,sy) to (mx,my) and (nx,ny).
// That is, indicates which offset triple intersects closer to (sx,sy)
// PRECONDITION: There exist single points at which the offsets at 'st', 'mt' and 'nt' intersect.
template<class FT>
Uncertain<Comparison_result>
compare_offset_lines_isec_dist_to_pointC2 ( Triedge<FT> const& s
                                          , Triedge<FT> const& m
                                          , Triedge<FT> const& n
                                          )
{
  Uncertain<Comparison_result> rResult = Uncertain<Comparison_result>::indeterminate();

  SortedTriedge<FT> s_sorted = collinear_sort(s);

  if ( !s_sorted.is_indeterminate() )
  {
    CGAL_assertion ( s_sorted.collinear_count() < 3 ) ;
    
    rResult = compare_offset_lines_isec_dist_to_pointC2(construct_offset_lines_isecC2(s_sorted),m,n);
  }
  
  return rResult ;
}



// Given a triple of oriented lines in _normalized_ implicit form: (e0,e1,e2) such that their offsets
// at somne distance intersects in a point (x,y); and another triple of oriented lines in _normalized_ implicit form:
// (el,ec,er); returns true if the intersection point (x,y) is inside the offset zone of 'ec' w.r.t 'el' and 'er';
// that is, the locus of points to the left of 'ec', to the right of the bisector (el,ec) and to the left of the bisector
// (ec,er).
// PRECONDITIONS:
//   There exist a single point at which the offset lines for e0,e1,e2 at 't' intersect.
//   'ec' must be one of (e0,e1,e2); that is, (x,y) must be exactly over the offseted 'ec' at time 't'.
//   (el,ec) and (ec,er) must be oblique; i.e, there must exist non-degenerate vertices at the intersections (el,ec) and (ec,er).
//
template<class FT>
Uncertain<bool>
is_offset_lines_isec_inside_offset_zoneC2 ( Triedge<FT> const& event, Triedge<FT> const& zone )
{
  Uncertain<bool> r = Uncertain<bool>::indeterminate();

  SortedTriedge<FT> e_sorted = collinear_sort(event);
  
  Uncertain<bool> degenerate_zone = is_offset_zone_degenerate(zone);
  
  if ( !(e_sorted.is_indeterminate() || is_indeterminate(degenerate_zone)) )
  {
    CGAL_assertion ( e_sorted.collinear_count() < 3 ) ;
    
    if ( !degenerate_zone )
    {
      Line<FT> zl = compute_normalized_line_ceoffC2(zone.e0()) ;
      Line<FT> zc = compute_normalized_line_ceoffC2(zone.e1()) ;
      Line<FT> zr = compute_normalized_line_ceoffC2(zone.e2()) ;
  
      // Construct intersection point (x,y)
      Vertex<FT> i = construct_offset_lines_isecC2(e_sorted);
  
      // Calculate scaled (signed) distance from (x,y) to 'zc'
      FT sdc = zc.a() * i.x() + zc.b() * i.y() + zc.c() ;
  
      CGAL_SSTRAITS_TRACE("\nsdc=" << sdc ) ;
  
      // NOTE:
      //   if (x,y) is not on the positive side of 'ec' it isn't on it's offset zone.
      //   Also, if (x,y) is over 'ec' (its signed distance to ec is not certainly positive) then by definition is not on its _offset_
      //   zone either.
      Uncertain<bool> cok = CGAL_NTS certified_is_positive(sdc);
      if ( ! CGAL_NTS is_indeterminate(cok) )
      {
        if ( cok == true )
        {
          CGAL_SSTRAITS_TRACE("\nright side of ec." ) ;
  
          // Calculate scaled (signed) distances from (x,y) to 'el' and 'er'
          FT sdl = zl.a() * i.x() + zl.b() * i.y() + zl.c() ;
          FT sdr = zr.a() * i.x() + zr.b() * i.y() + zr.c() ;
  
          CGAL_SSTRAITS_TRACE("\nsdl=" << sdl ) ;
          CGAL_SSTRAITS_TRACE("\nsdr=" << sdr ) ;
  
          // Determine if the vertices (el,ec) and (ec,er) are reflex.
          Uncertain<bool> lcx = CGAL_NTS certified_is_smaller(zl.a()*zc.b(),zc.a()*zl.b());
          Uncertain<bool> crx = CGAL_NTS certified_is_smaller(zc.a()*zr.b(),zr.a()*zc.b());
  
          if ( ! CGAL_NTS is_indeterminate(lcx) && ! CGAL_NTS is_indeterminate(crx) )
          {
            CGAL_SSTRAITS_TRACE("\n(el,ec) reflex:" << lcx ) ;
            CGAL_SSTRAITS_TRACE("\n(ec,er) reflex:" << crx ) ;
  
            // Is (x,y) to the right|left of the bisectors (el,ec) and (ec,er)?
            //  It depends on whether the vertex ((el,ec) and (ec,er)) is relfex or not.
            //  If it is reflex, then (x,y) is to the right|left of the bisector if sdl|sdr <= sdc; otherwise, if sdc <= sdl|srd
  
            Uncertain<bool> lok = lcx ? CGAL_NTS certified_is_smaller_or_equal(sdl,sdc)
                                      : CGAL_NTS certified_is_smaller_or_equal(sdc,sdl) ;
  
            Uncertain<bool> rok = crx ? CGAL_NTS certified_is_smaller_or_equal(sdr,sdc)
                                      : CGAL_NTS certified_is_smaller_or_equal(sdc,sdr) ;
  
            CGAL_SSTRAITS_TRACE("\nlok:" << lok) ;
            CGAL_SSTRAITS_TRACE("\nrok:" << rok) ;
  
            r = CGAL_NTS logical_and(lok , rok) ;
          }
          else
          {
            CGAL_SSTRAITS_TRACE("\nUnable to reliably determine side-of-line." ) ;
          }
        }
        else
        {
          CGAL_SSTRAITS_TRACE("\nWRONG side of ec." ) ;
          r = make_uncertain(false);
        }
      }
      else
      {
        CGAL_SSTRAITS_TRACE("\nUnable to reliably determine side-of-line." ) ;
      }
    }
    else
    {
      CGAL_SSTRAITS_TRACE("\nDegenerate offset zone. Edge collapsed." ) ;
      r = make_uncertain(false);
    }
  }
  else
  {
    CGAL_SSTRAITS_TRACE("\nUnable to determine collinearity of event triedge or parallelity of zone triedge." ) ;
  }

  return r ;
}

template<class FT>
Uncertain<bool> are_events_simultaneousC2 ( Triedge<FT> const& l, Triedge<FT> const& r )
{
  Uncertain<bool> rResult = Uncertain<bool>::indeterminate();

  SortedTriedge<FT> l_sorted = collinear_sort(l);
  SortedTriedge<FT> r_sorted = collinear_sort(r);

  if ( ! ( l_sorted.is_indeterminate() || r_sorted.is_indeterminate() ) )
  {
    CGAL_assertion ( l_sorted.collinear_count() < 3 ) ;
    CGAL_assertion ( r_sorted.collinear_count() < 3 ) ;
    
    Quotient<FT> lt = compute_offset_lines_isec_timeC2(l_sorted).to_quotient();
    Quotient<FT> rt = compute_offset_lines_isec_timeC2(r_sorted).to_quotient();

    CGAL_assertion ( CGAL_NTS certified_is_positive(lt) ) ;
    CGAL_assertion ( CGAL_NTS certified_is_positive(rt) ) ;

    Uncertain<bool> equal_times = CGAL_NTS certified_is_equal(lt,rt);
    if ( ! CGAL_NTS is_indeterminate(equal_times) )
    {
      if ( equal_times == true )
      {
        Vertex<FT> li = construct_offset_lines_isecC2(l_sorted);
        Vertex<FT> ri = construct_offset_lines_isecC2(r_sorted);

        rResult = CGAL_NTS logical_and( CGAL_NTS certified_is_equal(li.x(),ri.x())
                                      , CGAL_NTS certified_is_equal(li.y(),ri.y())
                                      ) ;
      }
      else rResult = make_uncertain(false);
    }
  }

  return rResult;
}

} // namespace CGAL_SLS_i

CGAL_END_NAMESPACE

#endif // CGAL_STRAIGHT_SKELETON_PREDICATES_FTC2_H //
// EOF //

