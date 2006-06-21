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

template<class K>
Uncertain<bool> certified_collinearC2( Point_2<K> const& p
                                     , Point_2<K> const& q
                                     , Point_2<K> const& r 
                                     )
{
  return CGAL_NTS certified_is_equal( ( q.x() - p.x() ) * ( r.y() - p.y() )
                                    , ( r.x() - p.x() ) * ( q.y() - p.y() )
                                    );
}

template<class K>
Uncertain<bool> are_edges_collinearC2( Segment_2<K> const& e0, Segment_2<K> const& e1 )
{
  return CGAL_NTS logical_and(  certified_collinearC2(e0.source(),e0.target(),e1.source())
                             ,  certified_collinearC2(e0.source(),e0.target(),e1.target())
                             ) ;
}

template<class K>
Uncertain<bool> are_edges_parallelC2( Segment_2<K> const& e0, Segment_2<K> const& e1 )
{
  Uncertain<Sign> s = certified_sign_of_determinant2x2(e0.target().x() - e0.source().x()
                                                      ,e0.target().y() - e0.source().y()
                                                      ,e1.target().x() - e1.source().x()
                                                      ,e1.target().y() - e1.source().y()
                                                     ) ;
  return s == Uncertain<Sign>(ZERO);  
}

//
// Constructs a Sorted_triedge_2 which stores 3 edges (segments) such that
// if two of them are collinear, they are put first, as e0, e1.
// Stores also the number of collinear edges. which should be 0 or 2.
//
// If the collinearity test is indeterminate for any pair of edges the
// resulting sorted triedge is itself indeterminate 
// (encoded as a collinear count of -1)
//
template<class K>
Sorted_triedge_2<K> collinear_sort ( Triedge_2<K> const& triedge )
{
  int lCollinearCount = -1 ;
  
  int  idx0=0, idx1=1, idx2=2 ;

  Uncertain<bool> is_01 = are_edges_collinearC2(triedge.e0(),triedge.e1());
  if ( !CGAL_NTS is_indeterminate(is_01) )
  {
    Uncertain<bool> is_02 = are_edges_collinearC2(triedge.e0(),triedge.e2());
    if ( !CGAL_NTS is_indeterminate(is_02) )
    {
      Uncertain<bool> is_12 = are_edges_collinearC2(triedge.e1(),triedge.e2());
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

  return Sorted_triedge_2<K>(triedge.e(idx0),triedge.e(idx1),triedge.e(idx2),lCollinearCount);
}

// Returns true if the offset-zone for the triedge 'zone' is degenerate.
//
// Since an offset zone [e0->e1->e2] is the region inside the polygon where the three edges remain connected as they
// move inward, if e0 and/or e2 is parallel to e1 then the zone is degenerate (collapsed to the bisecting line of the parallel edges)
//
// (that is, the parallel zone edges, since they share a vertex, have already collide and can't move any further)
//
template<class K>
Uncertain<bool> is_offset_zone_degenerate ( Triedge_2<K> const& zone )
{
  Uncertain<bool> is_01 = are_edges_parallelC2(zone.e0(),zone.e1());
  Uncertain<bool> is_12 = are_edges_parallelC2(zone.e1(),zone.e2());
  
  return CGAL_NTS logical_or(is_01,is_12);                
}

// Given 3 oriented straight line segments: e0, e1, e2 
// returns true if there exist some positive offset distance 't' for which the
// leftward-offsets of their supporting lines intersect at a single point.
//
// NOTE: This function allows e0 and e1 to be collinear if they are equally oriented,
// or parallel if they have opposite orientation. This allows the algorithm to handle degenerate vertices 
// (formed by 3 collinear consecutive points) and mutually collapsing edge events.
template<class K>
Uncertain<bool> exist_offset_lines_isec2 ( Triedge_2<K> const& triedge )
{
  typedef typename K::FT FT ;
  
  typedef Sorted_triedge_2<K> Sorted_triedge_2 ;
  
  typedef Rational<FT>       Rational ;
  typedef optional<Rational> Optional_rational ;
  
  Uncertain<bool> rResult = Uncertain<bool>::indeterminate();

  Sorted_triedge_2 sorted = collinear_sort(triedge);

  if ( !sorted.is_indeterminate() ) // Couldn't determine collinearity
  {
    if ( sorted.collinear_count() < 3 ) // If the 3 edges are collinear thre is no event.
    {
      CGAL_STSKEL_TRAITS_TRACE( ( sorted.collinear_count() == 0 ? " normal edges" : " collinear edges" ) ) ;
  
      Optional_rational t = compute_offset_lines_isec_timeC2(sorted) ;
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
  }
  else
    CGAL_STSKEL_TRAITS_TRACE("\nEdges uncertainly collinear, event existance is indeterminate." ) ;

  return rResult ;
}

// Given 2 triples of oriented straight line segments: (m0,m1,m2) and (n0,n1,n2), such that
// for each triple there exists distances 'mt' and 'nt' for which the offsets lines (at mt and nt resp.),
// (m0',m1',m2') and (n0',n1',n2') intersect each in a single point; returns the relative order of mt w.r.t nt.
// That is, indicates which offset triple intersects first (closer to the source lines)
// PRECONDITION: There exist distances mt and nt for which each offset triple intersect at a single point.
template<class K>
Uncertain<Comparison_result> compare_offset_lines_isec_timesC2 ( Triedge_2<K> const& m, Triedge_2<K> const& n )
{
  typedef typename K::FT FT ;
  
  typedef Sorted_triedge_2<K> Sorted_triedge_2 ;
  
  typedef Rational<FT>       Rational ;
  typedef Quotient<FT>       Quotient ;
  typedef optional<Rational> Optional_rational ;
  
  Uncertain<Comparison_result> rResult = Uncertain<Comparison_result>::indeterminate();

  Sorted_triedge_2 m_sorted = collinear_sort(m);
  Sorted_triedge_2 n_sorted = collinear_sort(n);

  if ( ! ( m_sorted.is_indeterminate() || n_sorted.is_indeterminate() ) )
  {
    CGAL_assertion ( m_sorted.collinear_count() < 3 ) ;
    CGAL_assertion ( n_sorted.collinear_count() < 3 ) ;
    
    Optional_rational mt_ = compute_offset_lines_isec_timeC2(m_sorted);
    Optional_rational nt_ = compute_offset_lines_isec_timeC2(n_sorted);
    
    if ( mt_ && nt_ )
    {
      Quotient mt = mt_->to_quotient();
      Quotient nt = nt_->to_quotient();
     
      CGAL_assertion ( CGAL_NTS certified_is_positive(mt) ) ;
      CGAL_assertion ( CGAL_NTS certified_is_positive(nt) ) ;
  
      rResult = CGAL_NTS certified_compare(mt,nt);
    }
    
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
                                          , Triedge_2<K>           const& m
                                          , Triedge_2<K>           const& n 
                                          )
{
  typedef typename K::FT FT ;
  
  typedef Sorted_triedge_2<K> Sorted_triedge_2 ;
  
  typedef Rational<FT>       Rational ;
  typedef Quotient<FT>       Quotient ;
  typedef optional<Rational> Optional_rational ;
  
  Uncertain<Comparison_result> rResult = Uncertain<Comparison_result>::indeterminate();
  
  if ( p )
  {
    Sorted_triedge_2 m_sorted = collinear_sort(m);
    Sorted_triedge_2 n_sorted = collinear_sort(n);
  
    if ( ! ( m_sorted.is_indeterminate() || n_sorted.is_indeterminate() ) )
    {
      CGAL_assertion ( m_sorted.collinear_count() < 3 ) ;
      CGAL_assertion ( n_sorted.collinear_count() < 3 ) ;
      
      optional<FT> dm = compute_offset_lines_isec_dist_to_pointC2(p,m_sorted);
      optional<FT> dn = compute_offset_lines_isec_dist_to_pointC2(p,n_sorted);
  
      if ( dm && dn )
        rResult = CGAL_NTS certified_compare(*dm,*dn);
    }
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
compare_offset_lines_isec_dist_to_pointC2 ( Triedge_2<K> const& s
                                          , Triedge_2<K> const& m
                                          , Triedge_2<K> const& n
                                          )
{
  Uncertain<Comparison_result> rResult = Uncertain<Comparison_result>::indeterminate();

  Sorted_triedge_2<K> s_sorted = collinear_sort(s);

  if ( !s_sorted.is_indeterminate() )
  {
    CGAL_assertion ( s_sorted.collinear_count() < 3 ) ;
    
    rResult = compare_offset_lines_isec_dist_to_pointC2(construct_offset_lines_isecC2(s_sorted),m,n);
  }
  
  return rResult ;
}



// Given a triple of oriented straight line segments: (e0,e1,e2) such that their offsets
// at some distance intersects in a point (x,y), returns true if (x,y) is inside the passed offset zone.
// 
// An offset zone [z0->z1->z2] is a region where the offset edges z0',z1',z2' remain connected in the absence of topological changes.
// It is the area (bounded or not) to the left of z1, to the right of the angular bisector (z0,z1) 
// and to the left of the angular bisector (z1,z2).
// This area represents all the possible offset edges z1'.
// If the event involves z1' (that is, z1 is one of (e0,e1,e2)) then a neccesary condition for the event to exist is
// that the point of coallision hits z1' as bounded by the vertices shared with z0' and z2',
// and not just the line supporting z1'. 
// This condition is equivalent to the condition that (x,y) be in the offset zone [z0->z1->z2] as described above
// (which refers to the input edges and not the offset edges).
//
// This condition is neccesary but not sufficient becasue z1' might end up connected with edges other than z0' and z2' 
// after some further event, so this predicate is valid only in the context of an event at a known time 't' such
// that it is known that at time t, z1' is indeed connected to z0' and z2'.
//
// If z0 and/or z2 are parallel to z1 then the offset zone is "degenerate" and the event is conventionally considered
// NOT inside the zone, so this predicate returns false in that case.
//
// PRECONDITIONS:
//   There exist a single point at which the offset lines for e0,e1,e2 at 't' intersect.
//   'z1' must be one of (e0,e1,e2); that is, (x,y) must be exactly over the offseted 'z1' at time 't'.
//
template<class K>
Uncertain<bool>
is_offset_lines_isec_inside_offset_zoneC2 ( Triedge_2<K> const& event, Triedge_2<K> const& zone )
{
  typedef typename K::FT FT ;
  
  typedef Point_2<K> Point_2 ;
  typedef Line_2<K>  Line_2 ;
  
  typedef Sorted_triedge_2<K> Sorted_triedge_2 ;
  
  typedef optional<Point_2> Optional_point_2 ;
  typedef optional<Line_2>  Optional_line_2 ;
  
  Uncertain<bool> r = Uncertain<bool>::indeterminate();

  Sorted_triedge_2 e_sorted = collinear_sort(event);
  
  Uncertain<bool> degenerate_zone = is_offset_zone_degenerate(zone);
  
  if ( !(e_sorted.is_indeterminate() || is_indeterminate(degenerate_zone)) )
  {
    CGAL_assertion ( e_sorted.collinear_count() < 3 ) ;
    
    if ( !degenerate_zone )
    {
      Optional_line_2 zl = compute_normalized_line_ceoffC2(zone.e0()) ;
      Optional_line_2 zc = compute_normalized_line_ceoffC2(zone.e1()) ;
      Optional_line_2 zr = compute_normalized_line_ceoffC2(zone.e2()) ;
  
      // Construct intersection point (x,y)
      Optional_point_2 i = construct_offset_lines_isecC2(e_sorted);
      
      if ( zl && zc && zr && i)
      {
        // Calculate scaled (signed) distance from (x,y) to 'zc'
        FT sdc = zc->a() * i->x() + zc->b() * i->y() + zc->c() ;
    
        CGAL_STSKEL_TRAITS_TRACE("\nsdc=" << sdc ) ;
    
        // NOTE:
        //   if (x,y) is not on the positive side of 'ec' it isn't on it's offset zone.
        //   Also, if (x,y) is over 'ec' (its signed distance to ec is not certainly positive) then by definition is not on its _offset_
        //   zone either.
        Uncertain<bool> cok = CGAL_NTS is_finite(sdc) ? CGAL_NTS certified_is_positive(sdc) : Uncertain<bool>::indeterminate() ;
        if ( ! CGAL_NTS is_indeterminate(cok) )
        {
          if ( cok == true )
          {
            CGAL_STSKEL_TRAITS_TRACE("\nright side of ec." ) ;
    
            // Determine if the vertices (el,ec) and (ec,er) are reflex.
            Uncertain<bool> lcx = CGAL_NTS certified_is_smaller(zl->a()*zc->b(),zc->a()*zl->b());
            Uncertain<bool> crx = CGAL_NTS certified_is_smaller(zc->a()*zr->b(),zr->a()*zc->b());
    
            if ( ! CGAL_NTS is_indeterminate(lcx) && ! CGAL_NTS is_indeterminate(crx) )
            {
              CGAL_STSKEL_TRAITS_TRACE("\n(el,ec) reflex:" << lcx ) ;
              CGAL_STSKEL_TRAITS_TRACE("\n(ec,er) reflex:" << crx ) ;
    
              // Calculate scaled (signed) distances from (x,y) to 'el' and 'er'
              FT sdl = zl->a() * i->x() + zl->b() * i->y() + zl->c() ;
              FT sdr = zr->a() * i->x() + zr->b() * i->y() + zr->c() ;
      
              if ( CGAL_NTS is_finite(sdl) && CGAL_NTS is_finite(sdc) )
              {
                CGAL_STSKEL_TRAITS_TRACE("\nsdl=" << sdl ) ;
                CGAL_STSKEL_TRAITS_TRACE("\nsdr=" << sdr ) ;
        
                // Is (x,y) to the right|left of the bisectors (el,ec) and (ec,er)?
                //  It depends on whether the vertex ((el,ec) and (ec,er)) is relfex or not.
                //  If it is reflex, then (x,y) is to the right|left of the bisector if sdl|sdr <= sdc; otherwise, if sdc <= sdl|srd
      
                Uncertain<bool> lok = lcx ? CGAL_NTS certified_is_smaller_or_equal(sdl,sdc)
                                          : CGAL_NTS certified_is_smaller_or_equal(sdc,sdl) ;
      
                Uncertain<bool> rok = crx ? CGAL_NTS certified_is_smaller_or_equal(sdr,sdc)
                                          : CGAL_NTS certified_is_smaller_or_equal(sdc,sdr) ;
      
                CGAL_STSKEL_TRAITS_TRACE("\nlok:" << lok) ;
                CGAL_STSKEL_TRAITS_TRACE("\nrok:" << rok) ;
      
                r = CGAL_NTS logical_and(lok , rok) ;
              }
              else
              {
                CGAL_STSKEL_TRAITS_TRACE("\nOverflow detected." ) ;
              }
            }
            else
            {
              CGAL_STSKEL_TRAITS_TRACE("\nUnable to reliably determine side-of-line." ) ;
            }
          }
          else
          {
            CGAL_STSKEL_TRAITS_TRACE("\nWRONG side of ec." ) ;
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
    }
    else
    {
      CGAL_STSKEL_TRAITS_TRACE("\nDegenerate offset zone. Edge collapsed." ) ;
      r = make_uncertain(false);
    }
  }
  else
  {
    CGAL_STSKEL_TRAITS_TRACE("\nUnable to determine collinearity of event triedge or parallelity of zone triedge." ) ;
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
Uncertain<bool> are_events_simultaneousC2 ( Triedge_2<K> const& l, Triedge_2<K> const& r )
{
  typedef typename K::FT FT ;
  
  typedef Point_2<K> Point_2 ;
  typedef Line_2<K>  Line_2 ;
  
  typedef Sorted_triedge_2<K> Sorted_triedge_2 ;
  
  typedef Rational<FT> Rational ;
  typedef Quotient<FT> Quotient ;
  
  typedef optional<Rational> Optional_rational ;
  typedef optional<Point_2>  Optional_point_2 ;
  
  Uncertain<bool> rResult = Uncertain<bool>::indeterminate();

  Sorted_triedge_2 l_sorted = collinear_sort(l);
  Sorted_triedge_2 r_sorted = collinear_sort(r);

  if ( ! ( l_sorted.is_indeterminate() || r_sorted.is_indeterminate() ) )
  {
    CGAL_assertion ( l_sorted.collinear_count() < 3 ) ;
    CGAL_assertion ( r_sorted.collinear_count() < 3 ) ;
    
    Optional_rational lt_ = compute_offset_lines_isec_timeC2(l_sorted);
    Optional_rational rt_ = compute_offset_lines_isec_timeC2(r_sorted);
    
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
          Optional_point_2 li = construct_offset_lines_isecC2(l_sorted);
          Optional_point_2 ri = construct_offset_lines_isecC2(r_sorted);
  
          if ( li && ri )
            rResult = CGAL_NTS logical_and( CGAL_NTS certified_is_equal(li->x(),ri->x())
                                          , CGAL_NTS certified_is_equal(li->y(),ri->y())
                                          ) ;
        }
        else rResult = make_uncertain(false);
      }
    }
  }

  return rResult;
}

} // namespace CGAL_SS_i

CGAL_END_NAMESPACE

#endif // CGAL_STRAIGHT_SKELETON_PREDICATES_FTC2_H //
// EOF //

