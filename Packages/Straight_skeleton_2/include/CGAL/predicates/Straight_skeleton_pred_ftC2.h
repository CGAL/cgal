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
// file          : include/CGAL/predicates/Straight_skeleton_ftC2.h
// package       : Straight_skeleton_2 (1.1.0)
//
// author(s)     : Fernando Cacciola
// maintainer    : Fernando Cacciola <fernando_cacciola@hotmail>
// coordinator   : Fernando Cacciola <fernando_cacciola@hotmail>
//
// ============================================================================
#ifndef CGAL_STRAIGHT_SKELETON_PREDICATES_FTC2_H
#define CGAL_STRAIGHT_SKELETON_PREDICATES_FTC2_H 1

#include <CGAL/constructions/Straight_skeleton_cons_ftC2.h>
#include <CGAL/Uncertain.h>
#include <CGAL/certified_quotient_predicates.h>

CGAL_BEGIN_NAMESPACE

// Given 3 oriented lines in _normalized_ implicit form: l0, l1 and l2, returns true if there exist some
// positive offset distance 't' for which the leftward-offsets of the lines intersect at a single point.
template<class FT>
Uncertain<bool>
exist_offset_lines_isec2 ( tuple<FT,FT,FT> const& l0
                         , tuple<FT,FT,FT> const& l1
                         , tuple<FT,FT,FT> const& l2
                         )
{
  FT n,d;
  tie(n,d) = compute_offset_lines_isec_timeC2(l0,l1,l2);

  Uncertain<bool> d_is_zero = CGAL_NTS certified_is_zero(d) ;
  
  if ( is_indeterminate(d_is_zero) )
  {
    CGAL_SSTRAITS_TRACE("\nDenominator is probably zero (but not exactly), event existance is indeterminate." )
    return Uncertain<bool>::indeterminate();
  }
  else if ( d_is_zero )
  {
    CGAL_SSTRAITS_TRACE("\nDenominator exactly zero, Event doesn't exist." )
    return make_uncertain(false);
  }
  else  
  {
    FT t = n/d ;
    bool rExist = CGAL_NTS certified_is_finite(t) && CGAL_NTS certified_is_positive(t) ;

    CGAL_SSTRAITS_TRACE("\nEvent time: " << t << ". Event " << ( rExist ? "exist." : "doesn't exist." ) )
    return rExist ;
  }
}

// Given 2 triples of oriented lines in _normalized_ implicit form: (m0,m1,m2) and (n0,n1,n2), such that
// for each triple there exists distances 'mt' and 'nt' for which the offsets lines (at mt and nt resp.),
// (m0',m1',m2') and (n0',n1',n2') intersect each in a single point; returns the relative order of mt w.r.t nt.
// That is, indicates which offset triple intersects first (closer to the source lines)
// PRECONDITION: There exist distances mt and nt for which each offset triple intersect at a single point.
template<class FT>
Uncertain<Comparison_result>
compare_offset_lines_isec_timesC2 ( tuple<FT,FT,FT> const& m0
                                  , tuple<FT,FT,FT> const& m1
                                  , tuple<FT,FT,FT> const& m2
                                  , tuple<FT,FT,FT> const& n0
                                  , tuple<FT,FT,FT> const& n1
                                  , tuple<FT,FT,FT> const& n2
                                  )
{
  FT mn, md, nn, nd ;
  
  tie(mn,md) = compute_offset_lines_isec_timeC2(m0,m1,m2);
  tie(nn,nd) = compute_offset_lines_isec_timeC2(n0,n1,n2);
           
  typedef Quotient<FT> QFT ;
  QFT mt(mn,md);
  QFT nt(nn,nd);
   
  CGAL_assertion ( CGAL_NTS certified_is_positive(mt) ) ;  
  CGAL_assertion ( CGAL_NTS certified_is_positive(nt) ) ;  
  
  return CGAL_NTS certified_compare(mt,nt);                                                                
                                                                
}

// Given a point (px,py) and 2 triples of oriented lines in _normalized_ implicit form: (m0,m1,m2) and (n0,n1,n2),
// such that their offsets at distances 'mt' and 'nt' intersects in points (mx,my) and (nx,ny),
// returns the relative order order of the distances from (px,py) to (mx,my) and (nx,ny).
// That is, indicates which offset triple intersects closer to (px,py)
// PRECONDITION: There exist single points at which the offset line triples 'm' and 'n' at 'mt' and 'nt' intersect.
template<class FT>
Uncertain<Comparison_result>
compare_offset_lines_isec_sdist_to_pointC2 ( tuple<FT,FT>    const& p
                                           , tuple<FT,FT,FT> const& m0
                                           , tuple<FT,FT,FT> const& m1
                                           , tuple<FT,FT,FT> const& m2
                                           , tuple<FT,FT,FT> const& n0
                                           , tuple<FT,FT,FT> const& n1
                                           , tuple<FT,FT,FT> const& n2
                                         )
{
  FT dm = compute_offset_lines_isec_sdist_to_pointC2(p,m0,m1,m2);
  FT dn = compute_offset_lines_isec_sdist_to_pointC2(p,n0,n1,n2);
                                                               
  return CGAL_NTS certified_compare(dm,dn);                                                                
}

// Given 3 triples of oriented lines in _normalized_ implicit form: (s0,s1,s2), (m0,m1,m2) and (n0,n1,n2),
// such that their offsets at distances 'st', 'mt' and 'nt' intersects in points (sx,sy), (mx,my) and (nx,ny),
// returns the relative order order of the distances from (sx,sy) to (mx,my) and (nx,ny).
// That is, indicates which offset triple intersects closer to (sx,sy)
// PRECONDITION: There exist single points at which the offsets at 'st', 'mt' and 'nt' intersect.
template<class FT>
Uncertain<Comparison_result>
compare_offset_lines_isec_sdist_to_pointC2 ( tuple<FT,FT,FT> const& s0
                                           , tuple<FT,FT,FT> const& s1
                                           , tuple<FT,FT,FT> const& s2
                                           , tuple<FT,FT,FT> const& m0
                                           , tuple<FT,FT,FT> const& m1
                                           , tuple<FT,FT,FT> const& m2
                                           , tuple<FT,FT,FT> const& n0
                                           , tuple<FT,FT,FT> const& n1
                                           , tuple<FT,FT,FT> const& n2
                                           )
{
  return compare_offset_lines_isec_sdist_to_pointC2(construct_offset_lines_isecC2(s0,s1,s2),m0,m1,m2,n0,n1,n2); 
}

// Given a triple of oriented lines in _normalized_ implicit form: (e0,e1,e2) such that their offsets 
// at a distance 't intersects in a point (x,y); and another triple of oriented lines in _normalized_ implicit form: 
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
is_offset_lines_isec_inside_offset_zoneC2 ( tuple<FT,FT,FT> const& e0
                                          , tuple<FT,FT,FT> const& e1
                                          , tuple<FT,FT,FT> const& e2
                                          , tuple<FT,FT,FT> const& el
                                          , tuple<FT,FT,FT> const& ec
                                          , tuple<FT,FT,FT> const& er
                                          )
{
  Uncertain<bool> r = Uncertain<bool>::indeterminate();
  
  FT x, y, ela, elb, elc, eca, ecb, ecc, era, erb, erc ;
  
  tie(ela,elb,elc) = el ;       
  tie(eca,ecb,ecc) = ec ;       
  tie(era,erb,erc) = er ;       
  
  // Construct intersection point (x,y)
  tie(x,y) = construct_offset_lines_isecC2(e0,e1,e2);

  // Calculate scaled (signed) distance from (x,y) to 'ec'
  FT sdc = eca * x + ecb * y + ecc ;
  
  CGAL_SSTRAITS_TRACE("\nsdc=" << sdc ) ;
  
  // NOTE:
  //   if (x,y) is not on the positive side of 'ec' it isn't on it's offset zone.
  //   Also, if (x,y) is over 'ec' (its signed distance to ec is not certainly positive) then by definition is not on its _offset_
  //   zone either.
  Uncertain<bool> cok = CGAL_NTS certified_is_positive(sdc);
  if ( !is_indeterminate(cok) && !!cok )
  {
    CGAL_SSTRAITS_TRACE("\nright side of ec." ) ;
    
    // Calculate scaled (signed) distances from (x,y) to 'el' and 'er'
    FT sdl = ela * x + elb * y + elc ;
    FT sdr = era * x + erb * y + erc ;
  
    CGAL_SSTRAITS_TRACE("\nsdl=" << sdl ) ;
    CGAL_SSTRAITS_TRACE("\nsdr=" << sdr ) ;
    
    // Determine if the vertices (el,ec) and (ec,er) are reflex.
    Uncertain<bool> lcx = CGAL_NTS certified_is_smaller(ela*ecb,eca*elb);
    Uncertain<bool> crx = CGAL_NTS certified_is_smaller(eca*erb,era*ecb);
    if ( !is_indeterminate(lcx) && !is_indeterminate(crx) )
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
      
      r = lok && rok ;                                
    }
    
  }
  else
  {
    CGAL_SSTRAITS_TRACE("\nWRONG side of ec." ) ;
  }
  
  
  
  return r ;
}

CGAL_END_NAMESPACE

#endif // CGAL_STRAIGHT_SKELETON_PREDICATES_FTC2_H //
// EOF //
 
