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

#include <CGAL/constructions/Straight_skeleton_ftC2.h>
#include <CGAL/certified_numeric_predicates.h>

CGAL_BEGIN_NAMESPACE

namespace certified {

// Given 3 oriented lines (a0,b0,c0),(a1,b1,c1) and (a2,b2,c2),
// returns 'true' if there exist a positive offset distance at which the 
// leftward-offsetted lines intersect in a single point
//
template<class FT>
optional<bool> 
exist_single_point_offset_lines_isec ( FT const& a0
                                     , FT const& b0
                                     , FT const& c0
                                     , FT const& a1
                                     , FT const& b1
                                     , FT const& c1
                                     , FT const& a2
                                     , FT const& b2
                                     , FT const& c2
                                    )
{
  FT den = (-a2*b1)+(a2*b0)+(b2*a1)-(b2*a0)+(b1*a0)-(b0*a1);
  return ! CGAL_CERTIFIED_NTS is_zero(den) ;
}


//
// Given 2 triples of oriented lines in implicit form, (m0,m1,m2) and (n0,n1,n2), such that
// for each triple there exists distances "mt" and "nt" for which the offsets lines (at t),
// (m0',m1',m2') and (n0',n1',n2') intersect each in a single point; returns the order
// comparison between mt and nt.
// That is, indicates which offset triple intersects first (closer to the source lines)
// PRECONDITION: There exist distances mt and nt for which each offset triple intersect at a single 
// point.
template<class FT>
optional<Comparison_result>
compare_offset_lines_isec_times ( FT const& ma0
                                , FT const& mb0
                                , FT const& mc0
                                , FT const& ma1
                                , FT const& mb1
                                , FT const& mc1
                                , FT const& ma2
                                , FT const& mb2
                                , FT const& mc2
                                , FT const& na0
                                , FT const& nb0
                                , FT const& nc0
                                , FT const& na1
                                , FT const& nb1
                                , FT const& nc1
                                , FT const& na2
                                , FT const& nb2
                                , FT const& nc2
                   )
{
  typedef Quotient<FT> QFT;
  
  QFT mt = compute_offset_lines_isec_time(ma0,mb0,mc0,ma1,mb1,mc1,ma2,mb2,mc2);
  QFT nt = compute_offset_lines_isec_time(na0,nb0,nc0,na1,nb1,nc1,na2,nb2,nc2);
                                                               
  return CGAL_CERTIFIED_NTS compare(mt,nt);                                                                
                                                                
}

//
// Given 2 triples of oriented lines in implicit form, (m0,m1,m2) and (n0,n1,n2), such that
// for each triple there exists distances "mt" and "nt" for which the offsets lines (at t),
// (m0',m1',m2') and (n0',n1',n2') intersect each in a single point; returns the order
// comparison between mt and nt.
// That is, indicates which offset triple intersects first (closer to the source lines)
// PRECONDITION: There exist distances mt and nt for which each offset triple intersect at a single 
// point.
template<class FT>
optional<Comparison_result>
compare_offset_lines_isec_dist_to_point ( FT const& px
                                        , FT const& py
                                        , FT const& ma0
                                        , FT const& mb0
                                        , FT const& mc0
                                        , FT const& ma1
                                        , FT const& mb1
                                        , FT const& mc1
                                        , FT const& ma2
                                        , FT const& mb2
                                        , FT const& mc2
                                        , FT const& na0
                                        , FT const& nb0
                                        , FT const& nc0
                                        , FT const& na1
                                        , FT const& nb1
                                        , FT const& nc1
                                        , FT const& na2
                                        , FT const& nb2
                                        , FT const& nc2
                           )
{
  typedef Quotient<FT> QFT;
  
  QFT dm = compute_offset_lines_isec_sdist_to_point(px,py,ma0,mb0,mc0,ma1,mb1,mc1,ma2,mb2,mc2);
  QFT dn = compute_offset_lines_isec_sdist_to_point(px,py,na0,nb0,nc0,na1,nb1,nc1,na2,nb2,nc2);
                                                               
  return CGAL_CERTIFIED_NTS compare(dm,dn);                                                                
}

//
// Given 2 triples of oriented lines in implicit form, (m0,m1,m2) and (n0,n1,n2), such that
// for each triple there exists distances "mt" and "nt" for which the offsets lines (at t),
// (m0',m1',m2') and (n0',n1',n2') intersect each in a single point; returns the order
// comparison between mt and nt.
// That is, indicates which offset triple intersects first (closer to the source lines)
// PRECONDITION: There exist distances mt and nt for which each offset triple intersect at a single 
// point.
template<class FT>
optional<Comparison_result>
compare_offset_lines_isec_dist_to_point ( FT const& sa0
                                        , FT const& sb0
                                        , FT const& sc0
                                        , FT const& sa1
                                        , FT const& sb1
                                        , FT const& sc1
                                        , FT const& sa2
                                        , FT const& sb2
                                        , FT const& sc2
                                        , FT const& ma0
                                        , FT const& mb0
                                        , FT const& mc0
                                        , FT const& ma1
                                        , FT const& mb1
                                        , FT const& mc1
                                        , FT const& ma2
                                        , FT const& mb2
                                        , FT const& mc2
                                        , FT const& na0
                                        , FT const& nb0
                                        , FT const& nc0
                                        , FT const& na1
                                        , FT const& nb1
                                        , FT const& nc1
                                        , FT const& na2
                                        , FT const& nb2
                                        , FT const& nc2
                           )
{
  typedef Quotient<FT> QFT;

  QFT sx, sy ;
  
  construct_offset_lines_isec(sa0,sb0,sc0,sa1,sb1,sc1,sa2,sb2,sc2,sx,sy);
    
  QFT dm = compute_offset_lines_isec_sdist_to_point(sx,sy,ma0,mb0,mc0,ma1,mb1,mc1,ma2,mb2,mc2);
  QFT dn = compute_offset_lines_isec_sdist_to_point(sx,sy,na0,nb0,nc0,na1,nb1,nc1,na2,nb2,nc2);
                                                               
  return CGAL_CERTIFIED_NTS compare(dm,dn);                                                                
}

template<class FT>
optional<bool>
is_offset_lines_isec_inside_offset_zone ( FT const& a0
                                        , FT const& b0
                                        , FT const& c0
                                        , FT const& a1
                                        , FT const& b1
                                        , FT const& c1
                                        , FT const& a2
                                        , FT const& b2
                                        , FT const& c2
                                        , FT const& al
                                        , FT const& bl
                                        , FT const& ac
                                        , FT const& bc
                                        , FT const& ar
                                        , FT const& br
                                        )
{
  typedef Quotient<FT> QFT;

  QFT ix, iy ;
  
  construct_offset_lines_isec(a0,b0,c0,a1,b1,c1,a2,b2,c2,ix,iy);
                   
  QFT sdl = al * ix + bl * iy ;
  QFT sdc = ac * ix + bc * iy ;
  QFT sde = ar * ix + br * iy ;
  
  return    CGAL_CERTIFIED_NTS compare(sdl, sdc) == SMALLER
         && CGAL_CERTIFIED_NTS compare(sdr, sdc) == SMALLER ;
}

} // namespace certified

CGAL_END_NAMESPACE

#endif // CGAL_STRAIGHT_SKELETON_PREDICATES_FTC2_H //
// EOF //
 
