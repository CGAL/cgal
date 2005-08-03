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
// file          : include/CGAL/Straight_skeleton_builder_traits_2.h
// package       : Straight_skeleton_2 (1.1.0)
//
// author(s)     : Fernando Cacciola
// maintainer    : Fernando Cacciola <fernando_cacciola@hotmail>
// coordinator   : Fernando Cacciola <fernando_cacciola@hotmail>
//
// ============================================================================
#ifndef CGAL_STRAIGHT_SKELETON_BUILDER_TRAITS_2_H
#define CGAL_STRAIGHT_SKELETON_BUILDER_TRAITS_2_H 1

#include <algorithm>

#include <CGAL/constructions/kernel_ftC2.h>
#include <CGAL/constructions/Straight_skeleton_ftC2.h>
#include <CGAL/predicates/Straight_skeleton_ftC2.h>
#include <CGAL/certified_numeric_predicates.h>
#include <CGAL/Unfiltered_predicate_adaptor.h>
#include <CGAL/Filtered_predicate.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Straight_skeleton_aux.h>
#include <CGAL/Threetuple.h>

CGAL_BEGIN_NAMESPACE

template<class K>
struct Exist_event
{
  typedef typename K::FT FT ;
  
  typedef typename K::Point_2 Point_2 ;
    
  typedef std::pair<Point_2,Point_2> Point_2_Pair ;

  typedef Uncertain<bool> result_type ;
  typedef Arity_tag<3>    Arity ;
  
  Uncertain<bool> operator() ( Point_2_Pair const& aM
                             , Point_2_Pair const& aN
                             , Point_2_Pair const& aS
                             ) const
  {
    FT  ma, mb, mc
       ,na, nb, nc
       ,sa, sb, sc ;
   
    line_from_pointsC2(aM.first.x(),aM.first.y(),aM.second.x(),aM.second.y()
                      ,ma,mb,mc
                      );          
                      
    line_from_pointsC2(aN.first.x(),aN.first.y(),aN.second.x(),aN.second.y()
                      ,na,nb,nc
                      );          
    
    line_from_pointsC2(aS.first.x(),aS.first.y(),aS.second.x(),aS.second.y()
                      ,sa,sb,sc
                      );          
    
    return exist_single_point_offset_lines_isec(ma,mb,mc,na,nb,nc,sa,sb,sc) ;
  }                             
};

template<class K>
struct Compare_event_times
{
  typedef typename K::FT FT ;
  
  typedef typename K::Point_2 Point_2 ;
  
  typedef std::pair<Point_2,Point_2> Point_2_Pair ;
  
  typedef Uncertain<Comparison_result> result_type ;
  typedef Arity_tag<6>                 Arity ;
  
  Uncertain<Comparison_result> operator() ( Point_2_Pair const& aL0
                                          , Point_2_Pair const& aL1
                                          , Point_2_Pair const& aL2
                                          , Point_2_Pair const& aR0
                                          , Point_2_Pair const& aR1
                                          , Point_2_Pair const& aR2
                                          ) const
  {
    FT  l0a, l0b, l0c
       ,l1a, l1b, l1c
       ,l2a, l2b, l2c
       ,r0a, r0b, r0c
       ,r1a, r1b, r1c
       ,r2a, r2b, r2c ;
   
    line_from_pointsC2(aL0.first.x(),aL0.first.y(),aL0.second.x(),aL0.second.y()
                      ,l0a,l0b,l0c
                      );          
                      
    line_from_pointsC2(aL1.first.x(),aL1.first.y(),aL1.second.x(),aL1.second.y()
                      ,l1a,l1b,l1c
                      );          
    
    line_from_pointsC2(aL2.first.x(),aL2.first.y(),aL2.second.x(),aL2.second.y()
                      ,l2a,l2b,l2c
                      );          
    
    line_from_pointsC2(aR0.first.x(),aR0.first.y(),aR0.second.x(),aR0.second.y()
                      ,r0a,r0b,r0c
                      );          
                      
    line_from_pointsC2(aR1.first.x(),aR1.first.y(),aR1.second.x(),aR1.second.y()
                      ,r1a,r1b,r1c
                      );          
    
    line_from_pointsC2(aR2.first.x(),aR2.first.y(),aR2.second.x(),aR2.second.y()
                      ,r2a,r2b,r2c
                      );        
   
    return compare_offset_lines_isec_times(l0a,l0b,l0c
                                          ,l1a,l1b,l1c
                                          ,l2a,l2b,l2c                                      
                                          ,r0a,r0b,r0c
                                          ,r1a,r1b,r1c
                                          ,r2a,r2b,r2c
                                          ) ;
  }                             
};

template<class K>
struct Compare_event_distance_to_seed
{
  typedef typename K::FT FT ;
  
  typedef typename K::Point_2 Point_2 ;
  
  typedef std::pair<Point_2,Point_2> Point_2_Pair ;
  
  typedef Threetuple<Point_2_Pair> Line_triple ;
  
  typedef Uncertain<Comparison_result> result_type ;
  typedef Arity_tag<7>                 Arity ;
  
  Uncertain<Comparison_result> operator() ( Point_2     const& aP
                                          , Line_triple const& aL
                                          , Line_triple const& aR
                                         ) const
  {
    FT  l0a, l0b, l0c
       ,l1a, l1b, l1c
       ,l2a, l2b, l2c
       ,r0a, r0b, r0c
       ,r1a, r1b, r1c
       ,r2a, r2b, r2c ;
   
    line_from_pointsC2(aL0.first.x(),aL0.first.y(),aL0.second.x(),aL0.second.y()
                      ,l0a,l0b,l0c
                      );          
                      
    line_from_pointsC2(aL1.first.x(),aL1.first.y(),aL1.second.x(),aL1.second.y()
                      ,l1a,l1b,l1c
                      );          
    
    line_from_pointsC2(aL2.first.x(),aL2.first.y(),aL2.second.x(),aL2.second.y()
                      ,l2a,l2b,l2c
                      );          
    
    line_from_pointsC2(aR0.first.x(),aR0.first.y(),aR0.second.x(),aR0.second.y()
                      ,r0a,r0b,r0c
                      );          
                      
    line_from_pointsC2(aR1.first.x(),aR1.first.y(),aR1.second.x(),aR1.second.y()
                      ,r1a,r1b,r1c
                      );          
    
    line_from_pointsC2(aR2.first.x(),aR2.first.y(),aR2.second.x(),aR2.second.y()
                      ,r2a,r2b,r2c
                      );        
   
    return compare_offset_lines_isec_sdist_to_point(aP.x()
                                                   ,aP.y()
                                                   ,l0a,l0b,l0c
                                                   ,l1a,l1b,l1c
                                                   ,l2a,l2b,l2c
                                                   ,r0a,r0b,r0c
                                                   ,r1a,r1b,r1c
                                                   ,r2a,r2b,r2c
                                                   ) ;
  }                             
  
  Uncertain<Comparison_result> operator() ( Point_2_Pair const& aS0
                                          , Point_2_Pair const& aS1
                                          , Point_2_Pair const& aS2
                                          , Point_2_Pair const& aL0
                                          , Point_2_Pair const& aL1
                                          , Point_2_Pair const& aL2
                                          , Point_2_Pair const& aR0
                                          , Point_2_Pair const& aR1
                                          , Point_2_Pair const& aR2
                                          ) const
  {
    FT  s0a, s0b, s0c
       ,s1a, s1b, s1c
       ,s2a, s2b, s2c
       ,l0a, l0b, l0c
       ,l1a, l1b, l1c
       ,l2a, l2b, l2c
       ,r0a, r0b, r0c
       ,r1a, r1b, r1c
       ,r2a, r2b, r2c ;
   
    line_from_pointsC2(aS0.first.x(),aS0.first.y(),aS0.second.x(),aS0.second.y()
                      ,s0a,s0b,s0c
                      );          
                      
    line_from_pointsC2(aS1.first.x(),aS1.first.y(),aS1.second.x(),aS1.second.y()
                      ,s1a,s1b,s1c
                      );          
                      
    line_from_pointsC2(aS2.first.x(),aS2.first.y(),aS2.second.x(),aS2.second.y()
                      ,s2a,s2b,s2c
                      );          
                      
    line_from_pointsC2(aL0.first.x(),aL0.first.y(),aL0.second.x(),aL0.second.y()
                      ,l0a,l0b,l0c
                      );          
                      
    line_from_pointsC2(aL1.first.x(),aL1.first.y(),aL1.second.x(),aL1.second.y()
                      ,l1a,l1b,l1c
                      );          
    
    line_from_pointsC2(aL2.first.x(),aL2.first.y(),aL2.second.x(),aL2.second.y()
                      ,l2a,l2b,l2c
                      );          
    
    line_from_pointsC2(aR0.first.x(),aR0.first.y(),aR0.second.x(),aR0.second.y()
                      ,r0a,r0b,r0c
                      );          
                      
    line_from_pointsC2(aR1.first.x(),aR1.first.y(),aR1.second.x(),aR1.second.y()
                      ,r1a,r1b,r1c
                      );          
    
    line_from_pointsC2(aR2.first.x(),aR2.first.y(),aR2.second.x(),aR2.second.y()
                      ,r2a,r2b,r2c
                      );        
   
    return compare_offset_lines_isec_sdist_to_point(s0a,s0b,s0c
                                                   ,s1a,s1b,s1c
                                                   ,s2a,s2b,s2c
                                                   ,l0a,l0b,l0c
                                                   ,l1a,l1b,l1c
                                                   ,l2a,l2b,l2c
                                                   ,r0a,r0b,r0c
                                                   ,r1a,r1b,r1c
                                                   ,r2a,r2b,r2c
                                                   ) ;
  }                             
};

template<class K>
struct Is_event_inside_offset_zone
{
  typedef typename K::FT FT ;
  
  typedef typename K::Point_2 Point_2 ;
  
  typedef std::pair<Point_2,Point_2> Point_2_Pair ;
  
  typedef Uncertain<bool> result_type ;
  typedef Arity_tag<5>    Arity ;
  
  Uncertain<bool> operator() ( Point_2_Pair const& aL
                             , Point_2_Pair const& aR
                             , Point_2_Pair const& aE
                             , Point_2_Pair const& aE_Prev
                             , Point_2_Pair const& aE_Next
                             ) const
  {
    FT  la, lb, lc
       ,ra, rb, rc
       ,ea, eb, ec
       ,pa, pb, pc
       ,na, nb, nc ;
   
    line_from_pointsC2(aL.first.x(),aL.first.y(),aL.second.x(),aL.second.y()
                      ,la,lb,lc
                      );          
                      
    line_from_pointsC2(aR.first.x(),aR.first.y(),aR.second.x(),aR.second.y()
                      ,ra,rb,rc
                      );          
    
    line_from_pointsC2(aE.first.x(),aE.first.y(),aE.second.x(),aE.second.y()
                      ,ea,eb,ec
                      );          
    
    line_from_pointsC2(aE_Prev.first.x(),aE_Prev.first.y(),aE_Prev.second.x(),aE_Prev.second.y()
                      ,pa,pb,pc
                      );          
                      
    line_from_pointsC2(aE_Next.first.x(),aE_Next.first.y(),aE_Next.second.x(),aE_Next.second.y()
                      ,na,nb,nc
                      );          
    
    return is_offset_lines_isec_inside_offset_zone(la,lb,lc
                                                  ,ra,rb,rc
                                                  ,ea,eb,ec
                                                  ,pa,pb
                                                  ,ea,eb
                                                  ,na,nb
                                                  ) ;
  }                             
};

template<class K>
struct Construct_event
{
  typedef typename K::FT FT ;
  
  typedef typename K::Point_2 Point_2 ;
  
  typedef std::pair<Point_2,Point_2> Point_2_Pair ;
  
  typedef std::pair<Point_2,FT> result_type ;
  typedef Arity_tag<3>          Arity ;
  
  std::pair<Point_2,FT> operator() ( Point_2_Pair const& aM
                                   , Point_2_Pair const& aN
                                   , Point_2_Pair const& aS
                                   ) const
  {
    FT  ma, mb, mc
       ,na, nb, nc
       ,sa, sb, sc ;
   
    line_from_pointsC2(aM.first.x(),aM.first.y(),aM.second.x(),aM.second.y()
                      ,ma,mb,mc
                      );          
                      
    line_from_pointsC2(aN.first.x(),aN.first.y(),aN.second.x(),aN.second.y()
                      ,na,nb,nc
                      );          
    
    line_from_pointsC2(aS.first.x(),aS.first.y(),aS.second.x(),aS.second.y()
                      ,sa,sb,sc
                      );          
    
    typedef Quotient<FT> QFT ;
                      
    QFT qx,qy ;
      
    QFT qt = construct_offset_lines_isec(ma,mb,mc,na,nb,nc,sa,sb,sc,qx,qy) ;
    
    FT x = qx.numerator() / qx.denominator();
    FT y = qy.numerator() / qy.denominator();
    FT t = qt.numerator() / qt.denominator();
                                                   
    return std::make_pair( Point_2(x,y), t ) ;
  }                             
};

template<class K>
struct Straight_skeleton_builder_traits_2_functors
{
  typedef Exist_event                   <K> Exist_event ;
  typedef Compare_event_times           <K> Compare_event_times ;
  typedef Compare_event_distance_to_seed<K> Compare_event_distance_to_seed ;
  typedef Is_event_inside_offset_zone   <K> Is_event_inside_offset_zone ;
  typedef Construct_event               <K> Construct_event ;
                                  
} ;
template<class K>
class Straight_skeleton_builder_traits_2_base
{
public:
  
  typedef typename K::FT      FT ;
  typedef typename K::Point_2 Point_2 ;
  
  typedef typename K::Left_turn_2 Left_turn_2 ;
  
  template<class F> F get() const { return F(); }
} ;

template<class K>
class Straight_skeleton_builder_traits_2 : public Straight_skeleton_builder_traits_2_base<K>
{
  typedef Straight_skeleton_builder_traits_2_functors<K> Unfiltering ;
  
public:
  
  typedef Unfiltered_predicate_adaptor<typename Unfiltering::Exist_event> 
    Exist_event ;
  
  typedef Unfiltered_predicate_adaptor<typename Unfiltering::Compare_event_times> 
    Compare_event_times ;
    
  typedef Unfiltered_predicate_adaptor<typename Unfiltering::Compare_event_distance_to_seed> 
    Compare_event_distance_to_seed ;  
    
  typedef Unfiltered_predicate_adaptor<typename Unfiltering::Is_event_inside_offset_zone> 
    Is_event_inside_offset_zone ;
  
  typedef typename Unfiltering::Construct_event Construct_event ;
  
} ;

template<>
class Straight_skeleton_builder_traits_2<Exact_predicates_inexact_constructions_kernel>
  :
  public Straight_skeleton_builder_traits_2_base<Exact_predicates_inexact_constructions_kernel>
{

  typedef Exact_predicates_inexact_constructions_kernel K ;
  
  typedef Straight_skeleton_builder_traits_2_functors<K::EK> Exact ;
  typedef Straight_skeleton_builder_traits_2_functors<K::FK> Filtering ;
  typedef Straight_skeleton_builder_traits_2_functors<K>     Unfiltering ;
  
  typedef K::C2E C2E ;
  typedef K::C2F C2F ;
  
public:
  
  
  typedef Filtered_predicate<Exact::Exist_event, Filtering::Exist_event, C2E, C2F> 
    Exist_event ;
  
  typedef Filtered_predicate< Exact::Compare_event_times, Filtering::Compare_event_times, C2E, C2F> 
    Compare_event_times ;
    
  typedef Filtered_predicate< Exact    ::Compare_event_distance_to_seed
                            , Filtering::Compare_event_distance_to_seed
                            , C2E
                            , C2F
                            > 
    Compare_event_distance_to_seed ;  
    
    
  typedef Filtered_predicate< Exact    ::Is_event_inside_offset_zone
                            , Filtering::Is_event_inside_offset_zone 
                            , C2E
                            , C2F
                            > 
    Is_event_inside_offset_zone ;
  
  typedef Unfiltering::Construct_event Construct_event ;
  
  template<class F> F get() const { return F(); }
} ;

CGAL_END_NAMESPACE


#endif // CGAL_STRAIGHT_SKELETON_BUILDER_TRAITS_2_H //
// EOF //
