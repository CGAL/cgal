// Copyright (c) 2005-2008 Fernando Luis Cacciola Carballal. All rights reserved.
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
#ifndef CGAL_STRAIGHT_SKELETON_BUILDER_TRAITS_2_H
#define CGAL_STRAIGHT_SKELETON_BUILDER_TRAITS_2_H 1

#include <CGAL/license/Straight_skeleton_2.h>


#include <CGAL/Filtered_construction.h>
#include <CGAL/Straight_skeleton_2/Straight_skeleton_aux.h>
#include <CGAL/Straight_skeleton_2/Straight_skeleton_builder_traits_2_aux.h>
#include <CGAL/predicates/Straight_skeleton_pred_ftC2.h>
#include <CGAL/constructions/Straight_skeleton_cons_ftC2.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

namespace CGAL {

namespace CGAL_SS_i {


template<class K>
struct Construct_ss_trisegment_2 : Functor_base_2<K>
{
  typedef Functor_base_2<K> Base ;
  
  typedef typename Base::Segment_2        Segment_2 ;
  typedef typename Base::Trisegment_2     Trisegment_2 ;
  typedef typename Base::Trisegment_2_ptr Trisegment_2_ptr ;
  
  typedef Trisegment_2_ptr result_type ;
  

  result_type operator() () const { return cgal_make_optional( Trisegment_2::null() ) ; }
  
  result_type operator() ( Segment_2 const& aS0, Segment_2 const& aS1, Segment_2 const& aS2 ) const
  {
    return construct_trisegment(aS0,aS1,aS2);
  }
};

template<class K>
struct Do_ss_event_exist_2 : Functor_base_2<K>
{
  typedef Functor_base_2<K> Base ;

  typedef typename Base::FT               FT ;
  typedef typename Base::Trisegment_2_ptr Trisegment_2_ptr ;

  typedef Uncertain<bool> result_type ;

  Uncertain<bool> operator() ( Trisegment_2_ptr const& aTrisegment, boost::optional<FT> aMaxTime ) const
  {
    Uncertain<bool> rResult = exist_offset_lines_isec2(aTrisegment,aMaxTime) ;

    CGAL_STSKEL_ASSERT_PREDICATE_RESULT(rResult,K,"Exist_event",aTrisegment);

    return rResult ;
  }
};

template<class K>
struct Is_edge_facing_ss_node_2 : Functor_base_2<K>
{
  typedef Functor_base_2<K> Base ;

  typedef typename Base::Point_2          Point_2 ;
  typedef typename Base::Segment_2        Segment_2 ;
  typedef typename Base::Trisegment_2_ptr Trisegment_2_ptr ;

  typedef Uncertain<bool> result_type ;

  Uncertain<bool> operator() ( Point_2 const& aContourNode, Segment_2 const& aEdge ) const
  {
    return is_edge_facing_pointC2(cgal_make_optional(aContourNode),aEdge) ;
  }

  Uncertain<bool> operator() ( Trisegment_2_ptr const& aSkeletonNode, Segment_2 const& aEdge ) const
  {
    return is_edge_facing_offset_lines_isecC2(aSkeletonNode,aEdge) ;
  }
};

template<class K>
struct Compare_ss_event_times_2 : Functor_base_2<K>
{
  typedef Functor_base_2<K> Base ;

  typedef typename Base::Trisegment_2_ptr Trisegment_2_ptr ;

  typedef Uncertain<Comparison_result> result_type ;

  Uncertain<Comparison_result> operator() ( Trisegment_2_ptr const& aL, Trisegment_2_ptr const& aR ) const
  {
    Uncertain<Comparison_result> rResult = compare_offset_lines_isec_timesC2(aL,aR) ;

    CGAL_STSKEL_ASSERT_PREDICATE_RESULT(rResult,K,"Compare_event_times","L: " << aL << "\nR:" << aR );

    return rResult ;
  }
};

template<class K>
struct Oriented_side_of_event_point_wrt_bisector_2 : Functor_base_2<K>
{
  typedef Functor_base_2<K> Base ;

  typedef typename Base::Segment_2        Segment_2 ;
  typedef typename Base::Trisegment_2_ptr Trisegment_2_ptr ;
  typedef typename Base::Point_2          Point_2 ;

  typedef Uncertain<Oriented_side> result_type ;

  Uncertain<Oriented_side> operator() ( Trisegment_2_ptr const& aEvent
                                      , Segment_2        const& aE0
                                      , Segment_2        const& aE1
                                      , Trisegment_2_ptr const& aE01Event
                                      , bool                    aE0isPrimary 
                                      ) const
  {
    Uncertain<Oriented_side> rResult = oriented_side_of_event_point_wrt_bisectorC2(aEvent,aE0,aE1,aE01Event,aE0isPrimary) ;

    CGAL_STSKEL_ASSERT_PREDICATE_RESULT(rResult,K,"Oriented_side_of_event_point_wrt_bisector_2","Event=" << aEvent << " E0=" << aE0 << " E1=" << aE1 );

    return rResult ;
  }
};


template<class K>
struct Are_ss_events_simultaneous_2 : Functor_base_2<K>
{
  typedef Functor_base_2<K> Base ;

  typedef typename Base::Trisegment_2_ptr Trisegment_2_ptr ;

  typedef Uncertain<bool> result_type ;

  Uncertain<bool> operator() ( Trisegment_2_ptr const& aA, Trisegment_2_ptr const& aB ) const
  {
    Uncertain<bool> rResult = are_events_simultaneousC2(aA,aB);

    CGAL_STSKEL_ASSERT_PREDICATE_RESULT(rResult,K,"Are_events_simultaneous","A=" << aA << "\nB=" << aB);

    return rResult ;
  }
};

template<class K>
struct Are_ss_edges_collinear_2 : Functor_base_2<K>
{
  typedef Functor_base_2<K> Base ;

  typedef typename Base::Segment_2 Segment_2 ;

  typedef Uncertain<bool> result_type ;

  Uncertain<bool> operator() ( Segment_2 const& aA, Segment_2 const& aB ) const
  {
    Uncertain<bool> rResult = are_edges_collinearC2(aA,aB);

    CGAL_STSKEL_ASSERT_PREDICATE_RESULT(rResult,K,"Are_ss_edges_collinear","A=" << aA << "\nB=" << aB);

    return rResult ;
  }
};

template<class K>
struct Are_ss_edges_parallel_2 : Functor_base_2<K>
{
  typedef Functor_base_2<K> Base ;

  typedef typename Base::Segment_2 Segment_2 ;

  typedef Uncertain<bool> result_type ;

  Uncertain<bool> operator() ( Segment_2 const& aA, Segment_2 const& aB ) const
  {
    Uncertain<bool> rResult = are_edges_parallelC2(aA,aB);

    CGAL_STSKEL_ASSERT_PREDICATE_RESULT(rResult,K,"Are_ss_edges_parallel","A=" << aA << "\nB=" << aB);

    return rResult ;
  }
};


template<class K>
struct Construct_ss_event_time_and_point_2 : Functor_base_2<K>
{
  typedef Functor_base_2<K> Base ;

  typedef typename Base::FT               FT ;
  typedef typename Base::Point_2          Point_2 ;
  typedef typename Base::Segment_2        Segment_2 ;
  typedef typename Base::Trisegment_2_ptr Trisegment_2_ptr ;

  typedef boost::tuple<FT,Point_2> rtype ;
  
  typedef boost::optional<rtype> result_type ;
  

  result_type operator() ( Trisegment_2_ptr const& aTrisegment ) const
  {
    typename Has_inexact_constructions<K>::type has_inexact_constructions
    CGAL_SUNPRO_INITIALIZE( = typename Has_inexact_constructions<K>::type()) ;
    
    return calc(aTrisegment, has_inexact_constructions);
  }
 
  result_type calc ( Trisegment_2_ptr const& aTrisegment
                   , Tag_false               // kernel already has exact constructions
                   ) const
  {
    bool lOK = false ;
    
    FT      t(0) ;
    Point_2 i = ORIGIN ;

    optional< Rational<FT> > ot = compute_offset_lines_isec_timeC2(aTrisegment);

    if ( !!ot && certainly( CGAL_NTS certified_is_not_zero(ot->d()) ) )
    {
      t = ot->n() / ot->d();
      
      optional<Point_2> oi = construct_offset_lines_isecC2(aTrisegment);
      if ( oi )
      {
        i = *oi ;
        CGAL_stskel_intrinsic_test_assertion(!is_point_calculation_clearly_wrong(t,i,aTrisegment));
        lOK = true ;
      } 
    }
    
    CGAL_STSKEL_ASSERT_CONSTRUCTION_RESULT(lOK,K,"Construct_ss_event_time_and_point_2",aTrisegment);

    return cgal_make_optional(lOK,boost::make_tuple(t,i)) ;
  }
  
  result_type calc ( Trisegment_2_ptr const& aTrisegment
                   , Tag_true         // kernel does not provides exact constructions
                   ) const
  {
    bool lOK = false ;
    
    FT      t(0) ;
    Point_2 i = ORIGIN ;

    optional< Rational<FT> > ot = compute_offset_lines_isec_timeC2(aTrisegment);

    if ( !!ot && certainly( CGAL_NTS certified_is_not_zero(ot->d()) ) )
    {
      t = ot->n() / ot->d();
      
      optional<Point_2> oi = construct_offset_lines_isecC2(aTrisegment);
      if ( oi )
      {
        if ( is_point_calculation_clearly_wrong(t,*oi,aTrisegment)) 
        {
          typedef Exact_predicates_exact_constructions_kernel EK ;
    
          typedef Cartesian_converter<K,EK> BaseC2E;
          typedef Cartesian_converter<EK,K> BaseE2C;
    
          SS_converter<BaseC2E> C2E ;  
          SS_converter<BaseE2C> E2C ;  
    
          oi = E2C(construct_offset_lines_isecC2(C2E(aTrisegment)));
          
          if ( oi )
          {
            i   = *oi ;
            CGAL_stskel_intrinsic_test_assertion(!is_point_calculation_clearly_wrong(t,i,aTrisegment));
            lOK = true ;
          }
        }
        else
        {
          i = *oi ;
          lOK = true ;
        }
      }
      
    }
    
    CGAL_STSKEL_ASSERT_CONSTRUCTION_RESULT(lOK,K,"Construct_ss_event_time_and_point_2",aTrisegment);

    return cgal_make_optional(lOK,boost::make_tuple(t,i)) ;
  }
  
  bool is_point_calculation_clearly_wrong( FT const& t, Point_2 const& p, Trisegment_2_ptr const& aTrisegment ) const 
  {
    bool rR = false ;
    
    if ( is_possibly_inexact_time_clearly_not_zero(t) )
    {
      Segment_2 const& e0 = aTrisegment->e0() ; 
      Segment_2 const& e1 = aTrisegment->e1() ;
      Segment_2 const& e2 = aTrisegment->e2() ;
      
      Point_2 const& e0s = e0.source();
      Point_2 const& e0t = e0.target();
      
      Point_2 const& e1s = e1.source();
      Point_2 const& e1t = e1.target();
      
      Point_2 const& e2s = e2.source();
      Point_2 const& e2t = e2.target();
    
      FT const very_short(0.1);
      FT const very_short_squared = CGAL_NTS square(very_short);
      
      FT l0 = squared_distance(e0s,e0t) ;
      FT l1 = squared_distance(e1s,e1t) ;
      FT l2 = squared_distance(e2s,e2t) ;
      
      bool e0_is_not_very_short = l0 > very_short_squared ;
      bool e1_is_not_very_short = l1 > very_short_squared ;
      bool e2_is_not_very_short = l2 > very_short_squared ;
      
      FT d0 = squared_distance_from_point_to_lineC2(p.x(),p.y(),e0s.x(),e0s.y(),e0t.x(),e0t.y()).to_nt();
      FT d1 = squared_distance_from_point_to_lineC2(p.x(),p.y(),e1s.x(),e1s.y(),e1t.x(),e1t.y()).to_nt();
      FT d2 = squared_distance_from_point_to_lineC2(p.x(),p.y(),e2s.x(),e2s.y(),e2t.x(),e2t.y()).to_nt();
    
      FT tt = CGAL_NTS square(t);
      
      bool e0_is_clearly_wrong = e0_is_not_very_short && is_possibly_inexact_distance_clearly_not_equal_to(d0,tt) ;
      bool e1_is_clearly_wrong = e1_is_not_very_short && is_possibly_inexact_distance_clearly_not_equal_to(d1,tt) ;
      bool e2_is_clearly_wrong = e2_is_not_very_short && is_possibly_inexact_distance_clearly_not_equal_to(d2,tt) ;
    
      rR = e0_is_clearly_wrong || e1_is_clearly_wrong || e2_is_clearly_wrong ;
      
      CGAL_stskel_intrinsic_test_trace_if(rR
                                        , "\nSkeleton node point calculation is clearly wrong:"
                                          << "\ntime=" << t << " p=" << p2str(p) << " e0=" << s2str(e0) << " e1=" << s2str(e1) << " e2=" << s2str(e2)
                                          << "\nl0=" << inexact_sqrt(l0) << " l1=" << inexact_sqrt(l1) << " l2=" << inexact_sqrt(l2)
                                          << "\nd0=" << d0 << " d1=" << d1 << " d2=" << d2 << " tt=" << tt 
                                        ) ;
    }
    
    return rR ;
  }
};

} // namespace CGAL_SS_i

template<class K>
struct Straight_skeleton_builder_traits_2_functors
{
  typedef CGAL_SS_i::Do_ss_event_exist_2                        <K> Do_ss_event_exist_2 ;
  typedef CGAL_SS_i::Compare_ss_event_times_2                   <K> Compare_ss_event_times_2 ;
  typedef CGAL_SS_i::Is_edge_facing_ss_node_2                   <K> Is_edge_facing_ss_node_2 ;
  typedef CGAL_SS_i::Oriented_side_of_event_point_wrt_bisector_2<K> Oriented_side_of_event_point_wrt_bisector_2 ;
  typedef CGAL_SS_i::Are_ss_events_simultaneous_2               <K> Are_ss_events_simultaneous_2 ;
  typedef CGAL_SS_i::Are_ss_edges_parallel_2                    <K> Are_ss_edges_parallel_2 ;
  typedef CGAL_SS_i::Are_ss_edges_collinear_2                   <K> Are_ss_edges_collinear_2 ;
  typedef CGAL_SS_i::Construct_ss_event_time_and_point_2        <K> Construct_ss_event_time_and_point_2 ;
  typedef CGAL_SS_i::Construct_ss_trisegment_2                  <K> Construct_ss_trisegment_2 ;
} ;

template<class K>
struct Straight_skeleton_builder_traits_2_base
{
  typedef K Kernel ;
  
  typedef typename K::FT          FT ;
  typedef typename K::Point_2     Point_2 ;
  typedef typename K::Vector_2    Vector_2 ;
  typedef typename K::Direction_2 Direction_2 ;
  typedef typename K::Segment_2   Segment_2 ;
  
  typedef CGAL_SS_i::Trisegment_2<K> Trisegment_2 ;
  
  typedef typename Trisegment_2::Self_ptr Trisegment_2_ptr ;

  template<class F> F get( F const* = 0 ) const { return F(); }
  
} ;


template<class Is_filtered_kernel, class K> class Straight_skeleton_builder_traits_2_impl ;

template<class K>
class Straight_skeleton_builder_traits_2_impl<Tag_false,K> : public Straight_skeleton_builder_traits_2_base<K>
{
  typedef Straight_skeleton_builder_traits_2_functors<K> Unfiltering ;

public:

  typedef Unfiltered_predicate_adaptor<typename Unfiltering::Do_ss_event_exist_2>
    Do_ss_event_exist_2 ;

  typedef Unfiltered_predicate_adaptor<typename Unfiltering::Compare_ss_event_times_2>
    Compare_ss_event_times_2 ;

  typedef Unfiltered_predicate_adaptor<typename Unfiltering::Is_edge_facing_ss_node_2>
     Is_edge_facing_ss_node_2 ;

  typedef Unfiltered_predicate_adaptor<typename Unfiltering::Oriented_side_of_event_point_wrt_bisector_2>
    Oriented_side_of_event_point_wrt_bisector_2 ;

  typedef Unfiltered_predicate_adaptor<typename Unfiltering::Are_ss_events_simultaneous_2>
    Are_ss_events_simultaneous_2 ;
    
  typedef Unfiltered_predicate_adaptor<typename Unfiltering::Are_ss_edges_parallel_2>
    Are_ss_edges_parallel_2 ;
    
  typedef Unfiltered_predicate_adaptor<typename Unfiltering::Are_ss_edges_collinear_2>
    Are_ss_edges_collinear_2 ;
    
  typedef typename Unfiltering::Construct_ss_event_time_and_point_2 Construct_ss_event_time_and_point_2 ;
  typedef typename Unfiltering::Construct_ss_trisegment_2           Construct_ss_trisegment_2 ;
} ;

template<class K>
class Straight_skeleton_builder_traits_2_impl<Tag_true,K> : public Straight_skeleton_builder_traits_2_base<K>
{
  typedef typename K::Exact_kernel       EK ;
  typedef typename K::Approximate_kernel FK ;

  typedef Straight_skeleton_builder_traits_2_functors<EK> Exact ;
  typedef Straight_skeleton_builder_traits_2_functors<FK> Filtering ;
  typedef Straight_skeleton_builder_traits_2_functors<K>  Unfiltering ;

  typedef Cartesian_converter<K,EK> BaseC2E;
  typedef Cartesian_converter<K,FK> BaseC2F;
  typedef Cartesian_converter<EK,K> BaseE2C;
  typedef Cartesian_converter<FK,K> BaseF2C;
  typedef Cartesian_converter<K,K>  BaseC2C;
  
  typedef CGAL_SS_i::SS_converter<BaseC2E> C2E ;
  typedef CGAL_SS_i::SS_converter<BaseC2F> C2F ;
  typedef CGAL_SS_i::SS_converter<BaseE2C> E2C ;
  typedef CGAL_SS_i::SS_converter<BaseF2C> F2C ;
  typedef CGAL_SS_i::SS_converter<BaseC2C> C2C ;

public:

  typedef Filtered_predicate<typename Exact    ::Do_ss_event_exist_2
                            ,typename Filtering::Do_ss_event_exist_2
                            , C2E
                            , C2F
                            >
                            Do_ss_event_exist_2 ;

  typedef Filtered_predicate< typename Exact    ::Compare_ss_event_times_2
                            , typename Filtering::Compare_ss_event_times_2
                            , C2E
                            , C2F
                            >
                            Compare_ss_event_times_2 ;

  typedef Filtered_predicate< typename Exact    ::Is_edge_facing_ss_node_2
                            , typename Filtering::Is_edge_facing_ss_node_2
                            , C2E
                            , C2F
                            >
                            Is_edge_facing_ss_node_2 ;
  
  typedef Filtered_predicate< typename Exact    ::Oriented_side_of_event_point_wrt_bisector_2 
                            , typename Filtering::Oriented_side_of_event_point_wrt_bisector_2
                            , C2E
                            , C2F
                            >
                            Oriented_side_of_event_point_wrt_bisector_2 ;

  typedef Filtered_predicate< typename Exact    ::Are_ss_events_simultaneous_2
                            , typename Filtering::Are_ss_events_simultaneous_2
                            , C2E
                            , C2F
                            >
                            Are_ss_events_simultaneous_2 ;

  typedef Filtered_predicate< typename Exact    ::Are_ss_edges_parallel_2
                            , typename Filtering::Are_ss_edges_parallel_2
                            , C2E
                            , C2F
                            >
                            Are_ss_edges_parallel_2 ;
                            
  typedef Filtered_predicate< typename Exact    ::Are_ss_edges_collinear_2
                            , typename Filtering::Are_ss_edges_collinear_2
                            , C2E
                            , C2F
                            >
                            Are_ss_edges_collinear_2 ;
                            
  typedef CGAL_SS_i::Exceptionless_filtered_construction< typename Unfiltering::Construct_ss_event_time_and_point_2
                                                        , typename Exact      ::Construct_ss_event_time_and_point_2
                                                        , typename Unfiltering::Construct_ss_event_time_and_point_2 
                                                        , C2E
                                                        , C2C 
                                                        , E2C
                                                        , C2C 
                                                        >
                                                        Construct_ss_event_time_and_point_2 ;

  typedef CGAL_SS_i::Exceptionless_filtered_construction< typename Unfiltering::Construct_ss_trisegment_2 
                                                        , typename Exact      ::Construct_ss_trisegment_2 
                                                        , typename Unfiltering::Construct_ss_trisegment_2 
                                                        , C2E
                                                        , C2C 
                                                        , E2C
                                                        , C2C 
                                                        >
                                                        Construct_ss_trisegment_2 ;
} ;

template<class K>
class Straight_skeleton_builder_traits_2
  : public Straight_skeleton_builder_traits_2_impl<typename CGAL_SS_i::Is_filtering_kernel<K>::type, K>
{
} ;

CGAL_STRAIGHT_SKELETON_CREATE_FUNCTOR_ADAPTER(Do_ss_event_exist_2)
CGAL_STRAIGHT_SKELETON_CREATE_FUNCTOR_ADAPTER(Compare_ss_event_times_2)
CGAL_STRAIGHT_SKELETON_CREATE_FUNCTOR_ADAPTER(Is_edge_facing_ss_node_2)
CGAL_STRAIGHT_SKELETON_CREATE_FUNCTOR_ADAPTER(Oriented_side_of_event_point_wrt_bisector_2)
CGAL_STRAIGHT_SKELETON_CREATE_FUNCTOR_ADAPTER(Are_ss_events_simultaneous_2)
CGAL_STRAIGHT_SKELETON_CREATE_FUNCTOR_ADAPTER(Are_ss_edges_parallel_2)
CGAL_STRAIGHT_SKELETON_CREATE_FUNCTOR_ADAPTER(Are_ss_edges_collinear_2)
CGAL_STRAIGHT_SKELETON_CREATE_FUNCTOR_ADAPTER(Construct_ss_event_time_and_point_2)
CGAL_STRAIGHT_SKELETON_CREATE_FUNCTOR_ADAPTER(Construct_ss_trisegment_2)

} // end namespace CGAL

#endif // CGAL_STRAIGHT_SKELETON_BUILDER_TRAITS_2_H //
// EOF //
