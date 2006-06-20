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
#ifndef CGAL_STRAIGHT_SKELETON_BUILDER_TRAITS_2_H
#define CGAL_STRAIGHT_SKELETON_BUILDER_TRAITS_2_H 1

#include <CGAL/Straight_skeleton_2/Straight_skeleton_aux.h>
#include <CGAL/Straight_skeleton_2/Straight_skeleton_builder_traits_2_aux.h>
#include <CGAL/predicates/Straight_skeleton_pred_ftC2.h>
#include <CGAL/constructions/Straight_skeleton_cons_ftC2.h>

CGAL_BEGIN_NAMESPACE

namespace CGAL_SS_i {

template<class K>
struct Construct_ss_triedge_2 : Functor_base_2<K>
{
  typedef Functor_base_2<K> Base ;
  
  typedef typename Base::FT        FT ;
  typedef typename Base::Segment_2 Segment_2 ;
  typedef typename Base::Triedge_2 Triedge_2 ;

  typedef Triedge_2    result_type ;
  typedef Arity_tag<3> Arity ;

  Triedge_2 operator() ( Segment_2 const& aA, Segment_2 const& aB, Segment_2 const& aC ) const
  {
    return Triedge_2(aA,aB,aC);
  }
};

template<class K>
struct Do_ss_event_exist_2 : Functor_base_2<K>
{
  typedef Functor_base_2<K> Base ;

  typedef typename Base::Triedge_2 Triedge_2 ;

  typedef Uncertain<bool> result_type ;
  typedef Arity_tag<1>    Arity ;

  Uncertain<bool> operator() ( Triedge_2 const& aTriedge ) const
  {
    Uncertain<bool> rResult = exist_offset_lines_isec2(aTriedge) ;

    CGAL_SLS_ASSERT_PREDICATE_RESULT(rResult,K,"Exist_event",aTriedge);

    return rResult ;
  }
};

template<class K>
struct Compare_ss_event_distance_to_seed_2 : Functor_base_2<K>
{
  typedef Functor_base_2<K> Base ;

  typedef typename Base::Point_2   Point_2 ;
  typedef typename Base::Triedge_2 Triedge_2 ;

  typedef Uncertain<Comparison_result> result_type ;
  typedef Arity_tag<3>                 Arity ;

  Uncertain<Comparison_result> operator() ( Point_2   const& aP
                                          , Triedge_2 const& aL
                                          , Triedge_2 const& aR
                                          ) const
  {
    return compare_offset_lines_isec_dist_to_pointC2(cgal_make_optional(aP),aL,aR) ;
  }

  Uncertain<Comparison_result> operator() ( Triedge_2 const& aS
                                          , Triedge_2 const& aL
                                          , Triedge_2 const& aR
                                          ) const
  {
    return compare_offset_lines_isec_dist_to_pointC2(aS,aL,aR) ;
  }

};

template<class K>
struct Compare_ss_event_times_2 : Functor_base_2<K>
{
  typedef Functor_base_2<K> Base ;

  typedef typename Base::Triedge_2 Triedge_2 ;

  typedef Uncertain<Comparison_result> result_type ;
  typedef Arity_tag<2>                 Arity ;

  Uncertain<Comparison_result> operator() ( Triedge_2 const& aL, Triedge_2 const& aR ) const
  {
    Uncertain<Comparison_result> rResult = compare_offset_lines_isec_timesC2(aL,aR) ;

    CGAL_SLS_ASSERT_PREDICATE_RESULT(rResult,K,"Compare_event_times","L: " << aL << "\nR:" << aR );

    return rResult ;
  }
};

template<class K>
struct Is_ss_event_inside_offset_zone_2 : Functor_base_2<K>
{
  typedef Functor_base_2<K> Base ;

  typedef typename Base::Triedge_2 Triedge_2 ;

  typedef Uncertain<bool> result_type ;
  typedef Arity_tag<2>    Arity ;

  Uncertain<bool> operator() ( Triedge_2 const& aE, Triedge_2 const& aZ ) const
  {
    Uncertain<bool> rResult = is_offset_lines_isec_inside_offset_zoneC2(aE,aZ) ;

    CGAL_SLS_ASSERT_PREDICATE_RESULT(rResult,K,"Is_event_inside_offset_zone","E=" << aE << "\nZ=" << aZ);

    return rResult ;
  }
};

template<class K>
struct Are_ss_events_simultaneous_2 : Functor_base_2<K>
{
  typedef Functor_base_2<K> Base ;

  typedef typename Base::Triedge_2 Triedge_2 ;

  typedef Uncertain<bool> result_type ;
  typedef Arity_tag<2>    Arity ;

  Uncertain<bool> operator() ( Triedge_2 const& aA, Triedge_2 const& aB ) const
  {
    Uncertain<bool> rResult = are_events_simultaneousC2(aA,aB);

    CGAL_SLS_ASSERT_PREDICATE_RESULT(rResult,K,"Are_events_simultaneous","A=" << aA << "\nB=" << aB);

    return rResult ;
  }
};

template<class K>
struct Are_ss_edges_collinear_2 : Functor_base_2<K>
{
  typedef Functor_base_2<K> Base ;

  typedef typename Base::Segment_2 Segment_2 ;

  typedef Uncertain<bool> result_type ;
  typedef Arity_tag<2>    Arity ;

  Uncertain<bool> operator() ( Segment_2 const& aA, Segment_2 const& aB ) const
  {
    Uncertain<bool> rResult = are_edges_collinearC2(aA,aB);

    CGAL_SLS_ASSERT_PREDICATE_RESULT(rResult,K,"Are_ss_edges_collinear","A=" << aA << "\nB=" << aB);

    return rResult ;
  }
};

template<class K>
struct Are_ss_edges_parallel_2 : Functor_base_2<K>
{
  typedef Functor_base_2<K> Base ;

  typedef typename Base::Segment_2 Segment_2 ;

  typedef Uncertain<bool> result_type ;
  typedef Arity_tag<2>    Arity ;

  Uncertain<bool> operator() ( Segment_2 const& aA, Segment_2 const& aB ) const
  {
    Uncertain<bool> rResult = are_edges_parallelC2(aA,aB);

    CGAL_SLS_ASSERT_PREDICATE_RESULT(rResult,K,"Are_ss_edges_parallel","A=" << aA << "\nB=" << aB);

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
  typedef typename Base::Triedge_2        Triedge_2 ;
  typedef typename Base::Sorted_triedge_2 Sorted_triedge_2 ;

  typedef boost::tuple< boost::optional<FT>, boost::optional<Point_2> > result_type ;
  
  typedef Arity_tag<1>  Arity ;

  result_type operator() ( Triedge_2 const& triedge ) const
  {
    optional< Rational<FT> > qt ;

    optional<Point_2> i ;
    
    Sorted_triedge_2 sorted = collinear_sort(triedge);

    FT t(0.0);
    
    if ( !sorted.is_indeterminate() )
    {
      CGAL_assertion(sorted.collinear_count() < 3) ;
  
      qt = compute_offset_lines_isec_timeC2(sorted);
  
      i = construct_offset_lines_isecC2(sorted);
      
      if ( qt )
        t = qt->n() / qt->d() ;
    }
    
    return boost::make_tuple( cgal_make_optional(qt,t), i ) ;
  }
};

} // namespace CGAL_SS_i

template<class K>
struct Straight_skeleton_builder_traits_2_functors
{
  typedef CGAL_SS_i::Do_ss_event_exist_2                <K> Do_ss_event_exist_2 ;
  typedef CGAL_SS_i::Compare_ss_event_times_2           <K> Compare_ss_event_times_2 ;
  typedef CGAL_SS_i::Compare_ss_event_distance_to_seed_2<K> Compare_ss_event_distance_to_seed_2 ;
  typedef CGAL_SS_i::Is_ss_event_inside_offset_zone_2   <K> Is_ss_event_inside_offset_zone_2 ;
  typedef CGAL_SS_i::Are_ss_events_simultaneous_2       <K> Are_ss_events_simultaneous_2 ;
  typedef CGAL_SS_i::Are_ss_edges_parallel_2            <K> Are_ss_edges_parallel_2 ;
  typedef CGAL_SS_i::Are_ss_edges_collinear_2           <K> Are_ss_edges_collinear_2 ;
  typedef CGAL_SS_i::Construct_ss_event_time_and_point_2<K> Construct_ss_event_time_and_point_2 ;
  typedef CGAL_SS_i::Construct_ss_triedge_2             <K> Construct_ss_triedge_2 ;
} ;

template<class K>
struct Straight_skeleton_builder_traits_2_base
{
  typedef K Kernel ;
  
  typedef typename K::FT        FT ;
  typedef typename K::Point_2   Point_2 ;
  typedef typename K::Segment_2 Segment_2 ;
  
  typedef CGAL_SS_i::Triedge_2<K> Triedge_2 ;

  typedef typename K::Equal_2     Equal_2 ;
  typedef typename K::Left_turn_2 Left_turn_2 ;
  typedef typename K::Collinear_2 Collinear_2 ;

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

  typedef Unfiltered_predicate_adaptor<typename Unfiltering::Compare_ss_event_distance_to_seed_2>
    Compare_ss_event_distance_to_seed_2 ;
    
  typedef Unfiltered_predicate_adaptor<typename Unfiltering::Is_ss_event_inside_offset_zone_2>
    Is_ss_event_inside_offset_zone_2 ;

  typedef Unfiltered_predicate_adaptor<typename Unfiltering::Are_ss_events_simultaneous_2>
    Are_ss_events_simultaneous_2 ;
    
  typedef Unfiltered_predicate_adaptor<typename Unfiltering::Are_ss_edges_parallel_2>
    Are_ss_edges_parallel_2 ;
    
  typedef Unfiltered_predicate_adaptor<typename Unfiltering::Are_ss_edges_collinear_2>
    Are_ss_edges_collinear_2 ;

  typedef typename Unfiltering::Construct_ss_event_time_and_point_2 Construct_ss_event_time_and_point_2 ;

  typedef typename Unfiltering::Construct_ss_triedge_2              Construct_ss_triedge_2 ;
} ;

template<class K>
class Straight_skeleton_builder_traits_2_impl<Tag_true,K> : public Straight_skeleton_builder_traits_2_base<K>
{
  typedef typename K::EK EK ;
  typedef typename K::FK FK ;
  
  typedef Straight_skeleton_builder_traits_2_functors<EK> Exact ;
  typedef Straight_skeleton_builder_traits_2_functors<FK> Filtering ;
  typedef Straight_skeleton_builder_traits_2_functors<K>  Unfiltering ;

  typedef Cartesian_converter<K,EK> BaseC2E;
  typedef Cartesian_converter<K,FK> BaseC2F;
  
  //typedef typename K::C2E BaseC2E ;
  //typedef typename K::C2F BaseC2F ;

  typedef CGAL_SS_i::SS_converter<BaseC2E> C2E ;
  typedef CGAL_SS_i::SS_converter<BaseC2F> C2F ;

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

  typedef Filtered_predicate< typename Exact    ::Compare_ss_event_distance_to_seed_2
                            , typename Filtering::Compare_ss_event_distance_to_seed_2
                            , C2E
                            , C2F
                            >
                            Compare_ss_event_distance_to_seed_2 ;
    
  typedef Filtered_predicate< typename Exact    ::Is_ss_event_inside_offset_zone_2
                            , typename Filtering::Is_ss_event_inside_offset_zone_2
                            , C2E
                            , C2F
                            >
                            Is_ss_event_inside_offset_zone_2 ;

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
                            
  typedef typename Unfiltering::Construct_ss_event_time_and_point_2 Construct_ss_event_time_and_point_2 ;
  
  typedef typename Unfiltering::Construct_ss_triedge_2              Construct_ss_triedge_2 ;

} ;

template<class K>
class Straight_skeleton_builder_traits_2
  : public Straight_skeleton_builder_traits_2_impl<typename CGAL_SS_i::Is_filtering_kernel<K>::type, K>
{
} ;

CGAL_STRAIGHT_SKELETON_CREATE_FUNCTOR_ADAPTER(Do_ss_event_exist_2);
CGAL_STRAIGHT_SKELETON_CREATE_FUNCTOR_ADAPTER(Compare_ss_event_times_2);
CGAL_STRAIGHT_SKELETON_CREATE_FUNCTOR_ADAPTER(Compare_ss_event_distance_to_seed_2);
CGAL_STRAIGHT_SKELETON_CREATE_FUNCTOR_ADAPTER(Is_ss_event_inside_offset_zone_2);
CGAL_STRAIGHT_SKELETON_CREATE_FUNCTOR_ADAPTER(Are_ss_events_simultaneous_2);
CGAL_STRAIGHT_SKELETON_CREATE_FUNCTOR_ADAPTER(Are_ss_edges_parallel_2);
CGAL_STRAIGHT_SKELETON_CREATE_FUNCTOR_ADAPTER(Are_ss_edges_collinear_2);
CGAL_STRAIGHT_SKELETON_CREATE_FUNCTOR_ADAPTER(Construct_ss_event_time_and_point_2);
CGAL_STRAIGHT_SKELETON_CREATE_FUNCTOR_ADAPTER(Construct_ss_triedge_2);

CGAL_END_NAMESPACE

#endif // CGAL_STRAIGHT_SKELETON_BUILDER_TRAITS_2_H //
// EOF //
