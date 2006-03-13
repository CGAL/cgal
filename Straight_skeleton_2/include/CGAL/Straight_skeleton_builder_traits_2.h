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
#ifndef CGAL_STRAIGHT_SKELETON_BUILDER_TRAITS_2_H
#define CGAL_STRAIGHT_SKELETON_BUILDER_TRAITS_2_H 1

#include <CGAL/Straight_skeleton_aux.h>
#include <CGAL/Straight_skeleton_builder_traits_2_aux.h>
#include <CGAL/predicates/Straight_skeleton_pred_ftC2.h>
#include <CGAL/constructions/Straight_skeleton_cons_ftC2.h>

CGAL_BEGIN_NAMESPACE



namespace CGAL_SS_i {

template<class K>
struct Construct_ss_vertex_2
{
  typedef typename K::FT      FT ;
  typedef typename K::Point_2 Point_2 ;
  
  typedef Vertex<FT> Vertex ;

  typedef Vertex       result_type ;
  typedef Arity_tag<1> Arity ;

  Vertex operator() ( Point_2 const& aP ) const
  {
    K k ;
    
    typename K::Compute_x_2 getx = k.compute_x_2_object();
    typename K::Compute_y_2 gety = k.compute_y_2_object();
    
    return Vertex(getx(aP),gety(aP));      
  }
};

template<class K>
struct Construct_ss_edge_2
{
  typedef typename K::FT      FT ;
  typedef typename K::Point_2 Point_2 ;
  
  typedef Edge<FT> Edge ;

  typedef Edge         result_type ;
  typedef Arity_tag<2> Arity ;

  Edge operator() ( Point_2 const& aS, Point_2 const& aT ) const
  {
    Construct_ss_vertex_2<K> construct_vertex ;
    return Edge(construct_vertex(aS),construct_vertex(aT));      
  }
};

template<class K>
struct Construct_ss_triedge_2
{
  typedef typename K::FT FT ;
  
  typedef Edge   <FT> Edge    ;
  typedef Triedge<FT> Triedge ;

  typedef Triedge      result_type ;
  typedef Arity_tag<3> Arity ;

  Triedge operator() ( Edge const& aA, Edge const& aB, Edge const& aC ) const
  {
    return Triedge(aA,aB,aC);
  }
};

template<class K>
struct Exist_ss_event_2 : Functor_base_2<K>
{
  typedef Functor_base_2<K> Base ;

  typedef typename Base::Triedge Triedge ;

  typedef Uncertain<bool> result_type ;
  typedef Arity_tag<1>    Arity ;

  Uncertain<bool> operator() ( Triedge const& aTriedge ) const
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

  typedef typename Base::Point_2 Point_2 ;
  typedef typename Base::Triedge Triedge ;

  typedef Uncertain<Comparison_result> result_type ;
  typedef Arity_tag<3>                 Arity ;

  Uncertain<Comparison_result> operator() ( Point_2 const& aP
                                          , Triedge const& aL
                                          , Triedge const& aR
                                          ) const
  {
    Construct_ss_vertex_2<K> construct_vertex ;
    return compare_offset_lines_isec_dist_to_pointC2(construct_vertex(aP),aL,aR) ;
  }

  Uncertain<Comparison_result> operator() ( Triedge const& aS
                                          , Triedge const& aL
                                          , Triedge const& aR
                                          ) const
  {
    return compare_offset_lines_isec_dist_to_pointC2(aS,aL,aR) ;
  }

};

template<class K>
struct Compare_ss_event_times_2 : Functor_base_2<K>
{
  typedef Functor_base_2<K> Base ;

  typedef typename Base::Triedge Triedge ;

  typedef Uncertain<Comparison_result> result_type ;
  typedef Arity_tag<2>                 Arity ;

  Uncertain<Comparison_result> operator() ( Triedge const& aL, Triedge const& aR ) const
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

  typedef typename Base::Triedge Triedge ;

  typedef Uncertain<bool> result_type ;
  typedef Arity_tag<2>    Arity ;

  Uncertain<bool> operator() ( Triedge const& aE, Triedge const& aZ ) const
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

  typedef typename Base::Triedge Triedge ;

  typedef Uncertain<bool> result_type ;
  typedef Arity_tag<2>    Arity ;

  Uncertain<bool> operator() ( Triedge const& aA, Triedge const& aB ) const
  {
    Uncertain<bool> rResult = are_events_simultaneousC2(aA,aB);

    CGAL_SLS_ASSERT_PREDICATE_RESULT(rResult,K,"Are_events_simultaneous","A=" << aA << "\nB=" << aB);

    return rResult ;
  }
};


template<class K>
struct Construct_ss_event_time_and_point_2 : Functor_base_2<K>
{
  typedef Functor_base_2<K> Base ;

  typedef typename Base::FT            FT ;
  typedef typename Base::Point_2       Point_2 ;
  typedef typename Base::Vertex        Vertex ;
  typedef typename Base::Edge          Edge ;
  typedef typename Base::Triedge       Triedge ;
  typedef typename Base::SortedTriedge SortedTriedge ;

  typedef boost::tuple<FT,Point_2> result_type ;
  typedef Arity_tag<1>             Arity ;

  boost::tuple<FT,Point_2> operator() ( Triedge const& triedge ) const
  {
    SortedTriedge sorted = collinear_sort(triedge);

    CGAL_assertion(!sorted.is_indeterminate()) ;
    CGAL_assertion(sorted.collinear_count() < 3) ;

    Rational<FT> qt = compute_offset_lines_isec_timeC2(sorted);

    FT t = qt.n() / qt.d() ;

    Vertex i = construct_offset_lines_isecC2(sorted);

    return boost::make_tuple(t,Point_2(i.x(),i.y())) ;
  }
};

} // namespace CGAL_SS_i

template<class K>
struct Straight_skeleton_builder_traits_2_functors
{
  typedef CGAL_SS_i::Exist_ss_event_2                   <K> Exist_ss_event_2 ;
  typedef CGAL_SS_i::Compare_ss_event_times_2           <K> Compare_ss_event_times_2 ;
  typedef CGAL_SS_i::Compare_ss_event_distance_to_seed_2<K> Compare_ss_event_distance_to_seed_2 ;
  typedef CGAL_SS_i::Is_ss_event_inside_offset_zone_2   <K> Is_ss_event_inside_offset_zone_2 ;
  typedef CGAL_SS_i::Are_ss_events_simultaneous_2       <K> Are_ss_events_simultaneous_2 ;
  typedef CGAL_SS_i::Construct_ss_event_time_and_point_2<K> Construct_ss_event_time_and_point_2 ;
  typedef CGAL_SS_i::Construct_ss_edge_2                <K> Construct_ss_edge_2 ;
  typedef CGAL_SS_i::Construct_ss_triedge_2             <K> Construct_ss_triedge_2 ;
} ;

template<class K>
struct Straight_skeleton_builder_traits_2_base
{
  typedef K Kernel ;
  
  typedef typename K::FT      FT ;
  typedef typename K::Point_2 Point_2 ;

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

  typedef Unfiltered_predicate_adaptor<typename Unfiltering::Exist_ss_event_2>
    Exist_ss_event_2 ;

  typedef Unfiltered_predicate_adaptor<typename Unfiltering::Compare_ss_event_times_2>
    Compare_ss_event_times_2 ;

  typedef Unfiltered_predicate_adaptor<typename Unfiltering::Compare_ss_event_distance_to_seed_2>
    Compare_ss_event_distance_to_seed_2 ;
    
  typedef Unfiltered_predicate_adaptor<typename Unfiltering::Is_ss_event_inside_offset_zone_2>
    Is_ss_event_inside_offset_zone_2 ;

  typedef Unfiltered_predicate_adaptor<typename Unfiltering::Are_ss_events_simultaneous_2>
    Are_ss_events_simultaneous_2 ;

  typedef typename Unfiltering::Construct_ss_event_time_and_point_2 Construct_ss_event_time_and_point_2 ;
  typedef typename Unfiltering::Construct_ss_edge_2                 Construct_ss_edge_2 ;
  typedef typename Unfiltering::Construct_ss_triedge_2              Construct_ss_triedge_2 ;

} ;

template<class K>
class Straight_skeleton_builder_traits_2_impl<Tag_true,K> : public Straight_skeleton_builder_traits_2_base<K>
{
  typedef Straight_skeleton_builder_traits_2_functors<typename K::EK> Exact ;
  typedef Straight_skeleton_builder_traits_2_functors<typename K::FK> Filtering ;
  typedef Straight_skeleton_builder_traits_2_functors<K>              Unfiltering ;

  typedef typename K::C2E BaseC2E ;
  typedef typename K::C2F BaseC2F ;

  typedef CGAL_SS_i::Triedge_converter<BaseC2E> C2E ;
  typedef CGAL_SS_i::Triedge_converter<BaseC2F> C2F ;

public:

  typedef Filtered_predicate<typename Exact    ::Exist_ss_event_2
                            ,typename Filtering::Exist_ss_event_2
                            , C2E
                            , C2F
                            >
                            Exist_ss_event_2 ;

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

  typedef typename Unfiltering::Construct_ss_event_time_and_point_2 Construct_ss_event_time_and_point_2 ;
  typedef typename Unfiltering::Construct_ss_edge_2                 Construct_ss_edge_2 ;
  typedef typename Unfiltering::Construct_ss_triedge_2              Construct_ss_triedge_2 ;

} ;

template<class K>
class Straight_skeleton_builder_traits_2
  : public Straight_skeleton_builder_traits_2_impl<typename CGAL_SS_i::Is_filtering_kernel<K>::type, K >
{
} ;

CGAL_STRAIGHT_SKELETON_CREATE_FUNCTOR_ADAPTER(Exist_ss_event_2);
CGAL_STRAIGHT_SKELETON_CREATE_FUNCTOR_ADAPTER(Compare_ss_event_times_2);
CGAL_STRAIGHT_SKELETON_CREATE_FUNCTOR_ADAPTER(Compare_ss_event_distance_to_seed_2);
CGAL_STRAIGHT_SKELETON_CREATE_FUNCTOR_ADAPTER(Is_ss_event_inside_offset_zone_2);
CGAL_STRAIGHT_SKELETON_CREATE_FUNCTOR_ADAPTER(Are_ss_events_simultaneous_2);
CGAL_STRAIGHT_SKELETON_CREATE_FUNCTOR_ADAPTER(Construct_ss_event_time_and_point_2);
CGAL_STRAIGHT_SKELETON_CREATE_FUNCTOR_ADAPTER(Construct_ss_edge_2);
CGAL_STRAIGHT_SKELETON_CREATE_FUNCTOR_ADAPTER(Construct_ss_triedge_2);

CGAL_END_NAMESPACE

#endif // CGAL_STRAIGHT_SKELETON_BUILDER_TRAITS_2_H //
// EOF //
