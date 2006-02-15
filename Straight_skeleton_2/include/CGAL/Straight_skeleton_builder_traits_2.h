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

#include <CGAL/Straight_skeleton_builder_traits_2_aux.h>
#include <CGAL/predicates/Straight_skeleton_pred_ftC2.h>
#include <CGAL/constructions/Straight_skeleton_cons_ftC2.h>

CGAL_BEGIN_NAMESPACE

namespace CGAL_SLS_i {

template<class K>
struct Exist_sls_event_2 : Sls_functor_base_2<K>
{
  typedef Sls_functor_base_2<K> Base ;

  typedef typename Base::Triedge Triedge ;

  typedef Uncertain<bool> result_type ;
  typedef Arity_tag<1>    Arity ;

  Uncertain<bool> operator() ( Triedge const& aTriedge ) const
  {
    CGAL_SSTRAITS_TRACE("Exist Event:" << aTriedge);
    return exist_offset_lines_isec2(aTriedge) ;
  }
};


template<class K>
struct Compare_sls_event_times_2 : Sls_functor_base_2<K>
{
  typedef Sls_functor_base_2<K> Base ;

  typedef typename Base::Triedge Triedge ;

  typedef Uncertain<Comparison_result> result_type ;
  typedef Arity_tag<2>                 Arity ;

  Uncertain<Comparison_result> operator() ( Triedge const& aL, Triedge const& aR ) const
  {
    return compare_offset_lines_isec_timesC2(aL,aR) ;
  }
};

template<class K>
struct Compare_sls_event_distance_to_seed_2 : Sls_functor_base_2<K>
{
  typedef Sls_functor_base_2<K> Base ;

  typedef typename Base::Point_2 Point_2 ;
  typedef typename Base::Triedge Triedge ;

  typedef Uncertain<Comparison_result> result_type ;
  typedef Arity_tag<3>                 Arity ;

  Uncertain<Comparison_result> operator() ( Point_2 const& aP
                                          , Triedge const& aL
                                          , Triedge const& aR
                                          ) const
  {
    return compare_offset_lines_isec_sdist_to_pointC2(toVertex(aP),aL,aR) ;
  }

  Uncertain<Comparison_result> operator() ( Triedge const& aS
                                          , Triedge const& aL
                                          , Triedge const& aR
                                          ) const
  {
    return compare_offset_lines_isec_sdist_to_pointC2(aS,aL,aR) ;
  }

};

template<class K>
struct Is_sls_event_inside_offset_zone_2 : Sls_functor_base_2<K>
{
  typedef Sls_functor_base_2<K> Base ;

  typedef typename Base::Triedge Triedge ;

  typedef Uncertain<bool> result_type ;
  typedef Arity_tag<2>    Arity ;

  Uncertain<bool> operator() ( Triedge const& aE, Triedge const& aZ ) const
  {
    return is_offset_lines_isec_inside_offset_zoneC2(aE,aZ) ;
   }
};

template<class K>
struct Are_sls_events_simultaneous_2 : Sls_functor_base_2<K>
{
  typedef Sls_functor_base_2<K> Base ;

  typedef typename Base::Triedge Triedge ;

  typedef Uncertain<bool> result_type ;
  typedef Arity_tag<2>    Arity ;

  Uncertain<bool> operator() ( Triedge const& aA, Triedge const& aB ) const
  {
    return are_events_simultaneousC2(aA,aB);
   }
};


template<class K>
struct Construct_sls_event_time_and_point_2 : Sls_functor_base_2<K>
{
  typedef Sls_functor_base_2<K> Base ;

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

    CGAL_assertion(sorted.is_valid()) ;

    Rational<FT> qt = compute_offset_lines_isec_timeC2(sorted);

    FT t = qt.n() / qt.d() ;

    Vertex i = construct_offset_lines_isecC2(sorted);

    return boost::make_tuple(t,Point_2(i.x(),i.y())) ;
  }
};

} // namespace CGAL_SLS_i

template<class K>
struct Straight_skeleton_builder_traits_2_functors
{
  typedef CGAL_SLS_i::Exist_sls_event_2                   <K> Exist_sls_event_2 ;
  typedef CGAL_SLS_i::Compare_sls_event_times_2           <K> Compare_sls_event_times_2 ;
  typedef CGAL_SLS_i::Compare_sls_event_distance_to_seed_2<K> Compare_sls_event_distance_to_seed_2 ;
  typedef CGAL_SLS_i::Is_sls_event_inside_offset_zone_2   <K> Is_sls_event_inside_offset_zone_2 ;
  typedef CGAL_SLS_i::Are_sls_events_simultaneous_2       <K> Are_sls_events_simultaneous_2 ;
  typedef CGAL_SLS_i::Construct_sls_event_time_and_point_2<K> Construct_sls_event_time_and_point_2 ;
} ;

template<class K>
struct Straight_skeleton_builder_traits_2_base
{
  typedef typename K::FT      FT ;
  typedef typename K::Point_2 Point_2 ;

  typedef typename K::Left_turn_2 Left_turn_2 ;
  typedef typename K::Collinear_2 Collinear_2 ;

  template<class F> F get() const { return F(); }
} ;


template<class Is_filtered_kernel, class K> class Straight_skeleton_builder_traits_2_impl ;

template<class K>
class Straight_skeleton_builder_traits_2_impl<Tag_false,K> : public Straight_skeleton_builder_traits_2_base<K>
{
  typedef Straight_skeleton_builder_traits_2_functors<K> Unfiltering ;

public:

  typedef Unfiltered_predicate_adaptor<typename Unfiltering::Exist_sls_event_2>
    Exist_sls_event_2 ;

  typedef Unfiltered_predicate_adaptor<typename Unfiltering::Compare_sls_event_times_2>
    Compare_sls_event_times_2 ;

  typedef Unfiltered_predicate_adaptor<typename Unfiltering::Compare_sls_event_distance_to_seed_2>
    Compare_sls_event_distance_to_seed_2 ;

  typedef Unfiltered_predicate_adaptor<typename Unfiltering::Is_sls_event_inside_offset_zone_2>
    Is_sls_event_inside_offset_zone_2 ;

  typedef Unfiltered_predicate_adaptor<typename Unfiltering::Are_sls_events_simultaneous_2>
    Are_sls_events_simultaneous_2 ;

  typedef typename Unfiltering::Construct_sls_event_time_and_point_2 Construct_sls_event_time_and_point_2 ;

} ;

template<class K>
class Straight_skeleton_builder_traits_2_impl<Tag_true,K> : public Straight_skeleton_builder_traits_2_base<K>
{
  typedef Straight_skeleton_builder_traits_2_functors<typename K::EK> Exact ;
  typedef Straight_skeleton_builder_traits_2_functors<typename K::FK> Filtering ;
  typedef Straight_skeleton_builder_traits_2_functors<K>              Unfiltering ;

  typedef typename K::C2E BaseC2E ;
  typedef typename K::C2F BaseC2F ;

  typedef CGAL_SLS_i::Triedge_converter<BaseC2E> C2E ;
  typedef CGAL_SLS_i::Triedge_converter<BaseC2F> C2F ;

public:


  typedef Filtered_predicate<typename Exact    ::Exist_sls_event_2
                            ,typename Filtering::Exist_sls_event_2
                            , C2E
                            , C2F
                            >
                            Exist_sls_event_2 ;

  typedef Filtered_predicate< typename Exact    ::Compare_sls_event_times_2
                            , typename Filtering::Compare_sls_event_times_2
                            , C2E
                            , C2F
                            >
                            Compare_sls_event_times_2 ;

  typedef Filtered_predicate< typename Exact    ::Compare_sls_event_distance_to_seed_2
                            , typename Filtering::Compare_sls_event_distance_to_seed_2
                            , C2E
                            , C2F
                            >
    Compare_sls_event_distance_to_seed_2 ;


  typedef Filtered_predicate< typename Exact    ::Is_sls_event_inside_offset_zone_2
                            , typename Filtering::Is_sls_event_inside_offset_zone_2
                            , C2E
                            , C2F
                            >
    Is_sls_event_inside_offset_zone_2 ;

  typedef Filtered_predicate< typename Exact    ::Are_sls_events_simultaneous_2
                            , typename Filtering::Are_sls_events_simultaneous_2
                            , C2E
                            , C2F
                            >
    Are_sls_events_simultaneous_2 ;

  typedef typename Unfiltering::Construct_sls_event_time_and_point_2 Construct_sls_event_time_and_point_2 ;

  template<class F> F get() const { return F(); }
} ;

template<class K>
class Straight_skeleton_builder_traits_2
  : public Straight_skeleton_builder_traits_2_impl<typename CGAL_SLS_i::Is_filtering_kernel<K>::type, K >
{
} ;

SLS_CREATE_FUNCTOR_ADAPTER(Exist_sls_event_2);
SLS_CREATE_FUNCTOR_ADAPTER(Compare_sls_event_times_2);
SLS_CREATE_FUNCTOR_ADAPTER(Compare_sls_event_distance_to_seed_2);
SLS_CREATE_FUNCTOR_ADAPTER(Is_sls_event_inside_offset_zone_2);
SLS_CREATE_FUNCTOR_ADAPTER(Construct_sls_event_time_and_point_2);
SLS_CREATE_FUNCTOR_ADAPTER(Are_sls_events_simultaneous_2);

CGAL_END_NAMESPACE

#endif // CGAL_STRAIGHT_SKELETON_BUILDER_TRAITS_2_H //
// EOF //
