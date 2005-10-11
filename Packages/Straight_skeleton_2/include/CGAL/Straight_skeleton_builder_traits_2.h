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

#include <CGAL/Straight_skeleton_builder_traits_2_aux.h>
#include <CGAL/predicates/Straight_skeleton_pred_ftC2.h>
#include <CGAL/constructions/Straight_skeleton_cons_ftC2.h>

CGAL_BEGIN_NAMESPACE

namespace CGALi {

template<class K>
struct Exist_sls_event_2 : Sls_functor_base_2<K>
{
  typedef Sls_functor_base_2<K> Base ;

  typedef typename Base::FT          FT ;
  typedef typename Base::Edge_triple Edge_triple ;
  typedef typename Base::LineC2      LineC2 ;

  typedef Uncertain<bool> result_type ;
  typedef Arity_tag<1>    Arity ;

  Uncertain<bool> operator() ( Edge_triple const& aET ) const
  {
    CGAL_SSTRAITS_TRACE("Exist Event:" << aET);
    LineC2 l0,l1,l2;
    tie(l0,l1,l2) = toLineC2_triple(aET);

    return exist_offset_lines_isec2(l0,l1,l2) ;
  }
};


template<class K>
struct Compare_sls_event_times_2 : Sls_functor_base_2<K>
{
  typedef Sls_functor_base_2<K> Base ;

  typedef typename Base::FT          FT ;
  typedef typename Base::Edge_triple Edge_triple ;
  typedef typename Base::LineC2      LineC2 ;

  typedef Uncertain<Comparison_result> result_type ;
  typedef Arity_tag<2>                 Arity ;

  Uncertain<Comparison_result> operator() ( Edge_triple const& aL, Edge_triple const& aR ) const
  {
    LineC2 l0,l1,l2,r0,r1,r2;
    tie(l0,l1,l2) = toLineC2_triple(aL);
    tie(r0,r1,r2) = toLineC2_triple(aR);

    return compare_offset_lines_isec_timesC2(l0,l1,l2,r0,r1,r2) ;
  }
};

template<class K>
struct Compare_sls_event_distance_to_seed_2 : Sls_functor_base_2<K>
{
  typedef Sls_functor_base_2<K> Base ;

  typedef typename Base::Point_2     Point_2 ;
  typedef typename Base::Edge_triple Edge_triple ;
  typedef typename Base::LineC2      LineC2 ;

  typedef Uncertain<Comparison_result> result_type ;
  typedef Arity_tag<3>                 Arity ;

  Uncertain<Comparison_result> operator() ( Point_2     const& aP
                                          , Edge_triple const& aL
                                          , Edge_triple const& aR
                                         ) const
  {
    LineC2 l0,l1,l2,r0,r1,r2;
    tie(l0,l1,l2) = toLineC2_triple(aL);
    tie(r0,r1,r2) = toLineC2_triple(aR);

    return compare_offset_lines_isec_sdist_to_pointC2(toVertexC2(aP),l0,l1,l2,r0,r1,r2) ;
  }

  Uncertain<Comparison_result> operator() ( Edge_triple const& aS
                                          , Edge_triple const& aL
                                          , Edge_triple const& aR
                                          ) const
  {
    LineC2 s0,s1,s2,l0,l1,l2,r0,r1,r2;
    tie(s0,s1,s2) = toLineC2_triple(aS);
    tie(l0,l1,l2) = toLineC2_triple(aL);
    tie(r0,r1,r2) = toLineC2_triple(aR);

    return compare_offset_lines_isec_sdist_to_pointC2(s0,s1,s2,l0,l1,l2,r0,r1,r2) ;
  }

};

template<class K>
struct Is_sls_event_inside_offset_zone_2 : Sls_functor_base_2<K>
{
  typedef Sls_functor_base_2<K> Base ;

  typedef typename Base::Edge_triple Edge_triple ;
  typedef typename Base::LineC2      LineC2 ;

  typedef Uncertain<bool> result_type ;
  typedef Arity_tag<2>    Arity ;

  Uncertain<bool> operator() ( Edge_triple const& aE, Edge_triple const& aO ) const
  {
    LineC2 e0,e1,e2,lo,co,ro ;
    tie(e0,e1,e2) = toLineC2_triple(aE);
    tie(lo,co,ro) = toLineC2_triple(aO);

    return is_offset_lines_isec_inside_offset_zoneC2(e0,e1,e2,lo,co,ro) ;
   }
};

template<class K>
struct Construct_sls_event_time_and_point_2 : Sls_functor_base_2<K>
{
  typedef Sls_functor_base_2<K> Base ;

  typedef typename Base::FT          FT ;
  typedef typename Base::Point_2     Point_2 ;
  typedef typename Base::Edge_triple Edge_triple ;
  typedef typename Base::LineC2      LineC2 ;

  typedef tuple<FT,Point_2> result_type ;
  typedef Arity_tag<1>      Arity ;

  tuple<FT,Point_2> operator() ( Edge_triple const& aE ) const
  {
    LineC2 l0,l1,l2;
    tie(l0,l1,l2) = toLineC2_triple(aE);

    FT tn,td;
    tie(tn,td) = compute_offset_lines_isec_timeC2(l0,l1,l2);

    FT t = tn / td ;

    FT x,y ;
    tie(x,y) = construct_offset_lines_isecC2(l0,l1,l2);

    return make_tuple(t,Point_2(x,y)) ;
  }
};

} // namespace CGALi

template<class K>
struct Straight_skeleton_builder_traits_2_functors
{
  typedef CGALi::Exist_sls_event_2                   <K> Exist_sls_event_2 ;
  typedef CGALi::Compare_sls_event_times_2           <K> Compare_sls_event_times_2 ;
  typedef CGALi::Compare_sls_event_distance_to_seed_2<K> Compare_sls_event_distance_to_seed_2 ;
  typedef CGALi::Is_sls_event_inside_offset_zone_2   <K> Is_sls_event_inside_offset_zone_2 ;
  typedef CGALi::Construct_sls_event_time_and_point_2<K> Construct_sls_event_time_and_point_2 ;

} ;

template<class K>
struct Straight_skeleton_builder_traits_2_base
{
  typedef typename K::FT      FT ;
  typedef typename K::Point_2 Point_2 ;

  typedef typename K::Left_turn_2 Left_turn_2 ;

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

  typedef CGALi::Edge_triple_converter_2<BaseC2E> C2E ;
  typedef CGALi::Edge_triple_converter_2<BaseC2F> C2F ;

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

  typedef typename Unfiltering::Construct_sls_event_time_and_point_2 Construct_sls_event_time_and_point_2 ;

  template<class F> F get() const { return F(); }
} ;

template<class K>
class Straight_skeleton_builder_traits_2
  : public Straight_skeleton_builder_traits_2_impl<typename CGALi::Is_filtering_kernel<K>::type, K >
{
} ;

SLS_CREATE_FUNCTOR_ADAPTER(Exist_sls_event_2);
SLS_CREATE_FUNCTOR_ADAPTER(Compare_sls_event_times_2);
SLS_CREATE_FUNCTOR_ADAPTER(Compare_sls_event_distance_to_seed_2);
SLS_CREATE_FUNCTOR_ADAPTER(Is_sls_event_inside_offset_zone_2);
SLS_CREATE_FUNCTOR_ADAPTER(Construct_sls_event_time_and_point_2);

CGAL_END_NAMESPACE

#endif // CGAL_STRAIGHT_SKELETON_BUILDER_TRAITS_2_H //
// EOF //
