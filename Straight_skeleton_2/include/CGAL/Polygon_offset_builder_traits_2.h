// Copyright (c) 2005-2008 Fernando Luis Cacciola Carballal. All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Fernando Cacciola <fernando_cacciola@ciudad.com.ar>
#ifndef CGAL_POLYGON_OFFSET_BUILDER_TRAITS_2_H
#define CGAL_POLYGON_OFFSET_BUILDER_TRAITS_2_H 1

#include <CGAL/license/Straight_skeleton_2.h>

#include <CGAL/Straight_skeleton_2/Straight_skeleton_aux.h>
#include <CGAL/Straight_skeleton_2/Straight_skeleton_builder_traits_2_aux.h>
#include <CGAL/Straight_skeleton_builder_traits_2.h>
#include <CGAL/predicates/Polygon_offset_pred_ftC2.h>
#include <CGAL/constructions/Polygon_offset_cons_ftC2.h>

#include <optional>

namespace CGAL {

namespace CGAL_SS_i {

template<class K>
struct Compare_offset_against_event_time_2 : Functor_base_2<K>
{
  typedef Functor_base_2<K> Base ;

  typedef typename Base::FT               FT ;
  typedef typename Base::Trisegment_2_ptr Trisegment_2_ptr ;

  typedef Uncertain<Comparison_result> result_type ;

  Uncertain<Comparison_result> operator() ( FT const& aT, Trisegment_2_ptr const& aE ) const
  {
    return compare_offset_against_isec_timeC2(aT,aE) ;
  }
};

template<class K>
struct Construct_offset_point_2 : Functor_base_2<K>
{
  typedef Functor_base_2<K> Base ;

  typedef typename Base::FT               FT ;
  typedef typename Base::Point_2          Point_2 ;
  typedef typename Base::Segment_2_with_ID Segment_2_with_ID ;
  typedef typename Base::Trisegment_2_ptr Trisegment_2_ptr ;

  typedef std::optional<Point_2> result_type ;

  result_type operator() ( FT                const& aT
                         , Segment_2_with_ID const& aE0
                         , FT                const& aWeight0
                         , Segment_2_with_ID const& aE1
                         , FT                const& aWeight1
                         , Trisegment_2_ptr  const& aNode
                         ) const
  {
    No_caches<K> no_caches;

    result_type p = construct_offset_pointC2(aT,aE0,aWeight0,aE1,aWeight1,aNode, no_caches);

    return p ;
  }
};


} // namespace CGAL_SS_i

template<class K>
struct Polygon_offset_builder_traits_2_functors
{
  typedef CGAL_SS_i::Compare_offset_against_event_time_2<K> Compare_offset_against_event_time_2 ;
  typedef CGAL_SS_i::Compare_ss_event_times_2           <K> Compare_ss_event_times_2 ;
  typedef CGAL_SS_i::Construct_offset_point_2           <K> Construct_offset_point_2 ;
  typedef CGAL_SS_i::Construct_ss_event_time_and_point_2<K> Construct_ss_event_time_and_point_2 ;
} ;

template<class K>
struct Polygon_offset_builder_traits_2_base
{
  typedef K Kernel ;

  typedef typename K::FT        FT ;
  typedef typename K::Point_2   Point_2 ;

  typedef CGAL_SS_i::Segment_2_with_ID<K> Segment_2 ;
  typedef CGAL_SS_i::Segment_2_with_ID<K> Segment_2_with_ID ; // for BOOST_MPL_HAS_XXX_TRAIT_DEF
  typedef CGAL::Trisegment_2<K, Segment_2_with_ID> Trisegment_2 ;
  typedef typename Trisegment_2::Self_ptr Trisegment_2_ptr ;
} ;

template<class Is_filtered_kernel, class K>
class Polygon_offset_builder_traits_2_impl ;

template<class K>
class Polygon_offset_builder_traits_2_impl<Tag_false,K>
  : public Polygon_offset_builder_traits_2_base<K>
{
  typedef Polygon_offset_builder_traits_2_functors<K> Unfiltering ;

public:
  typedef Unfiltered_predicate_adaptor<typename Unfiltering::Compare_offset_against_event_time_2>
    Compare_offset_against_event_time_2 ;

  typedef Unfiltered_predicate_adaptor<typename Unfiltering::Compare_ss_event_times_2>
    Compare_ss_event_times_2 ;

  typedef typename Unfiltering::Construct_offset_point_2            Construct_offset_point_2 ;
  typedef typename Unfiltering::Construct_ss_event_time_and_point_2 Construct_ss_event_time_and_point_2 ;
} ;

template<class K>
class Polygon_offset_builder_traits_2_impl<Tag_true,K>
  : public Polygon_offset_builder_traits_2_base<K>
{
  typedef typename K::Exact_kernel       EK ;
  typedef typename K::Approximate_kernel FK ;

  typedef Polygon_offset_builder_traits_2_functors<EK> Exact ;
  typedef Polygon_offset_builder_traits_2_functors<FK> Filtering ;
  typedef Polygon_offset_builder_traits_2_functors<K>  Unfiltering ;

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
  typedef Filtered_predicate<typename Exact    ::Compare_offset_against_event_time_2
                            ,typename Filtering::Compare_offset_against_event_time_2
                            , C2E
                            , C2F
                            >
                            Compare_offset_against_event_time_2 ;

  typedef Filtered_predicate< typename Exact    ::Compare_ss_event_times_2
                            , typename Filtering::Compare_ss_event_times_2
                            , C2E
                            , C2F
                            >
                            Compare_ss_event_times_2 ;

  typedef CGAL_SS_i::Exceptionless_filtered_construction< typename Unfiltering::Construct_offset_point_2
                                                        , typename Exact      ::Construct_offset_point_2
                                                        , typename Filtering::Construct_offset_point_2
                                                        , C2E
                                                        , C2F
                                                        , E2C
                                                        , F2C
                                                        >
                                                        Construct_offset_point_2 ;

  typedef CGAL_SS_i::Exceptionless_filtered_construction< typename Unfiltering::Construct_ss_event_time_and_point_2
                                                        , typename Exact      ::Construct_ss_event_time_and_point_2
                                                        , typename Filtering::Construct_ss_event_time_and_point_2
                                                        , C2E
                                                        , C2F
                                                        , E2C
                                                        , F2C
                                                        >
                                                        Construct_ss_event_time_and_point_2 ;
} ;

template<class K>
class Polygon_offset_builder_traits_2
  : public Polygon_offset_builder_traits_2_impl<typename CGAL_SS_i::Is_filtering_kernel<K>::type,K>
{
  typedef Polygon_offset_builder_traits_2_impl<typename CGAL_SS_i::Is_filtering_kernel<K>::type,K> Base ;

public:
  typedef typename Base::Compare_offset_against_event_time_2 Compare_offset_against_event_time_2 ;
  typedef typename Base::Compare_ss_event_times_2 Compare_ss_event_times_2 ;
  typedef typename Base::Construct_offset_point_2 Construct_offset_point_2 ;
  typedef typename Base::Construct_ss_event_time_and_point_2 Construct_ss_event_time_and_point_2 ;

  Compare_offset_against_event_time_2 compare_offset_against_event_time_2_object() const
  { return Compare_offset_against_event_time_2() ; }
  Compare_ss_event_times_2 compare_ss_event_times_2_object() const
  { return Compare_ss_event_times_2() ; }
  Construct_offset_point_2 construct_offset_point_2_object() const
  { return Construct_offset_point_2() ; }
  Construct_ss_event_time_and_point_2 construct_ss_event_time_and_point_2_object() const
  { return Construct_ss_event_time_and_point_2() ; }
} ;

} // namespace CGAL

#endif // CGAL_POLYGON_OFFSET_BUILDER_TRAITS_2_H
