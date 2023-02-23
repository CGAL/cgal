// Copyright (c) 2005-2008 Fernando Luis Cacciola Carballal. All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Fernando Cacciola <fernando_cacciola@ciudad.com.ar>
//
#ifndef CGAL_STRAIGHT_SKELETON_BUILDER_TRAITS_2_H
#define CGAL_STRAIGHT_SKELETON_BUILDER_TRAITS_2_H 1

#include <CGAL/license/Straight_skeleton_2.h>

#include <CGAL/constructions/Straight_skeleton_cons_ftC2.h>
#include <CGAL/predicates/Straight_skeleton_pred_ftC2.h>
#include <CGAL/Straight_skeleton_2/Straight_skeleton_aux.h>
#include <CGAL/Straight_skeleton_2/Straight_skeleton_builder_traits_2_aux.h>
#include <CGAL/Trisegment_2.h>

#include <CGAL/Filtered_construction.h>
#include <CGAL/Uncertain.h>

#include <boost/optional/optional.hpp>
#include <boost/tuple/tuple.hpp>

#include <iostream>
#include <vector>

namespace CGAL {

namespace CGAL_SS_i {

template<class K>
struct Construct_ss_trisegment_2 : Functor_base_2<K>
{
  typedef Functor_base_2<K> Base ;

  typedef typename Base::Segment_2_with_ID        Segment_2_with_ID ;
  typedef typename Base::Trisegment_2_ptr Trisegment_2_ptr ;
  typedef typename Base::FT               FT ;

  typedef Trisegment_2_ptr result_type ;

  template <class Traits>
  Construct_ss_trisegment_2(const Traits& traits)
    : mNext_ID(traits.trisegment_ID())
  {}

  result_type operator() ( Segment_2_with_ID const& aS0,
                           FT const& aW0,
                           Segment_2_with_ID const& aS1,
                           FT const& aW1,
                           Segment_2_with_ID const& aS2,
                           FT const& aW2 ) const
  {
    return construct_trisegment(aS0,aW0,aS1,aW1,aS2,aW2,mNext_ID++) ;
  }

  std::size_t& mNext_ID ;
};

template<class K>
struct Do_ss_event_exist_2 : Functor_base_2<K>
{
  typedef Functor_base_2<K> Base ;

  typedef typename Base::FT               FT ;
  typedef typename Base::Trisegment_2_ptr Trisegment_2_ptr ;

  typedef Uncertain<bool> result_type ;

  Do_ss_event_exist_2(Time_cache<K>& aTime_cache, Coeff_cache<K>& aCoeff_cache)
    : mTime_cache(aTime_cache), mCoeff_cache(aCoeff_cache)
  {}

  Uncertain<bool> operator() ( Trisegment_2_ptr const& aTrisegment, boost::optional<FT> aMaxTime ) const
  {
    Uncertain<bool> rResult = exist_offset_lines_isec2(aTrisegment,aMaxTime,mTime_cache,mCoeff_cache) ;

    CGAL_STSKEL_ASSERT_PREDICATE_RESULT(rResult,K,"Exist_event",aTrisegment);

    return rResult ;
  }

private:
  Time_cache<K>& mTime_cache ;
  Coeff_cache<K>& mCoeff_cache ;
};

template<class K>
struct Is_edge_facing_ss_node_2 : Functor_base_2<K>
{
  typedef Functor_base_2<K> Base ;

  typedef typename Base::Point_2          Point_2 ;
  typedef typename Base::Segment_2_with_ID        Segment_2_with_ID ;
  typedef typename Base::Trisegment_2_ptr Trisegment_2_ptr ;

  typedef Uncertain<bool> result_type ;

  Is_edge_facing_ss_node_2(Coeff_cache<K>& aCoeff_cache)
    : mCoeff_cache(aCoeff_cache)
  { }

  Uncertain<bool> operator() ( Point_2 const& aContourNode, Segment_2_with_ID const& aEdge ) const
  {
    return is_edge_facing_pointC2(cgal_make_optional(aContourNode),aEdge) ;
  }

  Uncertain<bool> operator() ( Trisegment_2_ptr const& aSkeletonNode, Segment_2_with_ID const& aEdge ) const
  {
    return is_edge_facing_offset_lines_isecC2(aSkeletonNode,aEdge,mCoeff_cache) ;
  }

private:
  Coeff_cache<K>& mCoeff_cache ;
};

template<class K>
struct Compare_ss_event_times_2 : Functor_base_2<K>
{
  typedef Functor_base_2<K> Base ;

  typedef typename Base::Trisegment_2_ptr Trisegment_2_ptr ;

  typedef Uncertain<Comparison_result> result_type ;

  Compare_ss_event_times_2(Time_cache<K>& aTime_cache, Coeff_cache<K>& aCoeff_cache)
    : mTime_cache(aTime_cache), mCoeff_cache(aCoeff_cache)
  {}

  Uncertain<Comparison_result> operator() ( Trisegment_2_ptr const& aL, Trisegment_2_ptr const& aR ) const
  {
    Uncertain<Comparison_result> rResult = compare_offset_lines_isec_timesC2(aL,aR,mTime_cache,mCoeff_cache) ;
    CGAL_STSKEL_ASSERT_PREDICATE_RESULT(rResult,K,"Compare_event_times","L: " << aL << "\nR:" << aR );
    return rResult ;
  }

private:
  Time_cache<K>& mTime_cache ;
  Coeff_cache<K>& mCoeff_cache ;
};

template<class K>
struct Compare_ss_event_angles_2 : Functor_base_2<K>
{
  typedef Functor_base_2<K> Base ;

  typedef typename Base::Trisegment_2_ptr Trisegment_2_ptr ;
  typedef typename Base::FT FT ;
  typedef typename Base::Vector_2 Vector_2 ;

  typedef Uncertain<Comparison_result> result_type ;

  // Two triedges with common e0 and e1 (BV1 and BV2)
  Uncertain<Comparison_result> operator() ( Vector_2 const& aBV1, Vector_2 const& aBV2,
                                            Vector_2 const& aLV, Vector_2 const& aRV ) const
  {
    Uncertain<Comparison_result> rResult = compare_isec_anglesC2(aBV1,aBV2,aLV,aRV) ;

    CGAL_STSKEL_ASSERT_PREDICATE_RESULT(rResult,K,"Compare_event_angles","Dump, @tmp");

    return rResult ;
  }
};

template<class K>
struct Oriented_side_of_event_point_wrt_bisector_2 : Functor_base_2<K>
{
  typedef Functor_base_2<K> Base ;

  typedef typename Base::Segment_2_with_ID        Segment_2_with_ID ;
  typedef typename Base::Trisegment_2_ptr Trisegment_2_ptr ;
  typedef typename Base::Point_2          Point_2 ;
  typedef typename Base::FT               FT;

  typedef Uncertain<Oriented_side> result_type ;

  Oriented_side_of_event_point_wrt_bisector_2(Coeff_cache<K>& aCoeff_cache)
    : mCoeff_cache(aCoeff_cache)
  {}

  Uncertain<Oriented_side> operator() ( Trisegment_2_ptr const& aEvent
                                      , Segment_2_with_ID        const& aE0
                                      , FT               const& aW0
                                      , Segment_2_with_ID        const& aE1
                                      , FT               const& aW1
                                      , Trisegment_2_ptr const& aE01Event
                                      , bool                    aE0isPrimary
                                      ) const
  {
    Uncertain<Oriented_side> rResult = oriented_side_of_event_point_wrt_bisectorC2(aEvent,aE0,aW0,aE1,aW1,aE01Event,aE0isPrimary,mCoeff_cache) ;

    CGAL_STSKEL_ASSERT_PREDICATE_RESULT(rResult,K,"Oriented_side_of_event_point_wrt_bisector_2","Event=" << aEvent << " E0=" << aE0 << " E1=" << aE1 );

    return rResult ;
  }

private:
  Coeff_cache<K>& mCoeff_cache ;
};


template<class K>
struct Are_ss_events_simultaneous_2 : Functor_base_2<K>
{
  typedef Functor_base_2<K> Base ;

  typedef typename Base::Trisegment_2_ptr Trisegment_2_ptr ;

  typedef Uncertain<bool> result_type ;

  Are_ss_events_simultaneous_2(Time_cache<K>& aTime_cache, Coeff_cache<K>& aCoeff_cache)
    : mTime_cache(aTime_cache), mCoeff_cache(aCoeff_cache)
  {}

  Uncertain<bool> operator() ( Trisegment_2_ptr const& aA, Trisegment_2_ptr const& aB ) const
  {
    Uncertain<bool> rResult = are_events_simultaneousC2(aA,aB,mTime_cache,mCoeff_cache);

    CGAL_STSKEL_ASSERT_PREDICATE_RESULT(rResult,K,"Are_events_simultaneous","A=" << aA << "\nB=" << aB);

    return rResult ;
  }

private:
  Time_cache<K>& mTime_cache ;
  Coeff_cache<K>& mCoeff_cache ;
};

// Not actually in use
template<class K>
struct Are_ss_edges_collinear_2 : Functor_base_2<K>
{
  typedef Functor_base_2<K> Base ;

  typedef typename Base::Segment_2_with_ID Segment_2_with_ID ;

  typedef Uncertain<bool> result_type ;

  Uncertain<bool> operator() ( Segment_2_with_ID const& aA, Segment_2_with_ID const& aB ) const
  {
    Uncertain<bool> rResult = are_edges_collinearC2(aA,aB);

    CGAL_STSKEL_ASSERT_PREDICATE_RESULT(rResult,K,"Are_ss_edges_collinear","A=" << aA << "\nB=" << aB);

    return rResult ;
  }
};

// Not actually in use
template<class K>
struct Are_ss_edges_parallel_2 : Functor_base_2<K>
{
  typedef Functor_base_2<K> Base ;

  typedef typename Base::Segment_2_with_ID Segment_2_with_ID ;

  typedef Uncertain<bool> result_type ;

  Uncertain<bool> operator() ( Segment_2_with_ID const& aA, Segment_2_with_ID const& aB ) const
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
  typedef typename Base::Segment_2_with_ID        Segment_2_with_ID ;
  typedef typename Base::Trisegment_2_ptr Trisegment_2_ptr ;

  typedef boost::tuple<FT,Point_2> rtype ;

  typedef boost::optional<rtype> result_type ;

  Construct_ss_event_time_and_point_2(Time_cache<K>& aTime_cache, Coeff_cache<K>& aCoeff_cache)
    : mTime_cache(aTime_cache), mCoeff_cache(aCoeff_cache)
  {}

  result_type operator() ( Trisegment_2_ptr const& aTrisegment ) const
  {
    bool lOK = false ;

    FT      t(0) ;
    Point_2 i = ORIGIN ;

    boost::optional< Rational<FT> > ot = compute_offset_lines_isec_timeC2(aTrisegment,mTime_cache,mCoeff_cache);

    if ( !!ot && certainly( CGAL_NTS certified_is_not_zero(ot->d()) ) )
    {
      t = ot->n() / ot->d();

      boost::optional<Point_2> oi = construct_offset_lines_isecC2(aTrisegment,mCoeff_cache);
      if ( oi )
      {
        i = *oi ;
        lOK = true ;
      }
    }

    CGAL_STSKEL_ASSERT_CONSTRUCTION_RESULT(lOK,K,"Construct_ss_event_time_and_point_2",aTrisegment);

    return cgal_make_optional(lOK,boost::make_tuple(t,i)) ;
  }

private:
  Time_cache<K>& mTime_cache ;
  Coeff_cache<K>& mCoeff_cache ;
};

} // namespace CGAL_SS_i

template<class K>
struct Straight_skeleton_builder_traits_2_functors
{
  typedef CGAL_SS_i::Do_ss_event_exist_2                        <K> Do_ss_event_exist_2 ;
  typedef CGAL_SS_i::Compare_ss_event_times_2                   <K> Compare_ss_event_times_2 ;
  typedef CGAL_SS_i::Compare_ss_event_angles_2                  <K> Compare_ss_event_angles_2 ;
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

  typedef CGAL_SS_i::Segment_2_with_ID<K> Segment_2 ;
  typedef CGAL_SS_i::Segment_2_with_ID<K> Segment_2_with_ID ; // for BOOST_MPL_HAS_XXX_TRAIT_DEF
  typedef CGAL::Trisegment_2<K, Segment_2_with_ID> Trisegment_2 ;
  typedef typename Trisegment_2::Self_ptr Trisegment_2_ptr ;

  template<class F> F get( F const* = 0 ) const { return F(); }

} ;


template<class Is_filtered_kernel, class K>
class Straight_skeleton_builder_traits_2_impl ;

template<class K>
class Straight_skeleton_builder_traits_2_impl<Tag_false /*Is_filtered_kernel*/, K>
  : public Straight_skeleton_builder_traits_2_base<K>
{
  typedef Straight_skeleton_builder_traits_2_functors<K> Unfiltering ;
  typedef Straight_skeleton_builder_traits_2_base<K>     Base ;

public:

  struct Protector {};

  typedef Unfiltered_predicate_adaptor<typename Unfiltering::Do_ss_event_exist_2>
    Do_ss_event_exist_2 ;

  typedef Unfiltered_predicate_adaptor<typename Unfiltering::Compare_ss_event_times_2>
    Compare_ss_event_times_2 ;

  typedef Unfiltered_predicate_adaptor<typename Unfiltering::Compare_ss_event_angles_2>
    Compare_ss_event_angles_2 ;

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

  using Straight_skeleton_builder_traits_2_base<K>::get;

  Is_edge_facing_ss_node_2
  get(Is_edge_facing_ss_node_2 const* = 0 ) const
  {
    return Is_edge_facing_ss_node_2(mCoeff_cache);
  }

  Compare_ss_event_times_2
  get(Compare_ss_event_times_2 const* = 0 ) const
  {
    return Compare_ss_event_times_2(mTime_cache, mCoeff_cache);
  }


  Oriented_side_of_event_point_wrt_bisector_2
  get(Oriented_side_of_event_point_wrt_bisector_2 const* = 0 ) const
  {
    return Oriented_side_of_event_point_wrt_bisector_2(mCoeff_cache);
  }

  Do_ss_event_exist_2
  get(Do_ss_event_exist_2 const* = 0 ) const
  {
    return Do_ss_event_exist_2(mTime_cache, mCoeff_cache);
  }

  Are_ss_events_simultaneous_2
  get(Are_ss_events_simultaneous_2 const* = 0 ) const
  {
    return Are_ss_events_simultaneous_2(mTime_cache, mCoeff_cache);
  }

  Construct_ss_event_time_and_point_2
  get(Construct_ss_event_time_and_point_2 const* = 0 ) const
  {
    return Construct_ss_event_time_and_point_2(mTime_cache, mCoeff_cache);
  }

  Construct_ss_trisegment_2
  get( Construct_ss_trisegment_2 const* = 0 ) const
  {
    return Construct_ss_trisegment_2(*this);
  }

// ID functions
  std::size_t& trisegment_ID() const { return mTrisegment_ID ; }

  void reset_trisegment(std::size_t i) const
  {
    --mTrisegment_ID ;
    mTime_cache.Reset(i) ;
  }

// functions to initialize (and harmonize) and cache speeds
  void InitializeLineCoeffs ( CGAL_SS_i::Segment_2_with_ID<K> const& aBorderS )
  {
    CGAL_SS_i::compute_normalized_line_coeffC2( aBorderS, mCoeff_cache ) ;
  }

  void InitializeLineCoeffs ( std::size_t aID, std::size_t aOtherID )
  {
    if ( mCoeff_cache.Get( aOtherID ) )
      mCoeff_cache.Set( aID, CGAL_SS_i::cgal_make_optional(*(mCoeff_cache.Get(aOtherID))) ) ;
    else
      mCoeff_cache.Set( aID, boost::none ) ;
  }

  // functions and tag for filtering split events
  struct Filters_split_events_tag{};

  template <class EventPtr>
  bool CanSafelyIgnoreSplitEvent(const EventPtr& lEvent) const
  {
    // filter event
    if ( ! mFilteringBound )
      return false;

    typename Base::Trisegment_2_ptr tri = lEvent->trisegment() ;
    boost::optional<CGAL_SS_i::Rational<typename K::FT> > lOptTime =
        CGAL_SS_i::compute_offset_lines_isec_timeC2(tri, mTime_cache, mCoeff_cache);

    if ( lOptTime && lOptTime->to_nt() > *mFilteringBound )
    {
      // avoid filling the cache vectors with times of trisegments that will be removed
      reset_trisegment(tri->id());
      return true;
    }

    return false;
  }

  // @todo there shouldn't be any combinatorial structures such as vertices in the traits
  template <class Vertex_handle, class Halfedge_handle_vector_iterator>
  void ComputeFilteringBound(Vertex_handle aNode,
                             Halfedge_handle_vector_iterator contour_halfedges_begin,
                             Halfedge_handle_vector_iterator contour_halfedges_end) const
  {
    typedef typename K::FT FT;
    typedef typename K::Vector_2 Vector_2;
    typedef typename K::Segment_2 Segment_2;
    typedef typename K::Ray_2 Ray_2;
    typedef typename K::Line_2 Line_2;

    mFilteringBound = boost::none;

    // get the contour input segments on each side of the bisector spawned ataNode
    auto lHL = aNode->halfedge()->defining_contour_edge();
    auto lHR = ( aNode->is_contour() ) ? lHL->opposite()->prev()->opposite()
                                       : aNode->halfedge()->opposite()->defining_contour_edge() ;

    CGAL_SS_i::Segment_2_with_ID<K> lSL (lHL->opposite()->vertex()->point(),
                                         lHL->vertex()->point(),
                                         lHL->id());
    CGAL_SS_i::Segment_2_with_ID<K> lSR (lHR->opposite()->vertex()->point(),
                                         lHR->vertex()->point(),
                                         lHR->id());

    boost::optional< Line_2 > lL = CGAL_SS_i::compute_weighted_line_coeffC2(lSL, FT(1)/lHL->weight(), mCoeff_cache);
    boost::optional< Line_2 > lR = CGAL_SS_i::compute_weighted_line_coeffC2(lSR, FT(1)/lHR->weight(), mCoeff_cache);

    // @fixme below needs to use inverted weights like in degenerate time/point computations
    Vector_2 lVL(lL->a(), lL->b()) ;
    Vector_2 lVR(lR->a(), lR->b()) ;
    Vector_2 lVLR = lVL + lVR ;
    Ray_2 bisect_ray(aNode->point(), lVLR) ;

    // @todo this should use spatial searching
    for ( Halfedge_handle_vector_iterator h = contour_halfedges_begin; h != contour_halfedges_end; ++h )
    {
      CGAL_assertion((*h)->vertex()->is_contour() && (*h)->opposite()->vertex()->is_contour() );

      // @todo could be a line as long as we are in a convex area
      Segment_2 s_h((*h)->opposite()->vertex()->point(), (*h)->vertex()->point());

      // we use segments of the input polygon intersected by the bisector and such that
      // they are oriented such that the reflex vertex is on the left side of the segment
      auto orient = K().orientation_2_object()(s_h[0], s_h[1], aNode->point());
      if (!is_certain(orient) || orient != LEFT_TURN)
        continue;

      auto inter = K().do_intersect_2_object()(s_h, bisect_ray);
      if (!is_certain(inter) || !inter)
        continue;

      CGAL_SS_i::Segment_2_with_ID<K> lSh (s_h, (*h)->id());
      boost::optional< Line_2 > lh = CGAL_SS_i::compute_weighted_line_coeffC2(lSh, (*h)->weight(), mCoeff_cache);

      FT lBound = ( - lh->c() - lh->a()*aNode->point().x() - lh->b()*aNode->point().y() ) /
                    ( lh->a()*lVLR.x() + lh->b()*lVLR.y() ) - aNode->time() ;

      if( ! is_finite(lBound) || ! is_positive(lBound) )
        continue;

      if ( ! mFilteringBound || *mFilteringBound > lBound )
        mFilteringBound = lBound;
    }
  }

public:
  mutable std::size_t mTrisegment_ID = 0 ;
  mutable CGAL_SS_i::Time_cache<K> mTime_cache ;
  mutable CGAL_SS_i::Coeff_cache<K> mCoeff_cache ;
  mutable boost::optional< typename K::FT > mFilteringBound ;
} ;

template<class K>
class Straight_skeleton_builder_traits_2_impl<Tag_true /*Is_filtered_kernel*/, K>
  : public Straight_skeleton_builder_traits_2_base<K>
{
  typedef typename K::Exact_kernel       EK ;
  typedef typename K::Approximate_kernel FK ;

  typedef Straight_skeleton_builder_traits_2_impl<Tag_false, EK> Exact ;
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

  typedef typename FK::FT::Protector Protector;

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

  typedef Filtered_predicate< typename Exact    ::Compare_ss_event_angles_2
                            , typename Filtering::Compare_ss_event_angles_2
                            , C2E
                            , C2F
                            >
                            Compare_ss_event_angles_2 ;

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
                                                        , typename Filtering  ::Construct_ss_event_time_and_point_2
                                                        , C2E
                                                        , C2F
                                                        , E2C
                                                        , F2C
                                                        >
                                                        Construct_ss_event_time_and_point_2 ;

  typedef typename Unfiltering::Construct_ss_trisegment_2        Construct_ss_trisegment_2 ;

  using Straight_skeleton_builder_traits_2_base<K>::get;

// constructor of predicates using time caching
  Compare_ss_event_times_2
  get(Compare_ss_event_times_2 const* = 0 ) const
  {
    return Compare_ss_event_times_2( typename Exact::Compare_ss_event_times_2(
                                       mExact_traits.mTime_cache, mExact_traits.mCoeff_cache),
                                     typename Filtering::Compare_ss_event_times_2(
                                       mApproximate_traits.mTime_cache, mApproximate_traits.mCoeff_cache) );
  }

  Is_edge_facing_ss_node_2
  get(Is_edge_facing_ss_node_2 const* = 0 ) const
  {
    return Is_edge_facing_ss_node_2( typename Exact::Is_edge_facing_ss_node_2(mExact_traits.mCoeff_cache),
                                     typename Filtering::Is_edge_facing_ss_node_2(mApproximate_traits.mCoeff_cache) );
  }

  Oriented_side_of_event_point_wrt_bisector_2
  get(Oriented_side_of_event_point_wrt_bisector_2 const* = 0 ) const
  {
    return Oriented_side_of_event_point_wrt_bisector_2( typename Exact::Oriented_side_of_event_point_wrt_bisector_2(
                                                          mExact_traits.mCoeff_cache),
                                                        typename Filtering::Oriented_side_of_event_point_wrt_bisector_2(
                                                          mApproximate_traits.mCoeff_cache) );
  }

  Do_ss_event_exist_2
  get(Do_ss_event_exist_2 const* = 0 ) const
  {
    return Do_ss_event_exist_2( typename Exact::Do_ss_event_exist_2(
                                  mExact_traits.mTime_cache, mExact_traits.mCoeff_cache),
                                typename Filtering::Do_ss_event_exist_2(
                                  mApproximate_traits.mTime_cache, mApproximate_traits.mCoeff_cache) );
  }

  Are_ss_events_simultaneous_2
  get(Are_ss_events_simultaneous_2 const* = 0 ) const
  {
    return Are_ss_events_simultaneous_2( typename Exact::Are_ss_events_simultaneous_2(
                                           mExact_traits.mTime_cache, mExact_traits.mCoeff_cache),
                                         typename Filtering::Are_ss_events_simultaneous_2(
                                           mApproximate_traits.mTime_cache, mApproximate_traits.mCoeff_cache) );
  }

  Construct_ss_event_time_and_point_2
  get(Construct_ss_event_time_and_point_2 const* = 0 ) const
  {
    return Construct_ss_event_time_and_point_2( typename Exact::Construct_ss_event_time_and_point_2(
                                                  mExact_traits.mTime_cache, mExact_traits.mCoeff_cache),
                                                typename Filtering::Construct_ss_event_time_and_point_2(
                                                  mApproximate_traits.mTime_cache, mApproximate_traits.mCoeff_cache) );
  }

// constructor of trisegments using global id stored in the traits
  Construct_ss_trisegment_2
  get( Construct_ss_trisegment_2 const* = 0 ) const
  {
    return Construct_ss_trisegment_2(*this);
  }

// ID functions
  std::size_t& trisegment_ID() const { return mApproximate_traits.mTrisegment_ID ; }

  void reset_trisegment(std::size_t i) const
  {
    if ( i+1 == trisegment_ID() )
    {
      --trisegment_ID() ;
      mExact_traits.mTime_cache.Reset(i) ;
      mApproximate_traits.mTime_cache.Reset(i) ;
    }
  }

// functions to initialize (and harmonize) and cache speeds
  void InitializeLineCoeffs ( CGAL_SS_i::Segment_2_with_ID<K> const& aBorderS )
  {
    C2E lToExact ;
    C2F lToFiltered ;

    mApproximate_traits.InitializeLineCoeffs( lToFiltered(aBorderS) );
    mExact_traits.InitializeLineCoeffs( lToExact(aBorderS) );
  }

  // This overload copies the coefficients from the halfedge `aOtherID` and caches them for the halfedge `aID`
  void InitializeLineCoeffs ( std::size_t aID, std::size_t aOtherID )
  {
    mApproximate_traits.InitializeLineCoeffs(aID, aOtherID) ;
    mExact_traits.InitializeLineCoeffs(aID, aOtherID) ;
  }

  // functions and tag for filtering split events
  struct Filters_split_events_tag{};

  // The kernel is filtered, and only the filtered part is used in the check (to avoid computing
  // exact stuff and losing time in a check that is there to gain speed)
  template <class EventPtr>
  bool CanSafelyIgnoreSplitEvent(const EventPtr& lEvent) const
  {
    // filter event
    if ( ! mApproximate_traits.mFilteringBound )
      return false;

    typename FK::FT::Protector p;

    typedef CGAL::Trisegment_2<K, CGAL_SS_i::Segment_2_with_ID<K> > Source_trisegment_2 ;
    typedef typename Source_trisegment_2::Self_ptr Source_trisegment_2_ptr;
    typedef CGAL::Trisegment_2<FK, CGAL_SS_i::Segment_2_with_ID<FK> > Target_trisegment_2 ;
    typedef typename Target_trisegment_2::Self_ptr Target_trisegment_2_ptr;

    C2F to_FK ;
    Source_trisegment_2_ptr src_tri = lEvent->trisegment() ;
    CGAL_assertion( src_tri != nullptr ) ;
    Target_trisegment_2_ptr tri = to_FK.cvt_single_trisegment(src_tri) ;

    CGAL_postcondition( src_tri->collinearity() == tri->collinearity() ) ;
    CGAL_postcondition( src_tri->id() == tri->id() ) ;

    try
    {
      boost::optional<CGAL_SS_i::Rational<typename FK::FT> > lOptTime =
          CGAL_SS_i::compute_offset_lines_isec_timeC2(
            tri, mApproximate_traits.mTime_cache, mApproximate_traits.mCoeff_cache);

      if ( lOptTime && lOptTime->to_nt() > *mApproximate_traits.mFilteringBound )
      {
        // avoid filling the cache vectors with times of trisegments that will be removed
        reset_trisegment(tri->id());
        return true;
      }
    }
    catch(Uncertain_conversion_exception&)
    {}

    return false;
  }

  // @todo there shouldn't be any combinatorial structures such as vertices in the traits
  template <class Vertex_handle, class Halfedge_handle_vector_iterator>
  void ComputeFilteringBound(Vertex_handle aNode,
                             Halfedge_handle_vector_iterator contour_halfedges_begin,
                             Halfedge_handle_vector_iterator contour_halfedges_end) const
  {
    typedef typename FK::FT Target_FT;
    typedef typename FK::Point_2 Target_Point_2;
    typedef typename FK::Vector_2 Target_Vector_2;
    typedef typename FK::Segment_2 Target_Segment_2;
    typedef typename FK::Ray_2 Target_Ray_2;
    typedef typename FK::Line_2 Target_Line_2;

    mApproximate_traits.mFilteringBound = boost::none;

    typename FK::FT::Protector protector;

    C2F lToFiltered;

    // get the contour input segments on each side of the bisector spawned ataNode
    auto lHL = aNode->halfedge()->defining_contour_edge();
    auto lHR = ( aNode->is_contour() ) ? lHL->opposite()->prev()->opposite()
                                       : aNode->halfedge()->opposite()->defining_contour_edge() ;

    CGAL_SS_i::Segment_2_with_ID<FK> lSL (lToFiltered(lHL->opposite()->vertex()->point()),
                                          lToFiltered(lHL->vertex()->point()),
                                          lHL->id());
    CGAL_SS_i::Segment_2_with_ID<FK> lSR (lToFiltered(lHR->opposite()->vertex()->point()),
                                          lToFiltered(lHR->vertex()->point()),
                                          lHR->id());

    boost::optional<Target_Line_2> lL = CGAL_SS_i::compute_weighted_line_coeffC2(lSL, Target_FT(1) / lToFiltered(lHL->weight()), mApproximate_traits.mCoeff_cache);
    boost::optional<Target_Line_2> lR = CGAL_SS_i::compute_weighted_line_coeffC2(lSR, Target_FT(1) / lToFiltered(lHR->weight()), mApproximate_traits.mCoeff_cache);

    Target_Point_2 laP = lToFiltered(aNode->point());

    // @fixme below needs to use inverted weights like in degenerate time/point computations
    Target_Vector_2 lVL(lL->a(), lL->b()) ;
    Target_Vector_2 lVR(lR->a(), lR->b()) ;
    Target_Vector_2 lVLR = lVL + lVR ;
    Target_Ray_2 bisect_ray(laP, lVLR) ;

    // @todo this should use spatial searching
    for ( Halfedge_handle_vector_iterator h = contour_halfedges_begin; h != contour_halfedges_end; ++h )
    {
      try
      {
        CGAL_assertion((*h)->vertex()->is_contour() && (*h)->opposite()->vertex()->is_contour() );

        // @todo could be a line as long as we are in a convex area
        Target_Segment_2 s_h(lToFiltered((*h)->opposite()->vertex()->point()),
                             lToFiltered((*h)->vertex()->point()));

        // we use segments of the input polygon intersected by the bisector and such that
        // they are oriented such that the reflex vertex is on the left side of the segment
        Uncertain<Oriented_side> orient = FK().orientation_2_object()(s_h[0], s_h[1], laP);
        if (!is_certain(orient) || orient != LEFT_TURN)
          continue;

        Uncertain<bool> inter = FK().do_intersect_2_object()(s_h, bisect_ray);
        if (!is_certain(inter) || !inter)
          continue;

        CGAL_SS_i::Segment_2_with_ID<FK> lSh (s_h, (*h)->id());
        auto lh = CGAL_SS_i::compute_weighted_line_coeffC2(lSh, lToFiltered((*h)->weight()), mApproximate_traits.mCoeff_cache);

        // @fixme precision issues...? lToFiltered earlier?
        Target_FT lBound = (- lh->c() - lh->a()*laP.x() - lh->b()*laP.y()) /
                             ( lh->a()*lVLR.x() + lh->b()*lVLR.y() ) - lToFiltered(aNode->time()) ;

        if ( ! is_finite(lBound) || ! is_positive(lBound) )
          continue;

        if ( ! mApproximate_traits.mFilteringBound || *mApproximate_traits.mFilteringBound > lBound )
          mApproximate_traits.mFilteringBound = lBound ;
      }
      catch(CGAL::Uncertain_conversion_exception&)
      {}
    }
  }

public:
  // @todo as soon as an exact value we could refine the interval one. Not sure if it is worth it
  Straight_skeleton_builder_traits_2_impl<Tag_false, EK> mExact_traits ;

  // Below is only used for the cached variables not the functor types
  Straight_skeleton_builder_traits_2_impl<Tag_false, FK> mApproximate_traits ;
} ;

template<class K>
class Straight_skeleton_builder_traits_2
  : public Straight_skeleton_builder_traits_2_impl<typename CGAL_SS_i::Is_filtering_kernel<K>::type, K>
{} ;

CGAL_STRAIGHT_SKELETON_CREATE_FUNCTOR_ADAPTER(Do_ss_event_exist_2)
CGAL_STRAIGHT_SKELETON_CREATE_FUNCTOR_ADAPTER(Compare_ss_event_times_2)
CGAL_STRAIGHT_SKELETON_CREATE_FUNCTOR_ADAPTER(Compare_ss_event_angles_2)
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
