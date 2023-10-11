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
#include <CGAL/Straight_skeleton_2/Straight_skeleton_builder_traits_2_caches.h>
#include <CGAL/Straight_skeleton_2/Straight_skeleton_builder_traits_2_aux.h>
#include <CGAL/Trisegment_2.h>

#include <CGAL/Filtered_construction.h>
#include <CGAL/Uncertain.h>

#include <optional>

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

  Construct_ss_trisegment_2(std::size_t& aNext_ID, Caches<K>& aCaches)
    : mNext_ID(aNext_ID), mCaches(aCaches)
  {}

  result_type operator() ( Segment_2_with_ID const& aS0,
                           FT const& aW0,
                           Segment_2_with_ID const& aS1,
                           FT const& aW1,
                           Segment_2_with_ID const& aS2,
                           FT const& aW2 ) const
  {
    Trisegment_2_ptr rRes = construct_trisegment(aS0,aW0,aS1,aW1,aS2,aW2,mNext_ID,mCaches) ;
    if(rRes)
      ++mNext_ID;
    return rRes;
  }

private:
  std::size_t& mNext_ID;
  Caches<K>& mCaches;
};

template<class K>
struct Do_ss_event_exist_2 : Functor_base_2<K>
{
  typedef Functor_base_2<K> Base ;

  typedef typename Base::FT               FT ;
  typedef typename Base::Trisegment_2_ptr Trisegment_2_ptr ;

  typedef Uncertain<bool> result_type ;

  Do_ss_event_exist_2(Caches<K>& aCaches)
    : mCaches(aCaches)
  {}

  Uncertain<bool> operator() ( Trisegment_2_ptr const& aTrisegment, std::optional<FT> aMaxTime ) const
  {
    Uncertain<bool> rResult = exist_offset_lines_isec2(aTrisegment, aMaxTime, mCaches);

    CGAL_STSKEL_ASSERT_PREDICATE_RESULT(rResult,K,"Exist_event",aTrisegment);

    return rResult ;
  }

private:
  Caches<K>& mCaches;
};

template<class K>
struct Is_edge_facing_ss_node_2 : Functor_base_2<K>
{
  typedef Functor_base_2<K> Base ;

  typedef typename Base::Point_2          Point_2 ;
  typedef typename Base::Segment_2_with_ID        Segment_2_with_ID ;
  typedef typename Base::Trisegment_2_ptr Trisegment_2_ptr ;

  typedef Uncertain<bool> result_type ;

  Is_edge_facing_ss_node_2(Caches<K>& aCaches)
    : mCaches(aCaches)
  { }

  Uncertain<bool> operator() ( Point_2 const& aContourNode, Segment_2_with_ID const& aEdge ) const
  {
    return is_edge_facing_pointC2(cgal_make_optional(aContourNode),aEdge) ;
  }

  Uncertain<bool> operator() ( Trisegment_2_ptr const& aSkeletonNode, Segment_2_with_ID const& aEdge ) const
  {
    return is_edge_facing_offset_lines_isecC2(aSkeletonNode, aEdge, mCaches);
  }

private:
  Caches<K>& mCaches;
};

template<class K>
struct Compare_ss_event_times_2 : Functor_base_2<K>
{
  typedef Functor_base_2<K> Base ;

  typedef typename Base::Trisegment_2_ptr Trisegment_2_ptr ;

  typedef Uncertain<Comparison_result> result_type ;

  Compare_ss_event_times_2(Caches<K>& aCaches)
    : mCaches(aCaches)
  {}

  Uncertain<Comparison_result> operator() ( Trisegment_2_ptr const& aL, Trisegment_2_ptr const& aR ) const
  {
    Uncertain<Comparison_result> rResult = compare_offset_lines_isec_timesC2(aL, aR, mCaches);
    CGAL_STSKEL_ASSERT_PREDICATE_RESULT(rResult,K,"Compare_event_times","L: " << aL << "\nR:" << aR );
    return rResult ;
  }

private:
  Caches<K>& mCaches;
};

template<class K>
struct Compare_ss_event_angles_2 : Functor_base_2<K>
{
  typedef Functor_base_2<K> Base ;

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

  Oriented_side_of_event_point_wrt_bisector_2(Caches<K>& aCaches)
    : mCaches(aCaches)
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
    Uncertain<Oriented_side> rResult = oriented_side_of_event_point_wrt_bisectorC2(aEvent,aE0,aW0,aE1,aW1,aE01Event,aE0isPrimary,mCaches) ;

    CGAL_STSKEL_ASSERT_PREDICATE_RESULT(rResult,K,"Oriented_side_of_event_point_wrt_bisector_2","Event=" << aEvent << " E0=" << aE0 << " E1=" << aE1 );

    return rResult ;
  }

private:
  Caches<K>& mCaches;
};


template<class K>
struct Are_ss_events_simultaneous_2 : Functor_base_2<K>
{
  typedef Functor_base_2<K> Base ;

  typedef typename Base::Trisegment_2_ptr Trisegment_2_ptr ;

  typedef Uncertain<bool> result_type ;

  Are_ss_events_simultaneous_2(Caches<K>& aCaches)
    : mCaches(aCaches)
  {}

  Uncertain<bool> operator() ( Trisegment_2_ptr const& aA, Trisegment_2_ptr const& aB ) const
  {
    Uncertain<bool> rResult = are_events_simultaneousC2(aA,aB, mCaches);

    CGAL_STSKEL_ASSERT_PREDICATE_RESULT(rResult,K,"Are_events_simultaneous","A=" << aA << "\nB=" << aB);

    return rResult ;
  }

private:
  Caches<K>& mCaches;
};

template<class K>
struct Construct_ss_event_time_and_point_2 : Functor_base_2<K>
{
  typedef Functor_base_2<K> Base ;

  typedef typename Base::FT               FT ;
  typedef typename Base::Point_2          Point_2 ;
  typedef typename Base::Segment_2_with_ID        Segment_2_with_ID ;
  typedef typename Base::Trisegment_2_ptr Trisegment_2_ptr ;

  typedef std::tuple<FT,Point_2> rtype ;

  typedef std::optional<rtype> result_type ;

  Construct_ss_event_time_and_point_2(Caches<K>& aCaches)
    : mCaches(aCaches)
  {}

  result_type operator() ( Trisegment_2_ptr const& aTrisegment ) const
  {
    bool lOK = false ;

    FT      t(0) ;
    Point_2 i = ORIGIN ;

    std::optional< Rational<FT> > ot = compute_offset_lines_isec_timeC2(aTrisegment, mCaches);

    if ( !!ot && certainly( CGAL_NTS certified_is_not_zero(ot->d()) ) )
    {
      t = ot->n() / ot->d();

      std::optional<Point_2> oi = construct_offset_lines_isecC2(aTrisegment, mCaches);
      if ( oi )
      {
        i = *oi ;
        lOK = true ;
      }
    }

    CGAL_STSKEL_ASSERT_CONSTRUCTION_RESULT(lOK,K,"Construct_ss_event_time_and_point_2",aTrisegment);

    return cgal_make_optional(lOK,std::make_tuple(t,i)) ;
  }

private:
  Caches<K>& mCaches;
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

  typedef typename Unfiltering::Construct_ss_event_time_and_point_2 Construct_ss_event_time_and_point_2 ;
  typedef typename Unfiltering::Construct_ss_trisegment_2           Construct_ss_trisegment_2 ;

  Is_edge_facing_ss_node_2 is_edge_facing_ss_node_2_object() const
  {
    return Is_edge_facing_ss_node_2(mCaches);
  }

  Compare_ss_event_times_2 compare_ss_event_times_2_object() const
  {
    return Compare_ss_event_times_2(mCaches);
  }

  Compare_ss_event_angles_2 compare_ss_event_angles_2_object() const
  {
    return Compare_ss_event_angles_2();
  }

  Oriented_side_of_event_point_wrt_bisector_2 oriented_side_of_event_point_wrt_bisector_2_object() const
  {
    return Oriented_side_of_event_point_wrt_bisector_2(mCaches);
  }

  Do_ss_event_exist_2 do_ss_event_exist_2_object() const
  {
    return Do_ss_event_exist_2(mCaches);
  }

  Are_ss_events_simultaneous_2 are_ss_events_simultaneous_2_object() const
  {
    return Are_ss_events_simultaneous_2(mCaches);
  }

  Construct_ss_event_time_and_point_2 construct_ss_event_time_and_point_2_object() const
  {
    return Construct_ss_event_time_and_point_2(mCaches);
  }

  Construct_ss_trisegment_2 construct_ss_trisegment_2_object() const
  {
    return Construct_ss_trisegment_2(trisegment_ID(), mCaches);
  }

// ID functions
  std::size_t& trisegment_ID() const { return mTrisegment_ID ; }

  void reset_trisegment(std::size_t i) const
  {
    --mTrisegment_ID ;
    mCaches.Reset(i) ;
  }

// functions to initialize (and harmonize) and cache speeds
  void InitializeLineCoeffs ( CGAL_SS_i::Segment_2_with_ID<K> const& aBorderS )
  {
    CGAL_SS_i::compute_normalized_line_coeffC2(aBorderS, mCaches);
  }

  void InitializeLineCoeffs ( std::size_t aID, std::size_t aOtherID )
  {
    if(mCaches.mCoeff_cache.Get(aOtherID))
      mCaches.mCoeff_cache.Set(aID, CGAL_SS_i::cgal_make_optional(*(mCaches.mCoeff_cache.Get(aOtherID))));
    else
      mCaches.mCoeff_cache.Set(aID, std::nullopt);
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
    std::optional<CGAL_SS_i::Rational<typename K::FT> > lOptTime =
        CGAL_SS_i::compute_offset_lines_isec_timeC2(tri, mCaches);

    if ( lOptTime && lOptTime->to_nt() > *mFilteringBound )
    {
      CGAL_STSKEL_TRAITS_TRACE("Ignoring potential split event");

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
    typedef typename K::Point_2 Point_2;
    typedef typename K::Vector_2 Vector_2;
    typedef typename K::Segment_2 Segment_2;
    typedef typename K::Ray_2 Ray_2;
    typedef typename K::Line_2 Line_2;

    CGAL_STSKEL_TRAITS_TRACE("Computing filtering bound of V" << aNode->id() << " [" << typeid(FT).name() << "]" );

    mFilteringBound = std::nullopt;

    // No gain observed on norway while doing it for more than contour nodes
    if(!aNode->is_contour())
      return;

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

    std::optional< Line_2 > lL = CGAL_SS_i::compute_weighted_line_coeffC2(lSL, lHL->weight(), mCaches);
    std::optional< Line_2 > lR = CGAL_SS_i::compute_weighted_line_coeffC2(lSR, lHR->weight(), mCaches);

    Vector_2 lVL(  lL->b(), - lL->a()) ;
    Vector_2 lVR(- lR->b(),   lR->a()) ;
    Vector_2 lVLR = lVL + lVR ;
    const Point_2& laP = aNode->point();
    Ray_2 bisect_ray(laP, lVLR) ;

    // @todo this should use some kind of spatial searching
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

      // See the other function for the equations
      CGAL_SS_i::Segment_2_with_ID<K> lSh (s_h, (*h)->id());
      std::optional< Line_2 > lh = CGAL_SS_i::compute_normalized_line_coeffC2(lSh, mCaches);

      FT lLambda = - ( lh->a()*laP.x() + lh->b()*laP.y() + lh->c() ) /
                       ( lh->a()*lVLR.x() + lh->b()*lVLR.y() ) ;

      Point_2 lP = laP + lVLR;
      FT lBound = lLambda * ( lL->a()*lP.x() + lL->b()*lP.y() + lL->c() ) ;

      if(!is_finite(lBound) || !is_positive(lBound))
        continue;

      if(!mFilteringBound || *mFilteringBound > lBound)
        mFilteringBound = lBound;
    }

    if(mFilteringBound)
    {
      CGAL_STSKEL_TRAITS_TRACE("Filtering bound: " << *mFilteringBound);
    } else {
      CGAL_STSKEL_TRAITS_TRACE("Filtering bound: none");
    }
  }

public:
  mutable std::size_t mTrisegment_ID = 0 ;
  mutable CGAL_SS_i::Caches<K> mCaches ;
  mutable std::optional< typename K::FT > mFilteringBound ;
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

  typedef CGAL_SS_i::Exceptionless_filtered_construction< typename Unfiltering::Construct_ss_event_time_and_point_2
                                                        , typename Exact      ::Construct_ss_event_time_and_point_2
                                                        , typename Filtering  ::Construct_ss_event_time_and_point_2
                                                        , C2E
                                                        , C2F
                                                        , E2C
                                                        , F2C
                                                        >
                                                        Construct_ss_event_time_and_point_2 ;

  typedef CGAL_SS_i::Exceptionless_filtered_construction< typename Unfiltering::Construct_ss_trisegment_2
                                                        , typename Exact      ::Construct_ss_trisegment_2
                                                        , typename Filtering  ::Construct_ss_trisegment_2
                                                        , C2E
                                                        , C2F
                                                        , E2C
                                                        , F2C
                                                        >
                                                        Construct_ss_trisegment_2 ;

// constructor of predicates using time caching
  Compare_ss_event_times_2 compare_ss_event_times_2_object() const
  {
    return Compare_ss_event_times_2(typename Exact::Compare_ss_event_times_2(mExact_traits.mCaches),
                                    typename Filtering::Compare_ss_event_times_2(mApproximate_traits.mCaches));
  }

  Compare_ss_event_angles_2 compare_ss_event_angles_2_object() const
  {
    return Compare_ss_event_angles_2();
  }

  Is_edge_facing_ss_node_2 is_edge_facing_ss_node_2_object() const
  {
    return Is_edge_facing_ss_node_2(typename Exact::Is_edge_facing_ss_node_2(mExact_traits.mCaches),
                                    typename Filtering::Is_edge_facing_ss_node_2(mApproximate_traits.mCaches));
  }

  Oriented_side_of_event_point_wrt_bisector_2 oriented_side_of_event_point_wrt_bisector_2_object() const
  {
    return Oriented_side_of_event_point_wrt_bisector_2(typename Exact::Oriented_side_of_event_point_wrt_bisector_2(mExact_traits.mCaches),
                                                       typename Filtering::Oriented_side_of_event_point_wrt_bisector_2(mApproximate_traits.mCaches));
  }

  Do_ss_event_exist_2 do_ss_event_exist_2_object() const
  {
    return Do_ss_event_exist_2(typename Exact::Do_ss_event_exist_2(mExact_traits.mCaches),
                               typename Filtering::Do_ss_event_exist_2(mApproximate_traits.mCaches));
  }

  Are_ss_events_simultaneous_2 are_ss_events_simultaneous_2_object() const
  {
    return Are_ss_events_simultaneous_2(typename Exact::Are_ss_events_simultaneous_2(mExact_traits.mCaches),
                                        typename Filtering::Are_ss_events_simultaneous_2(mApproximate_traits.mCaches));
  }

  Construct_ss_event_time_and_point_2 construct_ss_event_time_and_point_2_object() const
  {
    return Construct_ss_event_time_and_point_2(typename Exact::Construct_ss_event_time_and_point_2(mExact_traits.mCaches),
                                               typename Filtering::Construct_ss_event_time_and_point_2(mApproximate_traits.mCaches) );
  }

// constructor of trisegments using global id stored in the traits
  Construct_ss_trisegment_2 construct_ss_trisegment_2_object() const
  {
    return Construct_ss_trisegment_2(typename Exact::Construct_ss_trisegment_2(trisegment_ID(), mExact_traits.mCaches),
                                     typename Filtering::Construct_ss_trisegment_2(trisegment_ID(), mApproximate_traits.mCaches));
  }

// ID functions
  std::size_t& trisegment_ID() const { return mApproximate_traits.mTrisegment_ID ; }

  void reset_trisegment(std::size_t i) const
  {
    if ( i+1 == trisegment_ID() )
    {
      --trisegment_ID() ;
      mApproximate_traits.mCaches.Reset(i) ;
      mExact_traits.mCaches.Reset(i) ;
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
      std::optional<CGAL_SS_i::Rational<typename FK::FT> > lOptTime =
          CGAL_SS_i::compute_offset_lines_isec_timeC2(tri, mApproximate_traits.mCaches);

      if ( lOptTime && lOptTime->to_nt() > *mApproximate_traits.mFilteringBound )
      {
        CGAL_STSKEL_TRAITS_TRACE("Ignoring potential split event");

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

    typedef decltype(aNode->halfedge()) Halfedge_handle;
    typedef CGAL_SS_i::Segment_2_with_ID<FK> Target_Segment_with_ID_2;

    CGAL_STSKEL_TRAITS_TRACE("Computing approximate filtering bound of V" << aNode->id() << " [" << typeid(Target_FT).name() << "]" );

    mApproximate_traits.mFilteringBound = std::nullopt;

    // No gain observed on norway while doing it for more than contour nodes
    if(!aNode->is_contour())
      return;

    typename FK::FT::Protector protector;

    C2F lToFiltered;

    // get the contour input segments on each side of the bisector spawned ataNode
    Halfedge_handle lHL = aNode->halfedge()->defining_contour_edge();
    Halfedge_handle lHR = ( aNode->is_contour() ) ? lHL->opposite()->prev()->opposite()
                                                  : aNode->halfedge()->opposite()->defining_contour_edge() ;

    Target_Segment_with_ID_2 lSL (lToFiltered(lHL->opposite()->vertex()->point()),
                                  lToFiltered(lHL->vertex()->point()),
                                  lHL->id());
    Target_Segment_with_ID_2 lSR (lToFiltered(lHR->opposite()->vertex()->point()),
                                  lToFiltered(lHR->vertex()->point()),
                                  lHR->id());

    std::optional<Target_Line_2> lL = CGAL_SS_i::compute_weighted_line_coeffC2(lSL, lToFiltered(lHL->weight()), mApproximate_traits.mCaches);
    std::optional<Target_Line_2> lR = CGAL_SS_i::compute_weighted_line_coeffC2(lSR, lToFiltered(lHR->weight()), mApproximate_traits.mCaches);

    Target_Point_2 laP = lToFiltered(aNode->point());

    // These are weighted direction of the lines supporting the contour segments.
    // Coefficients for SR are negated because aNode is at the "source" of SR.
    Target_Vector_2 lVL(  lL->b(), - lL->a()) ;
    Target_Vector_2 lVR(- lR->b(),   lR->a()) ;
    Target_Vector_2 lVLR = lVL + lVR ;
    Target_Ray_2 bisect_ray(laP, lVLR) ;

    // @todo this should use some kind of spatial searching
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

        Uncertain<bool> inter = FK().do_intersect_2_object()(bisect_ray, s_h);
        if (!is_certain(inter) || !inter)
          continue;

        // We want the time it takes to get from aNode to h along the primary bisector of aNode.
        //
        // Let d0 be the weighted direction of the defining contour edge to the left of aNode, that is (-b0, a0) / w0.
        // Let d1 be the weighted, opposite direction of the defining contour edge to the right of aNode, that is (b1, -a1) / w1.
        // The bisector has direction n0 + n1.
        // Let n2 be the unweighted normal of 'h', that is v3 = (a, b) if h is defined as la*x+lb*y+lc = 0
        //
        // Projecting aNode onto h orthogonally, we create a right triangle.
        // Let theta be the angle between the line orthogonal to h through aNode and the primary bisector of aNode.
        // Let T be the distance between aNode and the projection of aNode onto h.
        // Let H be the distance between aNode and the intersection of aNode and h along the primary bisector of aNode.
        //
        // We have:
        //   cos(theta) = (d0 + d1) * n2 / (|d0 + d1| * |n2|)  [note that |n2| = 1]
        // and on the other hand:
        //   cos(theta) = T / H
        // If we express H as lambda * |d0 + d1|, we get:
        //   lambda = - T / (d0 + d1) * n2

        Target_Segment_with_ID_2 lSh (s_h, (*h)->id());
        std::optional<Target_Line_2> lh = CGAL_SS_i::compute_normalized_line_coeffC2(lSh, mApproximate_traits.mCaches);

        Target_FT lLambda = - ( lh->a()*laP.x() + lh->b()*laP.y() + lh->c() ) /
                                ( lh->a()*lVLR.x() + lh->b()*lVLR.y() ) ;

        // Scale it because |d0 + d1| doesn't send aNode to the t=1 line, but to t=lL(aNode+d0+d1)=lR(aNode+d0+d1)
        Target_Point_2 lP = laP + lVLR;
        Target_FT lBound = lLambda * ( lL->a()*lP.x() + lL->b()*lP.y() + lL->c() ) ;

#if 0
        std::cout << "E" << (*h)->id() << " s_h = " << s_h << std::endl;
        std::cout << "left/right E" << lHL->id() << " E" << lHR->id() << std::endl;
        std::cout << "V" << aNode->id() << std::endl;
        std::cout << "lSL: " << lSL << std::endl;
        std::cout << "lSR: " << lSR << std::endl;
        std::cout << "lHL->weight(): " << lHL->weight() << std::endl;
        std::cout << "lHR->weight(): " << lHR->weight() << std::endl;
        std::cout << "lVL " << lVL << std::endl;
        std::cout << "lVR " << lVR << std::endl;
        std::cout << "lVLR " << lVLR << std::endl;
        std::cout << "lAP " << laP <<  std::endl;
        std::cout << "lAP + lVL " << laP + lVL <<  std::endl;
        std::cout << "lAP + lVR " << laP + lVR <<  std::endl;
        std::cout << "lAP+v " << laP + lVLR << std::endl;
        std::cout << "Inter pt: " << *ip << std::endl;

        std::optional<Target_Line_2> lh1 = CGAL_SS_i::compute_weighted_line_coeffC2(lSR, lToFiltered(lHR->weight()), mApproximate_traits.mCaches);

        std::cout << "lh0 check" << square(lh0->a()) + square(lh0->b()) << std::endl;
        std::cout << "lh1 check" << square(lh1->a()) + square(lh1->b()) << std::endl;
        std::cout << "l0 time at aNode: " << lh0->a()*laP.x() + lh0->b()*laP.y() + lh0->c() << std::endl;
        std::cout << "l1 time at aNode: " << lh1->a()*laP.x() + lh1->b()*laP.y() + lh1->c() << std::endl;
        std::cout << "l0 time at aNode + lVLR: " << lh0->a()*(laP + lVLR).x() + lh0->b()*(laP + lVLR).y() + lh0->c() << std::endl;
        std::cout << "l1 time at aNode + lVLR: " << lh1->a()*(laP + lVLR).x() + lh1->b()*(laP + lVLR).y() + lh1->c() << std::endl;

        auto ipp = FK().intersect_2_object()(s_h, bisect_ray);
        Target_Point_2* ip = std::get<Target_Point_2>(&*ipp);
        std::cout << "l0 time at inter pt: " << lh0->a()*ip->x() + lh0->b()*ip->y() + lh0->c() << std::endl;
        std::cout << "l1 time at inter pt: " << lh1->a()*ip->x() + lh1->b()*ip->y() + lh1->c() << std::endl;
        std::cout << "lh-> " << lh->a() << " " << lh->b() << " " << square(lh->a()) + square(lh->b()) << std::endl;
        std::cout << "lh time at inter pt: " << lh->a()*ip->x() + lh->b()*ip->y() + lh->c() << std::endl;
        std::cout << "lLambda: " << lLambda << std::endl;
        std::cout << "lBound: " << lBound << std::endl;
        CGAL_assertion(lBound == lh0->a()*ip->x() + lh0->b()*ip->y() + lh0->c());
#endif

        if(!is_finite(lBound) || !is_positive(lBound))
          continue;

        if(!mApproximate_traits.mFilteringBound || *mApproximate_traits.mFilteringBound > lBound)
          mApproximate_traits.mFilteringBound = lBound;
      }
      catch(CGAL::Uncertain_conversion_exception&)
      {}
    }

    if(mApproximate_traits.mFilteringBound)
    {
      CGAL_STSKEL_TRAITS_TRACE("Filtering bound: " << *mApproximate_traits.mFilteringBound);
    } else {
      CGAL_STSKEL_TRAITS_TRACE("Filtering bound: none");
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

} // end namespace CGAL

#endif // CGAL_STRAIGHT_SKELETON_BUILDER_TRAITS_2_H //
// EOF //
