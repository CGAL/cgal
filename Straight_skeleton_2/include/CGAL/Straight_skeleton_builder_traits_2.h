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

  typedef Trisegment_2_ptr result_type ;

  template <class Traits>
  Construct_ss_trisegment_2(const Traits& traits)
    : mNext_ID(traits.trisegment_ID())
  {}

  result_type operator() ( Segment_2_with_ID const& aS0, Segment_2_with_ID const& aS1, Segment_2_with_ID const& aS2 ) const
  {
    return construct_trisegment(aS0,aS1,aS2,mNext_ID++) ;
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

  typedef Uncertain<Oriented_side> result_type ;

  Oriented_side_of_event_point_wrt_bisector_2(Coeff_cache<K>& aCoeff_cache)
    : mCoeff_cache(aCoeff_cache)
  {}

  Uncertain<Oriented_side> operator() ( Trisegment_2_ptr const& aEvent
                                      , Segment_2_with_ID        const& aE0
                                      , Segment_2_with_ID        const& aE1
                                      , Trisegment_2_ptr const& aE01Event
                                      , bool                    aE0isPrimary
                                      ) const
  {
    Uncertain<Oriented_side> rResult = oriented_side_of_event_point_wrt_bisectorC2(aEvent,aE0,aE1,aE01Event,aE0isPrimary,mCoeff_cache) ;

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
        CGAL_stskel_intrinsic_test_assertion(!is_point_calculation_clearly_wrong(t,i,aTrisegment));
        lOK = true ;
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
      Segment_2_with_ID const& e0 = aTrisegment->e0() ;
      Segment_2_with_ID const& e1 = aTrisegment->e1() ;
      Segment_2_with_ID const& e2 = aTrisegment->e2() ;

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
    CGAL_SS_i::compute_normalized_line_ceoffC2 ( aBorderS, mCoeff_cache ) ;
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
  void ComputeFilteringBound(Vertex_handle lPrev, Vertex_handle aNode, Vertex_handle lNext,
                             Halfedge_handle_vector_iterator contour_halfedges_begin,
                             Halfedge_handle_vector_iterator contour_halfedges_end) const
  {
    mFilteringBound = boost::none;

    if ( ! aNode->is_contour() )
      return ;

    typedef typename K::Vector_2 Vector_2;
    typedef typename K::Segment_2 Segment_2;
    typedef typename K::Ray_2 Ray_2;
    typedef typename K::Line_2 Line_2;

    Segment_2 s1(lPrev->point(), aNode->point());
    Segment_2 s2(aNode->point(), lNext->point());

    // @todo? These are not input segments, but it might still worth caching (just gotta assign them an ID)
    boost::optional< Line_2 > l1 = CGAL_SS_i::compute_normalized_line_ceoffC2(s1);
    boost::optional< Line_2 > l2 = CGAL_SS_i::compute_normalized_line_ceoffC2(s2);

    Vector_2 lV1(l1->a(), l1->b()) ;
    Vector_2 lV2(l2->a(), l2->b()) ;
    Vector_2 lV12 = lV1 + lV2 ;
    Ray_2 bisect_ray(aNode->point(), lV12) ;

    // F2C to_input;
    // std::cout << "bisect " << aNode->point() << " " << aNode->point() + to_input(lV12) << "\n";

    for ( Halfedge_handle_vector_iterator i = contour_halfedges_begin; i != contour_halfedges_end; ++ i )
    {
      CGAL_assertion((*i)->vertex()->is_contour() && (*i)->opposite()->vertex()->is_contour() );
      Segment_2 s_h((*i)->opposite()->vertex()->point(), (*i)->vertex()->point());

      auto inter = K().do_intersect_2_object()(s_h, bisect_ray);
      auto orient = K().orientation_2_object()(s_h[0], s_h[1], aNode->point());

      // we use segments of the input polygon intersected by the bisector and such that
      // they are oriented such that the reflex vertex is on the left side of the segment
      if (!is_certain(inter) || !is_certain(orient) || !inter || orient != LEFT_TURN)
        continue;

      // Note that we don't need the normalization
      boost::optional< Line_2 > lh = CGAL_SS_i::compute_normalized_line_ceoffC2(s_h);

      typename K::FT lBound = ( - lh->c() - lh->a()*aNode->point().x() - lh->b()*aNode->point().y() ) /
                                ( lh->a()*lV12.x() + lh->b()*lV12.y() );

      if( ! is_finite(lBound) )
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
  void ComputeFilteringBound(Vertex_handle lPrev, Vertex_handle aNode, Vertex_handle lNext,
                             Halfedge_handle_vector_iterator contour_halfedges_begin,
                             Halfedge_handle_vector_iterator contour_halfedges_end) const
  {
    mApproximate_traits.mFilteringBound = boost::none;

    if ( ! aNode->is_contour() )
      return ;

    typedef typename FK::Point_2 Target_Point_2;
    typedef typename FK::Vector_2 Target_Vector_2;
    typedef typename FK::Segment_2 Target_Segment_2;
    typedef typename FK::Ray_2 Target_Ray_2;
    typedef typename FK::Line_2 Target_Line_2;

    typename FK::FT::Protector protector;

    C2F to_FK;

    Target_Point_2 laP = to_FK(aNode->point());
    Target_Segment_2 s1(to_FK(lPrev->point()), laP);
    Target_Segment_2 s2(laP, to_FK(lNext->point()));

    // @todo? These are not input segments, but it might still worth caching (just gotta assign them an ID)
    boost::optional< Target_Line_2 > l1 = CGAL_SS_i::compute_normalized_line_ceoffC2(s1);
    boost::optional< Target_Line_2 > l2 = CGAL_SS_i::compute_normalized_line_ceoffC2(s2);

    Target_Vector_2 lV1(l1->a(), l1->b()) ;
    Target_Vector_2 lV2(l2->a(), l2->b()) ;
    Target_Vector_2 lV12 = lV1 + lV2 ;
    Target_Ray_2 bisect_ray(laP, lV12) ;

    // F2C to_input;
    // std::cout << "bisect " << aNode->point() << " " << aNode->point() + to_input(lV12) << "\n";

    for ( Halfedge_handle_vector_iterator i = contour_halfedges_begin; i != contour_halfedges_end; ++ i )
    {
      try
      {
        CGAL_assertion((*i)->vertex()->is_contour() && (*i)->opposite()->vertex()->is_contour() );
        Target_Segment_2 s_h(to_FK((*i)->opposite()->vertex()->point()), to_FK((*i)->vertex()->point()));

        Uncertain<bool> inter = FK().do_intersect_2_object()(s_h, bisect_ray);
        Uncertain<Oriented_side> orient = FK().orientation_2_object()(s_h[0], s_h[1], laP);

        // we use segments of the input polygon intersected by the bisector and such that
        // they are oriented such that the reflex vertex is on the left side of the segment
        if (!is_certain(inter) || !is_certain(orient) || !inter || orient != LEFT_TURN)
          continue;

        // Note that we don't need the normalization
        boost::optional< Target_Line_2 > lh = CGAL_SS_i::compute_normalized_line_ceoffC2(s_h);

        typename FK::FT lBound = ( - lh->c() - lh->a()*laP.x() - lh->b()*laP.y() ) /
                                   ( lh->a()*lV12.x() + lh->b()*lV12.y() );

        if ( ! is_finite(lBound) )
          continue;

        if ( ! mApproximate_traits.mFilteringBound || *mApproximate_traits.mFilteringBound > lBound )
          mApproximate_traits.mFilteringBound = lBound ;
      }
      catch(CGAL::Uncertain_conversion_exception&)
      {}
    }
  }

public:
  // TODO: as soon as an exact value we could refine the interval one. Not sure if it is worth it
  Exact mExact_traits ;

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
