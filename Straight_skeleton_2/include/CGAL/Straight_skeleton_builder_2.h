// Copyright (c) 2006-2008 Fernando Luis Cacciola Carballal. All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//

// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Fernando Cacciola <fernando_cacciola@ciudad.com.ar>
//
#ifndef CGAL_STRAIGHT_SKELETON_BUILDER_2_H
#define CGAL_STRAIGHT_SKELETON_BUILDER_2_H 1

#include <CGAL/license/Straight_skeleton_2.h>

#include <CGAL/algorithm.h>
#include <CGAL/Straight_skeleton_2/Straight_skeleton_aux.h>
#include <CGAL/Straight_skeleton_2/Straight_skeleton_builder_events_2.h>
#include <CGAL/Straight_skeleton_2.h>
#include <CGAL/Straight_skeleton_builder_traits_2.h>
#include <CGAL/HalfedgeDS_const_decorator.h>
#include <CGAL/enum.h>

#include <boost/tuple/tuple.hpp>
#include <boost/intrusive_ptr.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/mpl/bool.hpp>

#include <algorithm>
#include <exception>
#include <iostream>
#include <map>
#include <list>
#include <queue>
#include <string>
#include <sstream>
#include <utility>

namespace CGAL {

template<class SSkel_>
struct Dummy_straight_skeleton_builder_2_visitor
{
  typedef SSkel_ SSkel ;

  typedef typename SSkel::Halfedge_const_handle Halfedge_const_handle ;
  typedef typename SSkel::Vertex_const_handle   Vertex_const_handle ;

  void on_contour_edge_entered ( Halfedge_const_handle const&  ) const {}

  void on_initialization_started( std::size_t /* size_of_vertices */ ) const {}

  void on_initial_events_collected( Vertex_const_handle const& , bool /* is_reflex */, bool /*is_degenerate*/ )  const  {}

  void on_edge_event_created( Vertex_const_handle const& /* lnode */
                            , Vertex_const_handle const& /* rnode */
                            )  const {}

  void on_split_event_created( Vertex_const_handle const&  )  const {}

  void on_pseudo_split_event_created( Vertex_const_handle const& /* lnode */
                                    , Vertex_const_handle const& /* rnode */
                                    )  const {}

  void on_initialization_finished() const {}

  void on_propagation_started() const {}

  void on_anihiliation_event_processed ( Vertex_const_handle const& /* node0 */
                                       , Vertex_const_handle const& /* node1 */
                                       )  const  {}


  void on_edge_event_processed( Vertex_const_handle const& /* lseed */
                              , Vertex_const_handle const& /* rseed */
                              , Vertex_const_handle const& /* node */
                              )  const {}

  void on_split_event_processed( Vertex_const_handle const& /* seed */
                               , Vertex_const_handle const& /* node0 */
                               , Vertex_const_handle const& /* node1 */
                               )  const {}

  void on_pseudo_split_event_processed( Vertex_const_handle const& /* lseed */
                                      , Vertex_const_handle const& /* rseed */
                                      , Vertex_const_handle const& /* node0 */
                                      , Vertex_const_handle const& /* node1 */
                                      )  const {}

  void on_vertex_processed( Vertex_const_handle const& )  const {}

  void on_propagation_finished() const {}

  void on_cleanup_started() const {}

  void on_cleanup_finished() const {}

  void on_algorithm_finished ( bool /* finished_ok */ ) const {}

  void on_error( char const* ) const {}
} ;

template<class Traits_, class SSkel_, class Visitor_ = Dummy_straight_skeleton_builder_2_visitor<SSkel_> >
class Straight_skeleton_builder_2
{
public:

  typedef Traits_  Traits ;
  typedef SSkel_   SSkel ;
  typedef Visitor_ Visitor ;

  typedef boost::shared_ptr<SSkel> SSkelPtr ;

private :

  typedef typename Traits::Kernel K ;

  typedef typename Traits::FT               FT ;
  typedef typename Traits::Point_2          Point_2 ;
  typedef typename Traits::Vector_2         Vector_2 ;
  typedef typename Traits::Direction_2      Direction_2 ;
  typedef typename Traits::Segment_2        Segment_2 ;
  typedef typename Traits::Trisegment_2     Trisegment_2 ;
  typedef typename Traits::Trisegment_2_ptr Trisegment_2_ptr ;

  typedef typename SSkel::Vertex   Vertex ;
  typedef typename SSkel::Halfedge Halfedge ;
  typedef typename SSkel::Face     Face ;

  typedef typename SSkel::Vertex_const_handle   Vertex_const_handle ;
  typedef typename SSkel::Halfedge_const_handle Halfedge_const_handle ;
  typedef typename SSkel::Face_const_handle     Face_const_handle ;

  typedef typename SSkel::Vertex_const_iterator   Vertex_const_iterator ;
  typedef typename SSkel::Halfedge_const_iterator Halfedge_const_iterator ;
  typedef typename SSkel::Face_const_iterator     Face_const_iterator ;

  typedef typename SSkel::Vertex_handle   Vertex_handle ;
  typedef typename SSkel::Halfedge_handle Halfedge_handle ;
  typedef typename SSkel::Face_handle     Face_handle ;

  typedef typename SSkel::Vertex_iterator   Vertex_iterator ;
  typedef typename SSkel::Halfedge_iterator Halfedge_iterator ;
  typedef typename SSkel::Face_iterator     Face_iterator ;

  typedef typename Vertex::Halfedge_around_vertex_circulator Halfedge_around_vertex_circulator ;

  typedef typename SSkel::size_type size_type ;

  typedef CGAL_SS_i::Triedge<Halfedge_handle> Triedge ;

  typedef CGAL_SS_i::Event_2             <SSkel,Traits> Event ;
  typedef CGAL_SS_i::Edge_event_2        <SSkel,Traits> EdgeEvent ;
  typedef CGAL_SS_i::Split_event_2       <SSkel,Traits> SplitEvent ;
  typedef CGAL_SS_i::Pseudo_split_event_2<SSkel,Traits> PseudoSplitEvent ;

  typedef boost::intrusive_ptr<Event> EventPtr ;

  typedef std::pair<Vertex_handle,Vertex_handle> Vertex_handle_pair ;

  typedef std::vector<EventPtr>           EventPtr_Vector ;
  typedef std::vector<Halfedge_handle>    Halfedge_handle_vector ;
  typedef std::vector<Vertex_handle>      Vertex_handle_vector ;
  typedef std::vector<Vertex_handle_pair> Vertex_handle_pair_vector ;

  typedef typename Halfedge_handle_vector   ::iterator Halfedge_handle_vector_iterator ;
  typedef typename Vertex_handle_vector     ::iterator Vertex_handle_vector_iterator ;
  typedef typename Vertex_handle_pair_vector::iterator Vertex_handle_pair_vector_iterator ;

  typedef typename EventPtr_Vector::const_iterator event_const_iterator ;

  typedef Straight_skeleton_builder_2<Traits,SSkel,Visitor> Self ;

  typedef typename Halfedge::Base_base HBase_base ;
  typedef typename Halfedge::Base      HBase ;
  typedef typename Vertex::Base        VBase ;
  typedef typename Face::Base          FBase ;

  enum Site { AT_SOURCE = -1 , INSIDE = 0, AT_TARGET = +1 } ;

  struct Multinode : public Ref_counted_base
  {
    Multinode ( Halfedge_handle b, Halfedge_handle e )
      :
       begin(b)
      ,end  (e)
      ,v    (b->vertex())
      ,size (0)
    {}

    Halfedge_handle        begin ;
    Halfedge_handle        end ;
    Vertex_handle          v ;
    std::size_t            size ;
    Halfedge_handle_vector bisectors_to_relink ;
    Halfedge_handle_vector bisectors_to_remove ;
    Vertex_handle_vector   nodes_to_remove ;
  } ;

  typedef boost::intrusive_ptr<Multinode> MultinodePtr ;

  struct MultinodeComparer
  {
    bool operator() ( MultinodePtr const& x, MultinodePtr const& y ) { return x->size > y->size ; }
  } ;

  typedef std::vector<MultinodePtr> MultinodeVector ;

  struct Halfedge_ID_compare : CGAL::cpp98::binary_function<bool,Halfedge_handle,Halfedge_handle>
  {
    bool operator() ( Halfedge_handle const& aA, Halfedge_handle const& aB ) const
    {
      return aA->id() < aB->id() ;
    }
  } ;

public:

  Straight_skeleton_builder_2 ( boost::optional<FT> aMaxTime = boost::none, Traits const& = Traits(), Visitor const& aVisitor = Visitor() ) ;

  SSkelPtr construct_skeleton(  bool aNull_if_failed = true ) ;

private :
  // Only used to debug
  double ComputeApproximateAngle ( Vector_2 const& u, Vector_2 const& v) const
  {
    typename K::Compute_scalar_product_2 sp = K().compute_scalar_product_2_object();

    double product = CGAL::sqrt(to_double(sp(u,u)) * to_double(sp(v,v)));

    if(product == FT(0))
      return 0;

    // cosine
    double dot = to_double(sp(u,v));
    double cosine = dot / product;

    if(cosine > 1.)
      cosine = 1.;
    if(cosine < -1.)
      cosine = -1.;

    return std::acos(cosine) * 180. / CGAL_PI;
  }
// Only used to debug
  double ComputeSupportsAngleSplit ( EventPtr const& aEvent )
  {
    SplitEvent& lEvent = dynamic_cast<SplitEvent&>(*aEvent) ;

    Vector_2 lV1 ( aEvent->point(), lEvent.seed0()->point() );
    Vector_2 lV2 ( aEvent->point(), aEvent->point() + CreateVector(aEvent->triedge().e2()) ) ;

    return ComputeApproximateAngle(lV1, lV2) ;
  }
// Only used to debug
  double ComputeSupportsAnglePseudoSplit ( EventPtr const& aEvent )
  {
    PseudoSplitEvent& lEvent = dynamic_cast<PseudoSplitEvent&>(*aEvent) ;

    Vector_2 lV1 ;
    Vector_2 lV2 ;

    if(lEvent.is_at_source_vertex())
    {
      // is_at_source_vertex <=> opposite node is seed0
      lV1 = Vector_2( aEvent->point(), lEvent.seed1()->point() );
      lV2 = Vector_2( aEvent->point(), aEvent->point() + CreateVector(aEvent->triedge().e2()) ) ;
    }
    else
    {
      lV1 = Vector_2( aEvent->point(), lEvent.seed0()->point() );
      lV2 = Vector_2( aEvent->point(), aEvent->point() - CreateVector(aEvent->triedge().e2()) ) ;
    }

    return ComputeApproximateAngle(lV1, lV2) ;
  }
// Only used to debug
  double ComputeSupportsAngle ( EventPtr const& aEvent )
  {
    if ( aEvent->type() == Event::cSplitEvent )
    {
      Halfedge_handle lOppEdge = aEvent->triedge().e2() ;
      Site lSite;
      Vertex_handle_pair lOpp = LookupOnSLAV(lOppEdge,aEvent,lSite);
      if ( handle_assigned(lOpp.first) )
      {
        EventPtr lPseudoSplitEvent = IsPseudoSplitEvent(aEvent,lOpp,lSite);
        if ( lPseudoSplitEvent )
          return ComputeSupportsAnglePseudoSplit ( lPseudoSplitEvent ) ;
        else
          return ComputeSupportsAngleSplit ( aEvent ) ;
      }
      else
      {
        return 2 * CGAL_PI; // event isn't valid anymore
      }
    }
    else
    {
      CGAL_assertion ( aEvent->type() == Event::cPseudoSplitEvent );
      return ComputeSupportsAngleSplit ( aEvent ) ;
    }
  }

  // Real stuff

  Comparison_result CompareEventsSupportAnglesSplitSplit ( EventPtr const& aA, EventPtr const& aB )
  {
    CGAL_precondition ( aA->triedge().e0() == aB->triedge().e0() && aA->triedge().e1() == aB->triedge().e1() ) ;
    return Compare_ss_event_angles_2(mTraits)( CreateVector(aA->triedge().e0()),
                                               CreateVector(aA->triedge().e1()),
                                               CreateVector(aA->triedge().e2()),
                                               CreateVector(aB->triedge().e2()) );
  }

  Comparison_result CompareEventsSupportAnglesSplitPseudoSplit ( EventPtr const& aA, EventPtr const& aB )
  {
    CGAL_precondition ( aA->triedge().e0() == aB->triedge().e0() && aA->triedge().e1() == aB->triedge().e1() ) ;

    PseudoSplitEvent& lPSEvent = dynamic_cast<PseudoSplitEvent&>(*aB) ;
    if(lPSEvent.is_at_source_vertex())
    {
      return Compare_ss_event_angles_2(mTraits)( CreateVector(aA->triedge().e0()),
                                                 CreateVector(aA->triedge().e1()),
                                                 CreateVector(aA->triedge().e2()),
                                                 CreateVector(aB->triedge().e2()) );
    }
    else
    {
      return Compare_ss_event_angles_2(mTraits)( CreateVector(aA->triedge().e0()),
                                                 CreateVector(aA->triedge().e1()),
                                                 CreateVector(aA->triedge().e2()),
                                                 K().construct_opposite_vector_2_object()( CreateVector(aB->triedge().e2())) );
    }
  }

  Comparison_result CompareEventsSupportAnglesPseudoSplitPseudoSplit ( EventPtr const& aA, EventPtr const& aB )
  {
    CGAL_precondition ( aA->triedge().e0() == aB->triedge().e0() && aA->triedge().e1() == aB->triedge().e1() ) ;

    PseudoSplitEvent& lPSEventA = dynamic_cast<PseudoSplitEvent&>(*aA) ;
    PseudoSplitEvent& lPSEventB = dynamic_cast<PseudoSplitEvent&>(*aB) ;

    if(lPSEventA.is_at_source_vertex())
    {
      if(lPSEventB.is_at_source_vertex())
      {
        return Compare_ss_event_angles_2(mTraits)( CreateVector(aA->triedge().e0()),
                                                   CreateVector(aA->triedge().e1()),
                                                   CreateVector(aA->triedge().e2()),
                                                   CreateVector(aB->triedge().e2()) );
      }
      else
      {
        return Compare_ss_event_angles_2(mTraits)( CreateVector(aA->triedge().e0()),
                                                   CreateVector(aA->triedge().e1()),
                                                   CreateVector(aA->triedge().e2()),
                                                   K().construct_opposite_vector_2_object()( CreateVector(aB->triedge().e2())) );
      }
    }
    else // aA is a Pseudo-split Event at the target
    {
      if(lPSEventB.is_at_source_vertex())
      {
        return Compare_ss_event_angles_2(mTraits)( CreateVector(aA->triedge().e0()),
                                                   CreateVector(aA->triedge().e1()),
                                                   K().construct_opposite_vector_2_object()( CreateVector(aA->triedge().e2()) ),
                                                   CreateVector(aB->triedge().e2()) );
      }
      else
      {
        return Compare_ss_event_angles_2(mTraits)( CreateVector(aA->triedge().e0()),
                                                   CreateVector(aA->triedge().e1()),
                                                   K().construct_opposite_vector_2_object()( CreateVector(aA->triedge().e2())),
                                                   K().construct_opposite_vector_2_object()( CreateVector(aB->triedge().e2())) );
      }
    }
  }

  Comparison_result CompareEventsSupportAnglesSplitX ( EventPtr const& aA, EventPtr const& aB )
  {
    if ( aB->type() == Event::cSplitEvent )
    {
      Halfedge_handle lOppEdge = aB->triedge().e2() ;
      Site lSite;
      Vertex_handle_pair lOpp = LookupOnSLAV(lOppEdge,aB,lSite);
      if ( handle_assigned(lOpp.first) )
      {
        EventPtr lPseudoSplitEvent = IsPseudoSplitEvent(aB,lOpp,lSite);
        if ( lPseudoSplitEvent )
          return CompareEventsSupportAnglesSplitPseudoSplit ( aA,lPseudoSplitEvent ) ;
        else
          return CompareEventsSupportAnglesSplitSplit ( aA,aB ) ;
      }
      else
      {
        // Event B does not exist, so give it the lower priority by returning SMALLER
        // (meaning, operator() == false and A has higher priority than B)
        return CGAL::SMALLER;
      }
    }
    else
    {
      CGAL_assertion ( aB->type() == Event::cPseudoSplitEvent );
      return CompareEventsSupportAnglesSplitPseudoSplit ( aA,aB ) ;
    }
  }

  Comparison_result CompareEventsSupportAnglesPseudoSplitX ( EventPtr const& aA, EventPtr const& aB )
  {
    if ( aB->type() == Event::cSplitEvent )
    {
      Halfedge_handle lOppEdge = aB->triedge().e2() ;
      Site lSite;
      Vertex_handle_pair lOpp = LookupOnSLAV(lOppEdge,aB,lSite);
      if ( handle_assigned(lOpp.first) )
      {
        EventPtr lPseudoSplitEvent = IsPseudoSplitEvent(aB,lOpp,lSite);
        if ( lPseudoSplitEvent )
        {
          return CompareEventsSupportAnglesPseudoSplitPseudoSplit ( aA,lPseudoSplitEvent ) ;
        }
        else
        {
          Comparison_result lRes = CompareEventsSupportAnglesSplitPseudoSplit ( aB,aA ) ;
          if ( lRes == LARGER )
            return SMALLER ;
          else if ( lRes == SMALLER )
            return LARGER ;
          else
            return EQUAL ;
        }
      }
      else
      {
        // Event B does not exist, so give it the lower priority by returning SMALLER
        // (meaning, operator() == false and A has higher priority than B)
        return CGAL::SMALLER;
      }
    }
    else
    {
      CGAL_assertion ( aB->type() == Event::cPseudoSplitEvent );
      return CompareEventsSupportAnglesPseudoSplitPseudoSplit ( aA,aB ) ;
    }
  }

  // @todo Should split events always have lower priority than split events?
  Comparison_result CompareEventsSupportAngles ( EventPtr const& aA, EventPtr const& aB )
  {
    CGAL_precondition ( aA->type() != Event::cEdgeEvent && aB->type() != Event::cEdgeEvent ) ;

    if(aA->triedge() == aB->triedge())
      return EQUAL;

    if ( aA->type() == Event::cSplitEvent )
    {
      Halfedge_handle lOppEdge = aA->triedge().e2() ;
      Site lSite;
      Vertex_handle_pair lOpp = LookupOnSLAV(lOppEdge,aA,lSite);
      if ( handle_assigned(lOpp.first) )
      {
        EventPtr lPseudoSplitEvent = IsPseudoSplitEvent(aA,lOpp,lSite);
        if ( lPseudoSplitEvent )
          return CompareEventsSupportAnglesPseudoSplitX ( lPseudoSplitEvent,aB ) ;
        else
          return CompareEventsSupportAnglesSplitX ( aA,aB ) ;
      }
      else
      {
        // Event A does not exist, so give it the lower priority by returning LARGER
        // (meaning, operator() == true and A has lower priority than B)
        //
        // Exception: event B also is a pseudo split that does not exist...
        if ( aB->type() == Event::cSplitEvent )
        {
          Halfedge_handle lOppEdgeB = aB->triedge().e2() ;
          Site lSiteB;
          Vertex_handle_pair lOppB = LookupOnSLAV(lOppEdgeB,aB,lSiteB);
          if ( !handle_assigned(lOppB.first) )
            return EQUAL;
        }

        return LARGER;
      }
    }
    else
    {
      CGAL_assertion ( aA->type() == Event::cPseudoSplitEvent );
      return CompareEventsSupportAnglesPseudoSplitX ( aA,aB ) ;
    }
  }

public:
  // Event compare for the main queue
  class Event_compare
    : public CGAL::cpp98::binary_function<bool,EventPtr,EventPtr>
  {
  public:
    Event_compare ( Self* aBuilder ) : mBuilder(aBuilder) {}

    bool operator() ( EventPtr const& aA, EventPtr const& aB ) const
    {
      return mBuilder->CompareEvents(aA,aB) == LARGER ;
    }

  private:
    Self* mBuilder ;
  } ;

  // Event compare for the local queues
  //
  // Special ordering for simultaneous split events (i.e. both (pseudo)splits + same time + same point)
  // to prevent impossible-to-untangle knots
  class Split_event_compare
    : public CGAL::cpp98::binary_function<bool,EventPtr,EventPtr>
  {
  public:
    Split_event_compare ( Self* aBuilder, Vertex_handle aV ) : mBuilder(aBuilder), mV(aV) {}

    bool operator() ( EventPtr const& aA, EventPtr const& aB ) const
    {
      CGAL_precondition( aA->type() != Event::cEdgeEvent || aB->type() != Event::cEdgeEvent ) ;

      if ( ! mBuilder->AreEventsSimultaneous(aA,aB) )
      {
        Comparison_result res = mBuilder->CompareEvents(aA,aB);
        if ( res == EQUAL )
          return aA < aB;
        return ( res == LARGER ) ;
      }

      // There are simultaneous events, we will need to refresh the queue before calling top()
      // see PopNextSplitEvent()
      mBuilder->GetVertexData(mV).mHasSimultaneousEvents = true ;

      // Priority queue comparison: `A` has higher priority than `B` if `operator()(A, B)` is `false`.
      // We want to give priority to smaller angles, so we must return `false` if the angle is smaller
      // i.e. `true` if the angle is larger
      Comparison_result res = mBuilder->CompareEventsSupportAngles(aA, aB);
      if ( res == EQUAL )
        return aA < aB;
      return ( res == LARGER ) ;
    }

  private:
    Self* mBuilder ;
    const Vertex_handle mV ;
  } ;

  typedef std::priority_queue<EventPtr,std::vector<EventPtr>,Split_event_compare> SplitPQ ;

  struct Vertex_data : public Ref_counted_base
  {
    Vertex_data ( Vertex_handle aVertex, Split_event_compare const& aComparer )
      :
        mVertex(aVertex)
      , mIsReflex(false)
      , mIsDegenerate(false)
      , mIsProcessed(false)
      , mIsExcluded(false)
      , mPrevInLAV(-1)
      , mNextInLAV(-1)
      , mNextSplitEventInMainPQ(false)
      , mHasSimultaneousEvents(false)
      , mSplitEvents(aComparer)
    {}

    Vertex_handle     mVertex ;
    bool              mIsReflex ;
    bool              mIsDegenerate ;
    bool              mIsProcessed ;
    bool              mIsExcluded ;
    int               mPrevInLAV ;
    int               mNextInLAV ;
    bool              mNextSplitEventInMainPQ ;
    bool              mHasSimultaneousEvents ;
    SplitPQ           mSplitEvents ;
    Triedge           mTriedge ; // Here, E0,E1 corresponds to the vertex (unlike *event* triedges)
    Trisegment_2_ptr  mTrisegment ; // Skeleton nodes cache the full trisegment tree that defines the originating event
  } ;

  typedef boost::intrusive_ptr<Vertex_data> Vertex_data_ptr ;

  typedef std::priority_queue<EventPtr,std::vector<EventPtr>,Event_compare> PQ ;

private :

  Halfedge_handle validate( Halfedge_handle aH ) const ;
  Vertex_handle   validate( Vertex_handle   aH ) const ;

  std::list<Vertex_handle>& GetHalfedgeLAVList(Halfedge_handle aH)
  {
    return mLAVLists[aH->id()];
  }

  Halfedge_handle SSkelEdgesPushBack(const Halfedge& aH1, const Halfedge& aH2)
  {
    mLAVLists.resize(aH2.id()+1);
    return mSSkel->SSkel::Base::edges_push_back (aH1, aH2);
  }

  void InitVertexData( Vertex_handle aV )
  {
    mVertexData.push_back( Vertex_data_ptr( new Vertex_data(aV, Split_event_compare(this, aV) ) ) ) ;
  }

  Vertex_data const& GetVertexData( Vertex_const_handle aV ) const
  {
    CGAL_precondition( handle_assigned(aV) ) ;
    return *mVertexData[aV->id()];
  }

  Vertex_data&  GetVertexData( Vertex_const_handle aV )
  {
    CGAL_precondition( handle_assigned(aV) ) ;
    return *mVertexData[aV->id()];
  }

  Vertex_handle GetVertex ( int aIdx )
  {
    CGAL_precondition(aIdx>=0);
    return mVertexData[aIdx]->mVertex ;
  }

  // Null if aV is a contour or infinite node
  Trisegment_2_ptr const& GetTrisegment ( Vertex_handle aV ) const
  {
    return GetVertexData(aV).mTrisegment ;
  }

  void SetTrisegment ( Vertex_handle aV, Trisegment_2_ptr const& aTrisegment )
  {
    GetVertexData(aV).mTrisegment = aTrisegment ;
  }

  // Null if aV is a contour node
  Triedge const& GetVertexTriedge ( Vertex_handle aV ) const
  {
    return GetVertexData(aV).mTriedge ;
  }

  void SetVertexTriedge ( Vertex_handle aV, Triedge const& aTriedge )
  {
    GetVertexData(aV).mTriedge = aTriedge ;
    GetHalfedgeLAVList(aTriedge.e0()).push_back(aV);
  }

  void GLAV_push_back ( Vertex_handle /* aV */ )
  {}

  void GLAV_remove ( Vertex_handle aV )
  {
    GetHalfedgeLAVList(GetVertexData(aV).mTriedge.e0()).remove(aV);
  }

  typename K::Segment_2 CreateRawSegment ( Halfedge_const_handle aH ) const
  {
    Point_2 s = aH->opposite()->vertex()->point() ;
    Point_2 t = aH->vertex()->point() ;
    return K().construct_segment_2_object()(s,t);
  }

  template <typename GT> // actually just equal to 'Traits', but gotta template for SFINAE to work
  Segment_2 CreateSegment ( Halfedge_const_handle aH,
                            typename std::enable_if<
                              ! CGAL_SS_i::has_Segment_2_with_ID<GT>::value>::type* = nullptr ) const
  {
    return Segment_2(CreateRawSegment(aH)) ;
  }

  template <typename GT> // actually just equal to 'Traits', but gotta template for SFINAE to work
  Segment_2 CreateSegment ( Halfedge_const_handle aH,
                            typename std::enable_if<
                              CGAL_SS_i::has_Segment_2_with_ID<GT>::value>::type* = nullptr ) const
  {
    return Segment_2(CreateRawSegment(aH), aH->id());
  }

  Vector_2 CreateVector ( Halfedge_const_handle aH ) const
  {
    Point_2 s = aH->opposite()->vertex()->point() ;
    Point_2 t = aH->vertex()->point() ;
    return K().construct_vector_2_object()(s,t);
  }

  Direction_2 CreateDirection ( Halfedge_const_handle aH ) const
  {
    return K().construct_direction_2_object()( CreateVector(aH) );
  }

  Trisegment_2_ptr CreateTrisegment ( Triedge const& aTriedge ) const
  {
    CGAL_precondition(aTriedge.is_valid() && aTriedge.is_skeleton());

    Trisegment_2_ptr r = Construct_ss_trisegment_2(mTraits)(CreateSegment<Traits>(aTriedge.e0())
                                                           ,CreateSegment<Traits>(aTriedge.e1())
                                                           ,CreateSegment<Traits>(aTriedge.e2())
                                                           );

    CGAL_STSKEL_BUILDER_TRACE(5,"Trisegment for " << aTriedge << ":" << r ) ;

    CGAL_postcondition_msg((r!= Trisegment_2_ptr()), "Unable to determine edges collinearity");

    return r ;
  }

  Trisegment_2_ptr CreateTrisegment ( Triedge const& aTriedge, Vertex_handle aLSeed ) const
  {
    Trisegment_2_ptr r = CreateTrisegment( aTriedge ) ;
    r->set_child_l( GetTrisegment(aLSeed) ) ;
    return r ;
  }

  Trisegment_2_ptr CreateTrisegment ( Triedge const& aTriedge, Vertex_handle aLSeed, Vertex_handle aRSeed ) const
  {
    Trisegment_2_ptr r = CreateTrisegment( aTriedge ) ;
    r->set_child_l( GetTrisegment(aLSeed) ) ;
    r->set_child_r( GetTrisegment(aRSeed) ) ;
    return r ;
  }

  Vertex_handle GetPrevInLAV ( Vertex_handle aV )
  {
    return GetVertex ( GetVertexData(aV).mPrevInLAV ) ;
  }

  Vertex_handle GetNextInLAV ( Vertex_handle aV )
  {
    return GetVertex ( GetVertexData(aV).mNextInLAV ) ;
  }

  void SetPrevInLAV ( Vertex_handle aV, Vertex_handle aPrev )
  {
    GetVertexData(aV).mPrevInLAV = aPrev->id();
  }

  void SetNextInLAV ( Vertex_handle aV, Vertex_handle aPrev )
  {
    GetVertexData(aV).mNextInLAV = aPrev->id();
  }

  Halfedge_handle GetEdgeEndingAt ( Vertex_handle aV )
  {
    return GetVertexTriedge(aV).e0();
  }

  Halfedge_handle GetEdgeStartingAt ( Vertex_handle aV )
  {
    return GetEdgeEndingAt( GetNextInLAV(aV) ) ;
  }

  void Exclude ( Vertex_handle aV )
  {
    GetVertexData(aV).mIsExcluded = true ;
  }
  bool IsExcluded ( Vertex_const_handle aV ) const
  {
    return GetVertexData(aV).mIsExcluded ;
  }

  void SetIsReflex ( Vertex_handle aV )
  {
    GetVertexData(aV).mIsReflex = true ;
  }

  bool IsReflex ( Vertex_handle aV )
  {
    return GetVertexData(aV).mIsReflex ;
  }

  void SetIsDegenerate ( Vertex_handle aV )
  {
    GetVertexData(aV).mIsDegenerate = true ;
  }

  bool IsDegenerate ( Vertex_handle aV )
  {
    return GetVertexData(aV).mIsDegenerate ;
  }

  void SetIsProcessed ( Vertex_handle aV )
  {
    GetVertexData(aV).mIsProcessed = true ;

    mVisitor.on_vertex_processed(aV);
  }

  bool IsConvex ( Vertex_handle aV )
  {
    Vertex_data const& lData = GetVertexData(aV) ;
    return !lData.mIsReflex && !lData.mIsDegenerate ;
  }

  bool IsProcessed ( Vertex_handle aV )
  {
    return GetVertexData(aV).mIsProcessed ;
  }

  void AddSplitEvent ( Vertex_handle aV, EventPtr const& aEvent )
  {
    CGAL_STSKEL_BUILDER_TRACE(2, "V" << aV->id() << " PQ: " << *aEvent);
    GetVertexData(aV).mSplitEvents.push(aEvent);
  }

  EventPtr PopNextSplitEvent ( Vertex_handle aV )
  {
    EventPtr rEvent ;
    Vertex_data& lData = GetVertexData(aV) ;
    if ( !lData.mNextSplitEventInMainPQ )
    {
      SplitPQ& lPQ = lData.mSplitEvents ;
      if ( !lPQ.empty() )
      {
        // When there are simultaneous split events, we sort them to handle nearby pseudo split events
        // together as to avoid multiple fronts crossing each other without seeing
        // and creating an invalid SLS.
        //
        // Unfortunately, the way that this sorting is performed requires knowing whether an event
        // is a pseudo split, and requires looking into the LAV tree. More annoyingly, whether an event
        // exists or not (see calls to handle_assigned(lOpp) and such) can change, which will
        // change the position of the item in the split priority queue (as an invalid event is
        // at the bottom of the priority queue).
        //
        // Consequently, when a simultaneous event has been detected, we need to refresh the queue
        // to get the real top.
        //
        // In practice, things are usually fine with a roughly-correct priority queue because an
        // invalid event being popped will just be ignored, but MSVC 2015-Debug is a very zealous
        // compiler that checks that the full queue is sane when it pops() and will detect this kind
        // of unimportant inconsistencies...
        if ( lData.mHasSimultaneousEvents )
        {
          std::make_heap(const_cast<EventPtr*>(&lPQ.top()),
                         const_cast<EventPtr*>(&lPQ.top()) + lPQ.size(),
                         Split_event_compare(this, aV));
        }

        rEvent = lPQ.top();
        lPQ.pop();
        lData.mNextSplitEventInMainPQ = true ;
      }
    }
    return rEvent ;
  }

  void AllowNextSplitEvent ( Vertex_handle aV )
  {
    GetVertexData(aV).mNextSplitEventInMainPQ = false ;
  }

  void InsertEventInPQ( EventPtr aEvent ) ;

  EventPtr PopEventFromPQ() ;

  bool ExistEvent ( Trisegment_2_ptr const& aS )
  {
    return Do_ss_event_exist_2(mTraits)(aS,mMaxTime);
  }

  bool IsOppositeEdgeFacingTheSplitSeed( Vertex_handle aSeed, Halfedge_handle aOpposite ) const
  {
    if ( aSeed->is_skeleton() )
         return Is_edge_facing_ss_node_2(mTraits)( GetTrisegment(aSeed), CreateSegment<Traits>(aOpposite) ) ;
    else return Is_edge_facing_ss_node_2(mTraits)( aSeed->point()      , CreateSegment<Traits>(aOpposite) ) ;
  }

  Oriented_side EventPointOrientedSide( Event const&          aEvent
                                      , Halfedge_const_handle aE0
                                      , Halfedge_const_handle aE1
                                      , Vertex_handle         aV01
                                      , bool                  aE0isPrimary
                                      ) const
  {
    return Oriented_side_of_event_point_wrt_bisector_2(mTraits)( aEvent.trisegment()
                                                               , CreateSegment<Traits>(aE0)
                                                               , CreateSegment<Traits>(aE1)
                                                               , GetTrisegment(aV01) // Can be null
                                                               , aE0isPrimary
                                                               ) ;
  }

  Comparison_result CompareEvents ( Trisegment_2_ptr const& aA, Trisegment_2_ptr const& aB ) const
  {
    return Compare_ss_event_times_2(mTraits)(aA,aB) ;
  }

  Comparison_result CompareEvents ( EventPtr const& aA, EventPtr const& aB ) const
  {
    return aA->triedge() != aB->triedge() ? CompareEvents( aA->trisegment(), aB->trisegment() ) : EQUAL ;
  }

  Comparison_result CompareEvents( Trisegment_2_ptr const& aTrisegment, Vertex_handle aSeedNode ) const
  {
    return aSeedNode->is_skeleton() ? aSeedNode->has_infinite_time() ? SMALLER
                                                                     : CompareEvents( aTrisegment, GetTrisegment(aSeedNode) )
                                    : LARGER  ;
  }

  Comparison_result CompareEvents( Vertex_handle aLNode, Vertex_handle aRNode ) const
  {
    if ( aLNode->is_skeleton() )
    {
      if ( aRNode->is_skeleton() )
      {
        if ( aLNode->has_infinite_time() && aRNode->has_infinite_time() )
          return EQUAL ;
        else if ( aLNode->has_infinite_time() )
          return LARGER ;
        else if ( aRNode->has_infinite_time() )
          return SMALLER ;
        else
          return CompareEvents ( GetTrisegment(aLNode), GetTrisegment(aRNode) ) ;
      }
      else // left skeleton, right not skeleton
      {
        return LARGER ;
      }
    }
    else // left not skeleton
    {
      if ( aRNode->is_skeleton() )
        return SMALLER ;
      else // both not skeleton
        return EQUAL ;
    }
  }

  bool AreEventsSimultaneous( Trisegment_2_ptr const& x, Trisegment_2_ptr const& y ) const
  {
    return Are_ss_events_simultaneous_2(mTraits)(x,y) ;
  }

  bool AreEventsSimultaneous( EventPtr const& x, EventPtr const& y ) const
  {
    return AreEventsSimultaneous(x->trisegment(),y->trisegment());
  }

  bool AreContourNodesCoincident( Vertex_handle aX, Vertex_handle aY ) const
  {
    CGAL_precondition( aX->is_contour() );
    CGAL_precondition( aY->is_contour() );

    return CGAL::compare_xy(aX->point(),aY->point()) == EQUAL ;
  }

  bool AreSkeletonNodesCoincident( Vertex_handle aX, Vertex_handle aY ) const
  {
    CGAL_precondition(  aX->is_skeleton() );
    CGAL_precondition(  aY->is_skeleton() );
    CGAL_precondition( !aX->has_infinite_time() );
    CGAL_precondition( !aY->has_infinite_time() );

    return AreEventsSimultaneous( GetTrisegment(aX),  GetTrisegment(aY) ) ;
  }

  void SetBisectorSlope ( Halfedge_handle aBisector, Sign aSlope )
  {
    aBisector->HBase_base::set_slope(aSlope);
  }

  void SetBisectorSlope ( Vertex_handle aA, Vertex_handle aB )
  {
    Halfedge_handle lOBisector = aA->primary_bisector();
    Halfedge_handle lIBisector = lOBisector->opposite();

    CGAL_precondition( !aA->is_contour() || !aB->is_contour() ) ;

    if ( aA->is_contour() )
    {
      SetBisectorSlope(lOBisector,POSITIVE);
      SetBisectorSlope(lIBisector,NEGATIVE);
    }
    else if ( aB->is_contour())
    {
      SetBisectorSlope(lOBisector,NEGATIVE);
      SetBisectorSlope(lIBisector,POSITIVE);
    }
    else
    {
      if ( aA->has_infinite_time() )
      {
        CGAL_precondition( !aB->has_infinite_time());

        SetBisectorSlope(lOBisector,NEGATIVE);
        SetBisectorSlope(lIBisector,POSITIVE);
      }
      else if ( aB->has_infinite_time())
      {
        CGAL_precondition( !aA->has_infinite_time());

        SetBisectorSlope(lOBisector,NEGATIVE);
        SetBisectorSlope(lIBisector,POSITIVE);
      }
      else
      {
        CGAL_precondition( !aA->has_infinite_time());
        CGAL_precondition( !aB->has_infinite_time());

        Sign lSlope = CompareEvents(GetTrisegment(aB),GetTrisegment(aA));
        SetBisectorSlope(lOBisector,lSlope);
        SetBisectorSlope(lIBisector,opposite(lSlope));
      }
    }
  }

  boost::tuple<FT,Point_2> ConstructEventTimeAndPoint( Trisegment_2_ptr const& aS ) const
  {
    boost::optional< boost::tuple<FT,Point_2> > r = Construct_ss_event_time_and_point_2(mTraits)(aS);
    CGAL_postcondition_msg(!!r, "Unable to compute skeleton node coordinates");
    return *r ;
  }

  void SetEventTimeAndPoint( Event& aE )
  {
    FT      lTime ;
    Point_2 lP ;
    boost::tie(lTime,lP) = ConstructEventTimeAndPoint(aE.trisegment());

    aE.SetTimeAndPoint(lTime,lP);
  }


  void EraseBisector( Halfedge_handle aB )
  {
    CGAL_STSKEL_BUILDER_TRACE(1,"Dangling B" << aB->id() << " and B" << aB->opposite()->id() << " removed.");

    mSSkel->SSkel::Base::edges_erase(aB);
  }

#ifdef CGAL_STSKEL_TRACE_ON
  std::string wavefront2str( Vertex_handle v )
  {
    std::ostringstream ss ;

    ss << "N" << GetPrevInLAV(v)->id() << "->N" << v->id() << "->N" << GetNextInLAV(v)->id()
       << "  E" << GetVertexTriedge(v).e0()->id() << "->E" << GetVertexTriedge(v).e1()->id() ;

    return ss.str() ;

  }
#endif

  void Link( Halfedge_handle aH, Face_handle aF )
  {
    aH->HBase_base::set_face(aF);
  }

  void Link( Halfedge_handle aH, Vertex_handle aV )
  {
    aH->HBase_base::set_vertex(aV);
  }

  void Link( Vertex_handle aV, Halfedge_handle aH )
  {
    aV->VBase::set_halfedge(aH);
  }

  void CrossLinkFwd( Halfedge_handle aPrev, Halfedge_handle aNext )
  {
    aPrev->HBase_base::set_next(aNext);
    aNext->HBase_base::set_prev(aPrev);
  }

  void CrossLink( Halfedge_handle aH, Vertex_handle aV )
  {
    Link(aH,aV);
    Link(aV,aH);
  }

  Triedge GetCommonTriedge( Vertex_handle aA, Vertex_handle aB ) ;

  bool AreBisectorsCoincident ( Halfedge_const_handle aA, Halfedge_const_handle aB ) const ;

  EventPtr IsPseudoSplitEvent( EventPtr const& aEvent, Vertex_handle_pair aOpp, Site const& aSite ) ;

  void CollectSplitEvent( Vertex_handle aNode, Triedge const& aTriedge ) ;

  void CollectSplitEvents( Vertex_handle aNode, Triedge const& aPrevEventTriedge ) ;

  EventPtr FindEdgeEvent( Vertex_handle aLNode, Vertex_handle aRNode, Triedge const& aPrevEventTriedge  ) ;

  void HandleSimultaneousEdgeEvent( Vertex_handle aA, Vertex_handle aB ) ;

  void CollectNewEvents( Vertex_handle aNode, Triedge const& aPrevEventTriedge ) ;
  void UpdatePQ( Vertex_handle aV, Triedge const& aPrevEventTriedge ) ;
  void CreateInitialEvents();
  void CreateContourBisectors();
  void HarmonizeSpeeds(boost::mpl::bool_<false>) { }
  void HarmonizeSpeeds(boost::mpl::bool_<true>);

  // @fixme This function can create an inconsistency between the skeleton and the offset
  // when caching is used: this function will harmonize values in the cache used by the skeleton
  // builder, but not the values in the cache used in the offset builder (actually, no_cache is used
  // for the offset builder, at the moment). It is difficult to harmonize between caches as one might
  // wish to build the skeleton and the offset in completely different moments and independent functions...
  void HarmonizeSpeeds() {
    // Harmonize speed only if we have a segment type with id and a not exact square root
    return HarmonizeSpeeds(
      boost::mpl::bool_< CGAL_SS_i::has_Segment_2_with_ID<Traits>::value &&
                         ( !is_same_or_derived<Field_with_sqrt_tag,
                                               typename Algebraic_structure_traits<FT>::Algebraic_category
                          >::value || std::is_floating_point<FT>::value) >() );
  }
  void InitPhase();

  void SetupNewNode( Vertex_handle aNode );

  Vertex_handle_pair LookupOnSLAV ( Halfedge_handle aOBorder, EventPtr const& aEvent, Site& rSite ) ;

  Vertex_handle      ConstructEdgeEventNode         ( EdgeEvent&   aEvent ) ;
  Vertex_handle_pair ConstructSplitEventNodes       ( SplitEvent&  aEvent, Vertex_handle aOppR ) ;
  Vertex_handle_pair ConstructPseudoSplitEventNodes ( PseudoSplitEvent& aEvent ) ;

  bool IsValidEvent            ( EventPtr                aEvent ) ;
  bool IsValidEdgeEvent        ( EdgeEvent        const& aEvent ) ;
  bool IsValidSplitEvent       ( SplitEvent       const& aEvent ) ;
  bool IsValidPseudoSplitEvent ( PseudoSplitEvent const& aEvent ) ;

  void HandleEdgeEvent               ( EventPtr aEvent ) ;
  void HandleSplitEvent              ( EventPtr aEvent, Vertex_handle_pair aOpp ) ;
  void HandlePseudoSplitEvent        ( EventPtr aEvent ) ;
  void HandleSplitOrPseudoSplitEvent ( EventPtr aEvent ) ;

  void InsertNextSplitEventInPQ( Vertex_handle v ) ;
  void InsertNextSplitEventsInPQ() ;

  bool IsProcessed( EventPtr aEvent ) ;

  void Propagate();

  void EraseNode ( Vertex_handle aNode ) ;

  void MergeSplitNodes ( Vertex_handle_pair aSplitNodes ) ;

  void RelinkBisectorsAroundMultinode( Vertex_handle const& v0, Halfedge_handle_vector& aLinks ) ;

  void PreprocessMultinode( Multinode& aMN ) ;

  void ProcessMultinode( Multinode& aMN, Halfedge_handle_vector& rHalfedgesToRemove , Vertex_handle_vector& rNodesToRemove ) ;

  MultinodePtr CreateMultinode( Halfedge_handle begin, Halfedge_handle end ) ;

  // returns 'true' if something was merged
  bool MergeCoincidentNodes() ;

  bool FinishUp();

  bool Run();

private:

  //Input
  Traits  mTraits ;

  Visitor const& mVisitor ;

  std::vector<Vertex_data_ptr> mVertexData ;

  std::vector<std::list<Vertex_handle>> mLAVLists;

  Vertex_handle_vector   mReflexVertices ;
  Halfedge_handle_vector mDanglingBisectors ;
  Halfedge_handle_vector mContourHalfedges ;


  Vertex_handle_pair_vector mSplitNodes ;

  Event_compare mEventCompare ;

  int mVertexID ;
  int mEdgeID   ;
  int mFaceID   ;
  int mEventID  ;
  int mStepID   ;

  boost::optional<FT> mMaxTime ;

  PQ mPQ ;

  //Output
  SSkelPtr mSSkel ;

private :

  template<class InputPointIterator, class Converter>
  void enter_valid_contour ( InputPointIterator aBegin, InputPointIterator aEnd, Converter const& cvt )
  {
    CGAL_STSKEL_BUILDER_TRACE(0,"Inserting Connected Component of the Boundary....");

    Halfedge_handle lFirstCCWBorder ;
    Halfedge_handle lPrevCCWBorder ;
    Halfedge_handle lNextCWBorder ;
    Vertex_handle   lFirstVertex ;
    Vertex_handle   lPrevVertex ;

    InputPointIterator lCurr = aBegin ;
    // InputPointIterator lPrev = aBegin ;

    int c = 0 ;

    while ( lCurr != aEnd )
    {
      Halfedge_handle lCCWBorder = SSkelEdgesPushBack( Halfedge(mEdgeID),Halfedge(mEdgeID+1)  );
      Halfedge_handle lCWBorder = lCCWBorder->opposite();
      mEdgeID += 2 ;

      mContourHalfedges.push_back(lCCWBorder);

      Vertex_handle lVertex = mSSkel->SSkel::Base::vertices_push_back( Vertex(mVertexID++,cvt(*lCurr)) ) ;
      CGAL_STSKEL_BUILDER_TRACE(1,"Vertex: V" << lVertex->id() << " at " << lVertex->point() );
      InitVertexData(lVertex);

      Face_handle lFace = mSSkel->SSkel::Base::faces_push_back( Face(mFaceID++) ) ;

      ++ c ;

      lCCWBorder->HBase_base::set_face(lFace);
      lFace     ->FBase     ::set_halfedge(lCCWBorder);

      lVertex   ->VBase     ::set_halfedge(lCCWBorder);
      lCCWBorder->HBase_base::set_vertex(lVertex);

      if ( lCurr == aBegin )
      {
        lFirstVertex    = lVertex ;
        lFirstCCWBorder = lCCWBorder ;
      }
      else
      {
        SetPrevInLAV(lVertex    ,lPrevVertex);
        SetNextInLAV(lPrevVertex,lVertex    );

        SetVertexTriedge(lPrevVertex, Triedge(lPrevCCWBorder,lCCWBorder) ) ;

        lCWBorder->HBase_base::set_vertex(lPrevVertex);

        lCCWBorder    ->HBase_base::set_prev(lPrevCCWBorder);
        lPrevCCWBorder->HBase_base::set_next(lCCWBorder);

        lNextCWBorder->HBase_base::set_prev(lCWBorder);
        lCWBorder    ->HBase_base::set_next(lNextCWBorder);

        CGAL_STSKEL_BUILDER_TRACE(1,"CCW Border: E" << lCCWBorder->id() << ' ' << lPrevVertex->point() << " -> " << lVertex    ->point());
        CGAL_STSKEL_BUILDER_TRACE(1,"CW  Border: E" << lCWBorder ->id() << ' ' << lVertex    ->point() << " -> " << lPrevVertex->point() );

        mVisitor.on_contour_edge_entered(lCCWBorder);

      }

      // lPrev = lCurr ;

      ++ lCurr ;

      lPrevVertex    = lVertex ;
      lPrevCCWBorder = lCCWBorder ;
      lNextCWBorder  = lCWBorder ;
    }

    SetPrevInLAV(lFirstVertex,lPrevVertex );
    SetNextInLAV(lPrevVertex ,lFirstVertex);

    SetVertexTriedge( lPrevVertex, Triedge(lPrevCCWBorder,lFirstCCWBorder) ) ;

    lFirstCCWBorder->opposite()->HBase_base::set_vertex(lPrevVertex);

    lFirstCCWBorder->HBase_base::set_prev(lPrevCCWBorder);
    lPrevCCWBorder ->HBase_base::set_next(lFirstCCWBorder);

    lPrevCCWBorder ->opposite()->HBase_base::set_prev(lFirstCCWBorder->opposite());
    lFirstCCWBorder->opposite()->HBase_base::set_next(lPrevCCWBorder ->opposite());

    CGAL_precondition_msg(c >=3, "The contour must have at least 3 _distinct_ vertices" ) ;

    CGAL_STSKEL_BUILDER_TRACE(1
                             , "CCW Border: E" << lFirstCCWBorder->id()
                             << ' ' << lPrevVertex ->point() << " -> " << lFirstVertex->point() << '\n'
                             << "CW  Border: E" << lFirstCCWBorder->opposite()->id()
                             << ' ' << lFirstVertex->point() << " -> " << lPrevVertex ->point()
                             );

    mVisitor.on_contour_edge_entered(lFirstCCWBorder);
  }

// internal function to filter split event in case the traits is Filtered
  bool CanSafelyIgnoreSplitEventImpl(const EventPtr&, boost::mpl::bool_<false>) const
  {
    return false;
  }

  bool CanSafelyIgnoreSplitEventImpl(const EventPtr& lEvent, boost::mpl::bool_<true>) const
  {
    return mTraits.CanSafelyIgnoreSplitEvent(lEvent);
  }

  bool CanSafelyIgnoreSplitEvent(const EventPtr& lEvent) const
  {
    return CanSafelyIgnoreSplitEventImpl(lEvent, typename CGAL_SS_i::has_Filters_split_events_tag<Traits>::type());
  }

  void ComputeUpperBoundForValidSplitEventsImpl(Vertex_handle, Vertex_handle, Vertex_handle,
                                                Halfedge_handle_vector_iterator, Halfedge_handle_vector_iterator,
                                                boost::mpl::bool_<false>) const
  {
  }

  void ComputeUpperBoundForValidSplitEventsImpl(Vertex_handle lPrev, Vertex_handle aNode, Vertex_handle lNext,
                                                Halfedge_handle_vector_iterator contour_halfedges_begin,
                                                Halfedge_handle_vector_iterator contour_halfedges_end,
                                                boost::mpl::bool_<true>) const
  {
    return mTraits.ComputeFilteringBound(lPrev, aNode, lNext, contour_halfedges_begin, contour_halfedges_end);
  }

  void ComputeUpperBoundForValidSplitEvents(Vertex_handle lPrev, Vertex_handle aNode, Vertex_handle lNext,
                                            Halfedge_handle_vector_iterator contour_halfedges_begin,
                                            Halfedge_handle_vector_iterator contour_halfedges_end) const
  {
    return ComputeUpperBoundForValidSplitEventsImpl(lPrev, aNode, lNext, contour_halfedges_begin, contour_halfedges_end,
                                                    typename CGAL_SS_i::has_Filters_split_events_tag<Traits>::type());
  }

public:

  // This compares INPUT vertices so there is no need to filter it.
  struct AreVerticesEqual
  {
    template<class P>
    bool operator() ( P const&x, P const& y ) const
    {
      return CGAL::compare_xy(x,y) == EQUAL ;
    }
  } ;

  template<class InputPointIterator, class Converter>
  Straight_skeleton_builder_2& enter_contour ( InputPointIterator aBegin
                                             , InputPointIterator aEnd
                                             , Converter const&   cvt
                                             , bool               aCheckValidity = true
                                             )
  {
    if ( aCheckValidity )
    {
      typedef typename std::iterator_traits<InputPointIterator>::value_type Input_point;

      typedef std::vector<Input_point> Input_point_vector ;

      // Remove coincident consecutive vertices
      Input_point_vector lList;
      std::unique_copy(aBegin,aEnd,std::back_inserter(lList),AreVerticesEqual());

      while ( lList.size() > 0 && CGAL::compare_xy(lList.front(),lList.back()) == EQUAL )
        lList.pop_back();

      if ( lList.size() >= 3 )
      {
        enter_valid_contour(lList.begin(),lList.end(),cvt);
      }
      else
      {
        CGAL_STSKEL_BUILDER_TRACE(0,"Degenerate contour (less than 3 non-degenerate vertices).");
      }
    }
    else
    {
      enter_valid_contour(aBegin,aEnd,cvt);
    }

    return *this ;
  }

  template<class InputPointIterator>
  Straight_skeleton_builder_2& enter_contour ( InputPointIterator aBegin
                                             , InputPointIterator aEnd
                                             , bool               aCheckValidity = true
                                             )
  {
    return enter_contour(aBegin, aEnd, Cartesian_converter<K,K>(), aCheckValidity);
  }

} ;

} // end namespace CGAL

#include <CGAL/Straight_skeleton_2/Straight_skeleton_builder_2_impl.h>

#endif // CGAL_STRAIGHT_SKELETON_BUILDER_2_H //
// EOF //
