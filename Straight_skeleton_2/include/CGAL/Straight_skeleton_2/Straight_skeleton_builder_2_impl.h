// Copyright (c) 2005, 2006 Fernando Luis Cacciola Carballal. All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Fernando Cacciola <fernando_cacciola@ciudad.com.ar>
//
#ifndef CGAL_STRAIGHT_SKELETON_BUILDER_2_IMPL_H
#define CGAL_STRAIGHT_SKELETON_BUILDER_2_IMPL_H 1

#include <CGAL/license/Straight_skeleton_2.h>

#include <CGAL/number_type_config.h>
#include <CGAL/Real_timer.h>
#include <CGAL/Unique_hash_map.h>

#include <boost/bind.hpp>
#include <boost/utility.hpp>
#include <boost/version.hpp>
#if BOOST_VERSION == 106000
//ice_not is deprecated in boost 1.60 but used within adjacency_matrix.hpp
#include <boost/type_traits/detail/ice_not.hpp>
#endif
#include <boost/graph/adjacency_matrix.hpp>

#if defined(BOOST_MSVC)
#  pragma warning(push)
#  pragma warning(disable:4355) // complaint about using 'this' to initialize a member
#endif

#include <algorithm>
#include <iostream>
#include <sstream>

namespace CGAL {

template<class Gt, class Ss, class V>
Straight_skeleton_builder_2<Gt,Ss,V>::Straight_skeleton_builder_2 ( boost::optional<FT> aMaxTime, Traits const& aTraits, Visitor const& aVisitor )
  :
  mTraits(aTraits)
 ,mVisitor(aVisitor)
 ,mEventCompare(this)
 ,mVertexID(0)
 ,mEdgeID(0)
 ,mFaceID(0)
 ,mEventID(0)
 ,mStepID(0)
 ,mMaxTime(aMaxTime)
 ,mPQ(mEventCompare)
 ,mSSkel( new SSkel() )
{
}

template<class Gt, class Ss, class V>
typename Straight_skeleton_builder_2<Gt,Ss,V>::Halfedge_handle
Straight_skeleton_builder_2<Gt,Ss,V>::validate( Halfedge_handle aH ) const
{
  if ( !handle_assigned(aH) )
    throw std::runtime_error("Incomplete straight skeleton");
  return aH ;
}

template<class Gt, class Ss, class V>
typename Straight_skeleton_builder_2<Gt,Ss,V>::Vertex_handle
Straight_skeleton_builder_2<Gt,Ss,V>::validate( Vertex_handle aH ) const
{
  if ( !handle_assigned(aH) )
    throw std::runtime_error("Incomplete straight skeleton");
  return aH ;
}

template<class Gt, class Ss, class V>
void Straight_skeleton_builder_2<Gt,Ss,V>::InsertEventInPQ( EventPtr aEvent )
{
  mPQ.push(aEvent);
  CGAL_STSKEL_BUILDER_TRACE(4, "Enque: " << *aEvent);
}

template<class Gt, class Ss, class V>
typename Straight_skeleton_builder_2<Gt,Ss,V>::EventPtr
Straight_skeleton_builder_2<Gt,Ss,V>::PopEventFromPQ()
{
  EventPtr rR = mPQ.top(); mPQ.pop();
  return rR ;
}

// Tests whether there is an edge event between the 3 contour edges defining nodes 'aLnode' and 'aRNode'.
// If such event exits and is not in the past, it's returned. Otherwise the result is null.
//
template<class Gt, class Ss, class V>
typename Straight_skeleton_builder_2<Gt,Ss,V>::EventPtr
Straight_skeleton_builder_2<Gt,Ss,V>::FindEdgeEvent( Vertex_handle aLNode, Vertex_handle aRNode, Triedge const& aPrevEventTriedge  )
{
  EventPtr rResult ;

  CGAL_STSKEL_BUILDER_TRACE(4, "FindEdgeEvent(), Left/Right Nodes: N" << aLNode->id() << " N" << aRNode->id() ) ;

  Triedge lTriedge = GetVertexTriedge(aLNode) & GetVertexTriedge(aRNode) ;

  if ( lTriedge.is_valid() && lTriedge != aPrevEventTriedge )
  {
    Trisegment_2_ptr lTrisegment = CreateTrisegment(lTriedge,aLNode,aRNode);

    // The 02 collinearity configuration is problematic: 01 or 12 collinearity has a seed position
    // giving the point through which the bisector passes. However, for 02, it is not a given.
    //
    // If the seed exists, the information is passed to the traits as the "third" child of the trisegment.
    // Otherwise, ignore this as it should re-appear when the seed of 02 is created.
    //
    // Note that this is only for edge events; (pseudo-)split events are not concerned.
    if ( lTrisegment->collinearity() == TRISEGMENT_COLLINEARITY_02 )
    {
      // Check in the SLAV if the seed corresponding to 02 exists
      Vertex_handle lPrevNode = GetPrevInLAV(aLNode) ;
      CGAL_assertion( GetEdgeStartingAt(lPrevNode) == lTriedge.e0() ) ;

      if ( GetEdgeEndingAt(lPrevNode) == lTriedge.e2() )
      {
        // Note that this can be a contour node and in that case GetTrisegment is null and we get
        // the middle point, but in that case e2 and e0 are consecutive in the input
        // and the middle point is the common extremity and things are fine.
        lTrisegment->set_child_t( GetTrisegment(lPrevNode) ) ;
      }
      else
      {
        Orientation lOrientationS = CGAL::orientation( lTrisegment->e0().source(), lTrisegment->e0().target(), lTrisegment->e1().source() ) ;
        Orientation lOrientationT = CGAL::orientation( lTrisegment->e0().source(), lTrisegment->e0().target(), lTrisegment->e1().target() ) ;
        if ( lOrientationS != LEFT_TURN && lOrientationT != LEFT_TURN )
        {
          // Reasonning is: if the middle halfedge (e1) is "below" e0 and e2, then there is some
          // kind of concavity in between e0 and e2. This concavity will resolve itself and either:
          // - e0 and e2 will never meet, but in that case we would not be here
          // - e0 and e2 will meet. In that case, we can ignore all the details of the concavity
          //   and simply consider that in the end, all that matters is the e0, e2, next(e0),
          //   and prev(e2). In that case, we get two bisectors issued from e0 and e2, and one
          //   bisector issued from some seed S and splitting next(e0) and prev(e2). This can also
          //   be seen as two exterior bisectors and one interior bisector of a triangle
          //   target(e0) -- S - source(e2). It is a known result that these three bisectors
          //   meet in a single point. Thus, when we get here e0-e1-e2, we know that
          //   these will meet in a single, existing point, either the left or the right child (the oldest).

          if ( CompareEvents(aLNode, aRNode) == SMALLER )
            lTrisegment->set_child_t( GetTrisegment(aRNode) ) ;
          else
            lTrisegment->set_child_t( GetTrisegment(aLNode) ) ;
        }
        else
        {
          return rResult;
        }
      }
    }

    if ( ExistEvent(lTrisegment) )
    {
      Comparison_result lLNodeD = CompareEvents(lTrisegment,aLNode) ;
      Comparison_result lRNodeD = CompareEvents(lTrisegment,aRNode) ;

      if ( lLNodeD != SMALLER && lRNodeD != SMALLER )
      {
        rResult = EventPtr( new EdgeEvent( lTriedge, lTrisegment, aLNode, aRNode ) ) ;

        mVisitor.on_edge_event_created(aLNode, aRNode) ;

        CGAL_STSKEL_DEBUG_CODE( SetEventTimeAndPoint(*rResult) );
      }
      else
      {
        CGAL_STSKEL_BUILDER_TRACE(4, "Edge event: " << lTriedge << " is in the past. Compared to L=" << lLNodeD << " to R=" << lRNodeD ) ;
      }
    }
  }
  return rResult ;
}


template<class Gt, class Ss, class V>
typename Straight_skeleton_builder_2<Gt,Ss,V>::EventPtr
Straight_skeleton_builder_2<Gt,Ss,V>::IsPseudoSplitEvent( EventPtr const& aEvent, Vertex_handle_pair aOpp, Site const& aSite )
{
  EventPtr rPseudoSplitEvent ;

  if ( aSite != INSIDE )
  {
    SplitEvent& lEvent = dynamic_cast<SplitEvent&>(*aEvent) ;

    Triedge           const& lEventTriedge    = lEvent.triedge();
    Trisegment_2_ptr  const& lEventTrisegment = lEvent.trisegment();
    Vertex_handle            lSeedN           = lEvent.seed0();

    Vertex_handle lOppL = aOpp.first ;
    Vertex_handle lOppR = aOpp.second ;

    if ( aSite == AT_SOURCE )
    {
      Halfedge_handle lOppPrevBorder = GetVertexTriedge(lOppL).e0() ;

      if ( lEventTriedge.e0() != lOppPrevBorder && lEventTriedge.e1() != lOppPrevBorder )
      {
        rPseudoSplitEvent = EventPtr( new PseudoSplitEvent(lEventTriedge,lEventTrisegment,lOppL,lSeedN,true) ) ;

        CGAL_STSKEL_BUILDER_TRACE(1,"Pseudo-split-event found against " << v2str(*lOppL) ) ;

        mVisitor.on_pseudo_split_event_created(lOppL,lSeedN) ;
      }
    }
    else // aSite == AT_TARGET
    {
      Vertex_handle lOppNextN = GetNextInLAV(lOppR) ;

      Halfedge_handle lOppNextBorder = GetVertexTriedge(lOppNextN).e0() ;

      if ( lEventTriedge.e0() != lOppNextBorder && lEventTriedge.e1() != lOppNextBorder )
      {
        rPseudoSplitEvent = EventPtr( new PseudoSplitEvent(lEventTriedge, lEventTrisegment, lSeedN, lOppR,false) ) ;

        CGAL_STSKEL_BUILDER_TRACE(1,"Pseudo-split-event found against " << v2str(*lOppR) ) ;

        mVisitor.on_pseudo_split_event_created(lSeedN,lOppR) ;
      }
    }
  }

  if ( rPseudoSplitEvent )
    rPseudoSplitEvent->SetTimeAndPoint(aEvent->time(),aEvent->point());

  return rPseudoSplitEvent ;
}

// Tests whether there is a split event between the contour edges (aReflexLBorder,aReflexRBorder,aOppositeBorder).
// If such event exits and is not in the past, it's returned. Otherwise the result is null
// 'aReflexLBorder' and 'aReflexRBorder' are consecutive contour edges which 'aNode' as the vertex.
// 'aOppositeBorder' is some other edge in the polygon which, if the event exists, is split by the reflex wavefront.
//
// NOTE: 'aNode' can be a skeleton node (an interior split event produced by a previous vertex event). In that case,
// the 'reflex borders' are not consecutive in the input polygon but they are in the corresponding offset polygon that
// contains aNode as a vertex.
//
template<class Gt, class Ss, class V>
void Straight_skeleton_builder_2<Gt,Ss,V>::CollectSplitEvent( Vertex_handle aNode, Triedge const& aTriedge )
{
  CGAL_STSKEL_BUILDER_TRACE(3, "Collect SplitEvent for N" << aNode->id() << " triedge: " << aTriedge);

  if ( IsOppositeEdgeFacingTheSplitSeed(aNode,aTriedge.e2()) )
  {
    Trisegment_2_ptr lTrisegment = CreateTrisegment(aTriedge,aNode);

    if ( lTrisegment->collinearity() != TRISEGMENT_COLLINEARITY_02 && ExistEvent(lTrisegment) )
    {
      if ( CompareEvents(lTrisegment,aNode) != SMALLER )
      {
        EventPtr lEvent = EventPtr( new SplitEvent (aTriedge,lTrisegment,aNode) ) ;

        // filter split event
        if (CanSafelyIgnoreSplitEvent(lEvent))
          return;

        mVisitor.on_split_event_created(aNode) ;

        CGAL_STSKEL_DEBUG_CODE( SetEventTimeAndPoint(*lEvent) ) ;

        AddSplitEvent(aNode,lEvent);
      }
    }
  }
}

// Tests the reflex wavefront emerging from 'aNode' against the other contour edges in search for split events.
template<class Gt, class Ss, class V>
void Straight_skeleton_builder_2<Gt,Ss,V>::CollectSplitEvents( Vertex_handle aNode, Triedge const& aPrevEventTriedge )
{
  // lLBorder and lRBorder are the consecutive contour edges forming the reflex wavefront.
  Triedge const& lTriedge = GetVertexTriedge(aNode);

  Halfedge_handle lLBorder = lTriedge.e0();
  Halfedge_handle lRBorder = lTriedge.e1();

  CGAL_STSKEL_BUILDER_TRACE(3
                      ,"Finding SplitEvent for N" << aNode->id()
                      << " LBorder: E" << lLBorder->id() << " RBorder: E" << lRBorder->id()
                      );

  ComputeUpperBoundForValidSplitEvents(GetPrevInLAV(aNode), aNode, GetNextInLAV(aNode),
                                       mContourHalfedges.begin(), mContourHalfedges.end());

  for ( Halfedge_handle_vector_iterator i = mContourHalfedges.begin(); i != mContourHalfedges.end(); ++ i )
  {
    Halfedge_handle lOpposite = *i ;

    if ( lOpposite != lLBorder && lOpposite != lRBorder )
    {
      Triedge lEventTriedge(lLBorder, lRBorder, lOpposite);

      if ( lEventTriedge != aPrevEventTriedge )
      {
        CollectSplitEvent(aNode, lEventTriedge) ;
      }
    }
  }

#ifdef CGAL_STRAIGHT_SKELETON_ENABLE_TRACE
  std::cout << "  local queue size --> " << GetVertexData(aNode).mSplitEvents.size() << std::endl;
#endif
}


// Finds and enques all the new potential events produced by the vertex wavefront emerging from 'aNode' (which can be a reflex wavefront).
// This new events are simply stored in the priority queue, not processed.
template<class Gt, class Ss, class V>
void Straight_skeleton_builder_2<Gt,Ss,V>::CollectNewEvents( Vertex_handle aNode, Triedge const& aPrevEventTriedge )
{
  // A Straight Skeleton is the trace of the 'grassfire propagation' that corresponds to the inward move of all the vertices
  // of a polygon along their angular bisectors.
  // Since vertices are the common endpoints of contour edges, the propagation corresponds to contour edges moving inward,
  // shrinking and expanding as neccesasry to keep the vertices along the angular bisectors.
  // At each instant in time the current location of vertices (and edges) describe the current 'Offset polygon'
  // (with at time zero corresponds to the input polygon).
  //
  // An 'edge wavefront' is a moving contour edge.
  // A 'vertex wavefront' is the wavefront of two consecutive edge wavefronts (sharing a moving vertex).
  //
  // An 'Event' is the collision of 2 wavefronts.
  // Each event changes the topology of the shrinking polygon; that is, at the event, the current polygon differs from the
  // inmediately previous polygon in the number of vertices.
  //
  // If 2 vertex wavefronts sharing a common edge collide, the event is called an edge event. At the time of the event, the current
  // polygon doex not have the common edge anynmore, and the two vertices become one. This new 'skeleton' vertex generates a new
  // vertex wavefront which can further collide with other wavefronts, producing for instance, more edge events.
  //
  // If a refex vertex wavefront collide with an edge wavefront, the event is called a split event. At the time of the event, the current
  // polygon is split in two unconnected polygons, each one containing a portion of the edge hit and split by the reflex wavefront.
  //
  // If 2 reflex wavefronts collide each other, the event is called a vertex event. At the time of the event, the current polygon
  // is split in two unconnected polygons. Each one contains a different combination of the colliding reflex edges. That is, if the
  // wavefront (edgea,edgeb) collides with (edgec,edged), the two resulting polygons will contain (edgea,edgec) and (edgeb,edged).
  // Furthermore, one of the new vertices can be a reflex vertex generating a reflex wavefront which can further produce more split or
  // vertex events (or edge events of course)
  //
  // Each vertex wavefront (reflex or not) results in one and only one event from a set of possible events.
  // It can result in a edge event against the vertex wavefronts emerging from the adjacent vertices (in the current polygon, not
  // in the input polygon); or it can result in a split event (or vertex event) against any other wavefront in the rest of
  // current polygon.


  // Adjacent vertices in the current polygon containing aNode (called LAV)
  Vertex_handle lPrev = GetPrevInLAV(aNode) ;
  Vertex_handle lNext = GetNextInLAV(aNode) ;

  CGAL_STSKEL_BUILDER_TRACE
    ( 2
    , "Collecting new events generated by N" << aNode->id() << " at " << aNode->point() << " (Prev: N" << lPrev->id() << " Next: N"
       << lNext->id() << ")"
    ) ;

  if ( IsReflex(aNode) )
    CollectSplitEvents(aNode, aPrevEventTriedge) ;

  EventPtr lLEdgeEvent = FindEdgeEvent( lPrev , aNode, aPrevEventTriedge ) ;
  EventPtr lREdgeEvent = FindEdgeEvent( aNode , lNext, aPrevEventTriedge ) ;

  bool lAcceptL = !!lLEdgeEvent ;
  bool lAcceptR = !!lREdgeEvent ;

  if ( lAcceptL )
    InsertEventInPQ(lLEdgeEvent);

  if ( lAcceptR )
    InsertEventInPQ(lREdgeEvent);
}

// Handles the special case of two simultaneous edge events, that is, two edges
// collapsing along the line/point were they meet at the same time.
// This ocurrs when the bisector emerging from vertex 'aA' is defined by the same pair of
// contour edges as the bisector emerging from vertex 'aB' (but in opposite order).
//
template<class Gt, class Ss, class V>
void Straight_skeleton_builder_2<Gt,Ss,V>::HandleSimultaneousEdgeEvent( Vertex_handle aA, Vertex_handle aB )
{
  CGAL_STSKEL_BUILDER_TRACE ( 2, "Handling simultaneous EdgeEvent between N" << aA ->id() << " and N"  << aB ->id() ) ;

  mVisitor.on_anihiliation_event_processed(aA,aB) ;

  Halfedge_handle lOA = aA->primary_bisector() ;
  Halfedge_handle lOB = aB->primary_bisector() ;
  Halfedge_handle lIA = lOA->opposite();
  Halfedge_handle lIB = lOB->opposite();

  Vertex_handle lOAV = lOA->vertex() ;
  Vertex_handle lIAV = lIA->vertex() ;
  Vertex_handle lOBV = lOB->vertex() ;
  // Vertex_handle lIBV = lIB->vertex() ;

  CGAL_STSKEL_BUILDER_TRACE ( 2
                            ,    "OA: B" << lOA->id() << '\n'
                              << "IA: B" << lIA->id() << '\n'
                              << "OB: B" << lOB->id() << '\n'
                              << "IB: B" << lIB->id()
                            ) ;

  SetIsProcessed(aA) ;
  SetIsProcessed(aB) ;
  GLAV_remove(aA);
  GLAV_remove(aB);

  CGAL_STSKEL_BUILDER_TRACE ( 3, 'N' << aA->id() << " processed\nN" << aB->id() << " processed" ) ;

  Halfedge_handle lOA_Prev = lOA->prev() ;
  Halfedge_handle lIA_Next = lIA->next() ;

  Halfedge_handle lOB_Prev = lOB->prev() ;
  Halfedge_handle lIB_Next = lIB->next() ;
  (void) lOB_Prev; // may be unused
  (void) lIB_Next; // may be unused

  CGAL_STSKEL_BUILDER_TRACE ( 2
                            ,   "OA_Prev: B" << lOA_Prev->id() << '\n'
                              << "IA_Next: B" << lIA_Next->id() << '\n'
                              << "OB_Prev: B" << lOB_Prev->id() << '\n'
                              << "IB_Next: B" << lIB_Next->id()
                           ) ;

  CrossLinkFwd(lOB, lIA_Next );
  CrossLinkFwd(lOA_Prev, lIB );

  Link(lOB,aA);

  CGAL_STSKEL_BUILDER_TRACE ( 1, "B" << lOA->id() << " and B" << lIA->id() << " erased." ) ;
  mDanglingBisectors.push_back(lOA);

  //
  // The code above corrects the links for vertices aA/aB to the erased halfedges lOA and lIA.
  // However, any of these vertices (aA/aB) maybe one of the twin vertices of a split event.
  // If that's the case, the erased halfedge maybe be linked to a 'couple' of those vertices.
  // This situation is corrected below:


  if ( !lOAV->has_infinite_time() && lOAV != aA && lOAV != aB )
  {
    Link(lOAV,lIB);

    CGAL_STSKEL_BUILDER_TRACE ( 1, "N" << lOAV->id() << " has B" << lOA->id()
                              << " as it's halfedge. Replacing it with B" << lIB->id()
                              ) ;
  }
  if ( !lIAV->has_infinite_time() && lIAV != aA && lIAV != aB )
  {
    Link(lIAV,lOB);

    CGAL_STSKEL_BUILDER_TRACE ( 1, "N" << lIAV->id() << " has B" << lIA->id()
                              << " as it's halfedge. Replacing it with B" << lOB->id()
                              ) ;
  }

  CGAL_STSKEL_BUILDER_TRACE ( 2, "N" << aA->id() << " halfedge: B" << aA->halfedge()->id() ) ;
  CGAL_STSKEL_BUILDER_TRACE ( 2, "N" << aB->id() << " halfedge: B" << aB->halfedge()->id() ) ;

  SetBisectorSlope(aA,aB);

  CGAL_assertion( aA->primary_bisector() == lIB ) ;

  CGAL_STSKEL_BUILDER_TRACE ( 1, "Wavefront: E" << lIB->defining_contour_edge()->id() << " and E" << lIB->opposite()->defining_contour_edge()->id() << " annihilated each other." ) ;

  if ( lOAV->has_infinite_time() )
  {
    CGAL_STSKEL_BUILDER_TRACE ( 2, "Fictitious N" << lOAV->id() << " erased." ) ;
    EraseNode(lOAV);
  }

  if ( lOBV->has_infinite_time() )
  {
    CGAL_STSKEL_BUILDER_TRACE ( 2, "Fictitious N" << lOBV->id() << " erased." ) ;
    EraseNode(lOBV);
  }
}

// Returns true if the skeleton edges 'aA' and 'aB' are defined by the same pair of contour edges (but possibly in reverse order)
//
template<class Gt, class Ss, class V>
bool Straight_skeleton_builder_2<Gt,Ss,V>::AreBisectorsCoincident ( Halfedge_const_handle aA, Halfedge_const_handle aB  ) const
{
  CGAL_STSKEL_BUILDER_TRACE ( 3, "Testing for simultaneous EdgeEvents between B" << aA->id() << " and B" << aB->id() ) ;

  Halfedge_const_handle lA_LBorder = aA->defining_contour_edge();
  Halfedge_const_handle lA_RBorder = aA->opposite()->defining_contour_edge();
  Halfedge_const_handle lB_LBorder = aB->defining_contour_edge();
  Halfedge_const_handle lB_RBorder = aB->opposite()->defining_contour_edge();

  CGAL_STSKEL_BUILDER_TRACE ( 3, "aA = " << e2str(*aA)) ;
  CGAL_STSKEL_BUILDER_TRACE ( 3, "aB = " << e2str(*aB)) ;
  CGAL_STSKEL_BUILDER_TRACE ( 3, "lA_LBorder = " << e2str(*lA_LBorder)) ;
  CGAL_STSKEL_BUILDER_TRACE ( 3, "lA_RBorder = " << e2str(*lA_RBorder)) ;
  CGAL_STSKEL_BUILDER_TRACE ( 3, "lB_LBorder = " << e2str(*lB_LBorder)) ;
  CGAL_STSKEL_BUILDER_TRACE ( 3, "lB_RBorder = " << e2str(*lB_RBorder)) ;

  return    ( lA_LBorder == lB_LBorder && lA_RBorder == lB_RBorder )
         || ( lA_LBorder == lB_RBorder && lA_RBorder == lB_LBorder ) ;
}

template<class Gt, class Ss, class V>
void Straight_skeleton_builder_2<Gt,Ss,V>::UpdatePQ( Vertex_handle aNode, Triedge const& aPrevEventTriedge )
{
  Vertex_handle lPrev = GetPrevInLAV(aNode) ;
  Vertex_handle lNext = GetNextInLAV(aNode) ;

  CGAL_STSKEL_BUILDER_TRACE ( 3, "Updating PQ for N" << aNode->id() << " Prev N" << lPrev->id() << " Next N" << lNext->id() ) ;
  CGAL_STSKEL_BUILDER_TRACE ( 3, "Respective positions " << aNode->point() << " Prev " << lPrev->point() << " Next " << lNext->point() ) ;

  Halfedge_handle lOBisector_P = lPrev->primary_bisector() ;
  Halfedge_handle lOBisector_C = aNode->primary_bisector() ;
  Halfedge_handle lOBisector_N = lNext->primary_bisector() ;

  if ( AreBisectorsCoincident(lOBisector_C,lOBisector_P) )
    HandleSimultaneousEdgeEvent( aNode, lPrev ) ;
  else if ( AreBisectorsCoincident(lOBisector_C,lOBisector_N) )
    HandleSimultaneousEdgeEvent( aNode, lNext ) ;
  else
     CollectNewEvents(aNode,aPrevEventTriedge);
}
template<class Gt, class Ss, class V>
void Straight_skeleton_builder_2<Gt,Ss,V>::CreateInitialEvents()
{
  Triedge const cNull_triedge ;

  CGAL_STSKEL_BUILDER_TRACE(0, "Creating initial events...");
  for ( Vertex_iterator v = mSSkel->vertices_begin(); v != mSSkel->vertices_end(); ++ v )
  {
    if ( ! v->has_infinite_time() )
    {
      UpdatePQ(v,cNull_triedge);
      mVisitor.on_initial_events_collected(v,IsReflex(v),IsDegenerate(v)) ;
    }
  }
}


template<class Gt, class Ss, class V>
void Straight_skeleton_builder_2<Gt,Ss,V>::CreateContourBisectors()
{
  CGAL_STSKEL_BUILDER_TRACE(0, "Creating contour bisectors...");
  for ( Vertex_iterator v = mSSkel->vertices_begin(); v != mSSkel->vertices_end(); ++ v )
  {
    GLAV_push_back(static_cast<Vertex_handle>(v));
    Vertex_handle lPrev = GetPrevInLAV(v) ;
    Vertex_handle lNext = GetNextInLAV(v) ;

    Orientation lOrientation = CGAL::orientation(lPrev->point(),v->point(),lNext->point());
    if ( lOrientation == COLLINEAR )
    {
      SetIsDegenerate(v);
      CGAL_STSKEL_BUILDER_TRACE(1, "COLLINEAR vertex: N" << v->id() );
    }
    else if ( lOrientation == RIGHT_TURN )
    {
      mReflexVertices.push_back(v);
      SetIsReflex(v);
      CGAL_STSKEL_BUILDER_TRACE(1,"Reflex vertex: N" << v->id() );
    }

    Halfedge lOB(mEdgeID++), lIB(mEdgeID++);
    Halfedge_handle lOBisector = SSkelEdgesPushBack(lOB, lIB);
    Halfedge_handle lIBisector = lOBisector->opposite();
    lOBisector->HBase_base::set_face(v->halfedge()->face());
    lIBisector->HBase_base::set_face(v->halfedge()->next()->face());
    lIBisector->HBase_base::set_vertex(v);

    Halfedge_handle lIBorder = v->halfedge() ;
    Halfedge_handle lOBorder = v->halfedge()->next() ;
    lIBorder  ->HBase_base::set_next(lOBisector);
    lOBisector->HBase_base::set_prev(lIBorder);
    lOBorder  ->HBase_base::set_prev(lIBisector);
    lIBisector->HBase_base::set_next(lOBorder);
    CGAL_STSKEL_BUILDER_TRACE(3
                             ,"Adding Contour Bisector at " << v2str(*v)
                              << "\n B" << lOBisector->id()
                              << " (Out)\n B" << lIBisector->id() << " (In)"
                             ) ;
  }

  for( Face_iterator fit = mSSkel->SSkel::Base::faces_begin(); fit != mSSkel->SSkel::Base::faces_end(); ++fit)
  {
    Halfedge_handle lBorder    = fit->halfedge();
    Halfedge_handle lLBisector = lBorder->prev();
    Halfedge_handle lRBisector = lBorder->next();

    Vertex_handle lInfNode = mSSkel->SSkel::Base::vertices_push_back( Vertex( mVertexID++ ) ) ;
    InitVertexData(lInfNode);
    CGAL_assertion(lInfNode->has_null_point());

    lRBisector->HBase_base::set_next( lLBisector  );
    lLBisector->HBase_base::set_prev( lRBisector );

    lRBisector->HBase_base::set_vertex(lInfNode);

    lInfNode->VBase::set_halfedge(lRBisector);

    SetBisectorSlope(lRBisector,POSITIVE);
    SetBisectorSlope(lLBisector,NEGATIVE);

    CGAL_STSKEL_BUILDER_TRACE(3
                             ,"Closing face of " << e2str(*lBorder)
                             << " with a fictitious vertex. B" << lRBisector->id()
                             << "->N" << lInfNode->id()
                             << "->B" << lLBisector->id()
                             ) ;
  }
}

template<class Gt, class Ss, class V>
void Straight_skeleton_builder_2<Gt,Ss,V>::HarmonizeSpeeds(boost::mpl::bool_<true>)
{
  // Collinear input edges might not have the exact same speed if an inexact square root is used.
  // This might cause some inconsistencies in time, resulting in invalid skeletons. Therefore,
  // if the square root is not exact, we enforce that collinear input edges have the same speed,
  // by making them use the same line coefficients (which determines the speed of the front).
  //
  // That is achieved by creating a set of input edges, with two input edges being equal if they are collinear.
  // If a new input edge is not successfully inserted into the same, it takes the line coefficients
  // of the representative of this class.

  auto comparer = [&](Halfedge_handle lLH, Halfedge_handle lRH) -> bool
  {
    const Direction_2 lLD = CreateDirection(lLH) ;
    const Direction_2 lRD = CreateDirection(lRH) ;
    Comparison_result rRes = K().compare_angle_with_x_axis_2_object()(lLD, lRD) ;

    if ( rRes == EQUAL ) // parallel
    {
      if ( K().orientation_2_object()(lLH->vertex()->point(),
                                      lLH->opposite()->vertex()->point(),
                                      lRH->vertex()->point()) == EQUAL )
        return false; // collinear

      // parallel but not collinear, order arbitrarily (but consistently)
      return K().less_xy_2_object()(lLH->vertex()->point(), lRH->vertex()->point()) ;
    }
    else
    {
      // not parallel
      return ( rRes == SMALLER ) ;
    }
  } ;

  typedef std::set<Halfedge_handle, decltype(comparer)> Ordered_halfedges;
  Ordered_halfedges lOrdered_halfedges(comparer);

  typename CGAL_SS_i::Get_protector<Gt>::type protector;
  CGAL_USE(protector);

  for( Face_iterator fit = mSSkel->SSkel::Base::faces_begin(); fit != mSSkel->SSkel::Base::faces_end(); ++fit)
  {
    Halfedge_handle lBorder = fit->halfedge() ;
    Segment_2 lS = CreateSegment<Traits> ( lBorder ) ;

    std::pair<typename Ordered_halfedges::iterator, bool> rRes = lOrdered_halfedges.insert ( lBorder ) ;
    if ( ! rRes.second ) // some collinear edge is already in the set
      mTraits.InitializeLineCoeffs ( lBorder->id(), (*rRes.first)->id() );
    else
      mTraits.InitializeLineCoeffs ( lS );
  }
}

template<class Gt, class Ss, class V>
void Straight_skeleton_builder_2<Gt,Ss,V>::InitPhase()
{
  mVisitor.on_initialization_started(static_cast<int>(mSSkel->size_of_vertices()));
  CreateContourBisectors();
  HarmonizeSpeeds();
  CreateInitialEvents();
  mVisitor.on_initialization_finished();
}

template<class Gt, class Ss, class V>
typename Straight_skeleton_builder_2<Gt,Ss,V>::Vertex_handle
Straight_skeleton_builder_2<Gt,Ss,V>::ConstructEdgeEventNode( EdgeEvent& aEvent )
{
  CGAL_STSKEL_BUILDER_TRACE ( 2, "Creating EdgeEvent Node" ) ;

  Vertex_handle lLSeed = aEvent.seed0() ;
  Vertex_handle lRSeed = aEvent.seed1() ;

  Vertex_handle lNewNode = mSSkel->SSkel::Base::vertices_push_back( Vertex( mVertexID++, aEvent.point(), aEvent.time(), false, false) ) ;
  InitVertexData(lNewNode);

  GLAV_push_back(lNewNode);

  SetTrisegment(lNewNode,aEvent.trisegment());

  CGAL_STSKEL_BUILDER_TRACE
  ( 3
  ,    "LSeed: N" << lLSeed->id() << " processed\n"
    << "RSeed: N" << lRSeed->id() << " processed"
  ) ;

  SetIsProcessed(lLSeed) ;
  SetIsProcessed(lRSeed) ;
  GLAV_remove(lLSeed);
  GLAV_remove(lRSeed);

  Vertex_handle lLPrev = GetPrevInLAV(lLSeed) ;
  Vertex_handle lRNext = GetNextInLAV(lRSeed) ;

  SetPrevInLAV(lNewNode, lLPrev ) ;
  SetNextInLAV(lLPrev  , lNewNode  ) ;

  SetNextInLAV(lNewNode, lRNext ) ;
  SetPrevInLAV(lRNext  , lNewNode  ) ;

  CGAL_STSKEL_BUILDER_TRACE( 2, "New Node: N" << lNewNode->id() << " at " << lNewNode->point() << '\n'
                              << 'N' << lLSeed->id() << " removed from LAV\n"
                              << 'N' << lRSeed->id() << " removed from LAV\n"
                           );

  return lNewNode ;
}


template<class Gt, class Ss, class V>
typename Straight_skeleton_builder_2<Gt,Ss,V>::Vertex_handle_pair
Straight_skeleton_builder_2<Gt,Ss,V>::LookupOnSLAV ( Halfedge_handle aBorder, EventPtr const& aEvent, Site& rSite )
{
  Vertex_handle_pair rResult ;

  CGAL_STSKEL_DEBUG_CODE( bool lFound = false ; )

  // Vertex_handle lSeed = aEvent->seed0();

  CGAL_STSKEL_BUILDER_TRACE ( 3, "Looking up for E" << aBorder->id() << ". P=" << aEvent->point() ) ;

  for (Vertex_handle v : GetHalfedgeLAVList(aBorder))
  {
    Triedge const& lTriedge = GetVertexTriedge(v);

    Vertex_handle lPrevN = GetPrevInLAV(v);
    Vertex_handle lNextN = GetNextInLAV(v);

    if ( lTriedge.e0() == aBorder )
    {
      Halfedge_handle lPrevBorder = GetEdgeEndingAt(lPrevN) ;
      Halfedge_handle lNextBorder = GetEdgeEndingAt(lNextN) ;

      CGAL_STSKEL_DEBUG_CODE( lFound = true ; )

      CGAL_STSKEL_BUILDER_TRACE ( 3
                                , "Subedge found in SLAV: N" << lPrevN->id() << "->N" << v->id()
                                  << " (E" << lPrevBorder->id() << "->E" << aBorder->id() << "->E" << lNextBorder->id() << ")"
                                ) ;

      Oriented_side lLSide = EventPointOrientedSide(*aEvent, lPrevBorder, aBorder    , lPrevN, false ) ;
      Oriented_side lRSide = EventPointOrientedSide(*aEvent, aBorder    , lNextBorder, v     , true  ) ;

      if ( lLSide != ON_POSITIVE_SIDE && lRSide != ON_NEGATIVE_SIDE )
      {
        if ( lLSide != ON_ORIENTED_BOUNDARY || lRSide != ON_ORIENTED_BOUNDARY )
        {
          rSite = ( lLSide == ON_ORIENTED_BOUNDARY ? AT_SOURCE : ( lRSide == ON_ORIENTED_BOUNDARY ?  AT_TARGET : INSIDE ) ) ;

          rResult = std::make_pair(lPrevN,v) ;

          CGAL_STSKEL_BUILDER_TRACE ( 3, "Split point found at the "
                                    << ( rSite == AT_SOURCE ? "SOURCE vertex" : ( rSite == AT_TARGET ? "TARGET vertex" : "strict inside" ) )
                                    << " of the offset edge."
                                    ) ;
          break ;
        }
        else
        {
          CGAL_STSKEL_BUILDER_TRACE ( 3, "Opposite edge collapsed to a point" ) ;
        }
      }
    }
  }

#ifdef CGAL_STRAIGHT_SKELETON_ENABLE_TRACE
  if ( !handle_assigned(rResult.first) )
  {
    if ( !lFound )
    {
      CGAL_STSKEL_BUILDER_TRACE(1,"Split event is no longer valid. Opposite edge vanished.");
    }
    else
    {
      CGAL_STSKEL_BUILDER_TRACE(1,"Split event is no longer valid. Point not inside the opposite edge offset zone.");
    }
  }
#endif

  return rResult ;
}



template<class Gt, class Ss, class V>
typename Straight_skeleton_builder_2<Gt,Ss,V>::Vertex_handle_pair
Straight_skeleton_builder_2<Gt,Ss,V>::ConstructSplitEventNodes( SplitEvent& aEvent, Vertex_handle aOppR )
{
  Vertex_handle_pair rResult;

  CGAL_STSKEL_BUILDER_TRACE ( 2, "Creating SplitEvent Nodes" ) ;

  Vertex_handle lOppL = GetPrevInLAV(aOppR) ;

  Vertex_handle lNewNodeA = mSSkel->SSkel::Base::vertices_push_back( Vertex( mVertexID++, aEvent.point(), aEvent.time(), true, false ) ) ;
  Vertex_handle lNewNodeB = mSSkel->SSkel::Base::vertices_push_back( Vertex( mVertexID++, aEvent.point(), aEvent.time(), true, false ) ) ;

  InitVertexData(lNewNodeA);
  InitVertexData(lNewNodeB);
  SetTrisegment(lNewNodeA,aEvent.trisegment());
  SetTrisegment(lNewNodeB,aEvent.trisegment());

  GLAV_push_back(lNewNodeA);
  GLAV_push_back(lNewNodeB);

  Vertex_handle lSeed = aEvent.seed0() ;

  CGAL_STSKEL_BUILDER_TRACE ( 3, "Seed: N" << lSeed->id() << " processed" ) ;

  SetIsProcessed(lSeed) ;
  GLAV_remove(lSeed);

  CGAL_STSKEL_BUILDER_TRACE ( 2, 'N' << lNewNodeA->id() << " and N" << lNewNodeB->id() << " inserted into LAV." ) ;

  Vertex_handle lPrev = GetPrevInLAV(lSeed) ;
  Vertex_handle lNext = GetNextInLAV(lSeed) ;

  SetNextInLAV(lPrev    , lNewNodeA ) ;
  SetPrevInLAV(lNewNodeA, lPrev     ) ;

  SetNextInLAV(lNewNodeA, aOppR     ) ;
  SetPrevInLAV(aOppR    , lNewNodeA ) ;

  SetNextInLAV(lOppL    , lNewNodeB ) ;
  SetPrevInLAV(lNewNodeB, lOppL     ) ;

  SetNextInLAV(lNewNodeB, lNext     ) ;
  SetPrevInLAV(lNext    , lNewNodeB ) ;

  CGAL_STSKEL_BUILDER_TRACE( 2, 'N' << lSeed->id() << " removed from LAV" );

  rResult = std::make_pair(lNewNodeA,lNewNodeB);

  mSplitNodes.push_back(rResult);

  return rResult ;
}

template<class Gt, class Ss, class V>
typename Straight_skeleton_builder_2<Gt,Ss,V>::Vertex_handle_pair
Straight_skeleton_builder_2<Gt,Ss,V>::ConstructPseudoSplitEventNodes( PseudoSplitEvent& aEvent )
{
  Vertex_handle_pair rResult;

  CGAL_STSKEL_BUILDER_TRACE ( 2, "Creating PseudoSplitEvent Nodes" ) ;

  Vertex_handle lLSeed = aEvent.seed0() ;
  Vertex_handle lRSeed = aEvent.seed1() ;

  Vertex_handle lNewNodeA = mSSkel->SSkel::Base::vertices_push_back( Vertex( mVertexID++, aEvent.point(), aEvent.time(), true, false ) ) ;
  Vertex_handle lNewNodeB = mSSkel->SSkel::Base::vertices_push_back( Vertex( mVertexID++, aEvent.point(), aEvent.time(), true, false ) ) ;

  GLAV_push_back(lNewNodeA);
  GLAV_push_back(lNewNodeB);

  InitVertexData(lNewNodeA);
  InitVertexData(lNewNodeB);
  SetTrisegment(lNewNodeA,aEvent.trisegment());
  SetTrisegment(lNewNodeB,aEvent.trisegment());

  CGAL_STSKEL_BUILDER_TRACE
  (
   3
   ,   "LSeed: N" << lLSeed->id() << " processed\n"
    << "RSeed: N" << lRSeed->id() << " processed"
  ) ;


  SetIsProcessed(lLSeed) ;
  SetIsProcessed(lRSeed) ;
  GLAV_remove(lLSeed);
  GLAV_remove(lRSeed);

  Vertex_handle lLPrev = GetPrevInLAV(lLSeed) ;
  Vertex_handle lLNext = GetNextInLAV(lLSeed) ;
  Vertex_handle lRPrev = GetPrevInLAV(lRSeed) ;
  Vertex_handle lRNext = GetNextInLAV(lRSeed) ;

  SetPrevInLAV(lNewNodeA, lLPrev    ) ;
  SetNextInLAV(lLPrev   , lNewNodeA ) ;

  SetNextInLAV(lNewNodeA, lRNext    ) ;
  SetPrevInLAV(lRNext   , lNewNodeA ) ;

  SetPrevInLAV(lNewNodeB, lRPrev    ) ;
  SetNextInLAV(lRPrev   , lNewNodeB ) ;

  SetNextInLAV(lNewNodeB, lLNext    ) ;
  SetPrevInLAV(lLNext   , lNewNodeB ) ;

  CGAL_STSKEL_BUILDER_TRACE(2,   "NewNodeA: N" << lNewNodeA->id() << " at " << lNewNodeA->point() << '\n'
                              << "NewNodeB: N" << lNewNodeB->id() << " at " << lNewNodeB->point() << '\n'
                              << 'N' << lLSeed->id() << " removed from LAV\n"
                              << 'N' << lRSeed->id() << " removed from LAV"
                           );


  rResult = std::make_pair(lNewNodeA,lNewNodeB);

  mSplitNodes.push_back(rResult);

  return rResult ;
}

template<class Gt, class Ss, class V>
bool Straight_skeleton_builder_2<Gt,Ss,V>::IsProcessed( EventPtr aEvent )
{
  return IsProcessed(aEvent->seed0()) || IsProcessed(aEvent->seed1()) ;
}

template<class Gt, class Ss, class V>
bool Straight_skeleton_builder_2<Gt,Ss,V>::IsValidEvent( EventPtr aEvent )
{
  if ( IsProcessed(aEvent) )
    return false;

  SetEventTimeAndPoint(*aEvent) ;

  if ( aEvent->type() == Event::cEdgeEvent)
  {
    EdgeEvent& lEvent = dynamic_cast<EdgeEvent&>(*aEvent) ;
    return IsValidEdgeEvent(lEvent) ;
  }
  else if ( aEvent->type() == Event::cSplitEvent)
  {
    Halfedge_handle lOppEdge = aEvent->triedge().e2() ;
    Site lSite;
    Vertex_handle_pair lOpp = LookupOnSLAV(lOppEdge,aEvent,lSite);

    if ( handle_assigned(lOpp.first) )
    {
      EventPtr lPseudoSplitEvent = IsPseudoSplitEvent(aEvent,lOpp,lSite);
      if ( lPseudoSplitEvent )
      {
        PseudoSplitEvent& lEvent = dynamic_cast<PseudoSplitEvent&>(*lPseudoSplitEvent) ;
        return IsValidPseudoSplitEvent ( lEvent );
      }
      else
      {
        SplitEvent& lEvent = dynamic_cast<SplitEvent&>(*aEvent) ;
        return IsValidSplitEvent(lEvent);
      }
    }
    else
    {
      return false;
    }
  }
  else
  {
    CGAL_assertion( aEvent->type() == Event::cPseudoSplitEvent ) ;
    PseudoSplitEvent& lEvent = dynamic_cast<PseudoSplitEvent&>(*aEvent) ;
    return IsValidPseudoSplitEvent(lEvent);
  }
}

template<class Gt, class Ss, class V>
bool Straight_skeleton_builder_2<Gt,Ss,V>::IsValidEdgeEvent( EdgeEvent const& aEvent )
{
  bool rResult = false ;

  Vertex_handle lLSeed = aEvent.seed0() ;
  Vertex_handle lRSeed = aEvent.seed1() ;

  Vertex_handle lPrevLSeed = GetPrevInLAV(lLSeed) ;
  Vertex_handle lNextRSeed = GetNextInLAV(lRSeed) ;

  if ( lPrevLSeed != lNextRSeed )
  {
    Halfedge_handle lPrevE0 = GetEdgeEndingAt(lPrevLSeed) ;
    Halfedge_handle lE0     = aEvent.triedge().e0() ;
    Halfedge_handle lE2     = aEvent.triedge().e2() ;
    Halfedge_handle lNextE2 = GetEdgeStartingAt(lNextRSeed) ;

    CGAL_STSKEL_BUILDER_TRACE(3, "PrevLSeed=N" << lPrevLSeed->id() << " PrevE0=E" << lPrevE0->id() ) ;
    CGAL_STSKEL_BUILDER_TRACE(3, "NextRSeed=N" << lNextRSeed->id() << " NextE2=E" << lNextE2->id() ) ;

    Oriented_side lLSide = EventPointOrientedSide(aEvent, lPrevE0, lE0    , lPrevLSeed, false ) ;
    Oriented_side lRSide = EventPointOrientedSide(aEvent, lE2    , lNextE2, lNextRSeed, true  ) ;

    bool lLSideOK = ( lLSide != ON_POSITIVE_SIDE ) ;
    bool lRSideOK = ( lRSide != ON_NEGATIVE_SIDE ) ;

    CGAL_STSKEL_BUILDER_TRACE_IF( !lLSideOK
                                ,3
                                ,"Invalid edge event: " << aEvent.triedge() << " NewNode is before E" << lE0->id()
                                << " source N" << lPrevLSeed->id()
                                ) ;

    CGAL_STSKEL_BUILDER_TRACE_IF( !lRSideOK
                                ,3
                                ,"Invalid edge event: " << aEvent.triedge() << " NewNode is past E" << lE2->id()
                                 << " target N" << lNextRSeed->id()
                                ) ;

    rResult = lLSideOK && lRSideOK ;
  }
  else
  {
    // Triangle collapse. No need to test explicitely.
    rResult = true ;
  }
  return rResult ;
}

template<class Gt, class Ss, class V>
void Straight_skeleton_builder_2<Gt,Ss,V>::HandleEdgeEvent( EventPtr aEvent )
{
  EdgeEvent& lEvent = dynamic_cast<EdgeEvent&>(*aEvent) ;

  CGAL_STSKEL_BUILDER_TRACE( 2, "Edge event." );

  if ( IsValidEdgeEvent(lEvent) )
  {
    Vertex_handle lLSeed = lEvent.seed0() ;
    Vertex_handle lRSeed = lEvent.seed1() ;

    CGAL_STSKEL_BUILDER_TRACE( 3, "valid event." );

    CGAL_STSKEL_BUILDER_TRACE( 4, "TriEdge e0 = " << e2str(*lEvent.triedge().e0()) );
    CGAL_STSKEL_BUILDER_TRACE( 4, "TriEdge e1 = " << e2str(*lEvent.triedge().e1()) );
    CGAL_STSKEL_BUILDER_TRACE( 4, "TriEdge e2 = " << e2str(*lEvent.triedge().e2()) );

    Vertex_handle lNewNode = ConstructEdgeEventNode(lEvent);

    Halfedge_handle lLOBisector = lLSeed->primary_bisector() ;
    Halfedge_handle lROBisector = lRSeed->primary_bisector() ;
    Halfedge_handle lLIBisector = lLOBisector->opposite();
    Halfedge_handle lRIBisector = lROBisector->opposite();

    Vertex_handle lRIFicNode = lROBisector->vertex() ;
    Vertex_handle lLOFicNode = lLOBisector->vertex() ;

    CrossLink(lLOBisector,lNewNode);

    Link(lROBisector,lNewNode);

    CrossLinkFwd(lROBisector,lLIBisector) ;

    Halfedge_handle lDefiningBorderA = lNewNode->halfedge()->defining_contour_edge();
    Halfedge_handle lDefiningBorderB = lNewNode->halfedge()->opposite()->prev()->opposite()->defining_contour_edge();
    Halfedge_handle lDefiningBorderC = lNewNode->halfedge()->opposite()->prev()->defining_contour_edge();

    lNewNode->VBase::set_event_triedge( lEvent.triedge() ) ;

    Triedge lTri(lDefiningBorderA,lDefiningBorderB,lDefiningBorderC);

    SetVertexTriedge( lNewNode, lTri ) ;

    SetBisectorSlope(lLSeed,lNewNode);
    SetBisectorSlope(lRSeed,lNewNode);

    CGAL_STSKEL_BUILDER_TRACE( 1, "E" << e2str(*(lRSeed->halfedge()->defining_contour_edge())) << " collapsed." );
    CGAL_STSKEL_BUILDER_TRACE( 3, "fictitious node along collapsed face is N" << lRIFicNode->id()
                               << " between " << e2str(*lROBisector) << " and " << e2str(*lLIBisector) ) ;

    if ( lLOFicNode->has_infinite_time() )
    {
      CGAL_STSKEL_BUILDER_TRACE(3,"Creating new Edge Event's Bisector");

      Halfedge_handle lNOBisector = SSkelEdgesPushBack( Halfedge(mEdgeID),Halfedge(mEdgeID+1) );

      Halfedge_handle lNIBisector = lNOBisector->opposite();
      mEdgeID += 2 ;

      CrossLinkFwd(lNOBisector        , lLOBisector->next());
      CrossLinkFwd(lRIBisector->prev(), lNIBisector        );

      CrossLink(lNOBisector,lLOFicNode);

      SetBisectorSlope(lNOBisector,POSITIVE);
      SetBisectorSlope(lNIBisector,NEGATIVE);

      CrossLinkFwd(lNIBisector, lRIBisector);
      CrossLinkFwd(lLOBisector, lNOBisector);

      Link(lNOBisector,lLOBisector->face());
      Link(lNIBisector,lRIBisector->face());

      Link(lNIBisector,lNewNode);

      CGAL_STSKEL_BUILDER_TRACE( 2, newn2str("",lNewNode,lTri) ) ;
      CGAL_STSKEL_BUILDER_TRACE( 2, newb2str("O",lNOBisector) ) ;
      CGAL_STSKEL_BUILDER_TRACE( 2, newb2str("I",lNIBisector) ) ;

      CGAL_STSKEL_BUILDER_TRACE ( 2, "Fictitious N" << lRIFicNode->id() << " erased." ) ;
      EraseNode(lRIFicNode);

      SetupNewNode(lNewNode) ;

      UpdatePQ(lNewNode, lEvent.triedge());

      mVisitor.on_edge_event_processed(lLSeed,lRSeed,lNewNode) ;
    }
    else
    {
      CGAL_STSKEL_BUILDER_TRACE( 2, newn2str("",lNewNode,lTri)
                                    << ".\nThis is a multiple node (A node with these defining edges already exists in the LAV)"
                               );
    }

    CGAL_STSKEL_BUILDER_TRACE( 1, "Wavefront: " << wavefront2str(lNewNode) );
  }
}

template<class Gt, class Ss, class V>
bool Straight_skeleton_builder_2<Gt,Ss,V>::IsValidSplitEvent( SplitEvent const& /*aEvent*/ )
{
  return true ;
}

template<class Gt, class Ss, class V>
void Straight_skeleton_builder_2<Gt,Ss,V>::HandleSplitEvent( EventPtr aEvent, Vertex_handle_pair aOpp  )
{
  SplitEvent& lEvent = dynamic_cast<SplitEvent&>(*aEvent) ;

  CGAL_STSKEL_BUILDER_TRACE( 2, "Split event." );

  if ( IsValidSplitEvent(lEvent) )
  {
    CGAL_STSKEL_BUILDER_TRACE( 3, "valid event." );

    Vertex_handle lSeed = lEvent.seed0();

    Vertex_handle lOppL = aOpp.first ;
    Vertex_handle lOppR = aOpp.second ;

    CGAL_STSKEL_BUILDER_TRACE( 4, "TriEdge e0 = " << e2str(*lEvent.triedge().e0()) );
    CGAL_STSKEL_BUILDER_TRACE( 4, "TriEdge e1 = " << e2str(*lEvent.triedge().e1()) );
    CGAL_STSKEL_BUILDER_TRACE( 4, "TriEdge e2 = " << e2str(*lEvent.triedge().e2()) );
    CGAL_STSKEL_BUILDER_TRACE( 4, "lOppL = " << v2str(*lOppL) );
    CGAL_STSKEL_BUILDER_TRACE( 4, "lOppR = " << v2str(*lOppR) );

    Halfedge_handle lOppIBisector_L = lOppL->primary_bisector()->opposite();
    Halfedge_handle lOppOBisector_R = lOppR->primary_bisector();

    Vertex_handle lOppFicNode = lOppOBisector_R->vertex() ;
    (void)lOppFicNode; // variable may be unused

    CGAL_assertion(lOppOBisector_R->next() == lOppIBisector_L ) ;
    CGAL_assertion(lOppIBisector_L->prev() == lOppOBisector_R ) ;
    CGAL_assertion(lOppFicNode->has_infinite_time());

    CGAL_STSKEL_BUILDER_TRACE(2,"Splitted face: N" << lOppR->id()
                                << "->B" << lOppOBisector_R->id()
                                << "->N" << lOppFicNode->id()
                                << "->B" << lOppIBisector_L->id()
                                << "->N" << lOppL->id()
                             ) ;

    CGAL_STSKEL_BUILDER_TRACE(2,"fictitious node for right half of opposite edge: N" << lOppFicNode->id() ) ;

    Halfedge_handle lOppBorder = lEvent.triedge().e2() ;

    Vertex_handle lNewNode_L, lNewNode_R ;
    boost::tie(lNewNode_L,lNewNode_R) = ConstructSplitEventNodes(lEvent,lOppR);

    // Triedge lTriedge = aEvent->triedge();

    // Halfedge_handle lReflexLBorder = lTriedge.e0();
    // Halfedge_handle lReflexRBorder = lTriedge.e1();

    Halfedge_handle lNOBisector_L = SSkelEdgesPushBack( Halfedge(mEdgeID  ),Halfedge(mEdgeID+1) );
    Halfedge_handle lNOBisector_R = SSkelEdgesPushBack( Halfedge(mEdgeID+2),Halfedge(mEdgeID+3) );
    Halfedge_handle lNIBisector_L = lNOBisector_L->opposite();
    Halfedge_handle lNIBisector_R = lNOBisector_R->opposite();

    mEdgeID += 4 ;

    Halfedge_handle lXOBisector = lSeed->primary_bisector() ;
    Halfedge_handle lXIBisector = lXOBisector->opposite();

    Halfedge_handle lXONextBisector = lXOBisector->next();
    Halfedge_handle lXIPrevBisector = lXIBisector->prev();

    Vertex_handle lXOFicNode = lXOBisector->vertex() ;
    CGAL_assertion(lXOFicNode->has_infinite_time());

    CGAL_STSKEL_BUILDER_TRACE(2,"fictitious node for left reflex face: N" << lXOFicNode->id() ) ;
    CGAL_STSKEL_BUILDER_TRACE(2,"fictitious node for right reflex face: N" << lXIPrevBisector->vertex()->id() ) ;

    Link(lNewNode_L,lXOBisector);
    Link(lNewNode_R,lNIBisector_L) ;

    Link(lXOBisector,lNewNode_L);

    Link(lNOBisector_L,lXOBisector->face());
    Link(lNIBisector_L,lOppBorder ->face());
    Link(lNOBisector_R,lOppBorder ->face());
    Link(lNIBisector_R,lXIBisector->face());

    Link(lNIBisector_L,lNewNode_R);
    Link(lNIBisector_R,lNewNode_R);

    Link(lNOBisector_L,lXOFicNode);


    CrossLinkFwd(lXOBisector    ,lNOBisector_L);
    CrossLinkFwd(lNOBisector_L  ,lXONextBisector);
    CrossLinkFwd(lXIPrevBisector,lNIBisector_R);
    CrossLinkFwd(lNIBisector_R  ,lXIBisector);
    CrossLinkFwd(lOppOBisector_R,lNIBisector_L);
    CrossLinkFwd(lNIBisector_L  ,lNOBisector_R);
    CrossLinkFwd(lNOBisector_R  ,lOppIBisector_L);

    SetBisectorSlope(lSeed,lNewNode_L);

    Vertex_handle lNewFicNode = mSSkel->SSkel::Base::vertices_push_back( Vertex( mVertexID++ ) ) ;

    InitVertexData(lNewFicNode);
    CGAL_assertion(lNewFicNode->has_null_point());
    CrossLink(lNOBisector_R,lNewFicNode);

    SetBisectorSlope(lNOBisector_L,POSITIVE);
    SetBisectorSlope(lNIBisector_L,NEGATIVE);
    SetBisectorSlope(lNOBisector_R,POSITIVE);
    SetBisectorSlope(lNIBisector_R,NEGATIVE);

    CGAL_STSKEL_BUILDER_TRACE(2,"(New) fictitious node for left half of opposite edge: N" << lNewFicNode->id() ) ;

    Halfedge_handle lNewNode_L_DefiningBorderA = lNewNode_L->halfedge()->defining_contour_edge();
    Halfedge_handle lNewNode_L_DefiningBorderB = lNewNode_L->halfedge()->opposite()->prev()->opposite()->defining_contour_edge();
    Halfedge_handle lNewNode_L_DefiningBorderC = lNewNode_L->halfedge()->opposite()->prev()->defining_contour_edge();
    Halfedge_handle lNewNode_R_DefiningBorderA = lNewNode_R->halfedge()->defining_contour_edge();
    Halfedge_handle lNewNode_R_DefiningBorderB = lNewNode_R->halfedge()->opposite()->prev()->opposite()->defining_contour_edge();
    Halfedge_handle lNewNode_R_DefiningBorderC = lNewNode_R->halfedge()->opposite()->prev()->defining_contour_edge();

    lNewNode_L->VBase::set_event_triedge( lEvent.triedge() ) ;
    lNewNode_R->VBase::set_event_triedge( lEvent.triedge() ) ;

    Triedge lTriL( lNewNode_L_DefiningBorderA,lNewNode_L_DefiningBorderB,lNewNode_L_DefiningBorderC ) ;
    Triedge lTriR( lNewNode_R_DefiningBorderA,lNewNode_R_DefiningBorderB,lNewNode_R_DefiningBorderC ) ;

    SetVertexTriedge( lNewNode_L, lTriL ) ;
    SetVertexTriedge( lNewNode_R, lTriR ) ;

    CGAL_STSKEL_BUILDER_TRACE(2,   newn2str("L",lNewNode_L,lTriL) << std::endl
                                << newn2str("R",lNewNode_R,lTriR) << std::endl
                                << newb2str("OL",lNOBisector_L)   << std::endl
                                << newb2str("IL",lNIBisector_L)   << std::endl
                                << newb2str("OR",lNOBisector_R)   << std::endl
                                << newb2str("IR",lNIBisector_R)
                             ) ;

    CGAL_STSKEL_BUILDER_TRACE( 1, "Wavefronts:\n  " << wavefront2str(lNewNode_L) << "\n  " << wavefront2str(lNewNode_R) );

    SetupNewNode(lNewNode_L) ;
    SetupNewNode(lNewNode_R) ;

    UpdatePQ(lNewNode_L, lEvent.triedge());
    UpdatePQ(lNewNode_R, lEvent.triedge());

    mVisitor.on_split_event_processed(lSeed,lNewNode_L,lNewNode_R) ;
  }


}

template<class Gt, class Ss, class V>
void Straight_skeleton_builder_2<Gt,Ss,V>::SetupNewNode( Vertex_handle aNode )
{
  // In an edge-edge annihilation the current polygon becomes a two-node degenerate chain collapsed into a single point
  if ( GetPrevInLAV(aNode) != GetNextInLAV(aNode) )
  {
    Halfedge_handle lLE = GetEdgeEndingAt  (aNode);
    Halfedge_handle lRE = GetEdgeStartingAt(aNode);

    Vector_2 lLV = CreateVector(lLE);
    Vector_2 lRV = CreateVector(lRE);

    Orientation lOrientation = CGAL::orientation(lLV,lRV) ;
    if ( lOrientation == COLLINEAR )
    {
      SetIsDegenerate(aNode);
      CGAL_STSKEL_BUILDER_TRACE(1, "COLLINEAR *NEW* vertex: " << v2str(*aNode) << " (E" << lLE->id() << ",E" << lRE->id() << ")" ) ;
    }
    else if ( lOrientation == RIGHT_TURN )
    {
      mReflexVertices.push_back(aNode);
      SetIsReflex(aNode);
      CGAL_STSKEL_BUILDER_TRACE(1, "Reflex *NEW* vertex: N" << v2str(*aNode)  << " (E" << lLE->id() << ",E" << lRE->id() << ")" );
    }
  }
}

template<class Direction>
bool counterclockwise_at_or_in_between_2( Direction const& p, Direction const& q, Direction const& r )
{
  typedef typename Kernel_traits<Direction>::Kernel K ;

  return p == q || p == r || K().counterclockwise_in_between_2_object()(p,q,r) ;
}

template<class Gt, class Ss, class V>
bool Straight_skeleton_builder_2<Gt,Ss,V>::IsValidPseudoSplitEvent( PseudoSplitEvent const& aEvent )
{
  Vertex_handle lSeed0 = aEvent.seed0();
  Vertex_handle lSeed1 = aEvent.seed1();

  CGAL_STSKEL_BUILDER_TRACE(3, "Checking for tangleness..." );
  CGAL_STSKEL_BUILDER_TRACE(3, "lSeed0 = " << v2str(*lSeed0) );
  CGAL_STSKEL_BUILDER_TRACE(3, "lSeed1 = " << v2str(*lSeed1) );

  Halfedge_handle lEL0 = GetEdgeEndingAt  (lSeed0);
  Halfedge_handle lER0 = GetEdgeStartingAt(lSeed0);

  Halfedge_handle lEL1 = GetEdgeEndingAt  (lSeed1);
  Halfedge_handle lER1 = GetEdgeStartingAt(lSeed1);

  CGAL_STSKEL_BUILDER_TRACE(3, "lEL0 = " << e2str(*lEL0) );
  CGAL_STSKEL_BUILDER_TRACE(3, "lER0 = " << e2str(*lER0) );
  CGAL_STSKEL_BUILDER_TRACE(3, "lEL1 = " << e2str(*lEL1) );
  CGAL_STSKEL_BUILDER_TRACE(3, "lER1 = " << e2str(*lER1) );

  Direction_2 lDL0 = - CreateDirection(lEL0);
  Direction_2 lDL1 = - CreateDirection(lEL1);
  Direction_2 lDR0 =   CreateDirection(lER0);
  Direction_2 lDR1 =   CreateDirection(lER1);

  bool lV01Degenerate = (lDL0 == lDR1) ;
  bool lV10Degenerate = (lDL1 == lDR0) ;

  CGAL_STSKEL_BUILDER_TRACE(3, "Validating pseudo-split event. Resulting re-connection: "
                           << "\nE" << lEL0->id() << " [DL0:" << dir2str(lDL0) << "]"
                           << "->E" << lER1->id() << " [DR1:" << dir2str(lDR1) << "]" << ( lV01Degenerate ? " (degenerate)" : "" )
                           << "\nE" << lEL1->id() << " [DL1:" << dir2str(lDL1) << "]"
                           << "->E" << lER0->id() << " [DR0:" << dir2str(lDR0) << "]" << ( lV10Degenerate ? " (degenerate)" : "" )
                           ) ;

  bool lTangled ;

  if ( !lV01Degenerate )
  {
    bool lEL1V_Tangled = counterclockwise_at_or_in_between_2(lDL1,lDR1,lDL0);
    bool lER0V_Tangled = counterclockwise_at_or_in_between_2(lDR0,lDR1,lDL0);

    CGAL_STSKEL_BUILDER_TRACE(3, "lV01Degenerate not degenerate, CCW DL1,DR1,DL0 = " << lEL1V_Tangled );
    CGAL_STSKEL_BUILDER_TRACE(3, "lV01Degenerate not degenerate, CCW DR0,DR1,DL0 = " << lER0V_Tangled );

    lTangled = lEL1V_Tangled || lER0V_Tangled ;
  }
  else if ( !lV10Degenerate )
  {
    bool lEL0V_Tangled = counterclockwise_at_or_in_between_2(lDL0,lDR0,lDL1);
    bool lER1V_Tangled = counterclockwise_at_or_in_between_2(lDR1,lDR0,lDL1);

    CGAL_STSKEL_BUILDER_TRACE(3, "lV10Degenerate not degenerate, CCW DL0,DR0,DL1 = " << lEL0V_Tangled );
    CGAL_STSKEL_BUILDER_TRACE(3, "lV10Degenerate not degenerate, CCW DR1,DR0,DL1 = " << lER1V_Tangled );

    lTangled = lEL0V_Tangled || lER1V_Tangled ;
  }
  else
  {
    CGAL_STSKEL_BUILDER_TRACE(3, "Both degenerate, tangled = " << (lDL0 == lDL1) );

    lTangled = (lDL0 == lDL1) ;
  }

  CGAL_STSKEL_BUILDER_TRACE_IF(lTangled, 3, "Tangled profile. Pseudo-split event is invalid");

  return !lTangled ;
}

template<class Gt, class Ss, class V>
void Straight_skeleton_builder_2<Gt,Ss,V>::HandlePseudoSplitEvent( EventPtr aEvent )
{
  PseudoSplitEvent& lEvent = dynamic_cast<PseudoSplitEvent&>(*aEvent) ;

  CGAL_STSKEL_BUILDER_TRACE( 2, "Pseudo split event." );

  if ( IsValidPseudoSplitEvent(lEvent) )
  {
    CGAL_STSKEL_BUILDER_TRACE( 3, "valid event." );

    Vertex_handle lLSeed = lEvent.seed0() ;
    Vertex_handle lRSeed = lEvent.seed1() ;

    CGAL_STSKEL_BUILDER_TRACE( 4, "TriEdge e0 = " << e2str(*lEvent.triedge().e0()) );
    CGAL_STSKEL_BUILDER_TRACE( 4, "TriEdge e1 = " << e2str(*lEvent.triedge().e1()) );
    CGAL_STSKEL_BUILDER_TRACE( 4, "TriEdge e2 = " << e2str(*lEvent.triedge().e2()) );
    CGAL_STSKEL_BUILDER_TRACE( 4, "LSeed = " << v2str(*lLSeed) );
    CGAL_STSKEL_BUILDER_TRACE( 4, "LReed = " << v2str(*lRSeed) );

    Vertex_handle lNewNode_L, lNewNode_R ;
    boost::tie(lNewNode_L,lNewNode_R) = ConstructPseudoSplitEventNodes(lEvent);

    Halfedge_handle lNBisector_LO = SSkelEdgesPushBack( Halfedge(mEdgeID  ),Halfedge(mEdgeID+1) );
    Halfedge_handle lNBisector_RO = SSkelEdgesPushBack( Halfedge(mEdgeID+2),Halfedge(mEdgeID+3) );
    Halfedge_handle lNBisector_LI = lNBisector_LO->opposite();
    Halfedge_handle lNBisector_RI = lNBisector_RO->opposite();

    mEdgeID += 4 ;

    Halfedge_handle lSBisector_LO = lLSeed->primary_bisector() ;
    Halfedge_handle lSBisector_LI = lSBisector_LO->opposite();

    Halfedge_handle lSBisector_RO = lRSeed->primary_bisector() ;
    Halfedge_handle lSBisector_RI = lSBisector_RO->opposite();

    Halfedge_handle lSBisector_LO_Next = lSBisector_LO->next();
    Halfedge_handle lSBisector_RO_Next = lSBisector_RO->next();
    Halfedge_handle lSBisector_LI_Prev = lSBisector_LI->prev();
    Halfedge_handle lSBisector_RI_Prev = lSBisector_RI->prev();

    Vertex_handle lFicNod_SLO = lSBisector_LO->vertex();
    CGAL_assertion_code(Vertex_handle lFicNod_SLI = lSBisector_LI_Prev->vertex();) // unused
    Vertex_handle lFicNod_SRO = lSBisector_RO->vertex();
    CGAL_assertion_code(Vertex_handle lFicNod_SRI = lSBisector_RI_Prev->vertex();) // unused

    CGAL_assertion( lFicNod_SLO->has_infinite_time() ) ;
    CGAL_assertion( lFicNod_SLI->has_infinite_time() ) ;
    CGAL_assertion( lFicNod_SRO->has_infinite_time() ) ;
    CGAL_assertion( lFicNod_SRI->has_infinite_time() ) ;

    Link(lNBisector_LO, lSBisector_LO->face());
    Link(lNBisector_LI, lSBisector_RI->face());
    Link(lNBisector_RO, lSBisector_RO->face());
    Link(lNBisector_RI, lSBisector_LI->face());

    CrossLink(lSBisector_LO, lNewNode_L );
    CrossLink(lSBisector_RO, lNewNode_R );

    CrossLink(lNBisector_LO, lFicNod_SLO );
    CrossLink(lNBisector_RO, lFicNod_SRO );

    SetBisectorSlope(lNBisector_LO,POSITIVE);
    SetBisectorSlope(lNBisector_LI,NEGATIVE);
    SetBisectorSlope(lNBisector_RO,POSITIVE);
    SetBisectorSlope(lNBisector_RI,NEGATIVE);

    Link(lNBisector_LI, lNewNode_L);
    Link(lNBisector_RI, lNewNode_R);

    Link(lNewNode_L, lSBisector_LO);
    Link(lNewNode_R, lSBisector_RO);

    CrossLinkFwd(lSBisector_LO,lNBisector_LO);
    CrossLinkFwd(lNBisector_LO,lSBisector_LO_Next);
    CrossLinkFwd(lSBisector_LI_Prev,lNBisector_RI);
    CrossLinkFwd(lNBisector_RI,lSBisector_LI);
    CrossLinkFwd(lSBisector_RI_Prev,lNBisector_LI);
    CrossLinkFwd(lNBisector_LI,lSBisector_RI);
    CrossLinkFwd(lSBisector_RO,lNBisector_RO);
    CrossLinkFwd(lNBisector_RO,lSBisector_RO_Next);

    SetBisectorSlope(lLSeed,lNewNode_L);
    SetBisectorSlope(lRSeed,lNewNode_R);

    Halfedge_handle lNewNode_L_DefiningBorderA = lNewNode_L->halfedge()->defining_contour_edge();
    Halfedge_handle lNewNode_L_DefiningBorderB = lNewNode_L->halfedge()->next()->opposite()->defining_contour_edge();
    Halfedge_handle lNewNode_L_DefiningBorderC = lNewNode_L->halfedge()->opposite()->prev()->defining_contour_edge();
    Halfedge_handle lNewNode_R_DefiningBorderA = lNewNode_R->halfedge()->defining_contour_edge();
    Halfedge_handle lNewNode_R_DefiningBorderB = lNewNode_R->halfedge()->next()->opposite()->defining_contour_edge();
    Halfedge_handle lNewNode_R_DefiningBorderC = lNewNode_R->halfedge()->opposite()->prev()->defining_contour_edge();

    lNewNode_L->VBase::set_event_triedge( lEvent.triedge() ) ;
    lNewNode_R->VBase::set_event_triedge( lEvent.triedge() ) ;

    Triedge lTriL( lNewNode_L_DefiningBorderA, lNewNode_L_DefiningBorderB, lNewNode_L_DefiningBorderC ) ;
    Triedge lTriR( lNewNode_R_DefiningBorderA, lNewNode_R_DefiningBorderB, lNewNode_R_DefiningBorderC ) ;

    SetVertexTriedge( lNewNode_L, lTriL ) ;
    SetVertexTriedge( lNewNode_R, lTriR ) ;

    CGAL_STSKEL_BUILDER_TRACE(2,   newn2str("L",lNewNode_L,lTriL) << std::endl
                                << newn2str("R",lNewNode_R,lTriR) << std::endl
                                << newb2str("OL",lNBisector_LO)   << std::endl
                                << newb2str("IL",lNBisector_LI)   << std::endl
                                << newb2str("OR",lNBisector_RO)   << std::endl
                                << newb2str("IR",lNBisector_RI)
                             ) ;

    CGAL_STSKEL_BUILDER_TRACE( 1, "Wavefronts:\n  " << wavefront2str(lNewNode_L) << "\n  " << wavefront2str(lNewNode_R) );

    SetupNewNode(lNewNode_L) ;
    SetupNewNode(lNewNode_R) ;

    UpdatePQ(lNewNode_L, lEvent.triedge());
    UpdatePQ(lNewNode_R, lEvent.triedge());

    mVisitor.on_pseudo_split_event_processed(lLSeed,lRSeed,lNewNode_L,lNewNode_R) ;
  }
}

template<class Gt, class Ss, class V>
void Straight_skeleton_builder_2<Gt,Ss,V>::HandleSplitOrPseudoSplitEvent( EventPtr aEvent )
{
  Halfedge_handle lOppEdge = aEvent->triedge().e2() ;

  CGAL_STSKEL_BUILDER_TRACE( 2, "Split or Pseudo split event." );
  CGAL_STSKEL_BUILDER_TRACE( 3, "Opposite edge: " << e2str(*lOppEdge) );

  Site lSite;

  Vertex_handle_pair lOpp = LookupOnSLAV(lOppEdge,aEvent,lSite);

  if ( handle_assigned(lOpp.first) )
  {
    EventPtr lPseudoSplitEvent = IsPseudoSplitEvent(aEvent,lOpp,lSite);
    if ( lPseudoSplitEvent )
      HandlePseudoSplitEvent(lPseudoSplitEvent);
    else
      HandleSplitEvent(aEvent,lOpp);
  }
}

template<class Gt, class Ss, class V>
void Straight_skeleton_builder_2<Gt,Ss,V>::InsertNextSplitEventInPQ( Vertex_handle v )
{
  CGAL_STSKEL_BUILDER_TRACE(2,"Insert split event from N" << v->id() << ", "
                            << GetVertexData(v).mSplitEvents.size() << " potential splits");

  EventPtr lSplitEvent = PopNextSplitEvent(v);
  if ( !!lSplitEvent )
    InsertEventInPQ(lSplitEvent);
}

template<class Gt, class Ss, class V>
void Straight_skeleton_builder_2<Gt,Ss,V>::InsertNextSplitEventsInPQ()
{
  CGAL_STSKEL_BUILDER_TRACE(2,"Insert next split events...");
  for ( typename Vertex_handle_vector::iterator v = mReflexVertices.begin(), ev = mReflexVertices.end(); v != ev ; ++ v )
    if ( !IsProcessed(*v) )
      InsertNextSplitEventInPQ(*v);
}

template<class Gt, class Ss, class V>
void Straight_skeleton_builder_2<Gt,Ss,V>::Propagate()
{
  CGAL_STSKEL_BUILDER_TRACE(0,"Propagating events...");
  mVisitor.on_propagation_started();

  for (;;)
  {
    InsertNextSplitEventsInPQ();

    if ( !mPQ.empty() )
    {
#ifdef CGAL_SLS_PRINT_QUEUE_BEFORE_EACH_POP
      std::cout << "MAIN QUEUE -------------------------------------------------- " << std::endl;
      std::cout << "Queue size: " << mPQ.size() << std::endl;
      auto mpq = mPQ;
      while(!mpq.empty())
      {
        EventPtr event = mpq.top();
        mpq.pop();
        std::cout << *event << std::endl;
      }
      std::cout << "END MAIN QUEUE --------------------------------------------- " << std::endl;
#endif

      EventPtr lEvent = PopEventFromPQ();

      if ( lEvent->type() != Event::cEdgeEvent )
        AllowNextSplitEvent(lEvent->seed0());

      CGAL_STSKEL_BUILDER_TRACE (3,"\nTentative Event: " << *lEvent << " at ID: " << mStepID) ;

      if ( !IsProcessed(lEvent) )
      {
        CGAL_STSKEL_BUILDER_TRACE (1,"\nS" << mStepID << " Event: " << *lEvent ) ;

        SetEventTimeAndPoint(*lEvent) ;

        switch ( lEvent->type() )
        {
          case Event::cEdgeEvent       : HandleEdgeEvent              (lEvent) ; break ;
          case Event::cSplitEvent      : HandleSplitOrPseudoSplitEvent(lEvent) ; break ;
          case Event::cPseudoSplitEvent: HandlePseudoSplitEvent       (lEvent) ; break ;
        }

        ++ mStepID ;
      }
      else
      {
        CGAL_STSKEL_BUILDER_TRACE (3,"\nAlready processed") ;
      }
    }
    else break ;
  }

  mVisitor.on_propagation_finished();
}

template<class Gt, class Ss, class V>
void Straight_skeleton_builder_2<Gt,Ss,V>::MergeSplitNodes ( Vertex_handle_pair aSplitNodes )
{
  Vertex_handle lLNode, lRNode ;
  boost::tie(lLNode,lRNode)=aSplitNodes;

  Halfedge_handle lIBisectorL1 = lLNode->primary_bisector()->opposite();
  Halfedge_handle lIBisectorR1 = lRNode->primary_bisector()->opposite();
  Halfedge_handle lIBisectorL2 = lIBisectorL1->next()->opposite();
  Halfedge_handle lIBisectorR2 = lIBisectorR1->next()->opposite();

  CGAL_STSKEL_BUILDER_TRACE(2
                      ,"Merging SplitNodes: (L) N" << lLNode->id() << " and (R) N" << lRNode->id() << ".\n"
                       << "  LOut: B" << lLNode->primary_bisector()->id() << '\n'
                       << "  ROut: B" << lRNode->primary_bisector()->id() << '\n'
                       << "  LIn1: B" << lIBisectorL1->id() << '\n'
                       << "  RIn1: B" << lIBisectorR1->id() << '\n'
                       << "  LIn2: B" << lIBisectorL2->id() << '\n'
                       << "  RIn2: B" << lIBisectorR2->id()
                       );

  if ( lIBisectorL1->vertex() == lRNode )
    lIBisectorL1->HBase_base::set_vertex(lLNode);

  if ( lIBisectorR1->vertex() == lRNode )
    lIBisectorR1->HBase_base::set_vertex(lLNode);

  if ( lIBisectorL2->vertex() == lRNode )
    lIBisectorL2->HBase_base::set_vertex(lLNode);

  if ( lIBisectorR2->vertex() == lRNode )
    lIBisectorR2->HBase_base::set_vertex(lLNode);

  CGAL_STSKEL_BUILDER_TRACE(2
                      ,   "  N" << lRNode->id() << " removed.\n"
                       << "  LIn1 B" << lIBisectorL1->id() << " now linked to N" << lIBisectorL1->vertex()->id() << '\n'
                       << "  RIn1 B" << lIBisectorR1->id() << " now linked to N" << lIBisectorR1->vertex()->id() << '\n'
                       << "  LIn2 B" << lIBisectorL2->id() << " now linked to N" << lIBisectorL2->vertex()->id() << '\n'
                       << "  RIn2 B" << lIBisectorR2->id() << " now linked to N" << lIBisectorR2->vertex()->id()
                       );

  EraseNode(lRNode);
}


template<class Gt, class Ss, class V>
void Straight_skeleton_builder_2<Gt,Ss,V>::EraseNode ( Vertex_handle aNode )
{
  aNode->VBase::reset_id__internal__(-aNode->id());
  mSSkel->SSkel::Base::vertices_erase(aNode);
}

#ifdef CGAL_STRAIGHT_SKELETON_ENABLE_TRACE
template<class Halfedge_handle>
void TraceMultinode( char const* t, Halfedge_handle b, Halfedge_handle e )
{
  std::ostringstream ss ;
  ss.precision(19);
  ss << t ;

  ss << "before: B" << b->prev()->id() << " N" << b->prev()->vertex()->id() << " Pt: " << b->prev()->vertex()->point() << "\n" ;
  do
  {
    ss << "B" << b->id() << " N" << b->vertex()->id() << " Pt: " << b->vertex()->point() << "\n" ;
  }
  while ( b = b->next(), b != e ) ;

  ss << "after: B" << b->id() << " N" << b->vertex()->id() << " Pt: " << b->vertex()->point() << "\n" ;

  std::string s = ss.str();
  CGAL_STSKEL_BUILDER_TRACE(0, s);
}


template<class Point>
double angle_wrt_X ( Point const& a, Point const& b )
{
  double dx = to_double(b.x() - a.x() ) ;
  double dy = to_double(b.y() - a.y() ) ;
  double atan = std::atan2(dy,dx);
  double rad  = atan >= 0.0 ? atan : 2.0 * CGAL_PI + atan ;
  double deg  = rad * 180.0 / CGAL_PI;
  return deg ;
}

template<class Vertex_handle, class Halfedge_around_vertex_circulator>
void TraceFinalBisectors( Vertex_handle v, Halfedge_around_vertex_circulator cb )
{
  Halfedge_around_vertex_circulator c = cb ;

  do
  {
    double phi = angle_wrt_X((*c)->vertex()->point(),(*c)->opposite()->vertex()->point());

    CGAL_STSKEL_BUILDER_TRACE(2, "  N" << v->id() << " in=B" << (*c)->id()
                        << " E" << (*c)->defining_contour_edge()->id()
                        << " out=B" << (*c)->opposite()->id()
                        << " E" << (*c)->opposite()->defining_contour_edge()->id()
                        << " phi=" << phi
                        );

    ++ c ;
  }
  while( c != cb ) ;

}
#endif

template<class Vertex_handle, class Halfedge_around_vertex_circulator>
bool ValidateFinalBisectorsAfterMerge( Vertex_handle /* v */, Halfedge_around_vertex_circulator cb )
{
  bool rOK = true ;

  Halfedge_around_vertex_circulator c = cb ;

  do
  {
    if ( (*c)->defining_contour_edge() != (*c)->prev()->defining_contour_edge() )
      rOK = false ;

    ++ c ;
  }
  while( rOK && c != cb ) ;

  return rOK ;

}

template<class Gt, class Ss, class V>
void Straight_skeleton_builder_2<Gt,Ss,V>::RelinkBisectorsAroundMultinode( Vertex_handle const& v0, Halfedge_handle_vector& aLinks )
{
  CGAL_assertion( aLinks.size() > 0 ) ;

  CGAL_STSKEL_BUILDER_TRACE(4, "Relinking " << aLinks.size() << " bisectors around N" << v0->id() ) ;

  // Connect the bisectors with each other following the CCW ordering

  Halfedge_handle first_he = aLinks.front();
  Halfedge_handle prev_he  = first_he ;

  first_he->HBase_base::set_vertex(v0);

  for ( typename Halfedge_handle_vector::iterator i = std::next(aLinks.begin()), ei = aLinks.end(); i != ei ; ++ i )
  {
    Halfedge_handle he = *i ;

    he->HBase_base::set_vertex(v0);

    Halfedge_handle prev_he_opp = prev_he->opposite();

    he         ->HBase_base::set_next(prev_he_opp);
    prev_he_opp->HBase_base::set_prev(he);

    CGAL_STSKEL_BUILDER_TRACE(4, "Relinking B" << he->id() << "->B" << prev_he_opp->id() ) ;

    prev_he = he ;
  }

  Halfedge_handle prev_he_opp = prev_he->opposite();

  first_he   ->HBase_base::set_next(prev_he_opp);
  prev_he_opp->HBase_base::set_prev(first_he);

  CGAL_STSKEL_BUILDER_TRACE(4, "Relinking B" << first_he->id() << "->B" << prev_he_opp->id() ) ;

  // Reset the main halfedge for v0
  v0->VBase::set_halfedge(first_he) ;

  CGAL_STSKEL_DEBUG_CODE( TraceFinalBisectors(v0,v0->halfedge_around_vertex_begin()); )

  CGAL_postcondition( ValidateFinalBisectorsAfterMerge(v0,v0->halfedge_around_vertex_begin()) ) ;
}


template<class Gt, class Ss, class V>
void Straight_skeleton_builder_2<Gt,Ss,V>::PreprocessMultinode( Multinode& aMN )
{
  //
  // A Multinode is a run of coincident nodes along a face.
  // Its represented by a pair of halfedges describing a linear profile.
  // The first halfedge in the pair points to the first node in the multinode.
  // Each ->next() halfedge in the profile points to a subsequent node.
  // The second halfedge in the pair is past-the-end (it points to the first node around the face that IS NOT part of the multinode)
  //

  // Halfedge_handle oend = validate(aMN.end->opposite());

  CGAL_STSKEL_DEBUG_CODE( TraceMultinode("Preprocessing multinode:\n", aMN.begin,aMN.end) ) ;

  Halfedge_handle h = aMN.begin ;

  aMN.bisectors_to_relink.push_back(h);

  // Traverse the profile collecting:
  //  The nodes to be removed from the HDS (all but the first)
  //  The bisectors to be removed from the HDS (each bisector pointing to the next node in the multinode)
  //  The bisectors around each node that must be relinked to the first node (which will be kept in place of the multinode)
  do
  {
    ++ aMN.size ;
    Halfedge_handle nx = validate(h->next());
    if ( nx != aMN.end )
      aMN.bisectors_to_remove.push_back(nx);

    // Since each halfedge "h" in this lineal profile corresponds to a single face, all the bisectors around
    // each node which must be relinked are those found ccw between h and h->next()
    Halfedge_handle ccw = h ;
    Halfedge_handle ccw_end = validate(h->next()->opposite());
    for(;;)
    {
      ccw = validate(ccw->opposite()->prev()) ;
      if ( ccw != ccw_end )
        aMN.bisectors_to_relink.push_back(ccw);
      else
        break ;
    }

    if ( h != aMN.begin )
      aMN.nodes_to_remove.push_back(h->vertex());

    h = nx;
  }
  while ( h != aMN.end ) ;

  aMN.bisectors_to_relink.push_back(aMN.end->opposite());
}

//
// Replaces a run of coincident nodes with a single one by removing all but the first, removing node-to-node bisectors and
// relinking the other bisectors.
//
template<class Gt, class Ss, class V>
void Straight_skeleton_builder_2<Gt,Ss,V>::ProcessMultinode( Multinode&              aMN
                                                           , Halfedge_handle_vector& rBisectorsToRemove
                                                           , Vertex_handle_vector&   rNodesToRemove
                                                           )
{
  bool lDoNotProcess = false ;

  Halfedge_handle h = aMN.begin ;

  do
  {
    if ( h->vertex()->has_infinite_time() || IsExcluded(h->vertex()))
      lDoNotProcess = true ;
  }
  while ( h = h->next(), !lDoNotProcess && h != aMN.end ) ;

  if ( !lDoNotProcess )
  {
    CGAL_STSKEL_DEBUG_CODE( TraceMultinode("Processing multinode: ", aMN.begin,aMN.end) ) ;

    Halfedge_handle h = aMN.begin ;
    do
    {
      Exclude(h->vertex());
    }
    while ( h = h->next(), h != aMN.end ) ;

    std::copy(aMN.bisectors_to_remove.begin(), aMN.bisectors_to_remove.end(), std::back_inserter(rBisectorsToRemove));

    for( Vertex_handle_vector_iterator vi = aMN.nodes_to_remove.begin(), evi = aMN.nodes_to_remove.end() ; vi != evi ; ++ vi )
      rNodesToRemove.push_back(*vi) ;

    RelinkBisectorsAroundMultinode(aMN.v,aMN.bisectors_to_relink);
  }
}


template<class Gt, class Ss, class V>
typename Straight_skeleton_builder_2<Gt,Ss,V>::MultinodePtr
Straight_skeleton_builder_2<Gt,Ss,V>::CreateMultinode( Halfedge_handle begin, Halfedge_handle end )
{
  return MultinodePtr( new Multinode(begin,end) );
}


//
// Finds coincident skeleton nodes and merge them
//
// If moving edges Ei,Ej collide with moving edge Ek causing Ej to collapse, Ei and Ek becomes consecutive and a new
// polygon vertex (Ei,Ek) appears in the wavefront.
// If moving edges Ei,Ej collide with moving edge Ek causing Ek to be split in two halves, L(Ek) amd R(Ek) resp, two new
// polygon vertices appears in the wavefront; namely: (Ei,R(Ek)) and (L(Ek),Ej))
// If moving edge Ei,Ej collide with both Ek,El simultaneously causing the edges to cross-connect, two new vertices
// (Ei,Ek) and (El,Ej) appear in the wavefront.
//
// In all those 3 cases, each new polygon vertex is represented in the straight skeleton as a skeleton node.
// Every skeleton node is describing the collision of at least 3 edges (called the "defining edges" of the node)
// and it has at least 3 incident bisectors, each one pairing 2 out of the total number of defining edges.
//
// Any skeleton node has a degree of at least 3, but if more than 3 edges collide simultaneously, the corresponding
// skeleton node has a higher degree. (the degree of the node is exactly the number of colliding edges)
//
// However, the algorithm handles the coallison of 3 edges at a time so each skeleton node initially created
// has degree exactly 3 so this function which detects higher degree nodes and merge them into a single node
// of the proper degree is needed.
//
// Two skeleton nodes are "coincident" IFF they have 2 defining edges in common and each triedge of edges collide
// at the same time and point. IOW, 2 nodes are coincident if they represent the simultaneous
// coallison of exactly 4 edges (the union of 2 triedges with 2 common elements is a set of 4).
//
template<class Gt, class Ss, class V>
bool Straight_skeleton_builder_2<Gt,Ss,V>::MergeCoincidentNodes()
{
  //
  // NOTE: This code might be executed on a topologically incosistent HDS, thus the need to check
  // the structure along the way.
  //

  CGAL_STSKEL_BUILDER_TRACE(0, "Merging coincident nodes...");

  // ALGORITHM DESCRIPTION:
  //
  // While circulating the bisectors along the face for edge Ei we find all those edges E* which
  // are or become consecutive to Ei during the wavefront propagation. Each bisector along the face:
  // (Ei,Ea), (Ei,Eb), (Ei,Ec), etcc pairs Ei with such other edge.
  // Between one bisector (Ei,Ea) and the next (Ei,Eb) there is skeleton node which represents
  // the collision between the 3 edges (Ei,Ea,Eb).
  // It follows from the pairing that any skeleton node Ni, for example (Ei,Ea,Eb), neccesarily
  // shares two edges (Ei and Eb precisely) with any next skeleton node Ni+1 around the face.
  // That is, the triedge of defining edges that correspond to each skeleton node around the face follow this
  // sequence: (Ei,Ea,Eb), (Ei,Eb,Ec), (Ei,Ec,Ed), ...
  //
  // Any 2_ consecutive_ skeleton nodes around a face share 2 out of the 3 defining edges, which is one of the
  // neccesary conditions for "coincidence". Therefore, coincident nodes can only come as consecutive along a face
  //

  MultinodeVector lMultinodes ;

  for( Face_iterator fit = mSSkel->SSkel::Base::faces_begin(); fit != mSSkel->SSkel::Base::faces_end(); ++fit)
  {
    // 'h' is the first (CCW) skeleton halfedge.
    Halfedge_handle h = validate(validate(fit->halfedge())->next());

    CGAL_assertion ( h->is_bisector() ) ;

    // 'last' is the last (CCW) skeleton halfedge
    Halfedge_handle last = validate(fit->halfedge()->prev()) ;

    CGAL_assertion ( last->is_bisector() ) ;
    CGAL_assertion ( last->vertex()->is_contour() ) ;

    Halfedge_handle h0 = h ;
    Vertex_handle   v0 = validate(h0->vertex()) ;
    CGAL_assertion(handle_assigned(v0));

    if ( ! v0->has_infinite_time() )
    {
      CGAL_assertion ( v0->is_skeleton() ) ;

      h = validate(h->next()) ;

      while ( h != last )
      {
        Vertex_handle v = validate(h->vertex());

        if ( ! v->has_infinite_time() )
        {
          CGAL_assertion ( v->is_skeleton() ) ;

          if ( !AreSkeletonNodesCoincident(v0,v) )
          {
            if ( h0->next() != h )
              lMultinodes.push_back( CreateMultinode(h0,h) );

            v0 = v ;
            h0 = h ;
          }
        }

        h = validate(h->next());
      }

      if ( h0->next() != h )
        lMultinodes.push_back( CreateMultinode(h0,h) );
    }
  }

  if(lMultinodes.empty())
    return false;

  //
  // The merging loop removes all but one of the coincident skeleton nodes and the halfedges between them.
  // But it can't physically erase those from the HDS while looping, so the nodes/bisector to erase
  // are collected in these sequences are erased after the merging loop.
  //
  Halfedge_handle_vector lBisectorsToRemove ;
  Vertex_handle_vector   lNodesToRemove ;

  for ( typename MultinodeVector::iterator it = lMultinodes.begin(), eit = lMultinodes.end() ; it != eit ; ++ it )
    PreprocessMultinode(**it);

  std::sort(lMultinodes.begin(), lMultinodes.end(), MultinodeComparer());

  for ( typename MultinodeVector::iterator it = lMultinodes.begin(), eit = lMultinodes.end() ; it != eit ; ++ it )
    ProcessMultinode(**it,lBisectorsToRemove,lNodesToRemove);

  if(lBisectorsToRemove.empty())
    return false;

  for( Halfedge_handle_vector_iterator hi = lBisectorsToRemove.begin(), ehi = lBisectorsToRemove.end() ; hi != ehi ; ++ hi )
  {
    CGAL_STSKEL_BUILDER_TRACE(1, "B" << (*hi)->id() << " removed.");
    (*hi)->HBase_base::reset_id(-1);
    mSSkel->SSkel::Base::edges_erase(*hi);
  }

  for( Vertex_handle_vector_iterator vi = lNodesToRemove.begin(), evi = lNodesToRemove.end() ; vi != evi ; ++ vi )
    EraseNode(*vi);

  for( Vertex_iterator vit = mSSkel->SSkel::Base::vertices_begin(); vit != mSSkel->SSkel::Base::vertices_end(); ++vit)
    GetVertexData(vit).mIsExcluded = false;

  return true;
}

template<class Gt, class Ss, class V>
bool Straight_skeleton_builder_2<Gt,Ss,V>::FinishUp()
{
  CGAL_STSKEL_BUILDER_TRACE(0, "\n\nFinishing up...");

  mVisitor.on_cleanup_started();

  std::for_each( mSplitNodes.begin()
                ,mSplitNodes.end  ()
                ,boost::bind(&Straight_skeleton_builder_2<Gt,Ss,V>::MergeSplitNodes,this,_1)
               ) ;

  std::for_each( mDanglingBisectors.begin()
                ,mDanglingBisectors.end  ()
                ,boost::bind(&Straight_skeleton_builder_2<Gt,Ss,V>::EraseBisector,this,_1)
               ) ;

  // MergeCoincidentNodes() locks all extremities of halfedges that have a vertex involved in a multinode.
  // However, both extremities might have different (combinatorially and geometrically) vertices.
  // With a single pass, it would prevent one of the extremities from being properly simplified.
  // The simpliest is to just run it again as the skeleton structure is small compared to the rest
  // of the algorithm.
  for(;;)
  {
    if(!MergeCoincidentNodes())
      break;
  }

  mVisitor.on_cleanup_finished();

  // @todo if 'mMaxTime' is sufficiently large, it will be a full skeleton and should be validated as such
  if(mMaxTime) // might be a partial skeleton
    return mSSkel->is_valid(true);
  else
    return mSSkel->is_valid(false);
}

template<class Gt, class Ss, class V>
bool Straight_skeleton_builder_2<Gt,Ss,V>::Run()
{
  InitPhase();
  Propagate();
  return FinishUp();
}

template<class Gt, class Ss, class V>
typename Straight_skeleton_builder_2<Gt,Ss,V>::SSkelPtr Straight_skeleton_builder_2<Gt,Ss,V>::construct_skeleton( bool aNull_if_failed )
{
  bool ok = false ;

  try
  {
    ok = Run() ;
  }
  catch( std::exception const& e )
  {
    mVisitor.on_error ( e.what() ) ;
    CGAL_STSKEL_BUILDER_TRACE(0,"EXCEPTION THROWN (" << e.what() << ") during straight skeleton construction.");
  }
  catch(...)
  {
    mVisitor.on_error ( "Unhandled exception" ) ;
    CGAL_STSKEL_BUILDER_TRACE(0,"UNHANDLED EXCEPTION during straight skeleton construction.");
  }

  if ( !ok )
  {
    CGAL_STSKEL_BUILDER_TRACE(0,"Invalid result.");
    if ( aNull_if_failed )
      mSSkel = SSkelPtr() ;
  }

  mVisitor.on_algorithm_finished(ok);

  return mSSkel ;
}

} // end namespace CGAL

#if defined(BOOST_MSVC)
#  pragma warning(pop)
#endif

#endif // CGAL_STRAIGHT_SKELETON_BUILDER_2_IMPL_H //
// EOF //
