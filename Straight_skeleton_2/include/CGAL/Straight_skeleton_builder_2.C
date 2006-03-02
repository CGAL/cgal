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
#ifndef CGAL_STRAIGHT_SKELETON_BUILDER_2_C
#define CGAL_STRAIGHT_SKELETON_BUILDER_2_C 1

#include <boost/bind.hpp>
#include <boost/utility.hpp>
#include <CGAL/Unique_hash_map.h>

CGAL_BEGIN_NAMESPACE

namespace {

template<class Handle> inline bool handle_assigned ( Handle const& aH )
{
  Handle null ;
  return aH != null ;
}

}

template<class Gt, class SS>
Straight_skeleton_builder_2<Gt,SS>::Straight_skeleton_builder_2 ( Traits const& aTraits )
  :
  mTraits(aTraits)
 ,Left_turn(aTraits.get<typename Traits::Left_turn_2>())
 ,mEventCompare(*this)
 ,mVertexID(0)
 ,mEdgeID(0)
 ,mEventID(0)
 ,mStepID(0)
 ,mPQ(mEventCompare)
{
}

//
// This method returns the 3 distinct defining contour_edges of vertices aA and aB
// (As long as the vertices are proceesed in the right order there is 1 common defining contour_edge,
//  so there are 3 distinct contour_edges given these 2 vertices)
//
template<class Gt, class SS>
typename Straight_skeleton_builder_2<Gt,SS>::BorderTriple
Straight_skeleton_builder_2<Gt,SS>::GetDefiningBorders( Vertex_handle aA, Vertex_handle aB )
{
  Halfedge_handle lAL = GetDefiningBorderA(aA);
  Halfedge_handle lAR = GetDefiningBorderB(aA) ;
  Halfedge_handle lBL = GetDefiningBorderA(aB);
  Halfedge_handle lBR = GetDefiningBorderB(aB);

  return BorderTriple(lAL, lAR, ( lAL == lBL || lAR == lBL ) ? lBR : lBL ) ;
}

template<class Gt, class SS>
typename Straight_skeleton_builder_2<Gt,SS>::EventPtr
Straight_skeleton_builder_2<Gt,SS>::FindEdgeEvent( Vertex_handle aLNode, Vertex_handle aRNode )
{
  EventPtr rResult ;

  Halfedge_handle lBorderA, lBorderB, lBorderC ;
  boost::tie(lBorderA,lBorderB,lBorderC) = GetDefiningBorders(aLNode,aRNode);

  if ( lBorderA != lBorderB && lBorderB != lBorderC )
  {
    if ( ExistEvent(lBorderA,lBorderB,lBorderC) )
    {
      bool lAccepted = true ;

      if ( aLNode->is_skeleton() && IsNewEventInThePast(lBorderA,lBorderB,lBorderC,aLNode) )
      {
        CGAL_SSBUILDER_TRACE(1,"New edge event for Left seed N" << aLNode->id() << " is in the past. discarded." ) ;
        lAccepted = false ;
      }

      if ( aRNode->is_skeleton() && IsNewEventInThePast(lBorderA,lBorderB,lBorderC,aRNode) )
      {
        CGAL_SSBUILDER_TRACE(1,"New edge event for Right seed N" << aRNode->id() << " is in the past. discarded." ) ;
        lAccepted = false ;
      }

      if ( lAccepted )
      {
        rResult = EventPtr( new EdgeEvent( lBorderA, lBorderB, lBorderC, aLNode, aRNode ) ) ;
#ifdef CGAL_STRAIGHT_SKELETON_ENABLE_TRACE
        SetEventTimeAndPoint(*rResult);
#endif
      }
    }
  }
  return rResult ;
}


template<class Gt, class SS>
void Straight_skeleton_builder_2<Gt,SS>::CollectSplitEvent( Vertex_handle   aNode
                                                          , Halfedge_handle aReflexLBorder
                                                          , Halfedge_handle aReflexRBorder
                                                          , Halfedge_handle aOppositeBorder
                                                          )
{
  CGAL_SSBUILDER_TRACE(2,
                      "Computing potential split event between E" << aReflexLBorder->id()
                      << ",E" <<aReflexRBorder->id() << " and E" << aOppositeBorder->id()
                      );

  if ( ExistEvent(aReflexLBorder,aReflexRBorder,aOppositeBorder) )
  {
    bool lAccepted = true ;

    if ( aNode->is_skeleton() && IsNewEventInThePast(aReflexLBorder,aReflexRBorder,aOppositeBorder,aNode) )
    {
      CGAL_SSBUILDER_TRACE(1,"New split event for Seed N" << aNode->id() << " is in the past. discarded." ) ;
      lAccepted = false ;
    }

    if ( lAccepted )
    {
      EventPtr lEvent( new SplitEvent( aReflexLBorder,aReflexRBorder,aOppositeBorder,aNode,aOppositeBorder) ) ;
#ifdef CGAL_STRAIGHT_SKELETON_ENABLE_TRACE
      SetEventTimeAndPoint(*lEvent);
#endif
      mSplitEvents.push_back(lEvent);

      if ( IsReflex(aOppositeBorder->vertex()) )
        AddReflexSplit(aNode,lEvent);

      EnqueEvent(lEvent);
    }
  }
}

template<class Gt, class SS>
void Straight_skeleton_builder_2<Gt,SS>::CollectSplitEvents( Vertex_handle aNode )
{
  Vertex_handle lPrev = GetPrevInLAV(aNode) ;
  Vertex_handle lNext = GetNextInLAV(aNode) ;

  Halfedge_handle lLBorder = GetDefiningBorderA(aNode);
  Halfedge_handle lRBorder = GetDefiningBorderB(aNode) ;

  Halfedge_handle lLBorderP = lLBorder->opposite()->next()->opposite();
  Halfedge_handle lLBorderN = lLBorder->opposite()->prev()->opposite();
  Halfedge_handle lRBorderP = lRBorder->opposite()->next()->opposite();
  Halfedge_handle lRBorderN = lRBorder->opposite()->prev()->opposite();

  CGAL_SSBUILDER_TRACE(2
                      ,"Finding SplitEvent for N" << aNode->id()
                      << " LBorder: E" << lLBorder->id() << " RBorder: E" << lRBorder->id()
                      << " LBorderP: E" << lLBorderP->id() << " LBorderN: E" << lLBorderN->id()
                      << " RBorderP: E" << lRBorderP->id() << " RBorderN: E" << lRBorderN->id()
                      );

  for ( Halfedge_handle_vector_iterator i = mContourHalfedges.begin(); i != mContourHalfedges.end(); ++ i )
  {
    Halfedge_handle lOpposite = *i ;

    if (    lOpposite != lLBorder
         && lOpposite != lRBorder
         && lOpposite != lLBorderP
         && lOpposite != lLBorderN
         && lOpposite != lRBorderP
         && lOpposite != lRBorderN
       )
      CollectSplitEvent(aNode, lLBorder, lRBorder, lOpposite ) ;
  }
}


template<class Gt, class SS>
void Straight_skeleton_builder_2<Gt,SS>::CollectNewEvents( Vertex_handle aNode )
{
  Vertex_handle lPrev = GetPrevInLAV(aNode) ;
  Vertex_handle lNext = GetNextInLAV(aNode) ;

  CGAL_SSBUILDER_TRACE
    ( 1
    , "Collecting new events generated by N" << aNode->id() << " at " << aNode->point() << std::endl
      << "Prev: N" << lPrev->id() << " Next: N" << lNext->id()
    ) ;

  EventPtr lLEdgeEvent = FindEdgeEvent( lPrev , aNode ) ;
  EventPtr lREdgeEvent = FindEdgeEvent( aNode , lNext ) ;

  if ( IsReflex(aNode) )
    CollectSplitEvents(aNode) ;

  if ( lLEdgeEvent && !lREdgeEvent )
  {
    EnqueEvent(lLEdgeEvent);
  }
  else if ( lREdgeEvent && !lLEdgeEvent )
  {
    EnqueEvent(lREdgeEvent);
  }
  else if ( lLEdgeEvent && lREdgeEvent )
  {
    CGAL_SSBUILDER_TRACE(2
                        ,"Both Left and Right (candidate) EdgeEvents found:\n"
                         << "L: " << *lLEdgeEvent << '\n'
                         << "R: " << *lREdgeEvent
                        );

    bool lPickLeft = true ;

    Comparison_result lRel = CompareEvents(lLEdgeEvent,lREdgeEvent) ;
    if ( lRel == LARGER )
    {
      lPickLeft = false ;
    }
    else if ( lRel == EQUAL )
    {
      if ( CompareEventsDistanceToSeed(aNode,lLEdgeEvent,lREdgeEvent) == LARGER )
        lPickLeft = false ;
    }

    CGAL_SSBUILDER_TRACE(2,"Choosen: " << ( lPickLeft ? "Left" : "Right" ) ) ;

    if ( lPickLeft )
         EnqueEvent(lLEdgeEvent);
    else EnqueEvent(lREdgeEvent);
  }
}

//! Handles the special case of two simultaneous edge events, that is, two edges
//! collapsing along the line/point were they meet at the same time.\n
//! Algoritmically, this ocurrs when the bisector emerging from a vertex \a aA
//! is geometrically the same as the one emerging from \a aB.\n
//! Topologically, both bisectors (from \a aA and \a aB) are defined by the same pair of
//! edges.\n
//
template<class Gt, class SS>
void Straight_skeleton_builder_2<Gt,SS>::HandleSimultaneousEdgeEvent( Vertex_handle aA, Vertex_handle aB )
{
  CGAL_SSBUILDER_TRACE ( 2, "Handling simultaneous EdgeEvent between N" << aA ->id() << " and N"  << aB ->id() ) ;

  Halfedge_handle lOA = aA->primary_bisector() ;
  Halfedge_handle lOB = aB->primary_bisector() ;
  Halfedge_handle lIA = lOA->opposite();
  Halfedge_handle lIB = lOB->opposite();

  CGAL_SSBUILDER_TRACE ( 2
                       ,    "OA: B" << lOA->id() << '\n'
                         << "IA: B" << lIA->id() << '\n'
                         << "OB: B" << lOB->id() << '\n'
                         << "IB: B" << lIB->id()
                       ) ;

  SetIsProcessed(aA) ;
  SetIsProcessed(aB) ;
  mSLAV.remove(aA);
  mSLAV.remove(aB);

  CGAL_SSBUILDER_TRACE ( 2, 'N' << aA->id() << " processed\nN" << aB->id() << " processed" ) ;

  Halfedge_handle lOA_Prev = lOA->prev() ;
  Halfedge_handle lIA_Next = lIA->next() ;

  Halfedge_handle lOB_Prev = lOB->prev() ;
  Halfedge_handle lIB_Next = lIB->next() ;

  CGAL_SSBUILDER_TRACE ( 2
                       ,   "OA_Prev: B" << lOA_Prev->id() << '\n'
                         << "IA_Next: B" << lIA_Next->id() << '\n'
                         << "OB_Prev: B" << lOB_Prev->id() << '\n'
                         << "IB_Next: B" << lIB_Next->id()
                      ) ;

  lOB     ->HBase::set_next( lIA_Next );
  lIA_Next->HBase::set_prev( lOB      );
  lIB     ->HBase::set_prev( lOA_Prev );
  lOA_Prev->HBase::set_next( lIB      );

  lOB->HBase::set_vertex  (aA);

  CGAL_SSBUILDER_SHOW (1, DrawBisector(lOB) ) ;

  CGAL_SSBUILDER_TRACE ( 2, "B" << lOA->id() << " and B" << lIA->id() << " erased." ) ;
  mDanglingBisectors.push_back(lOA);

  //
  // The code above corrects the links for vertices aA/aB to the erased halfedges lOA and lIA.
  // However, any of these vertices (aA/aB) maybe one of the twin vertices of a split event.
  // If that's the case, the erased halfedge maybe be linked to a 'couple' of those vertices.
  // This situation is corrected below:

  if ( handle_assigned(lOA->vertex()) && lOA->vertex() != aA && lOA->vertex() != aB )
  {
    lOA->vertex()->VBase::set_halfedge(lIB);
    CGAL_SSBUILDER_TRACE ( 2, "N" << lOA->vertex()->id() << " has B" << lOA->id() << " as it's halfedge. Replacing it with B" << lIB->id() ) ;
  }
  if ( handle_assigned(lIA->vertex()) && lIA->vertex() != aA && lIA->vertex() != aB )
  {
    lIA->vertex()->VBase::set_halfedge(lOB);
    CGAL_SSBUILDER_TRACE ( 2, "N" << lIA->vertex()->id() << " has B" << lIA->id() << " as it's halfedge. Replacing it with B" << lOB->id() ) ;
  }

  CGAL_SSBUILDER_TRACE ( 2, "N" << aA->id() << " halfedge: B" << aA->halfedge()->id() ) ;
  CGAL_SSBUILDER_TRACE ( 2, "N" << aB->id() << " halfedge: B" << aB->halfedge()->id() ) ;

  CGAL_assertion( aA->primary_bisector() == lIB ) ;
}

template<class Gt, class SS>
bool Straight_skeleton_builder_2<Gt,SS>::AreBisectorsCoincident ( Halfedge_const_handle aA, Halfedge_const_handle aB  ) const
{
  CGAL_SSBUILDER_TRACE ( 2, "Testing for simultaneous EdgeEvents between B" << aA->id() << " and B" << aB->id() ) ;

  Halfedge_const_handle lA_LBorder = aA->defining_contour_edge();
  Halfedge_const_handle lA_RBorder = aA->opposite()->defining_contour_edge();
  Halfedge_const_handle lB_LBorder = aB->defining_contour_edge();
  Halfedge_const_handle lB_RBorder = aB->opposite()->defining_contour_edge();

  return    ( lA_LBorder == lB_LBorder && lA_RBorder == lB_RBorder )
         || ( lA_LBorder == lB_RBorder && lA_RBorder == lB_LBorder ) ;
}

template<class Gt, class SS>
void Straight_skeleton_builder_2<Gt,SS>::UpdatePQ( Vertex_handle aNode )
{
  Vertex_handle lPrev = GetPrevInLAV(aNode) ;
  Vertex_handle lNext = GetNextInLAV(aNode) ;

  CGAL_SSBUILDER_TRACE ( 2, "Updating PQ for N" << aNode->id() << " Prev N" << lPrev->id() << " Next N" << lNext->id() ) ;

  Halfedge_handle lOBisector_P = lPrev->primary_bisector() ;
  Halfedge_handle lOBisector_C = aNode->primary_bisector() ;
  Halfedge_handle lOBisector_N = lNext->primary_bisector() ;

  if ( AreBisectorsCoincident(lOBisector_C,lOBisector_P) )
    HandleSimultaneousEdgeEvent( aNode, lPrev ) ;
  else
  if ( AreBisectorsCoincident(lOBisector_C,lOBisector_N) )
    HandleSimultaneousEdgeEvent( aNode, lNext ) ;
  else
     CollectNewEvents(aNode);
}
template<class Gt, class SS>
void Straight_skeleton_builder_2<Gt,SS>::CreateInitialEvents()
{
  CGAL_SSBUILDER_TRACE(0, "Creating initial events...");
  for ( Vertex_iterator v = mSS.vertices_begin(); v != mSS.vertices_end(); ++ v )
    UpdatePQ(v);
}

template<class Gt, class SS>
typename Straight_skeleton_builder_2<Gt,SS>::EventPtr
Straight_skeleton_builder_2<Gt,SS>::FindVertexEvent( EventPtr aE0, Vertex_handle aOV )
{
  EventPtr rResult ;
  
  EventPtr_Vector lReflexSplits = GetReflexSplits(aOV) ;
  
  for ( event_iterator i = lReflexSplits.begin(), ei = lReflexSplits.end() ; i != ei ; ++ i )
  {
    EventPtr lE1 = *i ;
  
    CGAL_assertion(lE1->type() == Event::cSplitEvent ) ;
     
    if ( !lE1->is_excluded() && AreEventsSimultaneous(aE0,lE1) )
    {
      CGAL_SSBUILDER_TRACE(2, "Vertex Event found with\n" << *aE0 << "\n" << *lE1 );
  
      aE0->Exclude();
      lE1->Exclude();
  
      Halfedge_handle lBorderX[3], lBorderY[3];
  
      lBorderX[0] = aE0->border_a();
      lBorderX[1] = aE0->border_b();
      lBorderX[2] = aE0->border_c();
  
      lBorderY[0] = lE1->border_a();
      lBorderY[1] = lE1->border_b();
      lBorderY[2] = lE1->border_c();
  
      Halfedge_handle lDistinct1, lDistinct2, lEqual1, lEqual2 ;
  
      boost::tie(lDistinct1, lDistinct2, lEqual1, lEqual2) = SortTwoDistinctAndTwoEqual(lBorderX,lBorderY);
  
      CGAL_SSBUILDER_TRACE ( 2
                           , "Distinct1 E" << lDistinct1->id() << " Distinct2 E" << lDistinct2->id()
                           << " Equal1 E" << lEqual1->id() << " Equal2 E" << lEqual2->id()
                           ) ;
  
      if (    ExistEvent(lDistinct1, lDistinct2, lEqual1   )
           && ExistEvent(lEqual1   , lEqual2   , lDistinct1)
         )
      {
        Vertex_handle lSeedX = aE0->seed0();
        Vertex_handle lSeedY = lE1->seed0();
  
        rResult = EventPtr( new VertexEvent( lDistinct1, lDistinct2, lEqual1 ,lEqual2, lSeedX, lSeedY ) ) ;
  
  #ifdef CGAL_STRAIGHT_SKELETON_ENABLE_TRACE
        SetEventTimeAndPoint(*rResult);
  #endif
        CGAL_SSBUILDER_SHOW (1,  SS_IO_AUX::ScopedPointDrawing lDrawA(rResult->point(),CGAL::BLUE,"Event"); )
      }
    }
  }
  
  return rResult ;
}

template<class Gt, class SS>
typename Straight_skeleton_builder_2<Gt,SS>::EventPtr
Straight_skeleton_builder_2<Gt,SS>::FindVertexEvent( EventPtr aSplitEventPtr )
{
  EventPtr rResult ;
  
  SplitEvent& lSplitEvent = dynamic_cast<SplitEvent&>(*aSplitEventPtr) ;
  Halfedge_handle lOppBorder = lSplitEvent.opposite_border() ;
  Vertex_handle lOV1 = lOppBorder->vertex();
  if ( IsReflex(lOV1) )
    rResult = FindVertexEvent(aSplitEventPtr,lOV1);
  if ( !rResult )
  {
    Vertex_handle lOV2 = lOppBorder->opposite()->vertex();
    if ( IsReflex(lOV2) )
      rResult = FindVertexEvent(aSplitEventPtr,lOV2);
  }
  
  return rResult ;
}

/*
template<class Gt, class SS>
void Straight_skeleton_builder_2<Gt,SS>::FindVertexEvents()
{

int c = 0 ;
int k = 0 ;

  CGAL_SSBUILDER_TRACE(0, "Finding vertex events...");
  for ( Vertex_iterator v0 = mSS.vertices_begin(), ev0 = mSS.vertices_end() ; v0 != ev0 ; ++ v0 )
  {
    if ( IsReflex(v0) )
    {
      // These are all the split events against a reflex vertex.
      // That is, the event's opposite border is incident upon a reflex vertex which itself has split events,
      // some possibly also against a reflex vertex.

      EventPtr_Vector lV0ReflexSplits = GetReflexSplits(v0) ;

std::cout << "v" << v0->id() << " has " << lV0ReflexSplits.size() <<  " reflex splits" << std::endl ;

      for ( event_iterator e0 = lV0ReflexSplits.begin(), ee0 = lV0ReflexSplits.end() ; e0 != ee0 ; ++ e0 )
      {
++ k ;
        EventPtr lE0 = *e0 ;

        if ( !lE0->is_excluded() )
        {
          SplitEvent& lSE0 = dynamic_cast<SplitEvent&>(*lE0) ;

          Vertex_handle lV1 = lSE0.opposite_border()->vertex();

          CGAL_assertion(IsReflex(lV1));

          EventPtr_Vector lV1ReflexSplits = GetReflexSplits(lV1) ;

          for ( event_iterator e1 = lV1ReflexSplits.begin(), ee1 = lV1ReflexSplits.end() ; e1 != ee1 ; ++ e1 )
          {
++ c ;
            EventPtr lE1 = *e1 ;

            if ( !lE1->is_excluded() && AreEventsSimultaneous(lE0,lE1) )
            {
              CGAL_SSBUILDER_TRACE(2, "Vertex Event found with\n" << *lE0 << "\n" << *lE1 );

              lE0->Exclude();
              lE1->Exclude();

              Halfedge_handle lBorderX[3], lBorderY[3];

              lBorderX[0] = lE0->border_a();
              lBorderX[1] = lE0->border_b();
              lBorderX[2] = lE0->border_c();

              lBorderY[0] = lE1->border_a();
              lBorderY[1] = lE1->border_b();
              lBorderY[2] = lE1->border_c();

              Halfedge_handle lDistinct1, lDistinct2, lEqual1, lEqual2 ;

              boost::tie(lDistinct1, lDistinct2, lEqual1, lEqual2) = SortTwoDistinctAndTwoEqual(lBorderX,lBorderY);

              CGAL_SSBUILDER_TRACE ( 2
                                   , "Distinct1 E" << lDistinct1->id() << " Distinct2 E" << lDistinct2->id()
                                   << " Equal1 E" << lEqual1->id() << " Equal2 E" << lEqual2->id()
                                   ) ;

              if (    ExistEvent(lDistinct1, lDistinct2, lEqual1   )
                   && ExistEvent(lEqual1   , lEqual2   , lDistinct1)
                 )
              {
                Vertex_handle lSeedX = lE0->seed0();
                Vertex_handle lSeedY = lE1->seed0();

                EventPtr lVertexEvent = EventPtr( new VertexEvent( lDistinct1, lDistinct2, lEqual1 ,lEqual2, lSeedX, lSeedY ) ) ;

    #ifdef CGAL_STRAIGHT_SKELETON_ENABLE_TRACE
                SetEventTimeAndPoint(*lVertexEvent);
    #endif
                CGAL_SSBUILDER_SHOW (1,  SS_IO_AUX::ScopedPointDrawing lDrawA(lVertexEvent->point(),CGAL::BLUE,"Event"); )

                EnqueEvent(lVertexEvent);

                break ;
              }
            }
          }
        }

        if ( lE0->is_excluded() )
          break ;
      }
    }
  }

std::cout << "total reflex splits:" << k << ". Total innner iterations: " << c << std::endl ;
}
*/

template<class Gt, class SS>
void Straight_skeleton_builder_2<Gt,SS>::CreateContourBisectors()
{
  CGAL_SSBUILDER_TRACE(0, "Creating contour bisectors...");
  for ( Vertex_iterator v = mSS.vertices_begin(); v != mSS.vertices_end(); ++ v )
  {
    // NOTE: Bisectors are always contructed with no geometric embedding.
    Halfedge lOB(mEdgeID++), lIB(mEdgeID++);
    Halfedge_handle lOBisector = mSS.SBase::edges_push_back (lOB, lIB);
    Halfedge_handle lIBisector = lOBisector->opposite();
    lOBisector->HBase::set_face(v->halfedge()->face());
    lIBisector->HBase::set_face(v->halfedge()->next()->face());
    lIBisector->HBase::set_vertex(v);

    Halfedge_handle lIBorder = v->halfedge() ;
    Halfedge_handle lOBorder = v->halfedge()->next() ;
    lIBorder  ->HBase::set_next(lOBisector);
    lOBisector->HBase::set_prev(lIBorder);
    lOBorder  ->HBase::set_prev(lIBisector);
    lIBisector->HBase::set_next(lOBorder);
    CGAL_SSBUILDER_TRACE(2
                        ,"Adding Contour Bisector at N:" << v->id() << "\n B" << lOBisector->id()
                        << " (Out)\n B" << lIBisector->id() << " (In)"
                        ) ;
  }
}

template<class Gt, class SS>
void Straight_skeleton_builder_2<Gt,SS>::InitPhase()
{
  CreateContourBisectors();
  CreateInitialEvents();
  //FindVertexEvents();
}

template<class Gt, class SS>
typename Straight_skeleton_builder_2<Gt,SS>::Vertex_handle
Straight_skeleton_builder_2<Gt,SS>::ConstructEdgeEventNode( EdgeEvent& aEvent )
{
  CGAL_SSBUILDER_TRACE ( 2, "Creating EdgeEvent Node" ) ;

  Vertex_handle lLSeed = aEvent.seed0() ;
  Vertex_handle lRSeed = aEvent.seed1() ;

  SetEventTimeAndPoint(aEvent);

  Vertex_handle lNewNode = mSS.SBase::vertices_push_back( Vertex( mVertexID++, aEvent.point(), aEvent.time()) ) ;
  mSLAV.push_back(lNewNode);

  mWrappedVertices.push_back( VertexWrapper(lNewNode) ) ;

  Halfedge_handle lLOBisector = lLSeed->primary_bisector();
  Halfedge_handle lROBisector = lRSeed->primary_bisector();
  Halfedge_handle lLIBisector = lLOBisector->opposite();
  Halfedge_handle lRIBisector = lROBisector->opposite();

  lNewNode->VBase::set_halfedge(lLOBisector);
  lLOBisector->HBase::set_vertex(lNewNode);
  lROBisector->HBase::set_vertex(lNewNode);

  lLIBisector->HBase::set_prev( lROBisector ) ;
  lROBisector->HBase::set_next( lLIBisector ) ;

  CGAL_SSBUILDER_SHOW(1, DrawBisector(lLOBisector); DrawBisector(lROBisector); ) ;

  CGAL_SSBUILDER_TRACE
  ( 2
  ,    "LSeed: N" << lLSeed->id() << " proccesed\n"
    << "RSeed: N" << lRSeed->id() << " proccesed"
  ) ;

  SetIsProcessed(lLSeed) ;
  SetIsProcessed(lRSeed) ;
  mSLAV.remove(lLSeed);
  mSLAV.remove(lRSeed);

  Vertex_handle lLPrev = GetPrevInLAV(lLSeed) ;
  Vertex_handle lRNext = GetNextInLAV(lRSeed) ;

  SetPrevInLAV(lNewNode, lLPrev ) ;
  SetNextInLAV(lLPrev  , lNewNode  ) ;

  SetNextInLAV(lNewNode, lRNext ) ;
  SetPrevInLAV(lRNext  , lNewNode  ) ;

  CGAL_SSBUILDER_TRACE
  ( 2
  ,    "LO: B" << lLOBisector->id() << " LI: B" << lLIBisector->id() << " RO: B" << lROBisector->id() << " RI: B" << lRIBisector->id() << '\n'
    << "New Node: N" << lNewNode->id() << " at " << lNewNode->point() << '\n'
    << "New Links: B" << lROBisector->id() << "->B" << lLIBisector->id() << '\n'
    << 'N' << lNewNode->id() << " inserted into LAV: N" << lLPrev->id() << "->N" << lNewNode->id() << "->N" << lRNext->id() << std::endl
    << 'N' << lLSeed->id() << " removed from LAV\n"
    << 'N' << lRSeed->id() << " removed from LAV"
  );


  return lNewNode ;
}


template<class Gt, class SS>
typename Straight_skeleton_builder_2<Gt,SS>::Vertex_handle
Straight_skeleton_builder_2<Gt,SS>::LookupOnSLAV ( Halfedge_handle aBorder, Event const& aEvent )
{
  Vertex_handle rResult ;

  CGAL_SSBUILDER_TRACE ( 2, "Looking up for E" << aBorder->id() << " on SLAV. P=" << aEvent.point() ) ;

  CGAL_SSBUILDER_SHOW ( 2
                      , SS_IO_AUX::ScopedSegmentDrawing draw_(aBorder->vertex()->point()
                                                             ,aBorder->opposite()->vertex()->point()
                                                             ,CGAL::YELLOW,"OppBorder") ;
                      )

  for ( typename std::list<Vertex_handle>::const_iterator vi = mSLAV.begin(); vi != mSLAV.end(); ++ vi )
  {
    Vertex_handle v = *vi;

    if (  handle_assigned(GetPrevInLAV(v))
       && handle_assigned(GetNextInLAV(v))
       && GetDefiningBorderA(v) == aBorder
       )
    {
      Vertex_handle lPrev = GetPrevInLAV(v);
      Halfedge_handle lPrevBorder = GetDefiningBorderA(lPrev);
      Halfedge_handle lNextBorder = GetDefiningBorderB(v);
      if ( IsEventInsideOffsetZone( aEvent.border_a()
                                  , aEvent.border_b()
                                  , aBorder
                                  , lPrevBorder
                                  , lNextBorder
                                 )
        )
      {
        rResult = v ;
        CGAL_SSBUILDER_TRACE ( 2, 'N' << rResult->id() << " found" ) ;
        break ;
      }
    }
  }

  return rResult ;
}


template<class Gt, class SS>
typename Straight_skeleton_builder_2<Gt,SS>::Vertex_handle_pair
Straight_skeleton_builder_2<Gt,SS>::ConstructSplitEventNodes( SplitEvent& aEvent, Vertex_handle aOppR )
{
  Vertex_handle_pair rResult;

  CGAL_SSBUILDER_TRACE ( 2, "Creating SplitEvent Nodes" ) ;

  Vertex_handle lOppL = GetPrevInLAV(aOppR) ;

  SetEventTimeAndPoint(aEvent);

  Vertex_handle lNodeA = mSS.SBase::vertices_push_back( Vertex( mVertexID++, aEvent.point(), aEvent.time() ) ) ;
  Vertex_handle lNodeB = mSS.SBase::vertices_push_back( Vertex( mVertexID++, aEvent.point(), aEvent.time() ) ) ;

  mSLAV.push_back(lNodeA);
  mSLAV.push_back(lNodeB);
  mWrappedVertices.push_back( VertexWrapper(lNodeA) ) ;
  mWrappedVertices.push_back( VertexWrapper(lNodeB) ) ;

  Vertex_handle lSeed = aEvent.seed0() ;

  Halfedge_handle lXOutBisector = lSeed->primary_bisector() ;
  Halfedge_handle lXInBisector  = lXOutBisector->opposite();

  lNodeA->VBase::set_halfedge(lXOutBisector);
  // lNodeB hafledge is set outside with the New In Bisector to the Right.

  lXOutBisector->HBase::set_vertex(lNodeA);

  CGAL_SSBUILDER_SHOW(1, DrawBisector(lXOutBisector); ) ;

  CGAL_SSBUILDER_TRACE ( 2, "Seed: N" << lSeed->id() << " proccesed" ) ;

  SetIsProcessed(lSeed) ;
  mSLAV.remove(lSeed);

  CGAL_SSBUILDER_TRACE ( 2, 'N' << lNodeA->id() << " and N" << lNodeB->id() << " inserted into LAV." ) ;

  Vertex_handle lPrev = GetPrevInLAV(lSeed) ;
  Vertex_handle lNext = GetNextInLAV(lSeed) ;

  SetNextInLAV(lPrev , lNodeA ) ;
  SetPrevInLAV(lNodeA, lPrev  ) ;

  SetNextInLAV(lNodeA, aOppR  ) ;
  SetPrevInLAV(aOppR , lNodeA ) ;

  SetNextInLAV(lOppL , lNodeB ) ;
  SetPrevInLAV(lNodeB, lOppL  ) ;

  SetNextInLAV(lNodeB, lNext  ) ;
  SetPrevInLAV(lNext , lNodeB ) ;

  CGAL_SSBUILDER_TRACE
  (
   2
   ,   "Updated LAV: N" << lPrev->id() << "->N" << lNodeA->id() << "->N" << aOppR->id() << std::endl
    << "Updated LAV: N" << lOppL->id() << "->N" << lNodeB->id() << "->N" << lNext->id() << std::endl
    << 'N' << lSeed->id() << " removed from LAV"
  );

  rResult = std::make_pair(lNodeA,lNodeB);

  mSplitNodes.push_back(rResult);

  return rResult ;
}

template<class Gt, class SS>
typename Straight_skeleton_builder_2<Gt,SS>::Vertex_handle_pair
Straight_skeleton_builder_2<Gt,SS>::ConstructVertexEventNodes( VertexEvent& aEvent )
{
  Vertex_handle_pair rResult;

  CGAL_SSBUILDER_TRACE ( 2, "Creating VertexEvent Nodes" ) ;

  Vertex_handle lLSeed = aEvent.seed0() ;
  Vertex_handle lRSeed = aEvent.seed1() ;

  SetEventTimeAndPoint(aEvent);

  Vertex_handle lNewNodeA = mSS.SBase::vertices_push_back( Vertex( mVertexID++, aEvent.point(), aEvent.time()) ) ;
  Vertex_handle lNewNodeB = mSS.SBase::vertices_push_back( Vertex( mVertexID++, aEvent.point(), aEvent.time()) ) ;

  mSLAV.push_back(lNewNodeA);
  mSLAV.push_back(lNewNodeB);

  mWrappedVertices.push_back( VertexWrapper(lNewNodeA) ) ;
  mWrappedVertices.push_back( VertexWrapper(lNewNodeB) ) ;

  Halfedge_handle lLOBisector = lLSeed->primary_bisector();
  Halfedge_handle lROBisector = lRSeed->primary_bisector();
  Halfedge_handle lLIBisector = lLOBisector->opposite();
  Halfedge_handle lRIBisector = lROBisector->opposite();

  lNewNodeA->VBase::set_halfedge(lLOBisector);
  lNewNodeB->VBase::set_halfedge(lROBisector);
  lLOBisector->HBase::set_vertex(lNewNodeA);
  lROBisector->HBase::set_vertex(lNewNodeB);

  lLIBisector->HBase::set_prev( lROBisector ) ;
  lROBisector->HBase::set_next( lLIBisector ) ;

  lLOBisector->HBase::set_next( lRIBisector ) ;
  lRIBisector->HBase::set_prev( lLOBisector ) ;

  CGAL_SSBUILDER_SHOW(1, DrawBisector(lLOBisector); DrawBisector(lROBisector); ) ;

  CGAL_SSBUILDER_TRACE
  (
   2
   ,   "LSeed: N" << lLSeed->id() << " proccesed\n"
    << "RSeed: N" << lRSeed->id() << " proccesed"
  ) ;

  SetIsProcessed(lLSeed) ;
  SetIsProcessed(lRSeed) ;
  mSLAV.remove(lLSeed);
  mSLAV.remove(lRSeed);

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


  CGAL_SSBUILDER_TRACE
  (
  2
  ,    "LO: B" << lLOBisector->id() << " LI: B" << lLIBisector->id() << " RO: B" << lROBisector->id() << " RI: B" << lRIBisector->id() << '\n'
    << "NewNodeA: N" << lNewNodeA->id() << " at " << lNewNodeA->point() << '\n'
    << "NewNodeB: N" << lNewNodeB->id() << " at " << lNewNodeB->point() << '\n'
    << "New Links: B" << lROBisector->id() << "->B" << lLIBisector->id() << '\n'
    << 'N' << lNewNodeA->id() << " and N" << lNewNodeB->id() << " inserted into LAV:\n"
    << 'N' << lLPrev->id() << "->N" << lNewNodeA->id() << "->N" << lRNext->id() << '\n'
    << 'N' << lRPrev->id() << "->N" << lNewNodeB->id() << "->N" << lLNext->id() << '\n'
    << 'N' << lLSeed->id() << " removed from LAV\n"
    << 'N' << lRSeed->id() << " removed from LAV"
  );

  rResult = std::make_pair(lNewNodeA,lNewNodeB);

  mSplitNodes.push_back(rResult);

  return rResult ;
}

template<class Gt, class SS>
bool Straight_skeleton_builder_2<Gt,SS>::IsProcessed( EventPtr aEvent )
{
  return IsProcessed(aEvent->seed0()) || IsProcessed(aEvent->seed1()) ;
}

template<class Gt, class SS>
void Straight_skeleton_builder_2<Gt,SS>::HandleEdgeEvent( EventPtr aEvent )
{
  EdgeEvent& lEvent = dynamic_cast<EdgeEvent&>(*aEvent) ;

  Vertex_handle lLSeed = lEvent.seed0() ;
  Vertex_handle lRSeed = lEvent.seed1() ;

  Vertex_handle lNewNode = ConstructEdgeEventNode(lEvent);

  Halfedge_handle lLOBisector = lLSeed->primary_bisector() ;
  Halfedge_handle lROBisector = lRSeed->primary_bisector() ;
  Halfedge_handle lLIBisector = lLOBisector->opposite();
  Halfedge_handle lRIBisector = lROBisector->opposite();

  if ( !handle_assigned(lLOBisector->next()) && !handle_assigned(lRIBisector->prev()) )
  {
    CGAL_SSBUILDER_TRACE(2,"Creating new Edge Event's Bisector");

    Halfedge_handle lNOBisector = mSS.SBase::edges_push_back ( Halfedge(mEdgeID),Halfedge(mEdgeID+1) );

    Halfedge_handle lNIBisector = lNOBisector->opposite();
    mEdgeID += 2 ;

    lRIBisector->HBase::set_prev(lNIBisector);
    lNIBisector->HBase::set_next(lRIBisector);

    lNOBisector->HBase::set_face(lLOBisector->HBase::face());
    lNIBisector->HBase::set_face(lRIBisector->HBase::face());
    lNIBisector->HBase::set_vertex(lNewNode);

    lLOBisector->HBase::set_next(lNOBisector);
    lNOBisector->HBase::set_prev(lLOBisector);

    Halfedge_handle lDefiningBorderA = lNewNode->halfedge()->face()->halfedge();
    Halfedge_handle lDefiningBorderB = lNewNode->halfedge()->opposite()->prev()->opposite()->face()->halfedge();
    Halfedge_handle lDefiningBorderC = lNewNode->halfedge()->opposite()->prev()->face()->halfedge();

    SetDefiningBorderA(lNewNode, lDefiningBorderA) ;
    SetDefiningBorderB(lNewNode, lDefiningBorderB) ;
    SetDefiningBorderC(lNewNode, lDefiningBorderC) ;

    CGAL_SSBUILDER_TRACE
      ( 2
      , "NewNode N" << lNewNode->id() << " at " << lNewNode->point() << " defining borders: E"
        << lDefiningBorderA->id()
        << ",E" << lDefiningBorderB->id()
        << ",E" << lDefiningBorderC->id() << '\n'
        << "New Bisectors:\nB" << lNOBisector->id() << " [E" << lNOBisector->defining_contour_edge()->id()
        << ",E" << lNOBisector->opposite()->defining_contour_edge()->id()
        << "] (Out: Prev: B" << lNOBisector->prev()->id() << ")\nB"
        << lNIBisector->id() << " [E" << lNIBisector->defining_contour_edge()->id()
        << ",E" << lNIBisector->opposite()->defining_contour_edge()->id()
        << "] (In: Next: B" << lNIBisector->next()->id() << ")\n"
        << "N" << lNewNode->id() << " halfedge: " << lNewNode->halfedge()->id()
        << " primary bisector: B" << lNewNode->primary_bisector()->id()
      ) ;

    UpdatePQ(lNewNode);
  }
  else
  {
    CGAL_SSBUILDER_TRACE(3, "N" << lNewNode->id() << " is a multiple node. Bisector not created");
  }
}



template<class Gt, class SS>
void Straight_skeleton_builder_2<Gt,SS>::HandleSplitEvent( EventPtr aEvent, Vertex_handle aOppR )
{
  SplitEvent& lEvent = dynamic_cast<SplitEvent&>(*aEvent) ;

  Vertex_handle lSeed = lEvent.seed0();

  Vertex_handle lNewNode_L, lNewNode_R ;
  boost::tie(lNewNode_L,lNewNode_R) = ConstructSplitEventNodes(lEvent,aOppR);

  Halfedge_handle lOppBorder = lEvent.opposite_border();

  Halfedge_handle lReflexLBorder = GetDefiningBorderA(lSeed);
  Halfedge_handle lReflexRBorder = GetDefiningBorderB(lSeed);

  Halfedge_handle lNOBisector_L = mSS.SBase::edges_push_back ( Halfedge(mEdgeID++),Halfedge(mEdgeID++) );
  Halfedge_handle lNOBisector_R = mSS.SBase::edges_push_back ( Halfedge(mEdgeID++),Halfedge(mEdgeID++) );
  Halfedge_handle lNIBisector_L = lNOBisector_L->opposite();
  Halfedge_handle lNIBisector_R = lNOBisector_R->opposite();

  lNewNode_R->VBase::set_halfedge(lNIBisector_L) ;

  Halfedge_handle lXOBisector = lSeed->primary_bisector() ;
  Halfedge_handle lXIBisector = lXOBisector->opposite();

  lNOBisector_L->HBase::set_face(lXOBisector->HBase::face());
  lNIBisector_L->HBase::set_face(lOppBorder ->HBase::face());
  lNOBisector_R->HBase::set_face(lOppBorder ->HBase::face());
  lNIBisector_R->HBase::set_face(lXIBisector->HBase::face());

  lNIBisector_L->HBase::set_vertex(lNewNode_R);
  lNIBisector_R->HBase::set_vertex(lNewNode_R);

  lXOBisector  ->HBase::set_next(lNOBisector_L);
  lNOBisector_L->HBase::set_prev(lXOBisector);

  lXIBisector  ->HBase::set_prev(lNIBisector_R);
  lNIBisector_R->HBase::set_next(lXIBisector);

  lNIBisector_L->HBase::set_next(lNOBisector_R);
  lNOBisector_R->HBase::set_prev(lNIBisector_L);

  Halfedge_handle lNewNodeLDefiningBorderA = lNewNode_L->halfedge()->face()->halfedge();
  Halfedge_handle lNewNodeLDefiningBorderB = lNewNode_L->halfedge()->opposite()->prev()->opposite()->face()->halfedge();
  Halfedge_handle lNewNodeLDefiningBorderC = lNewNode_L->halfedge()->opposite()->prev()->face()->halfedge();
  Halfedge_handle lNewNodeRDefiningBorderA = lNewNode_R->halfedge()->face()->halfedge();
  Halfedge_handle lNewNodeRDefiningBorderB = lNewNode_R->halfedge()->opposite()->prev()->opposite()->face()->halfedge();
  Halfedge_handle lNewNodeRDefiningBorderC = lNewNode_R->halfedge()->opposite()->prev()->face()->halfedge();

  SetDefiningBorderA(lNewNode_L, lNewNodeLDefiningBorderA) ;
  SetDefiningBorderB(lNewNode_L, lNewNodeLDefiningBorderB) ;
  SetDefiningBorderC(lNewNode_L, lNewNodeLDefiningBorderC) ;

  SetDefiningBorderA(lNewNode_R, lNewNodeRDefiningBorderA) ;
  SetDefiningBorderB(lNewNode_R, lNewNodeRDefiningBorderB) ;
  SetDefiningBorderC(lNewNode_R, lNewNodeRDefiningBorderC) ;

  CGAL_SSBUILDER_TRACE
    (
     2
    ,    "New Node L: N" << lNewNode_L->id() << " at " << lNewNode_L->point()
      << " defining borders: E" << lNewNodeLDefiningBorderA->id()
      << ",E" << lNewNodeLDefiningBorderB->id()
      << ",E" << lNewNodeLDefiningBorderC->id() << '\n'
      << "New Node R: N" << lNewNode_R->id() << " at " << lNewNode_R->point()
      << " defining borders: E" << lNewNodeRDefiningBorderA->id()
      << ",E" << lNewNodeLDefiningBorderB->id() << '\n'
      << ",E" << lNewNodeRDefiningBorderC->id() << '\n'
      << "New Bisector OL:\nB" << lNOBisector_L->id()
      << "[E"  << lNOBisector_L            ->defining_contour_edge()->id()
      << ",E" << lNOBisector_L->opposite()->defining_contour_edge()->id() << "]"
      << " (Out: Prev: B" << lNOBisector_L->prev()->id() << ")\n"
      << "New Bisector IL:\nB" << lNIBisector_L->id()
      << "[E" << lNIBisector_L            ->defining_contour_edge()->id()
      << ",E" << lNIBisector_L->opposite()->defining_contour_edge()->id() << "]"
      << " (In: Next: B" << lNIBisector_L->next()->id() << ")\n"
      << "New Bisector OR:\nB" << lNOBisector_R->id()
      << "[E" << lNOBisector_R            ->defining_contour_edge()->id()
      << ",E" << lNOBisector_R->opposite()->defining_contour_edge()->id() << "]"
      << " (Out: Prev: B" << lNOBisector_R->prev()->id() << ")\n"
      << "New Bisector IR:\nB" << lNIBisector_R->id()
      << "[E" << lNIBisector_R            ->defining_contour_edge()->id()
      << ",E" << lNIBisector_R->opposite()->defining_contour_edge()->id()
      << "] (In: Next: B" << lNIBisector_R->next()->id() << ")\n"
      << "N" << lNewNode_L->id() << " halfedge: " << lNewNode_L->halfedge()->id()
      << " primary bisector: B" << lNewNode_L->primary_bisector()->id() << '\n'
      << "N" << lNewNode_R->id() << " halfedge: " << lNewNode_R->halfedge()->id()
      << " primary bisector: B" << lNewNode_R->primary_bisector()->id()
    ) ;

  UpdatePQ(lNewNode_L);
  UpdatePQ(lNewNode_R);
}

template<class Gt, class SS>
bool Straight_skeleton_builder_2<Gt,SS>::SetupVertexEventNode( Vertex_handle   aNode
                                                             , Halfedge_handle aDefiningBorderA
                                                             , Halfedge_handle aDefiningBorderB
                                                             )
{
  bool rR = false ;

  Point_2 p = aDefiningBorderA->opposite()->vertex()->point() ;
  Point_2 q = aDefiningBorderA->opposite()->prev()->vertex()->point() ;
  Point_2 r = aDefiningBorderB->opposite()->prev()->vertex()->point() ;

  bool lCollinear = Collinear(p,q,r) ;

  if ( lCollinear || !Left_turn(p,q,r) )
  {
    rR = true ;
    SetIsReflex(aNode);
    CGAL_SSBUILDER_TRACE(2, ( lCollinear ? "COLLINEAR ":"Reflex " ) << "*NEW* vertex: N" << aNode->id() );
  }

  return rR ;
}

template<class Gt, class SS>
void Straight_skeleton_builder_2<Gt,SS>::HandleVertexEvent( EventPtr aEvent )
{
  VertexEvent& lEvent = dynamic_cast<VertexEvent&>(*aEvent) ;

  Vertex_handle lLSeed = lEvent.seed0() ;
  Vertex_handle lRSeed = lEvent.seed1() ;

  Vertex_handle lNewNode_L, lNewNode_R ;
  boost::tie(lNewNode_L,lNewNode_R) = ConstructVertexEventNodes(lEvent);

  Halfedge_handle lNBisector_LO = mSS.SBase::edges_push_back ( Halfedge(mEdgeID++),Halfedge(mEdgeID++) );
  Halfedge_handle lNBisector_RO = mSS.SBase::edges_push_back ( Halfedge(mEdgeID++),Halfedge(mEdgeID++) );
  Halfedge_handle lNBisector_LI = lNBisector_LO->opposite();
  Halfedge_handle lNBisector_RI = lNBisector_RO->opposite();

  Halfedge_handle lSBisector_LO = lLSeed->primary_bisector() ;
  Halfedge_handle lSBisector_LI = lSBisector_LO->opposite();

  Halfedge_handle lSBisector_RO = lRSeed->primary_bisector() ;
  Halfedge_handle lSBisector_RI = lSBisector_RO->opposite();

  lNBisector_LO->HBase::set_face(lSBisector_LO->HBase::face());
  lNBisector_LI->HBase::set_face(lSBisector_RI->HBase::face());
  lNBisector_RO->HBase::set_face(lSBisector_RO->HBase::face());
  lNBisector_RI->HBase::set_face(lSBisector_LI->HBase::face());

  lNBisector_LI->HBase::set_vertex(lNewNode_L);
  lNBisector_RI->HBase::set_vertex(lNewNode_R);

  lSBisector_LO->HBase::set_next(lNBisector_LO);
  lNBisector_LO->HBase::set_prev(lSBisector_LO);

  lSBisector_LI->HBase::set_prev(lNBisector_RI);
  lNBisector_RI->HBase::set_next(lSBisector_LI);

  lSBisector_RI->HBase::set_prev(lNBisector_LI);
  lNBisector_LI->HBase::set_next(lSBisector_RI);

  lSBisector_RO->HBase::set_next(lNBisector_RO);
  lNBisector_RO->HBase::set_prev(lSBisector_RO);

  lNewNode_L->VBase::set_halfedge(lSBisector_LO);
  lNewNode_R->VBase::set_halfedge(lSBisector_RO);

  Halfedge_handle lNewNodeLDefiningBorderA = lNewNode_L->halfedge()->face()->halfedge();
  Halfedge_handle lNewNodeLDefiningBorderB = lNewNode_L->halfedge()->next()->opposite()->face()->halfedge();
  Halfedge_handle lNewNodeLDefiningBorderC = lNewNode_L->halfedge()->opposite()->prev()->face()->halfedge();
  Halfedge_handle lNewNodeRDefiningBorderA = lNewNode_R->halfedge()->face()->halfedge();
  Halfedge_handle lNewNodeRDefiningBorderB = lNewNode_R->halfedge()->next()->opposite()->face()->halfedge();
  Halfedge_handle lNewNodeRDefiningBorderC = lNewNode_R->halfedge()->opposite()->prev()->face()->halfedge();

  SetDefiningBorderA(lNewNode_L, lNewNodeLDefiningBorderA) ;
  SetDefiningBorderB(lNewNode_L, lNewNodeLDefiningBorderB) ;
  SetDefiningBorderC(lNewNode_L, lNewNodeLDefiningBorderC) ;

  SetDefiningBorderA(lNewNode_R, lNewNodeRDefiningBorderA) ;
  SetDefiningBorderB(lNewNode_R, lNewNodeRDefiningBorderB) ;
  SetDefiningBorderC(lNewNode_R, lNewNodeRDefiningBorderC) ;

  CGAL_SSBUILDER_TRACE
    (2
    ,    "New Node L: N" << lNewNode_L->id() << " at " << lNewNode_L->point()
      << " defining borders: E" << lNewNodeLDefiningBorderA->id()
      << ",E" << lNewNodeLDefiningBorderB->id()
      << ",E" << lNewNodeLDefiningBorderC->id() << '\n'
      << "New Node R: N" << lNewNode_R->id() << " at " << lNewNode_R->point()
      << " defining borders: E" << lNewNodeRDefiningBorderA->id()
      << ",E" << lNewNodeRDefiningBorderB->id()
      << ",E" << lNewNodeRDefiningBorderC->id() << '\n'
      << "New Bisector LO: B" << lNBisector_LO->id()
      << "[E"  << lNBisector_LO            ->defining_contour_edge()->id()
      << ",E" << lNBisector_LO->opposite()->defining_contour_edge()->id() << "]\n"
      << "New Bisector LI: B" << lNBisector_LI->id()
      << "[E" << lNBisector_LI            ->defining_contour_edge()->id()
      << ",E" << lNBisector_LI->opposite()->defining_contour_edge()->id() << "]\n"
      << "New Bisector RO: B" << lNBisector_RO->id()
      << "[E" << lNBisector_RO            ->defining_contour_edge()->id()
      << ",E" << lNBisector_RO->opposite()->defining_contour_edge()->id() << "]\n"
      << "New Bisector RI: B" << lNBisector_RI->id()
      << "[E" << lNBisector_RI            ->defining_contour_edge()->id()
      << ",E" << lNBisector_RI->opposite()->defining_contour_edge()->id() << "]\n"
      << "N" << lNewNode_L->id() << " halfedge: " << lNewNode_L->halfedge()->id()
      << " primary bisector: B" << lNewNode_L->primary_bisector()->id() << '\n'
      << "N" << lNewNode_R->id() << " halfedge: " << lNewNode_R->halfedge()->id()
      << " primary bisector: B" << lNewNode_R->primary_bisector()->id()
    ) ;

  bool lNodeLIsNonConvex = SetupVertexEventNode(lNewNode_L,lNewNodeLDefiningBorderA,lNewNodeLDefiningBorderB) ;
  if ( !lNodeLIsNonConvex )
    SetupVertexEventNode(lNewNode_R,lNewNodeRDefiningBorderA,lNewNodeRDefiningBorderB) ;

  UpdatePQ(lNewNode_L);
  UpdatePQ(lNewNode_R);
}

template<class Gt, class SS>
void Straight_skeleton_builder_2<Gt,SS>::HandlePotentialSplitEvent( EventPtr aEvent )
{
  SplitEvent& lEvent = dynamic_cast<SplitEvent&>(*aEvent);
  
  Halfedge_handle lOppBorder = lEvent.opposite_border() ;

  Vertex_handle lOppVertex = LookupOnSLAV(lOppBorder,lEvent);

  if ( handle_assigned(lOppVertex) )
  {
    EventPtr lVertexEvent = FindVertexEvent(aEvent);
    if ( !lVertexEvent )
         HandleSplitEvent (aEvent,lOppVertex);
    else HandleVertexEvent(lVertexEvent); 
  }
  else 
  {
    CGAL_SSBUILDER_TRACE(3,"Split event is no longer valid. Opposite edge vanished.");
  }
}

template<class Gt, class SS>
void Straight_skeleton_builder_2<Gt,SS>::Propagate()
{
  CGAL_SSBUILDER_TRACE(0,"Propagating events...");

  while ( !mPQ.empty() )
  {
    EventPtr lEvent = PopEventFromPQ();

    if ( !lEvent->is_excluded() && !IsProcessed(lEvent) )
    {
      CGAL_SSBUILDER_TRACE (0,"\nStep: " << mStepID << " Event: " << *lEvent ) ;
      CGAL_SSBUILDER_SHOW (1, SS_IO_AUX::ScopedPointDrawing lDraw(lEvent->point(),CGAL::BLUE,"Event"); )

      switch ( lEvent->type() )
      {
        case Event::cEdgeEvent : HandleEdgeEvent          (lEvent) ;  break ;
        case Event::cSplitEvent: HandlePotentialSplitEvent(lEvent) ;  break ;
      }

      ++ mStepID ;
    }
  }

}

template<class Gt, class SS>
void Straight_skeleton_builder_2<Gt,SS>::MergeSplitNodes ( Vertex_handle_pair aSplitNodes )
{
  Vertex_handle lLNode, lRNode ;
  boost::tie(lLNode,lRNode)=aSplitNodes;

  Halfedge_handle lIBisectorL1 = lLNode->primary_bisector()->opposite();
  Halfedge_handle lIBisectorR1 = lRNode->primary_bisector()->opposite();
  Halfedge_handle lIBisectorL2 = lIBisectorL1->next()->opposite();
  Halfedge_handle lIBisectorR2 = lIBisectorR1->next()->opposite();

  Exclude(lRNode);

  if ( lIBisectorL1->vertex() == lRNode )
    lIBisectorL1->HBase::set_vertex(lLNode);

  if ( lIBisectorR1->vertex() == lRNode )
    lIBisectorR1->HBase::set_vertex(lLNode);

  if ( lIBisectorL2->vertex() == lRNode )
    lIBisectorL2->HBase::set_vertex(lLNode);

  if ( lIBisectorR2->vertex() == lRNode )
    lIBisectorR2->HBase::set_vertex(lLNode);

  CGAL_SSBUILDER_TRACE(2
                      ,"SplitNodes: N" << lLNode->id() << " and N" << lRNode->id() << " merged.\n"
                       << ". N" << lRNode->id() << " excluded.\n"
                       << 'B' << lIBisectorL1->id() << " now linked to N" << lIBisectorL1->vertex()->id() << '\n'
                       << 'B' << lIBisectorR1->id() << " now linked to N" << lIBisectorR1->vertex()->id() << '\n'
                       << 'B' << lIBisectorL2->id() << " now linked to N" << lIBisectorL2->vertex()->id() << '\n'
                       << 'B' << lIBisectorR2->id() << " now linked to N" << lIBisectorR2->vertex()->id()
                       );

  mSS.SBase::vertices_erase(lRNode);
}

template<class Gt, class SS>
void Straight_skeleton_builder_2<Gt,SS>::FinishUp()
{
  std::for_each( mSplitNodes.begin()
                ,mSplitNodes.end  ()
                ,boost::bind(&Straight_skeleton_builder_2<Gt,SS>::MergeSplitNodes,this,_1)
               ) ;

  std::for_each( mDanglingBisectors.begin()
                ,mDanglingBisectors.end  ()
                ,boost::bind(&Straight_skeleton_builder_2<Gt,SS>::EraseBisector,this,_1)
               ) ;

}

template<class Gt, class SS>
void Straight_skeleton_builder_2<Gt,SS>::Run()
{
  InitPhase();
  Propagate();
  FinishUp ();
}

template<class Gt, class SS>
typename Straight_skeleton_builder_2<Gt,SS>::Ssds Straight_skeleton_builder_2<Gt,SS>::construct_skeleton()
{
  Run();
  return mSS ;
}
CGAL_END_NAMESPACE

#endif // CGAL_STRAIGHT_SKELETON_BUILDER_2_C //
// EOF //
