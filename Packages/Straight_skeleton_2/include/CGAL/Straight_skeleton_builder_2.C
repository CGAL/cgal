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
// file          : include/CGAL/Straight_skeleton_builder_2.c
// package       : Straight_skeleton_2 (1.1.0)
//
// author(s)     : Fernando Cacciola
// maintainer    : Fernando Cacciola <fernando_cacciola@hotmail>
// coordinator   : Fernando Cacciola <fernando_cacciola@hotmail>
//
// ============================================================================
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
  tie(lBorderA,lBorderB,lBorderC) = GetDefiningBorders(aLNode,aRNode);

  if ( lBorderA != lBorderB && lBorderB != lBorderC )
  {
    if ( ExistEvent(lBorderA,lBorderB,lBorderC) )
    {
      rResult = EventPtr( new EdgeEvent( lBorderA, lBorderB, lBorderC, aLNode, aRNode ) ) ;
#ifdef CGAL_STRAIGHT_SKELETON_ENABLE_TRACE
      SetEventTimeAndPoint(*rResult);
#endif
    }
  }
  return rResult ;
}


template<class Gt, class SS>
void Straight_skeleton_builder_2<Gt,SS>::CollectSplitEvent( Vertex_handle    aNode
                                                           ,Halfedge_handle  aReflexLBorder
                                                           ,Halfedge_handle  aReflexRBorder
                                                           ,Halfedge_handle  aOppositeBorder
                                                           ,EventPtr_Vector& aCandidates
                                                          )
{
  CGAL_SSBUILDER_TRACE("Computing potential split event between E" << aReflexLBorder->id()
                      << ",E" <<aReflexRBorder->id() << " and E" << aOppositeBorder->id()
                      );

  if ( ExistEvent(aReflexLBorder,aReflexRBorder,aOppositeBorder) )
  {
    Halfedge_handle lPrevOppBorder = GetPrevInCCB(aOppositeBorder);
    Halfedge_handle lNextOppBorder = GetNextInCCB(aOppositeBorder);

    if ( IsEventInsideOffsetZone(aReflexLBorder, aReflexRBorder, aOppositeBorder, lPrevOppBorder, lNextOppBorder) )
    {
       EventPtr lEvent( new SplitEvent( aReflexLBorder,aReflexRBorder,aOppositeBorder,aNode,aOppositeBorder) ) ;
#ifdef CGAL_STRAIGHT_SKELETON_ENABLE_TRACE
       SetEventTimeAndPoint(*lEvent);
#endif
       aCandidates.push_back(lEvent) ;
     }
  }
}


template<class Gt, class SS>
typename Straight_skeleton_builder_2<Gt,SS>::EventPtr
Straight_skeleton_builder_2<Gt,SS>::FindSplitEvent( Vertex_handle aNode )
{
  EventPtr rEvent ;

  Halfedge_handle lLBorder1 = GetDefiningBorderA(aNode);
  Halfedge_handle lRBorder1 = GetDefiningBorderB(aNode) ;

  Halfedge_handle lLBorder2, lRBorder2;

  if ( aNode->is_skeleton() )
  {
    Vertex_handle lPrev = GetPrevInLAV(aNode) ;
    Vertex_handle lNext = GetNextInLAV(aNode) ;

    Halfedge_handle lPrevBorder1 = GetDefiningBorderA(lPrev);
    Halfedge_handle lPrevBorder2 = GetDefiningBorderB(lPrev) ;
    Halfedge_handle lNextBorder1 = GetDefiningBorderA(lNext);
    Halfedge_handle lNextBorder2 = GetDefiningBorderB(lNext) ;

    lLBorder2 = ( lPrevBorder1 == lLBorder1 ? lPrevBorder2 : lPrevBorder1 ) ;
    lRBorder2 = ( lNextBorder1 == lRBorder1 ? lNextBorder2 : lNextBorder1 ) ;
  }
  else
  {
    lLBorder2 = lLBorder1 ;
    lRBorder2 = lRBorder1 ;
  }

  CGAL_SSBUILDER_TRACE("Finding SplitEvent for N" << aNode->id()
                      << " LBorder1: E" << lLBorder1->id() << " RBorder1: E" << lRBorder1->id()
                      << " LBorder2: E" << lLBorder2->id() << " RBorder2: E" << lRBorder2->id()
                      );

  EventPtr_Vector lCandidates ;

  for( typename Halfedge_handle_vector::const_iterator it = mBorderHalfedges.begin()
       , eit = mBorderHalfedges.end()
       ; it != eit
       ; ++ it
     )
  {
    Halfedge_handle lOppBorder = *it ;

    if (    lOppBorder != lLBorder1
         && lOppBorder != lLBorder2
         && lOppBorder != lRBorder1
         && lOppBorder != lRBorder2
       )
    {
      CollectSplitEvent(aNode, lLBorder1, lRBorder1, *it, lCandidates ) ;
    }
  }


  if ( !lCandidates.empty() )
  {
    CGAL_SSBUILDER_TRACE( (int)lCandidates.size() << " candidate Split Isecs computed.");

    rEvent = *std::min_element( lCandidates.begin(), lCandidates.end(), mEventCompare) ;

    CGAL_SSBUILDER_TRACE("NewSplitEvent found: " << *rEvent);

    return rEvent ;
  }

  return rEvent ;
}


template<class Gt, class SS>
void Straight_skeleton_builder_2<Gt,SS>::CollectNewEvents( Vertex_handle aNode, EventPtr_Vector& aEvents )
{
  Vertex_handle lPrev = GetPrevInLAV(aNode) ;
  Vertex_handle lNext = GetNextInLAV(aNode) ;

  CGAL_SSBUILDER_TRACE
    ( "Collecting new events generated by N" << aNode->id() << " at " << aNode->point() << std::endl
      << "Prev: N" << lPrev->id() << " Next: N" << lNext->id()
    ) ;

  EventPtr lLEdgeEvent = FindEdgeEvent( lPrev , aNode ) ;
  EventPtr lREdgeEvent = FindEdgeEvent( aNode , lNext ) ;
  EventPtr lSplitEvent = IsReflex(aNode) ? FindSplitEvent(aNode) : EventPtr() ;

  if ( lLEdgeEvent && !lREdgeEvent )
  {
    aEvents.push_back(lLEdgeEvent) ;
  }
  else if ( lREdgeEvent && !lLEdgeEvent )
  {
    aEvents.push_back(lREdgeEvent) ;
  }
  else if ( lLEdgeEvent && lREdgeEvent )
  {
    CGAL_SSBUILDER_TRACE("Both Left and Right (candidate) EdgeEvents found:\n"
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

    CGAL_SSBUILDER_TRACE("Choosen: " << ( lPickLeft ? "Left" : "Right" ) ) ;

    if ( lPickLeft )
         aEvents.push_back(lLEdgeEvent) ;
    else aEvents.push_back(lREdgeEvent) ;
  }

  if ( lSplitEvent )
  {
    CGAL_SSBUILDER_TRACE("Candidate Split Event found: " << *lSplitEvent ) ;
    aEvents.push_back(lSplitEvent) ;
  }
}

template<class Gt, class SS>
typename Straight_skeleton_builder_2<Gt,SS>::EventPtr
Straight_skeleton_builder_2<Gt,SS>::ChooseBestNewEvent( Vertex_handle aNode, EventPtr_Vector& aEvents )
{
  size_type lC = aEvents.size();

  if ( lC == 0 )
    return EventPtr();
  else
  if ( lC == 1 )
    return aEvents.front();
  else
    return *std::min_element(aEvents.begin(),aEvents.end(),mEventCompare) ;
}


template<class Gt, class SS>
typename Straight_skeleton_builder_2<Gt,SS>::EventPtr
Straight_skeleton_builder_2<Gt,SS>::FindBestNewEvent( Vertex_handle aNode )
{
  EventPtr_Vector lEvents;
  CollectNewEvents(aNode,lEvents);
  return ChooseBestNewEvent(aNode,lEvents);
}

template<class Gt, class SS>
void Straight_skeleton_builder_2<Gt,SS>::AddNewEvent( Vertex_handle aNode )
{
  if ( EventPtr lEvent = FindBestNewEvent(aNode) )
  {
    CGAL_SSBUILDER_TRACE ( "New Event: " << *lEvent ) ;

    CGAL_SSBUILDER_SHOW ( SS_IO_AUX::ScopedPointDrawing lDraw(lEvent->point(),CGAL::BLUE,"Event"); )

    if ( lEvent->type() == Event::cSplitEvent )
      mSplitEvents.push_back(lEvent);

    mPQ.push_back(lEvent);
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
  CGAL_SSBUILDER_TRACE ( "Handling simultaneous EdgeEvent between N" << aA ->id() << " and N"  << aB ->id() ) ;

  Halfedge_handle lOA = aA->primary_bisector() ;
  Halfedge_handle lOB = aB->primary_bisector() ;
  Halfedge_handle lIA = lOA->opposite();
  Halfedge_handle lIB = lOB->opposite();

  CGAL_SSBUILDER_TRACE (    "OA: B" << lOA->id() << '\n'
                         << "IA: B" << lIA->id() << '\n'
                         << "OB: B" << lOB->id() << '\n'
                         << "IB: B" << lIB->id()
                       ) ;

  SetIsProcessed(aA) ;
  SetIsProcessed(aB) ;
  CGAL_SSBUILDER_TRACE ( 'N' << aA->id() << " processed\nN" << aB->id() << " processed" ) ;

  Halfedge_handle lOA_Prev = lOA->prev() ;
  Halfedge_handle lIA_Next = lIA->next() ;

  Halfedge_handle lOB_Prev = lOB->prev() ;
  Halfedge_handle lIB_Next = lIB->next() ;

  CGAL_SSBUILDER_TRACE (    "OA_Prev: B" << lOA_Prev->id() << '\n'
                         << "IA_Next: B" << lIA_Next->id() << '\n'
                         << "OB_Prev: B" << lOB_Prev->id() << '\n'
                         << "IB_Next: B" << lIB_Next->id()
                      ) ;

  lOB     ->HBase::set_next( lIA_Next );
  lIA_Next->HBase::set_prev( lOB      );
  lIB     ->HBase::set_prev( lOA_Prev );
  lOA_Prev->HBase::set_next( lIB      );

  lOB->HBase::set_vertex(aA);

  CGAL_SSBUILDER_SHOW ( DrawBisector(lOB) ) ;

  CGAL_SSBUILDER_TRACE ( "B" << lOA->id() << " and B" << lIA->id() << " erased.\n" ) ;
  mDanglingBisectors.push_back(lOA);

  CGAL_assertion( aA->primary_bisector() == lIB ) ;
}

template<class Gt, class SS>
bool Straight_skeleton_builder_2<Gt,SS>::AreBisectorsCoincident ( Halfedge_const_handle aA, Halfedge_const_handle aB  ) const
{
  CGAL_SSBUILDER_TRACE ( "Testing for simultaneous EdgeEvents between B" << aA->id() << " and B" << aB->id() ) ;

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

  CGAL_SSBUILDER_TRACE ( "Updating PQ for N" << aNode->id() << " Prev N" << lPrev->id() << " Next N" << lNext->id() ) ;

  Halfedge_handle lOBisector_P = lPrev->primary_bisector() ;
  Halfedge_handle lOBisector_C = aNode->primary_bisector() ;
  Halfedge_handle lOBisector_N = lNext->primary_bisector() ;

  if ( AreBisectorsCoincident(lOBisector_C,lOBisector_P) )
    HandleSimultaneousEdgeEvent( aNode, lPrev ) ;
  else
  if ( AreBisectorsCoincident(lOBisector_C,lOBisector_N) )
    HandleSimultaneousEdgeEvent( aNode, lNext ) ;
  else
     AddNewEvent(aNode);
}
template<class Gt, class SS>
void Straight_skeleton_builder_2<Gt,SS>::CreateInitialEvents()
{
  CGAL_SSBUILDER_TRACE("Creating initial events...");
  for ( Vertex_iterator v = mSS.vertices_begin(); v != mSS.vertices_end(); ++ v )
    UpdatePQ(v);
}

template<class Gt, class SS>
void Straight_skeleton_builder_2<Gt,SS>::FindVertexEvents()
{
  // true if an event was merged with another
  Unique_hash_map<Event*,bool> lMap(false,mSplitEvents.size()) ;

  CGAL_SSBUILDER_TRACE("Finding vertex events...");
  for ( event_iterator itx = mSplitEvents.begin(); itx != mSplitEvents.end(); ++ itx )
  {
    EventPtr x = *itx ;

    for ( event_iterator ity = boost::next(itx); ity != mSplitEvents.end(); ++ ity )
    {
      EventPtr y = *ity ;

      if ( !lMap[boost::get_pointer(y)] ) // y was not already merged
      {
        if ( AreEventsSimultaneous(x,y) )
        {
          SplitEvent& lX = dynamic_cast<SplitEvent&>(*x) ;
          SplitEvent& lY = dynamic_cast<SplitEvent&>(*y) ;

          CGAL_SSBUILDER_TRACE("Vertex Event found with\n" << *x << "\n" << *y);
          lMap[boost::get_pointer(x)]=true;
          lMap[boost::get_pointer(y)]=true;

          Vertex_handle lSeedX = lX.seed();
          Vertex_handle lSeedY = lY.seed();

          Halfedge_handle lBorderX[3], lBorderY[3];

          lBorderX[0] = lX.border_a();
          lBorderX[1] = lX.border_b();
          lBorderX[2] = lX.border_c();

          lBorderY[0] = lY.border_a();
          lBorderY[1] = lY.border_b();
          lBorderY[2] = lY.border_c();

          Halfedge_handle lDistinct1, lDistinct2, lEqual1, lEqual2 ;

          tie(lDistinct1, lDistinct2, lEqual1, lEqual2) = SortTwoDistinctAndTwoEqual(lBorderX,lBorderY);

          CGAL_SSBUILDER_TRACE ( "Distinct1 E" << lDistinct1->id() << " Distinct2 E" << lDistinct2->id()
                               << " Equal1 E" << lEqual1->id() << " Equal2 E" << lEqual2->id()
                               ) ;

          CGAL_assertion( ExistEvent(lDistinct1, lDistinct2, lEqual1   ) ) ;
          CGAL_assertion( ExistEvent(lEqual1   , lEqual2   , lDistinct1) ) ;

          EventPtr lVertexEvent = EventPtr( new VertexEvent( lDistinct1, lDistinct2, lEqual1 ,lEqual2, lSeedX, lSeedY ) ) ;

#ifdef CGAL_STRAIGHT_SKELETON_ENABLE_TRACE
          SetEventTimeAndPoint(*lVertexEvent);
#endif
          CGAL_SSBUILDER_TRACE ( "VertexEvent: " << *lVertexEvent ) ;

          CGAL_SSBUILDER_SHOW ( SS_IO_AUX::ScopedPointDrawing lDrawA(lVertexEvent->point(),CGAL::BLUE,"Event"); )

          mPQ.push_back(lVertexEvent);
          break ;
        }
      }
    }
  }

  for ( event_iterator itx = mSplitEvents.begin(); itx != mSplitEvents.end(); ++ itx )
  {
    EventPtr x = *itx ;
    if ( lMap[boost::get_pointer(x)] )
      mPQ.remove(x);
  }
}

template<class Gt, class SS>
void Straight_skeleton_builder_2<Gt,SS>::CreateContourBisectors()
{
  CGAL_SSBUILDER_TRACE("Creating contour bisectors...");
  for ( Vertex_iterator v = mSS.vertices_begin(); v != mSS.vertices_end(); ++ v )
  {
    // NOTE: Bisectors are always contructed with no geometric embedding.
    Halfedge lOB(mEdgeID++), lIB(mEdgeID++);
    Halfedge_handle lOBisector = mSS.SBase::edges_push_back (lOB, lIB);
    Halfedge_handle lIBisector = lOBisector->opposite();
    mWrappedHalfedges.push_back( HalfedgeWrapper(lOBisector) ) ;
    mWrappedHalfedges.push_back( HalfedgeWrapper(lIBisector) ) ;
    lOBisector->HBase::set_face(v->halfedge()->face());
    lIBisector->HBase::set_face(v->halfedge()->next()->face());
    lIBisector->HBase::set_vertex(v);

    Halfedge_handle lIBorder = v->halfedge() ;
    Halfedge_handle lOBorder = v->halfedge()->next() ;
    lIBorder  ->HBase::set_next(lOBisector);
    lOBisector->HBase::set_prev(lIBorder);
    lOBorder  ->HBase::set_prev(lIBisector);
    lIBisector->HBase::set_next(lOBorder);
    CGAL_SSBUILDER_TRACE("Adding Contour Bisector at N:" << v->id() << "\n B" << lOBisector->id()
                        << " (Out)\n B" << lIBisector->id() << " (In)"
                        ) ;
  }
}

template<class Gt, class SS>
void Straight_skeleton_builder_2<Gt,SS>::InitPhase()
{
  CreateContourBisectors();
  CreateInitialEvents();
  FindVertexEvents();
}

template<class Gt, class SS>
typename Straight_skeleton_builder_2<Gt,SS>::Vertex_handle
Straight_skeleton_builder_2<Gt,SS>::ConstructEdgeEventNode( EdgeEvent& aEvent )
{
  CGAL_SSBUILDER_TRACE ( "Creating EdgeEvent Node" ) ;

  Vertex_handle lLSeed = aEvent.left_seed () ;
  Vertex_handle lRSeed = aEvent.right_seed() ;

  SetEventTimeAndPoint(aEvent);

  Vertex_handle lNewNode = mSS.SBase::vertices_push_back( Vertex( mVertexID++, aEvent.point(), aEvent.time()) ) ;

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

  CGAL_SSBUILDER_SHOW( DrawBisector(lLOBisector); DrawBisector(lROBisector); ) ;

  CGAL_SSBUILDER_TRACE
  (
       "LSeed: N" << lLSeed->id() << " proccesed\n"
    << "RSeed: N" << lRSeed->id() << " proccesed"
  ) ;

  SetIsProcessed(lLSeed) ;
  SetIsProcessed(lRSeed) ;

  Vertex_handle lLPrev = GetPrevInLAV(lLSeed) ;
  Vertex_handle lRNext = GetNextInLAV(lRSeed) ;

  SetPrevInLAV(lNewNode, lLPrev ) ;
  SetNextInLAV(lLPrev  , lNewNode  ) ;

  SetNextInLAV(lNewNode, lRNext ) ;
  SetPrevInLAV(lRNext  , lNewNode  ) ;

  CGAL_SSBUILDER_TRACE
  (    "LO: B" << lLOBisector->id() << " LI: B" << lLIBisector->id() << " RO: B" << lROBisector->id() << " RI: B" << lRIBisector->id() << '\n'
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

  CGAL_SSBUILDER_TRACE ( "Looking up for E" << aBorder->id() << " on SLAV. P=" << aEvent.point() ) ;

  CGAL_SSBUILDER_SHOW_AUX ( SS_IO_AUX::ScopedSegmentDrawing draw_(aBorder->vertex()->point()
                                                                 ,aBorder->opposite()->vertex()->point()
                                                                 ,CGAL::YELLOW,"OppBorder") ;
                          )

  for ( Vertex_iterator vi = mSS.vertices_begin(); vi != mSS.vertices_end(); ++ vi )
  {
    Vertex_handle v = static_cast<Vertex_handle>(vi);

    if (     handle_assigned(GetPrevInLAV(v))
         &&  handle_assigned(GetNextInLAV(v))
         &&  GetDefiningBorderA(v) == aBorder
         && !IsProcessed(v)
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
        CGAL_SSBUILDER_TRACE ( 'N' << rResult->id() << " found" ) ;
        break ;
      }
    }
  }

  return rResult ;
}

template<class Gt, class SS>
typename Straight_skeleton_builder_2<Gt,SS>::Vertex_handle_pair
Straight_skeleton_builder_2<Gt,SS>::ConstructSplitEventNodes( SplitEvent& aEvent )
{
  Vertex_handle_pair rResult;

  CGAL_SSBUILDER_TRACE ( "Creating SplitEvent Nodes" ) ;

  Halfedge_handle lOppBorder = aEvent.opposite_border() ;

  Vertex_handle lOppR = LookupOnSLAV(lOppBorder,aEvent);

  if ( handle_assigned(lOppR) )
  {
    Vertex_handle lOppL = GetPrevInLAV(lOppR) ;

    CGAL_SSBUILDER_TRACE ( "Opposite E" << lOppBorder->id() << " is between N" << lOppL->id() << " and N" << lOppR->id() ) ;

    SetEventTimeAndPoint(aEvent);

    Vertex_handle lNodeA = mSS.SBase::vertices_push_back( Vertex( mVertexID++, aEvent.point(), aEvent.time() ) ) ;
    Vertex_handle lNodeB = mSS.SBase::vertices_push_back( Vertex( mVertexID++, aEvent.point(), aEvent.time() ) ) ;

    mWrappedVertices.push_back( VertexWrapper(lNodeA) ) ;
    mWrappedVertices.push_back( VertexWrapper(lNodeB) ) ;

    Halfedge_handle lXOutBisector = aEvent.seed()->primary_bisector() ;
    Halfedge_handle lXInBisector  = lXOutBisector->opposite();

    lNodeA->VBase::set_halfedge(lXOutBisector);
    // lNodeB hafledge is set outside with the New In Bisector to the Right.

    lXOutBisector->HBase::set_vertex(lNodeA);

    Vertex_handle lSeed = aEvent.seed() ;

    CGAL_SSBUILDER_SHOW( DrawBisector(lXOutBisector); ) ;

    CGAL_SSBUILDER_TRACE ( "Seed: N" << lSeed->id() << " proccesed" ) ;

    SetIsProcessed(lSeed) ;

    CGAL_SSBUILDER_TRACE ( 'N' << lNodeA->id() << " and N" << lNodeB->id() << " inserted into LAV." ) ;

    Vertex_handle lPrev = GetPrevInLAV(lSeed) ;
    Vertex_handle lNext = GetNextInLAV(lSeed) ;

    SetNextInLAV(lPrev , lNodeA ) ;
    SetPrevInLAV(lNodeA, lPrev  ) ;

    SetNextInLAV(lNodeA, lOppR  ) ;
    SetPrevInLAV(lOppR , lNodeA ) ;

    SetNextInLAV(lOppL , lNodeB ) ;
    SetPrevInLAV(lNodeB, lOppL  ) ;

    SetNextInLAV(lNodeB, lNext  ) ;
    SetPrevInLAV(lNext , lNodeB ) ;

    CGAL_SSBUILDER_TRACE
    (
         "Updated LAV: N" << lPrev->id() << "->N" << lNodeA->id() << "->N" << lOppR->id() << std::endl
      << "Updated LAV: N" << lOppL->id() << "->N" << lNodeB->id() << "->N" << lNext->id() << std::endl
      << 'N' << lSeed->id() << " removed from LAV"
    );

    rResult = std::make_pair(lNodeA,lNodeB);

    mSplitNodes.push_back(rResult);
  }

  return rResult ;
}

template<class Gt, class SS>
typename Straight_skeleton_builder_2<Gt,SS>::Vertex_handle_pair
Straight_skeleton_builder_2<Gt,SS>::ConstructVertexEventNodes( VertexEvent& aEvent )
{
  Vertex_handle_pair rResult;

  CGAL_SSBUILDER_TRACE ( "Creating VertexEvent Nodes" ) ;

  Vertex_handle lLSeed = aEvent.left_seed () ;
  Vertex_handle lRSeed = aEvent.right_seed() ;

  SetEventTimeAndPoint(aEvent);

  Vertex_handle lNewNodeA = mSS.SBase::vertices_push_back( Vertex( mVertexID++, aEvent.point(), aEvent.time()) ) ;
  Vertex_handle lNewNodeB = mSS.SBase::vertices_push_back( Vertex( mVertexID++, aEvent.point(), aEvent.time()) ) ;

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

  CGAL_SSBUILDER_SHOW( DrawBisector(lLOBisector); DrawBisector(lROBisector); ) ;

  CGAL_SSBUILDER_TRACE
  (
       "LSeed: N" << lLSeed->id() << " proccesed\n"
    << "RSeed: N" << lRSeed->id() << " proccesed"
  ) ;

  SetIsProcessed(lLSeed) ;
  SetIsProcessed(lRSeed) ;

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
  (    "LO: B" << lLOBisector->id() << " LI: B" << lLIBisector->id() << " RO: B" << lROBisector->id() << " RI: B" << lRIBisector->id() << '\n'
    << "NewNodeA: N" << lNewNodeA->id() << " at " << lNewNodeA->point() << '\n'
    << "NewNodeB: N" << lNewNodeB->id() << " at " << lNewNodeB->point() << '\n'
    << "New Links: B" << lROBisector->id() << "->B" << lLIBisector->id() << '\n'
    << 'N' << lNewNodeA->id() << " and N" << lNewNodeB->id() << " inserted into LAV:\n"
    << 'N' << lLPrev->id() << "->N" << lNewNodeA->id() << "->N" << lRNext->id() << '\n'
    << 'N' << lRPrev->id() << "->N" << lNewNodeB->id() << "->N" << lLNext->id() << '\n'
    << 'N' << lLSeed->id() << " removed from LAV\n"
    << 'N' << lRSeed->id() << " removed from LAV"
  );

  SetIsReflex(lNewNodeB);

  rResult = std::make_pair(lNewNodeA,lNewNodeB);

  mSplitNodes.push_back(rResult);

  return rResult ;
}

template<class Gt, class SS>
void Straight_skeleton_builder_2<Gt,SS>::HandleEdgeEvent( EdgeEvent& aEvent )
{
  if ( !IsProcessed(aEvent.left_seed()) && !IsProcessed(aEvent.right_seed()) )
  {
    Vertex_handle lNewNode = ConstructEdgeEventNode(aEvent);

    Halfedge_handle lLOBisector = aEvent.left_seed ()->primary_bisector() ;
    Halfedge_handle lROBisector = aEvent.right_seed()->primary_bisector() ;
    Halfedge_handle lLIBisector = lLOBisector->opposite();
    Halfedge_handle lRIBisector = lROBisector->opposite();

    if ( !handle_assigned(lLOBisector->next()) && !handle_assigned(lRIBisector->prev()) )
    {
      CGAL_SSBUILDER_TRACE("Creating new Edge Event's Bisector");

      Halfedge_handle lNOBisector = mSS.SBase::edges_push_back ( Halfedge(mEdgeID),Halfedge(mEdgeID+1) );

      Halfedge_handle lNIBisector = lNOBisector->opposite();
      mEdgeID += 2 ;

      mWrappedHalfedges.push_back( HalfedgeWrapper(lNOBisector) ) ;
      mWrappedHalfedges.push_back( HalfedgeWrapper(lNIBisector) ) ;

      lRIBisector->HBase::set_prev(lNIBisector);
      lNIBisector->HBase::set_next(lRIBisector);

      lNOBisector->HBase::set_face(lLOBisector->HBase::face());
      lNIBisector->HBase::set_face(lRIBisector->HBase::face());
      lNIBisector->HBase::set_vertex(lNewNode);

      lLOBisector->HBase::set_next(lNOBisector);
      lNOBisector->HBase::set_prev(lLOBisector);

      CGAL_SSBUILDER_TRACE
        ( "New Bisectors:\nB" << lNOBisector->id() << " [E" << lNOBisector->defining_contour_edge()->id()
          << ",E" << lNOBisector->opposite()->defining_contour_edge()->id()
          << "] (Out: Prev: B" << lNOBisector->prev()->id() << ")\nB"
          << lNIBisector->id() << " [E" << lNIBisector->defining_contour_edge()->id()
          << ",E" << lNIBisector->opposite()->defining_contour_edge()->id()
          << "] (In: Next: B" << lNIBisector->next()->id() << ")"
        ) ;

      UpdatePQ(lNewNode);
    }
    else
    {
      CGAL_SSBUILDER_TRACE("N" << lNewNode->id() << " is a multiple node. Bisector not created");
    }
  }
}

template<class Gt, class SS>
void Straight_skeleton_builder_2<Gt,SS>::HandleSplitEvent( SplitEvent& aEvent )
{
  if ( !IsProcessed(aEvent.seed()) )
  {
    Vertex_handle lNewNode_L, lNewNode_R ;
    boost::tie(lNewNode_L,lNewNode_R) = ConstructSplitEventNodes(aEvent);

    if ( handle_assigned(lNewNode_L) && handle_assigned(lNewNode_R) )
    {
      Vertex_handle lSeed = aEvent.seed();

      Halfedge_handle lOppBorder = aEvent.opposite_border();

      Halfedge_handle lReflexLBorder = GetDefiningBorderA(lSeed);
      Halfedge_handle lReflexRBorder = GetDefiningBorderB(lSeed);

      Halfedge_handle lNOBisector_L = mSS.SBase::edges_push_back ( Halfedge(mEdgeID++),Halfedge(mEdgeID++) );
      Halfedge_handle lNOBisector_R = mSS.SBase::edges_push_back ( Halfedge(mEdgeID++),Halfedge(mEdgeID++) );
      Halfedge_handle lNIBisector_L = lNOBisector_L->opposite();
      Halfedge_handle lNIBisector_R = lNOBisector_R->opposite();

      mWrappedHalfedges.push_back( HalfedgeWrapper(lNOBisector_L) ) ;
      mWrappedHalfedges.push_back( HalfedgeWrapper(lNOBisector_R) ) ;
      mWrappedHalfedges.push_back( HalfedgeWrapper(lNIBisector_L) ) ;
      mWrappedHalfedges.push_back( HalfedgeWrapper(lNIBisector_R) ) ;

      lNewNode_R->VBase::set_halfedge(lNIBisector_L) ;

      Halfedge_handle lXOBisector = lSeed->primary_bisector() ;
      Halfedge_handle lXIBisector = lXOBisector->opposite();

      lNOBisector_L->HBase::set_face(lXOBisector->HBase::face());
      lNIBisector_L->HBase::set_face(lOppBorder ->HBase::face());
      lNOBisector_R->HBase::set_face(lOppBorder ->HBase::face());
      lNIBisector_R->HBase::set_face(lXIBisector->HBase::face());

      lNIBisector_L->HBase::set_vertex(lNewNode_L);
      lNIBisector_R->HBase::set_vertex(lNewNode_R);

      lXOBisector  ->HBase::set_next(lNOBisector_L);
      lNOBisector_L->HBase::set_prev(lXOBisector);

      lXIBisector  ->HBase::set_prev(lNIBisector_R);
      lNIBisector_R->HBase::set_next(lXIBisector);

      lNIBisector_L->HBase::set_next(lNOBisector_R);
      lNOBisector_R->HBase::set_prev(lNIBisector_L);

      CGAL_SSBUILDER_TRACE
        (    "New Node L: N" << lNewNode_L->id() << " at " << lNewNode_L->point() << "\n"
          << "New Node R: N" << lNewNode_R->id() << " at " << lNewNode_R->point() << "\n"
          << "New Bisector OL:\nB" << lNOBisector_L->id()
          << "E"  << lNOBisector_L            ->defining_contour_edge()->id()
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
          << "] (In: Next: B" << lNIBisector_R->next()->id() << ")"
        ) ;

      UpdatePQ(lNewNode_L);
      UpdatePQ(lNewNode_R);
    }
  }
}

template<class Gt, class SS>
void Straight_skeleton_builder_2<Gt,SS>::HandleVertexEvent( VertexEvent& aEvent )
{
  Vertex_handle lLSeed = aEvent.left_seed() ;
  Vertex_handle lRSeed = aEvent.right_seed() ;

  if ( !IsProcessed(lLSeed) && !IsProcessed(lRSeed) )
  {
    Vertex_handle lNewNode_L, lNewNode_R ;
    boost::tie(lNewNode_L,lNewNode_R) = ConstructVertexEventNodes(aEvent);

    Halfedge_handle lNBisector_LO = mSS.SBase::edges_push_back ( Halfedge(mEdgeID++),Halfedge(mEdgeID++) );
    Halfedge_handle lNBisector_RO = mSS.SBase::edges_push_back ( Halfedge(mEdgeID++),Halfedge(mEdgeID++) );
    Halfedge_handle lNBisector_LI = lNBisector_LO->opposite();
    Halfedge_handle lNBisector_RI = lNBisector_RO->opposite();

    mWrappedHalfedges.push_back( HalfedgeWrapper(lNBisector_LO) ) ;
    mWrappedHalfedges.push_back( HalfedgeWrapper(lNBisector_RO) ) ;
    mWrappedHalfedges.push_back( HalfedgeWrapper(lNBisector_LI) ) ;
    mWrappedHalfedges.push_back( HalfedgeWrapper(lNBisector_RI) ) ;

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

    lSBisector_RO->HBase::set_prev(lNBisector_RO);
    lNBisector_RO->HBase::set_prev(lSBisector_RO);

    lNewNode_L->VBase::set_halfedge(lNBisector_LI);
    lNewNode_R->VBase::set_halfedge(lNBisector_RI);

    CGAL_SSBUILDER_TRACE
      (    "New Node L: N" << lNewNode_L->id() << " at " << lNewNode_L->point() << "\n"
        << "New Node R: N" << lNewNode_R->id() << " at " << lNewNode_R->point() << "\n"
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
        << ",E" << lNBisector_RI->opposite()->defining_contour_edge()->id() << "]"
      ) ;


    UpdatePQ(lNewNode_L);
    UpdatePQ(lNewNode_R);

  }
}

template<class Gt, class SS>
void Straight_skeleton_builder_2<Gt,SS>::Propagate()
{
  CGAL_SSBUILDER_TRACE("Propagating events...");

  while ( !mPQ.empty() )
  {
    EventPtr lEvent = PopEventFromPQ();
    CGAL_SSBUILDER_TRACE ("\nStep: " << mStepID << " Event: " << *lEvent ) ;
    CGAL_SSBUILDER_SHOW ( SS_IO_AUX::ScopedPointDrawing lDraw(lEvent->point(),CGAL::BLUE,"Event"); )

    switch ( lEvent->type() )
    {
      case Event::cEdgeEvent :  HandleEdgeEvent   ( dynamic_cast<EdgeEvent&>  (*lEvent) ) ;  break ;
      case Event::cSplitEvent:  HandleSplitEvent  ( dynamic_cast<SplitEvent&> (*lEvent) ) ;  break ;
      case Event::cVertexEvent: HandleVertexEvent ( dynamic_cast<VertexEvent&>(*lEvent) ) ;  break ;
    }

    ++ mStepID ;
  }

}

template<class Gt, class SS>
void Straight_skeleton_builder_2<Gt,SS>::MergeSplitNodes ( Vertex_handle_pair aSplitNodes )
{
  Vertex_handle lLNode, lRNode ;
  tie(lLNode,lRNode)=aSplitNodes;
  Halfedge_handle lIBisector = lRNode->primary_bisector()->opposite();

  Exclude(lRNode);

  lIBisector->HBase::set_vertex(lLNode);

  CGAL_SSBUILDER_TRACE("SplitNodes: N" << lLNode->id() << " and N" << lRNode->id() << " merged.\n"
                       << 'B' << lIBisector->id() << " now linked to N" << lLNode->id()
                       << ". N" << lRNode->id() << " excluded."
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
/*
template<class Gt, class SS>
template<class InputPointIterator>
Straight_skeleton_builder_2<Gt,SS>&
Straight_skeleton_builder_2<Gt,SS>::enter_contour ( InputPointIterator aBegin, InputPointIterator aEnd  )
{
  CGAL_SSBUILDER_TRACE("Inserting Connected Component of the Boundary....");

  if ( std::distance(aBegin,aEnd) >= 3 )
  {
    Halfedge_handle lFirstCCWBorder ;
    Halfedge_handle lPrevCCWBorder ;
    Halfedge_handle lNextCWBorder ;
    Vertex_handle   lFirstVertex ;
    Vertex_handle   lPrevVertex ;
    int             lFirstBorderIdx ;
    int             lPrevBorderIdx ;

    InputPointIterator lCurr = aBegin ;
    while ( lCurr != aEnd )
    {
      Halfedge_handle lCCWBorder = mSS.SBase::edges_push_back ( Halfedge(mEdgeID++),Halfedge(mEdgeID++)  );
      Halfedge_handle lCWBorder = lCCWBorder->opposite();

      Vertex_handle lVertex = mSS.SBase::vertices_push_back( Vertex(mVertexID++,*lCurr) ) ;
      CGAL_SSBUILDER_TRACE("Vertex: V" << lVertex->id() << " at " << lVertex->point() );
      mWrappedVertices.push_back( VertexWrapper(lVertex) ) ;

      mWrappedHalfedges.push_back( HalfedgeWrapper(lCCWBorder) ) ;
      mWrappedHalfedges.push_back( HalfedgeWrapper(lCWBorder ) ) ;
      Face_handle lFace = mSS.SBase::faces_push_back( Face() ) ;

      int lBorderIdx = (int)mBorderHalfedges.size() ;
      mBorderHalfedges.push_back(lCCWBorder);

      lCCWBorder->HBase::set_face    (lFace);
      lFace     ->FBase::set_halfedge(lCCWBorder);

      lVertex   ->VBase::set_halfedge(lCCWBorder);
      lCCWBorder->HBase::set_vertex  (lVertex);

      if ( lCurr == aBegin )
      {
        lFirstVertex    = lVertex ;
        lFirstCCWBorder = lCCWBorder ;
        lFirstBorderIdx = lBorderIdx ;
      }
      else
      {
        SetPrevInLAV(lVertex    ,lPrevVertex);
        SetNextInLAV(lPrevVertex,lVertex    );

        lCWBorder->HBase::set_vertex(lPrevVertex);

        lCCWBorder    ->HBase::set_prev(lPrevCCWBorder);
        lPrevCCWBorder->HBase::set_next(lCCWBorder);

        lNextCWBorder->HBase::set_prev(lCWBorder);
        lCWBorder    ->HBase::set_next(lNextCWBorder);

        SetPrevInCCB(lCCWBorder    ,lPrevBorderIdx);
        SetNextInCCB(lPrevCCWBorder,lBorderIdx);

        CGAL_SSBUILDER_TRACE("CCW Border: E" << lCCWBorder->id() << ' ' << lPrevVertex->point() << " -> " << lVertex    ->point());
        CGAL_SSBUILDER_TRACE("CW  Border: E" << lCWBorder ->id() << ' ' << lVertex    ->point() << " -> " << lPrevVertex->point() );

        CGAL_SSBUILDER_SHOW_AUX
        (
          SS_IO_AUX::ScopedSegmentDrawing draw_(lPrevVertex->point(),lVertex->point(), CGAL::RED, "Border" ) ;
          draw_.Release();
        )
      }

      ++ lCurr ;

      lPrevVertex    = lVertex ;
      lPrevCCWBorder = lCCWBorder ;
      lNextCWBorder  = lCWBorder ;
      lPrevBorderIdx = lBorderIdx ;
    }

    SetPrevInLAV(lFirstVertex,lPrevVertex );
    SetNextInLAV(lPrevVertex ,lFirstVertex);

    lFirstCCWBorder->opposite()->HBase::set_vertex(lPrevVertex);

    CGAL_SSBUILDER_SHOW_AUX
    (
      SS_IO_AUX::ScopedSegmentDrawing draw_(lPrevVertex->point(),lFirstVertex->point(), CGAL::RED, "Border" ) ;
      draw_.Release();
    )

    lFirstCCWBorder->HBase::set_prev(lPrevCCWBorder);
    lPrevCCWBorder ->HBase::set_next(lFirstCCWBorder);

    lPrevCCWBorder ->opposite()->HBase::set_prev(lFirstCCWBorder->opposite());
    lFirstCCWBorder->opposite()->HBase::set_next(lPrevCCWBorder ->opposite());

    SetPrevInCCB(lFirstCCWBorder, lPrevBorderIdx);
    SetNextInCCB(lPrevCCWBorder , lFirstBorderIdx);

    CGAL_SSBUILDER_TRACE("CCW Border: E" << lFirstCCWBorder->id()
                        << ' ' << lPrevVertex ->point() << " -> " << lFirstVertex->point() << '\n'
                        << "CW  Border: E" << lFirstCCWBorder->opposite()->id()
                        << ' ' << lFirstVertex->point() << " -> " << lPrevVertex ->point()
                        );
  }

  for ( Vertex_iterator v = mSS.SBase::vertices_begin(); v != mSS.SBase::vertices_end(); ++ v )
  {
    bool lIsReflex = !Left_turn( GetPrevInLAV(v)->point(),v->point(),GetNextInLAV(v)->point() );
    if ( lIsReflex )
    {
      SetIsReflex(v);
      CGAL_SSBUILDER_TRACE("Reflex vertex: N" << v->id() );
      mReflexVertices.push_back(v);
    }
  }

  return *this ;
}
*/
CGAL_END_NAMESPACE

#endif // CGAL_STRAIGHT_SKELETON_BUILDER_2_C //
// EOF //
