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

#include <boost/tuple/tuple.hpp>
#include <boost/bind.hpp>


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
 ,mEventCompare(aTraits)
 ,mVertexID(0)
 ,mEdgeID(0)
 ,mEventID(0)
 ,mStepID(0)
 ,Left_turn                     (aTraits.template get<typename Traits::Left_turn_2>())
 ,Exist_event                   (aTraits.template get<typename Traits::Exist_event>())
 ,Is_event_inside_offset_zone   (aTraits.template get<typename Traits::Is_event_inside_offset_zone>())
 ,Compare_event_times           (aTraits.template get<typename Traits::Compare_event_times>())
 ,Compare_event_distance_to_seed(aTraits.template get<typename Traits::Compare_event_distance_to_seed>())
 ,Construct_event               (aTraits.template get<typename Traits::Construct_event>())
{
}

//
// This method returns the 3 distinct defining borders of vertices aA and aB
// (As long as the vertices are proceesed in the right order there is 1 common defining border,
//  so there are 3 distinct borders given these 2 vertices)
//
template<class Gt, class SS>
typename Straight_skeleton_builder_2<Gt,SS>::BorderTriple
Straight_skeleton_builder_2<Gt,SS>::GetDefiningBorders( Vertex_handle aA
                                                       ,Vertex_handle aB
                                                      )
{
  Halfedge_handle lAL = GetDefiningBorder0(aA);
  Halfedge_handle lAR = GetDefiningBorder2(aA);
  Halfedge_handle lBL = GetDefiningBorder0(aB);
  Halfedge_handle lBR = GetDefiningBorder2(aB);

  return BorderTriple(lAL, lAR, ( lAL == lBL || lAR == lBL ) ? lBR : lBL ) ;
}

template<class Gt, class SS>
typename Straight_skeleton_builder_2<Gt,SS>::EventPtr
Straight_skeleton_builder_2<Gt,SS>::FindEdgeEvent( Vertex_handle aLNode
                                                  ,Vertex_handle aRNode
                                                 )
{
  EventPtr rResult ;

  Halfedge_handle lBorderA, lBorderB, lBorderC ;
  boost::tie(lBorderA,lBorderB,lBorderC) = GetDefiningBorders(aLNode,aRNode);

  if ( lBorderA != lBorderB && lBorderB != lBorderC )
  {
    if ( ExistEvent(lBorderA,lBorderB,lBorderC)  )
      rResult = EventPtr( new EdgeEvent( lBorderA, lBorderB, lBorderC, aLNode, aRNode ) ) ;
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
  CGAL_SSBUILDER_TRACE("Computing potential split event between E" << aReflexLBorder->id() << ",E" <<aReflexRBorder->id() << " and E" << aOppositeBorder->id() );
  
  if (    Halfedge_const_handle(aOppositeBorder) != aReflexLBorder
       && Halfedge_const_handle(aOppositeBorder) != aReflexRBorder
     )
  {
    if ( ExistEvent(aReflexLBorder,aReflexRBorder,aOppositeBorder) )
    {
      Halfedge_handle lPrevOppBorder = GetPrevInCCB(aOppositeBorder);
      Halfedge_handle lNextOppBorder = GetNextInCCB(aOppositeBorder);

      if ( IsEventInsideOffsetZone( aReflexLBorder 
                                  , aReflexRBorder 
                                  , aOppositeBorder
                                  , lPrevOppBorder 
                                  , lNextOppBorder 
                                  )
         )
       {
         aCandidates.push_back( EventPtr( new SplitEvent( aReflexLBorder
                                                         ,aReflexRBorder
                                                         ,aOppositeBorder
                                                         ,aNode
                                                         ,aOppositeBorder
                                                        )
                                        )
                              ) ;
       }
    }
  }
}


template<class Gt, class SS>
typename Straight_skeleton_builder_2<Gt,SS>::EventPtr
Straight_skeleton_builder_2<Gt,SS>::FindSplitEvent( Vertex_handle aNode )
{
  EventPtr rEvent ;

  Halfedge_handle lLBorder = GetDefiningBorder0(aNode);
  Halfedge_handle lRBorder = GetDefiningBorder2(aNode);

  CGAL_SSBUILDER_TRACE("Finding SplitEvent for N" << aNode->id()
                       << " LBorder: E" << lLBorder->id() << " RBorder: E" << lRBorder->id()
                      );

  EventPtr_Vector lCandidates ;
  
  for( typename HalfedgeCollection::const_iterator it = mBorderHalfedges.begin()
       , eit = mBorderHalfedges.end()
       ; it != eit
       ; ++ it
     )
    CollectSplitEvent(aNode, lLBorder, lRBorder, *it, lCandidates ) ;
  
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
void Straight_skeleton_builder_2<Gt,SS>
  ::CollectNewEvents( Vertex_handle aNode, EventPtr_Vector& aEvents )
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
Straight_skeleton_builder_2<Gt,SS>::ChooseBestNewEvent( Vertex_handle    aNode
                                                       ,EventPtr_Vector& aEvents
                                                      )
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
void Straight_skeleton_builder_2<Gt,SS>::HandleSimultaneousEdgeEvent( Vertex_handle aA
                                                                    , Vertex_handle aB
                                                                    )
{
  Halfedge_handle lOA = aA->primary_bisector() ;
  Halfedge_handle lOB = aB->primary_bisector() ;
  Halfedge_handle lIA = lOA->opposite();
  Halfedge_handle lIB = lOB->opposite();

  CGAL_SSBUILDER_TRACE ( "Handling simultaneous EdgeEvent "
                         << "A: N"  << aA ->id() << '\n'
                         << "B: N"  << aB ->id() << '\n'
                         << "OA: B" << lOA->id() << '\n'
                         << "IA: B" << lIA->id() << '\n'
                         << "OB: B" << lOB->id() << '\n'
                         << "IB: B" << lIB->id()
                       ) ;

  CGAL_SSBUILDER_TRACE ( 'N' << aA->id() << " processed\nN" << aB->id() << " processed" ) ;
  SetIsProcessed(aA) ;
  SetIsProcessed(aB) ;

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

  Exclude(lOA);
  Exclude(lIA);
}

template<class Gt, class SS>
bool Straight_skeleton_builder_2<Gt,SS>::AreBisectorsCoincident ( Halfedge_const_handle aA
                                                                 ,Halfedge_const_handle aB
                                                                ) const
{
  Halfedge_const_handle lA_LBorder = aA->defining_border();
  Halfedge_const_handle lA_RBorder = aA->opposite()->defining_border();
  Halfedge_const_handle lB_LBorder = aB->defining_border();
  Halfedge_const_handle lB_RBorder = aB->opposite()->defining_border();

  return    ( lA_LBorder == lB_LBorder && lA_RBorder == lB_RBorder )
         || ( lA_LBorder == lB_RBorder && lA_RBorder == lB_LBorder ) ;
}

template<class Gt, class SS>
void Straight_skeleton_builder_2<Gt,SS>::UpdatePQ( Vertex_handle aNode )
{
  Vertex_handle lPrev = GetPrevInLAV(aNode) ;
  Vertex_handle lNext = GetNextInLAV(aNode) ;

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
    CGAL_SSBUILDER_TRACE("Adding Contour Bisector at N:" << v->id() << "\n B" << lOBisector->id() << " (Out)\n B" << lIBisector->id() << " (In)" ) ;
  }
}

template<class Gt, class SS>
void Straight_skeleton_builder_2<Gt,SS>::InitPhase()
{
  CreateContourBisectors();
  CreateInitialEvents();
}

template<class Gt, class SS>
typename Straight_skeleton_builder_2<Gt,SS>::Vertex_handle
Straight_skeleton_builder_2<Gt,SS>::ConstructEdgeEventNode( EdgeEvent& aEvent )
{
  CGAL_SSBUILDER_TRACE ( "Creating EdgeEvent Node" ) ;

  Vertex_handle lLSeed = aEvent.left_seed () ;
  Vertex_handle lRSeed = aEvent.right_seed() ;

  Point_2 lP ; FT lTime ;
  boost::tie(lP,lTime) = ConstructEventPointAndTime(aEvent);  
  
  Vertex_handle lNewNode = mSS.SBase::vertices_push_back( Vertex( mVertexID++, lP, lTime) ) ;
  
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

  CGAL_SSBUILDER_SHOW_AUX ( SS_IO_AUX::ScopedSegmentDrawing draw_(aBorder->vertex()->point(),aBorder->opposite()->vertex()->point(),CGAL::YELLOW,"OppBorder") ; )

  for ( Vertex_iterator vi = mSS.vertices_begin(); vi != mSS.vertices_end(); ++ vi )
  {
    Vertex_handle v = static_cast<Vertex_handle>(vi);

    CGAL_SSBUILDER_TRACE (    "PrevInLAV exist: " << handle_assigned(GetPrevInLAV(v))
                          << " NextInLAV exist: " << handle_assigned(GetNextInLAV(v))
                          << " HasBorder: " << ( GetDefiningBorder0(v) == aBorder )
                          << " NotProcesed: " << !IsProcessed(v)
                        ) ;

    if (     handle_assigned(GetPrevInLAV(v))
         &&  handle_assigned(GetNextInLAV(v))
         &&  GetDefiningBorder0(v) == aBorder
         && !IsProcessed(v)
       )
    {
      Vertex_handle lPrev = GetPrevInLAV(v);
      Halfedge_handle lPrevBorder = GetDefiningBorder0(lPrev);
      Halfedge_handle lNextBorder = GetDefiningBorder2(v);
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

    Point_2 lP ; FT lTime ;
    boost::tie(lP,lTime) = ConstructEventPointAndTime(aEvent);
    
    Vertex_handle lNodeA = mSS.SBase::vertices_push_back( Vertex( mVertexID++, lP, lTime ) ) ;
    Vertex_handle lNodeB = mSS.SBase::vertices_push_back( Vertex( mVertexID++, lP, lTime ) ) ;

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
      
      Halfedge_handle lNOBisector = mSS.SBase::edges_push_back ( Halfedge(mEdgeID)
                                                                ,Halfedge(mEdgeID+1)
                                                               );
                                                               
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
        ( "New Bisectors:\nB" << lNOBisector->id() << " [E" << lNOBisector->defining_border()->id()
          << ",E" << lNOBisector->opposite()->defining_border()->id()
          << "] (Out: Prev: B" << lNOBisector->prev()->id() << ")\nB"
          << lNIBisector->id() << " [E" << lNIBisector->defining_border()->id()
          << ",E" << lNIBisector->opposite()->defining_border()->id()
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

      Halfedge_handle lReflexLBorder = GetDefiningBorder0(lSeed);
      Halfedge_handle lReflexRBorder = GetDefiningBorder2(lSeed);

      Halfedge_handle lNOBisector_L = mSS.SBase::edges_push_back ( Halfedge(mEdgeID++)
                                                                  ,Halfedge(mEdgeID++)
                                                                  );
      Halfedge_handle lNOBisector_R = mSS.SBase::edges_push_back ( Halfedge(mEdgeID++)
                                                                  ,Halfedge(mEdgeID++)

                                                                 );
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
          << "E"  << lNOBisector_L            ->defining_border()->id() 
          << ",E" << lNOBisector_L->opposite()->defining_border()->id() << "]"
          << " (Out: Prev: B" << lNOBisector_L->prev()->id() << ")\n" 
          << "New Bisector IL:\nB" << lNIBisector_L->id() 
          << "[E" << lNIBisector_L            ->defining_border()->id()
          << ",E" << lNIBisector_L->opposite()->defining_border()->id() << "]"
          << " (In: Next: B" << lNIBisector_L->next()->id() << ")\n" 
          << "New Bisector OR:\nB" << lNOBisector_R->id()
          << "[E" << lNOBisector_R            ->defining_border()->id() 
          << ",E" << lNOBisector_R->opposite()->defining_border()->id() << "]"
          << " (Out: Prev: B" << lNOBisector_R->prev()->id() << ")\n" 
          << "New Bisector IR:\nB" << lNIBisector_R->id() 
          << "[E" << lNIBisector_R            ->defining_border()->id()
          << ",E" << lNIBisector_R->opposite()->defining_border()->id() 
          << "] (In: Next: B" << lNIBisector_R->next()->id() << ")" 
        ) ;

      UpdatePQ(lNewNode_L);
      UpdatePQ(lNewNode_R);
    }
  }
}

template<class Gt, class SS>
void Straight_skeleton_builder_2<Gt,SS>::HandleMultiSplitEvent( SplitEvent* aEvent )
{
  do
  {
    HandleSplitEvent(*aEvent);
    aEvent = aEvent->next();
  }
  while ( aEvent ) ;
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
      case Event::cEdgeEvent :
             HandleEdgeEvent( static_cast<EdgeEvent&> (*lEvent) ) ;
             break ;

      case Event::cSplitEvent:
             HandleMultiSplitEvent ( static_cast<SplitEvent*>(&(*lEvent)) ) ;
             break ;
    }

    ++ mStepID ;
  }

}

template<class Gt, class SS>
void Straight_skeleton_builder_2<Gt,SS>::MergeSplitNodes ( Vertex_handle_pair aSplitNodes )
{
  Vertex_handle lLNode, lRNode ;
  boost::tie(lLNode,lRNode)=aSplitNodes;
  Halfedge_handle lIBisector = lRNode->primary_bisector()->opposite();

  Exclude(lRNode);

  lIBisector->HBase::set_vertex(lLNode);

  CGAL_SSBUILDER_TRACE("SplitNodes: N" << lLNode->id() << " and N" << lRNode->id() << " merged.\n"
                       << 'B' << lIBisector->id() << " now linked to N" << lLNode->id()
                       << ". N" << lRNode->id() << " excluded."
                       ); 

  //mSS.Base::vertices_erase(lRNode);
}

template<class Gt, class SS>
void Straight_skeleton_builder_2<Gt,SS>::FinishUp()
{
  std::for_each( mSplitNodes.begin()
                ,mSplitNodes.end  ()
                ,boost::bind(&Straight_skeleton_builder_2<Gt,SS>::MergeSplitNodes,this,_1)
               ) ;

/*
  mSS.Base::edges_erase( std::remove_if( mSS.Base::halfedges_begin()
                                        ,mSS.Base::halfedges_end  ()
                                        //,boost::bind(&Halfedge::is_excluded,_1)
                                        ,is_excluded()
                                       )
                         ,mSS.Base::halfedges_end()
                        ) ;
*/
}

template<class Gt, class SS>
void Straight_skeleton_builder_2<Gt,SS>::Run()
{
  InitPhase();
  Propagate();
  FinishUp ();
}

template<class Gt, class SS>
typename Straight_skeleton_builder_2<Gt,SS>::Ssds Straight_skeleton_builder_2<Gt,SS>::proceed()
{
  Run();
  return mSS ;
}

template<class Gt, class SS>
template<class InputPointIterator>
Straight_skeleton_builder_2<Gt,SS>&
Straight_skeleton_builder_2<Gt,SS>::insert_CCB (  InputPointIterator aBegin
                                                , InputPointIterator aEnd
                                               )
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
      Halfedge_handle lCCWBorder = mSS.SBase::edges_push_back ( Halfedge(mEdgeID++)
                                                               ,Halfedge(mEdgeID++)
                                                             );
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

        CGAL_SSBUILDER_TRACE("CCW Border: E" << lCCWBorder->id() << lPrevVertex->point() << " -> " << lVertex->point());
        CGAL_SSBUILDER_TRACE("CW  Border: E" << lCWBorder->id() << lVertex->point() << " -> " << lPrevVertex->point() );

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
  }

  for ( Vertex_iterator v = mSS.SBase::vertices_begin(); v != mSS.SBase::vertices_end(); ++ v )
  {
    bool lIsReflex = !Left_turn( GetPrevInLAV(v)->point()
                                ,v              ->point()
                                ,GetNextInLAV(v)->point()
                               );
    if ( lIsReflex )
    {
      SetIsReflex(v);
      CGAL_SSBUILDER_TRACE("Reflex vertex: N" << v->id() );
      mReflexVertices.push_back(v);
    }
  }

  return *this ;
}

CGAL_END_NAMESPACE

#endif // CGAL_STRAIGHT_SKELETON_BUILDER_2_C //
// EOF //
