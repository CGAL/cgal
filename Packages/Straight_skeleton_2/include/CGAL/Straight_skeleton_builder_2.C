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

#ifdef CGAL_STRAIGHT_SKELETON_ENABLE_TRACE
#  include<string>
#  include<iostream>
#  include<sstream>
#  define CGAL_SSBUILDER_TRACE(m) \
     { \
       std::ostringstream ss ; \
       ss << m << std::ends ; \
       std::string s = ss.str(); \
       Straight_skeleton_external_trace(s); \
     }  
#else
#  define CGAL_SSBUILDER_TRACE(m)
#endif

#ifdef CGAL_STRAIGHT_SKELETON_ENABLE_SHOW_AUX
#  define CGAL_SSBUILDER_SHOW_AUX(code) code
#else
#  define CGAL_SSBUILDER_SHOW_AUX(code)
#endif

#ifdef CGAL_STRAIGHT_SKELETON_ENABLE_SHOW
#  define CGAL_SSBUILDER_SHOW(code) code
#else
#  define CGAL_SSBUILDER_SHOW(code)
#endif

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
  Halfedge_handle lAL = GetLeftDefiningBorder (aA);  
  Halfedge_handle lAR = GetRightDefiningBorder(aA);  
  Halfedge_handle lBL = GetLeftDefiningBorder (aB);  
  Halfedge_handle lBR = GetRightDefiningBorder(aB);  
  
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
  
  OptionalEventData lEventData = mTraits.compute_event( lBorderA->segment()
                                                       ,lBorderB->segment()
                                                       ,lBorderC->segment()
                                                      ) ;
                                                      
  if ( !!lEventData )
    rResult = EventPtr( new EdgeEvent( lBorderA
                                      ,lBorderB
                                      ,lBorderC
                                      ,lEventData->first
                                      ,lEventData->second
                                      ,aLNode
                                      ,aRNode
                                     ) 
                      ) ;
    
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
  if (    Halfedge_const_handle(aOppositeBorder) != aReflexLBorder 
       && Halfedge_const_handle(aOppositeBorder) != aReflexRBorder 
     )
  {
    OptionalEventData lEventData = mTraits.compute_event( aReflexLBorder ->segment()
                                                         ,aReflexRBorder ->segment()
                                                         ,aOppositeBorder->segment()
                                                        ) ;
    if ( !!lEventData )                                                        
    {
      Halfedge_handle lPrevOppBorder = GetPrevInCCB(aOppositeBorder);
      Halfedge_handle lNextOppBorder = GetNextInCCB(aOppositeBorder);

      if ( mTraits.is_event_inside_bounded_offset_zone( aReflexLBorder ->segment()
                                                       ,aReflexRBorder ->segment()
                                                       ,aOppositeBorder->segment()
                                                       ,aOppositeBorder->segment()
                                                       ,lPrevOppBorder ->segment()
                                                       ,lNextOppBorder ->segment()
                                                      ) 
         )                                             
       {
         aCandidates.push_back( EventPtr( new SplitEvent( aReflexLBorder
                                                         ,aReflexRBorder
                                                         ,aOppositeBorder
                                                         ,lEventData->first
                                                         ,lEventData->second
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
  
  Halfedge_handle lLBorder = GetLeftDefiningBorder (aNode);
  Halfedge_handle lRBorder = GetRightDefiningBorder(aNode);
  
  CGAL_SSBUILDER_TRACE("Finding SplitEvent for N" << aNode->id() 
                       << " LBorder=" << lLBorder->id() << " RBorder=" << lRBorder->id()
                      );
  
  EventPtr_Vector lCandidates ;
  for_each( mBorderHalfedges.begin()
           ,mBorderHalfedges.end  ()
           ,boost::bind(CollectSplitEvent,
                        this
                        ,aNode
                        ,lLBorder
                        ,lRBorder
                        ,_1
                        ,boost::ref(lCandidates)
                       )
          ) ; 
  
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
      << "Prev: " << lPrev->id() << " Next:" << lNext->id()
    ) ; 

  EventPtr lLEdgeEvent = FindEdgeEvent( lPrev , aNode ) ;  
  EventPtr lREdgeEvent = FindEdgeEvent( aNode , lNext ) ;
  EventPtr lSplitEvent = IsReflex(aNode) ? FindSplitEvent(aNode) : EventPtr() ;
  
  if ( lLEdgeEvent )
    aEvents.push_back(lLEdgeEvent) ;
  if ( lREdgeEvent )
    aEvents.push_back(lREdgeEvent) ;
  if ( lSplitEvent )
    aEvents.push_back(lSplitEvent) ;
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
    
    CGAL_SSBUILDER_SHOW ( SS_IO_AUX::ScopedPointDrawing lDraw(lEvent->point(),2,"Event"); )
    
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
  Halfedge_handle lOA = GetPrimaryBisector(aA) ;
  Halfedge_handle lOB = GetPrimaryBisector(aB) ;
  Halfedge_handle lIA = lOA->opposite();
  Halfedge_handle lIB = lOB->opposite();

  CGAL_SSBUILDER_TRACE ( "Handling simultaneous EdgeEvent " 
                         << "A: "  << aA ->id() << '\n'
                         << "B: "  << aB ->id() << '\n'
                         << "OA: " << lOA->id() << '\n'
                         << "IA: " << lIA->id() << '\n'
                         << "OB: " << lOB->id() << '\n'
                         << "IB: " << lIB->id() 
                       ) ;

  CGAL_SSBUILDER_TRACE ( 'N' << aA->id() << " processed\nN" << aB->id() << " processed" ) ;
  SetIsProcessed(aA) ;
  SetIsProcessed(aB) ;

  Halfedge_handle lOA_Prev = lOA->prev() ;
  Halfedge_handle lIA_Next = lIA->next() ;

  Halfedge_handle lOB_Prev = lOB->prev() ;
  Halfedge_handle lIB_Next = lIB->next() ;
  
  CGAL_SSBUILDER_TRACE (    "OA_Prev: " << lOA_Prev->id() << '\n'
                         << "IA_Next: " << lIA_Next->id() << '\n'
                         << "OB_Prev: " << lOB_Prev->id() << '\n'
                         << "IB_Next: " << lIB_Next->id()
                      ) ;
                   
  lOB     ->Base::set_next( lIA_Next );
  lIA_Next->Base::set_prev( lOB      );
  lIB     ->Base::set_prev( lOA_Prev );
  lOA_Prev->Base::set_next( lIB      );

  lOB->Base::set_vertex(aA);
  
  lOB->Base::set_segment( Segment_2(aB->point(),aA->point()) ) ;
  lIB->Base::set_segment( Segment_2(aA->point(),aB->point()) ) ;

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

  Halfedge_handle lOBisector_P = GetPrimaryBisector(lPrev) ;
  Halfedge_handle lOBisector_C = GetPrimaryBisector(aNode) ;
  Halfedge_handle lOBisector_N = GetPrimaryBisector(lNext) ;
  
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
    Halfedge_handle lOBisector = mSS.Base::edges_push_back (lOB, lIB);
    Halfedge_handle lIBisector = lOBisector->opposite();
    mWrappedHalfedges.push_back( HalfedgeWrapper(lOBisector) ) ;
    mWrappedHalfedges.push_back( HalfedgeWrapper(lIBisector) ) ;
    lOBisector->Base::set_face(v->halfedge()->face());
    lIBisector->Base::set_face(v->halfedge()->next()->face());
    lIBisector->Base::set_vertex(v);
    
    Halfedge_handle lIBorder = v->halfedge() ;
    Halfedge_handle lOBorder = v->halfedge()->next() ;
    lIBorder  ->Base::set_next(lOBisector);
    lOBisector->Base::set_prev(lIBorder);
    lOBorder  ->Base::set_prev(lIBisector);
    lIBisector->Base::set_next(lOBorder);
    CGAL_SSBUILDER_TRACE("Adding Contour Bisector:\n E" << lOBisector->id() << " (Out)\nE" << lIBisector->id() << " (In)" ) ;
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

  Vertex_handle lNode = mSS.Base::vertices_push_back( Vertex( mVertexID++
                                                             ,aEvent.point()
                                                             ,aEvent.time()
                                                            )
                                                    ) ;         
  mWrappedVertices.push_back( VertexWrapper(lNode) ) ;
  
  Halfedge_handle lLBisector = GetPrimaryBisector(lLSeed);
  Halfedge_handle lRBisector = GetPrimaryBisector(lRSeed);
  
  lNode->Base::set_halfedge(lLBisector);
  lLBisector->Base::set_vertex(lNode);
  lRBisector->Base::set_vertex(lNode);
  
  lLBisector            ->Base::set_segment( Segment_2(lLSeed->point(), lNode ->point()) ) ;
  lLBisector->opposite()->Base::set_segment( Segment_2(lNode ->point(), lLSeed->point()) ) ;
  lRBisector            ->Base::set_segment( Segment_2(lRSeed->point(), lNode ->point()) ) ;
  lRBisector->opposite()->Base::set_segment( Segment_2(lNode ->point(), lRSeed->point()) ) ;
  
  lLBisector->opposite()->Base::set_prev( lRBisector ) ;
  lRBisector            ->Base::set_next( lLBisector->opposite() ) ;

  CGAL_SSBUILDER_SHOW( DrawBisector(lLBisector); DrawBisector(lRBisector); ) ;
  
  CGAL_SSBUILDER_TRACE 
  (
       "LSeed: " << lLSeed->id() << " proccesed\n"  
    << "RSeed: " << lRSeed->id() << " proccesed"
  ) ;
  
  SetIsProcessed(lLSeed) ;
  SetIsProcessed(lRSeed) ;

  CGAL_SSBUILDER_TRACE ( 'N' << lNode->id() << " inserted into LAV." ) ;
  
  Vertex_handle lLPrev = GetPrevInLAV(lLSeed) ;
  Vertex_handle lRNext = GetNextInLAV(lRSeed) ;

  SetPrevInLAV(lNode , lLPrev ) ;
  SetNextInLAV(lLPrev, lNode  ) ;

  SetNextInLAV(lNode , lRNext ) ;
  SetPrevInLAV(lRNext, lNode  ) ;

  CGAL_SSBUILDER_TRACE 
  (
    "Updated LAV: " << lLPrev->id() << "->" << lNode->id() << "->" << lRNext->id() << std::endl
    << lLSeed->id() << " removed from LAV\n" 
    << lRSeed->id() << " removed from LAV" 
  );

  return lNode ;
}

template<class Gt, class SS>
typename Straight_skeleton_builder_2<Gt,SS>::Vertex_handle
Straight_skeleton_builder_2<Gt,SS>::LookupOnSLAV ( Halfedge_handle aBorder, Event const& aEvent )
{
  Vertex_handle rResult ;

  CGAL_SSBUILDER_TRACE ( "Looking up for E" << aBorder->id() << " on SLAV. P=" << aEvent.point() ) ;

  CGAL_SSBUILDER_SHOW_AUX ( SS_IO_AUX::ScopedSegmentDrawing draw_(aBorder->segment(),4,"OppBorder") ; ) 
   
  for ( Vertex_iterator vi = mSS.vertices_begin(); vi != mSS.vertices_end(); ++ vi )
  {
    Vertex_handle v = static_cast<Vertex_handle>(vi);
    
CGAL_SSBUILDER_TRACE (    "PrevInLAV exist: " << handle_assigned(GetPrevInLAV(v))
                      << " NextInLAV exist: " << handle_assigned(GetNextInLAV(v))
                      << " HasBorder: " << ( GetLeftDefiningBorder(v) == aBorder )
                      << " NotProcesed: " << !IsProcessed(v) 
                     ) ;
                       
    if (     handle_assigned(GetPrevInLAV(v))
         &&  handle_assigned(GetNextInLAV(v))
         &&  GetLeftDefiningBorder(v) == aBorder
         && !IsProcessed(v) 
       )
    {
      Vertex_handle lPrev = GetPrevInLAV(v);
      Halfedge_handle lPrevBorder = GetLeftDefiningBorder (lPrev);  
      Halfedge_handle lNextBorder = GetRightDefiningBorder(v);  
      if ( mTraits.is_event_inside_bounded_offset_zone
            ( aEvent.border_a()->segment() 
             ,aEvent.border_b()->segment()
             ,aEvent.border_c()->segment()
             ,aBorder          ->segment()
             ,lPrevBorder      ->segment()
             ,lNextBorder      ->segment()
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
    
    Vertex_handle lNodeA = mSS.Base::vertices_push_back( Vertex( mVertexID++
                                                                ,aEvent.point()
                                                                ,aEvent.time ()
                                                               ) 
                                                       ) ;
    Vertex_handle lNodeB = mSS.Base::vertices_push_back( Vertex( mVertexID++
                                                                ,aEvent.point()
                                                                ,aEvent.time ()
                                                               ) 
                                                       ) ;
    
    mWrappedVertices.push_back( VertexWrapper(lNodeA) ) ;
    mWrappedVertices.push_back( VertexWrapper(lNodeB) ) ;
    
    Halfedge_handle lXOutBisector = GetPrimaryBisector(aEvent.seed()) ;
    Halfedge_handle lXInBisector  = lXOutBisector->opposite();
    
    lNodeA->Base::set_halfedge(lXOutBisector);
    // lNodeB hafledge is set outside with the New In Bisector to the Right.
    
    lXOutBisector->Base::set_vertex(lNodeA);

    Vertex_handle lSeed = aEvent.seed() ;
    
    lXOutBisector->Base::set_segment( Segment_2(lSeed->point(),aEvent.point()) ) ;
    lXInBisector ->Base::set_segment( Segment_2(aEvent.point(),lSeed->point()) ) ;
    
    CGAL_SSBUILDER_SHOW( DrawBisector(lXOutBisector); ) ;
    
    CGAL_SSBUILDER_TRACE ( "Seed: " << lSeed->id() << " proccesed" ) ;
    
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
         "Updated LAV: " << lPrev->id() << "->" << lNodeA->id() << "->" << lOppR->id() << std::endl
      << "Updated LAV: " << lOppL->id() << "->" << lNodeB->id() << "->" << lNext->id() << std::endl
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
  if ( !IsProcessed(aEvent.left_seed()) || !IsProcessed(aEvent.right_seed()) ) 
  {
    Vertex_handle lNewNode = ConstructEdgeEventNode(aEvent);

    Halfedge_handle lNOBisector = mSS.Base::edges_push_back ( Halfedge(mEdgeID++)
                                                             ,Halfedge(mEdgeID++)
                                                            );
    Halfedge_handle lNIBisector = lNOBisector->opposite();
    
    mWrappedHalfedges.push_back( HalfedgeWrapper(lNOBisector) ) ;
    mWrappedHalfedges.push_back( HalfedgeWrapper(lNIBisector) ) ;
    
    Halfedge_handle lLOBisector = GetPrimaryBisector(aEvent.left_seed ()) ;
    Halfedge_handle lROBisector = GetPrimaryBisector(aEvent.right_seed());
    Halfedge_handle lLIBisector = lLOBisector->opposite();
    Halfedge_handle lRIBisector = lROBisector->opposite();

    lNOBisector->Base::set_face(lLOBisector->Base::face());
    lNIBisector->Base::set_face(lRIBisector->Base::face());
    lNIBisector->Base::set_vertex(lNewNode);
    
    lLOBisector->Base::set_next(lNOBisector);
    lNOBisector->Base::set_prev(lLOBisector);
    
    lRIBisector->Base::set_prev(lNIBisector);
    lNIBisector->Base::set_next(lRIBisector);
    
    CGAL_SSBUILDER_TRACE
      (     "LO: " << lLOBisector->id() << " LI: " << lLIBisector->id() 
        << " RO: " << lROBisector->id() << " RI: " << lRIBisector->id()  << '\n'
        << "New Node: " << lNewNode->id() << " at " << lNewNode->point() << '\n' 
        << "New Bisectors:\n" 
        << 'E' << lNOBisector->id() << " [" << lNOBisector->defining_border()->id()
        << ',' << lNOBisector->opposite()->defining_border()->id()
        << "] (Out)\n" 
        << 'E' << lNIBisector->id() << " [" << lNIBisector->defining_border()->id()
        << ',' << lNIBisector->opposite()->defining_border()->id()
        << "] (In)"
      ) ;
    
    UpdatePQ(lNewNode);
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
      
      Halfedge_handle lReflexLBorder = GetLeftDefiningBorder (lSeed);
      Halfedge_handle lReflexRBorder = GetRightDefiningBorder(lSeed);
      
      Halfedge_handle lNOBisector_L = mSS.Base::edges_push_back ( Halfedge(mEdgeID++)
                                                                 ,Halfedge(mEdgeID++)
                                                                 );
      Halfedge_handle lNOBisector_R = mSS.Base::edges_push_back ( Halfedge(mEdgeID++)
                                                                 ,Halfedge(mEdgeID++)
                                                                );
      Halfedge_handle lNIBisector_L = lNOBisector_L->opposite();
      Halfedge_handle lNIBisector_R = lNOBisector_R->opposite();
      
      mWrappedHalfedges.push_back( HalfedgeWrapper(lNOBisector_L) ) ;
      mWrappedHalfedges.push_back( HalfedgeWrapper(lNOBisector_R) ) ;
      mWrappedHalfedges.push_back( HalfedgeWrapper(lNIBisector_L) ) ;
      mWrappedHalfedges.push_back( HalfedgeWrapper(lNIBisector_R) ) ;
      
      lNewNode_R->Base::set_halfedge(lNIBisector_L) ;
      
      Halfedge_handle lXOBisector = GetPrimaryBisector(lSeed) ;
      Halfedge_handle lXIBisector = lXOBisector->opposite();
      
      lNOBisector_L->Base::set_face(lXOBisector->Base::face());
      lNIBisector_L->Base::set_face(lOppBorder ->Base::face());
      lNOBisector_R->Base::set_face(lOppBorder ->Base::face());
      lNIBisector_R->Base::set_face(lXIBisector->Base::face());
      
      lNIBisector_L->Base::set_vertex(lNewNode_L);
      lNIBisector_R->Base::set_vertex(lNewNode_R);
      
      lXOBisector  ->Base::set_next(lNOBisector_L);
      lNOBisector_L->Base::set_prev(lXOBisector);
      
      lXIBisector  ->Base::set_prev(lNIBisector_R);
      lNIBisector_R->Base::set_next(lXIBisector);
      
      lNIBisector_L->Base::set_next(lNOBisector_R);
      lNOBisector_R->Base::set_prev(lNIBisector_L);
      
      CGAL_SSBUILDER_TRACE
        (    "New Node L: " << lNewNode_L->id() << " at " << lNewNode_L->point() << '\n' 
          << "New Node R: " << lNewNode_R->id() << " at " << lNewNode_R->point() << '\n' 
          << "New Bisector OL:\n" 
          << lNOBisector_L->id() 
          << '[' << lNOBisector_L            ->defining_border()->id() << ','
          <<        lNOBisector_L->opposite()->defining_border()->id() 
          << ']'
          << " (Out)\n" 
          << "New Bisector IL:\n" 
          << lNIBisector_L->id() 
          << '[' << lNIBisector_L            ->defining_border()->id() << ','
          <<        lNIBisector_L->opposite()->defining_border()->id() 
          << ']'
          << " (In)\n" 
          << "New Bisector OR:\n" 
          << lNOBisector_R->id()
          << '[' << lNOBisector_R            ->defining_border()->id() << ','
          <<        lNOBisector_R->opposite()->defining_border()->id() 
          << ']'
          << " (Out)\n"
          << "New Bisector IR:\n" 
          << lNIBisector_R->id() 
          << '[' << lNIBisector_R            ->defining_border()->id() << ','
          <<        lNIBisector_R->opposite()->defining_border()->id() 
          << ']'
          << " (In)"
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
    CGAL_SSBUILDER_SHOW_TRAITS ( SS_IO_AUX::ScopedPointDrawing lDraw(lEvent->point(),2,"Event"); )

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
  Halfedge_handle lIBisector = GetPrimaryBisector(lRNode)->opposite();
 
  Exclude(lRNode);
  
  lIBisector->Base::set_vertex(lLNode);
  
  CGAL_SSBUILDER_TRACE("SplitNodes " << lLNode->id() << " and " << lRNode->id() << " merged.\n"
                       << 'E' << lIBisector->id() << " now linked to N" << lLNode->id()
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
  
  if ( aBegin != aEnd )
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
      Vertex_handle lVertex = mSS.Base::vertices_push_back( Vertex(mVertexID++,*lCurr) ) ; 
      
      mWrappedVertices.push_back( VertexWrapper(lVertex) ) ;
      
      Halfedge_handle lCCWBorder = mSS.Base::edges_push_back ( Halfedge(mEdgeID++)
                                                              ,Halfedge(mEdgeID++)
                                                             );
      Halfedge_handle lCWBorder = lCCWBorder->opposite();
      
      mWrappedHalfedges.push_back( HalfedgeWrapper(lCCWBorder) ) ;
      mWrappedHalfedges.push_back( HalfedgeWrapper(lCWBorder ) ) ;
      
      Face_handle lFace = mSS.Base::faces_push_back( Face() ) ;
      
      int lBorderIdx = (int)mBorderHalfedges.size() ;
      mBorderHalfedges.push_back(lCCWBorder);
      
      lCCWBorder->Base::set_face    (lFace); 
      lFace     ->Base::set_halfedge(lCCWBorder);
      
      lVertex   ->Base::set_halfedge(lCCWBorder);
      lCCWBorder->Base::set_vertex  (lVertex);
      
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
        
        lCWBorder->Base::set_vertex(lPrevVertex);
                
        lCCWBorder    ->Base::set_prev(lPrevCCWBorder);
        lPrevCCWBorder->Base::set_next(lCCWBorder);

        lNextCWBorder->Base::set_prev(lCWBorder);
        lCWBorder    ->Base::set_next(lNextCWBorder);

        SetPrevInCCB(lCCWBorder    ,lPrevBorderIdx);
        SetNextInCCB(lPrevCCWBorder,lBorderIdx);
        
        Segment_2 lCCWSeg(lPrevVertex->point(),lVertex    ->point());
        Segment_2 lCWSeg (lVertex    ->point(),lPrevVertex->point());
        
        lCCWBorder->Base::set_segment(lCCWSeg);
        lCWBorder ->Base::set_segment(lCWSeg );
        
        CGAL_SSBUILDER_SHOW_AUX
        ( 
          SS_IO_AUX::ScopedSegmentDrawing draw_(lCCWSeg, 1, "Border" ) ; 
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
    
    lFirstCCWBorder->opposite()->Base::set_vertex(lPrevVertex);

    Segment_2 lCCWSeg(lPrevVertex ->point(), lFirstVertex->point()) ;
    Segment_2 lCWSeg (lFirstVertex->point(), lPrevVertex ->point()) ; 
    
    lFirstCCWBorder                  ->Base::set_segment( lCCWSeg ) ;
    lFirstCCWBorder->Base::opposite()->Base::set_segment( lCWSeg ) ;
    
    CGAL_SSBUILDER_SHOW_AUX
    ( 
      SS_IO_AUX::ScopedSegmentDrawing draw_(lCCWSeg, 1, "Border" ) ; 
      draw_.Release(); 
    ) 
    
    lFirstCCWBorder->Base::set_prev(lPrevCCWBorder);
    lPrevCCWBorder ->Base::set_next(lFirstCCWBorder);

    lPrevCCWBorder ->opposite()->Base::set_prev(lFirstCCWBorder->opposite());
    lFirstCCWBorder->opposite()->Base::set_next(lPrevCCWBorder ->opposite());
    
    SetPrevInCCB(lFirstCCWBorder, lPrevBorderIdx);
    SetNextInCCB(lPrevCCWBorder , lFirstBorderIdx);
  }

  for ( Vertex_iterator v = mSS.Base::vertices_begin(); v != mSS.Base::vertices_end(); ++ v )
  {
    typename Traits::Rep::Left_turn_2 left_turn ;
                  
    bool lIsReflex = !left_turn( GetPrevInLAV(v)->point()
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
