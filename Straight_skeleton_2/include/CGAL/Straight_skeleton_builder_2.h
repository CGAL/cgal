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
#ifndef CGAL_STRAIGHT_SKELETON_BUILDER_2_H
#define CGAL_STRAIGHT_SKELETON_BUILDER_2_H 1

#include <list>
#include <queue>

#include <boost/tuple/tuple.hpp>
#include <boost/intrusive_ptr.hpp>

#include <CGAL/Straight_skeleton_aux.h>
#include <CGAL/Straight_skeleton_2.h>
#include <CGAL/Straight_skeleton_builder_traits_2.h>
#include <CGAL/Straight_skeleton_builder_events_2.h>
#include <CGAL/enum.h>

CGAL_BEGIN_NAMESPACE

template<class Traits_, class SSkel_>
class Straight_skeleton_builder_2
{
public:

  typedef Traits_ Traits ;
  typedef SSkel_  SSkel ;

  typedef typename SSkel::Traits::Segment_2 Segment_2 ;

private :

  typedef typename Traits::FT      FT ;
  typedef typename Traits::Point_2 Point_2 ;

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

  typedef typename SSkel::size_type size_type ;

  typedef Straight_skeleton_builder_event_2<SSkel>        Event ;
  typedef Straight_skeleton_builder_edge_event_2<SSkel>   EdgeEvent ;
  typedef Straight_skeleton_builder_split_event_2<SSkel>  SplitEvent ;
  typedef Straight_skeleton_builder_vertex_event_2<SSkel> VertexEvent ;

  typedef boost::intrusive_ptr<Event> EventPtr ;

  typedef std::vector<EventPtr>        EventPtr_Vector ;
  typedef std::vector<Halfedge_handle> Halfedge_handle_vector ;

  typedef typename Halfedge_handle_vector::iterator Halfedge_handle_vector_iterator ;
  typedef typename EventPtr_Vector       ::iterator event_iterator ;

  typedef boost::tuple<Halfedge_handle, Halfedge_handle, Halfedge_handle> BorderTriple ;

  typedef CGAL_SLS_i::Vertex <FT> iVertex ;
  typedef CGAL_SLS_i::Edge   <FT> iEdge ;
  typedef CGAL_SLS_i::Triedge<FT> iTriedge ;

  typedef typename Halfedge::Base HBase;
  typedef typename Vertex::Base   VBase;
  typedef typename Face::Base     FBase;
  typedef typename SSkel::Base    SBase ;

  typedef Straight_skeleton_builder_2<Traits,SSkel> Self ;

public:

  Straight_skeleton_builder_2 ( Traits const& = Traits() ) ;

  // NOTE: The following public method is implemented here in this header file to support some broken compilers.
  // But it's right at the end of the class declaration becuause it needs all of the class.
  //
  // template<class InputPointIterator> Straight_skeleton_builder_2& enter_contour ( InputPointIterator aBegin, InputPointIterator aEnd ) ;


  SSkel construct_skeleton() ;

private :

  class Event_compare : public std::binary_function<bool,EventPtr,EventPtr>
  {
  public:

    Event_compare ( Self const& aBuilder ) : mBuilder(aBuilder) {}

    bool operator() ( EventPtr const& aA, EventPtr const& aB ) const
    {
      return mBuilder.CompareEvents(aA,aB) == LARGER ;
    }

  private:

    Self const& mBuilder ;
  } ;


  typedef std::priority_queue<EventPtr,std::vector<EventPtr>,Event_compare> PQ ;

  typedef std::pair<Vertex_handle,Vertex_handle> Vertex_handle_pair ;

  typedef std::vector<Vertex_handle_pair> SplitNodesVector ;

  struct VertexWrapper
  {
    VertexWrapper( Vertex_handle aVertex )
        :
        mVertex(aVertex)
      , mIsReflex(false)
      , mIsProcessed(false)
      , mIsExcluded(false)
      , mPrev(-1)
      , mNext(-1)
    {}

    Vertex_handle   mVertex ;
    bool            mIsReflex ;
    bool            mIsProcessed ;
    bool            mIsExcluded ;
    int             mPrev ;
    int             mNext ;
    Halfedge_handle mDefiningBorderA ;
    Halfedge_handle mDefiningBorderB ;
    Halfedge_handle mDefiningBorderC ;
    EventPtr_Vector mReflexSplits ; // For fast vertex-event discovery.
  } ;

private :

  inline Halfedge_handle GetDefiningBorderA ( Vertex_handle aV ) const
  {
    return mWrappedVertices[aV->id()].mDefiningBorderA ;
  }
  inline Halfedge_handle GetDefiningBorderB ( Vertex_handle aV ) const
  {
    return mWrappedVertices[aV->id()].mDefiningBorderB ;
  }
  inline Halfedge_handle GetDefiningBorderC ( Vertex_handle aV ) const
  {
    return mWrappedVertices[aV->id()].mDefiningBorderC ;
  }
  inline void SetDefiningBorderA ( Vertex_handle aV, Halfedge_handle aH )
  {
    mWrappedVertices[aV->id()].mDefiningBorderA = aH ;
  }
  inline void SetDefiningBorderB ( Vertex_handle aV, Halfedge_handle aH )
  {
    mWrappedVertices[aV->id()].mDefiningBorderB = aH ;
  }
  inline void SetDefiningBorderC ( Vertex_handle aV, Halfedge_handle aH )
  {
    mWrappedVertices[aV->id()].mDefiningBorderC = aH ;
  }

  static inline iEdge CreateEdge ( Halfedge_const_handle aH )
  {
    Point_2 s = aH->opposite()->vertex()->point() ;
    Point_2 t = aH->vertex()->point() ;
    return iEdge( iVertex(s.x(),s.y()), iVertex(t.x(),t.y()) );
  }

  static inline iTriedge CreateTriedge ( Halfedge_const_handle aE0
                                       , Halfedge_const_handle aE1
                                       , Halfedge_const_handle aE2
                                       )
  {
    return iTriedge(CreateEdge(aE0),CreateEdge(aE1),CreateEdge(aE2));
  }

  static inline iTriedge CreateTriedge ( BorderTriple const& aTriple )
  {
    return iTriedge(CreateEdge(aTriple.get<0>()),CreateEdge(aTriple.get<1>()),CreateEdge(aTriple.get<2>()));
  }

  Vertex_handle GetVertex ( int aIdx )
  {
    return mWrappedVertices[aIdx].mVertex ;
  }

  Vertex_handle GetPrevInLAV ( Vertex_handle aV )
  {
    return GetVertex ( mWrappedVertices[aV->id()].mPrev ) ;
  }

  Vertex_handle GetNextInLAV ( Vertex_handle aV )
  {
    return GetVertex ( mWrappedVertices[aV->id()].mNext ) ;
  }

  Vertex_handle GetNextInLAV_NonSkeleton ( Vertex_handle aV )
  {
    Vertex_handle lNext = GetNextInLAV(aV);
    while ( lNext->is_skeleton() )
      lNext = GetNextInLAV(lNext);
    return lNext;
  }

  void SetPrevInLAV ( Vertex_handle aV, Vertex_handle aPrev )
  {
    mWrappedVertices[aV->id()].mPrev = aPrev->id();
  }

  void SetNextInLAV ( Vertex_handle aV, Vertex_handle aPrev )
  {
    mWrappedVertices[aV->id()].mNext = aPrev->id();
  }

  BorderTriple GetSkeletonVertexDefiningBorders( Vertex_handle aVertex ) const
  {
    CGAL_precondition(aVertex->is_skeleton() ) ;

    return boost::make_tuple( GetDefiningBorderA(aVertex)
                            , GetDefiningBorderB(aVertex)
                            , GetDefiningBorderC(aVertex)
                            ) ;
  }

  void Exclude ( Vertex_handle aVertex )
  {
    mWrappedVertices[aVertex->id()].mIsExcluded = true ;
  }

  void SetIsReflex ( Vertex_handle aVertex )
  {
    mWrappedVertices[aVertex->id()].mIsReflex = true ;
  }

  bool IsReflex ( Vertex_handle aVertex )
  {
    return mWrappedVertices[aVertex->id()].mIsReflex ;
  }

  void SetIsProcessed ( Vertex_handle aVertex )
  {
    mWrappedVertices[aVertex->id()].mIsProcessed = true ;
  }

  bool IsProcessed ( Vertex_handle aVertex )
  {
    return mWrappedVertices[aVertex->id()].mIsProcessed ;
  }

  void AddReflexSplit ( Vertex_handle aSeed, EventPtr aReflexSplit )
  {
    return mWrappedVertices[aSeed->id()].mReflexSplits.push_back(aReflexSplit) ;
  }
  EventPtr_Vector const& GetReflexSplits ( Vertex_handle aSeed )
  {
    return mWrappedVertices[aSeed->id()].mReflexSplits  ;
  }

  void EnqueEvent( EventPtr aEvent )
  {
    mPQ.push(aEvent);
    CGAL_SSBUILDER_TRACE(0, *aEvent);
  }

  EventPtr PopEventFromPQ()
  {
    EventPtr rR = mPQ.top();
    mPQ.pop();
    return rR ;
  }

  // Returns 1 aE is in the set (aA,aB,aC), 0 otherwise
  int CountInCommon( Halfedge_handle aE, Halfedge_handle aA, Halfedge_handle aB, Halfedge_handle aC ) const
  {
    return aE == aA || aE == aB || aE == aC ? 1 : 0 ;
  }

  // Returns the number of common halfedges in the sets (aXA,aXB,aXC) and (aYA,aYB,aYC)
  int CountInCommon( Halfedge_handle aXA, Halfedge_handle aXB, Halfedge_handle aXC
                   , Halfedge_handle aYA, Halfedge_handle aYB, Halfedge_handle aYC
                   ) const
  {
    return   CountInCommon(aXA,aYA,aYB,aYC)
           + CountInCommon(aXB,aYA,aYB,aYC)
           + CountInCommon(aXC,aYA,aYB,aYC) ;
  }

  // Returns true if the intersection of the sets (aXA,aXB,aXC) and (aYA,aYB,aYC) has size exactly 2
  // (that is, both sets have 2 elements in common)
  bool HaveTwoInCommon( Halfedge_handle aXA, Halfedge_handle aXB, Halfedge_handle aXC
                      , Halfedge_handle aYA, Halfedge_handle aYB, Halfedge_handle aYC
                      ) const
  {
    return CountInCommon(aXA,aXB,aXC,aYA,aYB,aYC) == 2 ;
  }

  // Returns true if the sets of halfedges (aXA,aXB,aXC) and (aYA,aYB,aYC) are equivalent
  // (one is a permutation of the other)
  bool AreTheSameTriple( Halfedge_handle aXA, Halfedge_handle aXB, Halfedge_handle aXC
                       , Halfedge_handle aYA, Halfedge_handle aYB, Halfedge_handle aYC
                       ) const
  {
    return CountInCommon(aXA,aXB,aXC,aYA,aYB,aYC) == 3 ;
  }

  // Returns the 0-base index of the one element from (aX[3]) NOT IN (aY[3])
  // NOTE: This function shall be called only when it is known that such an element exists
  // as 2 is returned by default without proper testing. That is, this function is for vertex-event analysis only.
  int GetUnique( Halfedge_handle aX[], Halfedge_handle aY[] ) const
  {
    return CountInCommon(aX[0],aY[0],aY[1],aY[2]) == 0 ? 0
             : CountInCommon(aX[1],aY[0],aY[1],aY[2]) == 0 ? 1
               : 2 ;
  }

  // Sorts the elements in the sets aX[2] and aY[3] returing (D0,D1,E0,E1)
  // where D0,D1 are unique elements in aX and aY respectively and E0,E1 are elements in common.
  // NOTE: This function shall only be called when it is known that thet sets aX and aY can indeed be sorted this way.
  // That is, this function is for vertex-event analysis only.
  boost::tuple<Halfedge_handle,Halfedge_handle,Halfedge_handle,Halfedge_handle>
    SortTwoDistinctAndTwoEqual( Halfedge_handle aX[], Halfedge_handle aY[] ) const
  {
     int lUniqueX = GetUnique(aX,aY) ;
     int lUniqueY = GetUnique(aY,aX) ;
     int lCommon1 = ( lUniqueX + 1 ) % 3 ;
     int lCommon2 = ( lUniqueX + 2 ) % 3 ;
     return boost::make_tuple(aX[lUniqueX],aY[lUniqueY],aX[lCommon1],aX[lCommon2]);
  }


  bool ExistEvent ( Halfedge_const_handle aE0, Halfedge_const_handle aE1, Halfedge_const_handle aE2 ) const
  {
    return Exist_sls_event_2<Traits>(mTraits)()(CreateTriedge(aE0, aE1, aE2));
  }

  bool IsEventInsideOffsetZone( Halfedge_const_handle aReflexL
                              , Halfedge_const_handle aReflexR
                              , Halfedge_const_handle aOpposite
                              , Halfedge_const_handle aOppositePrev
                              , Halfedge_const_handle aOppositeNext
                              ) const
  {
    return Is_sls_event_inside_offset_zone_2<Traits>(mTraits)()( CreateTriedge(aReflexL     , aReflexR, aOpposite)
                                                               , CreateTriedge(aOppositePrev,aOpposite, aOppositeNext)
                                                               ) ;
  }

  Comparison_result CompareEvents ( iTriedge const& aA, iTriedge const& aB ) const
  {
    return Compare_sls_event_times_2<Traits>(mTraits)()(aA,aB) ;
  }

  Comparison_result CompareEvents ( EventPtr const& aA, EventPtr const& aB ) const
  {
    if ( !AreTheSameTriple( aA->border_a(), aA->border_b(), aA->border_c()
                          , aB->border_a(), aB->border_b(), aB->border_c()
                          )
        )
    {
      return CompareEvents( CreateTriedge(aA->border_a(), aA->border_b(), aA->border_c())
                          , CreateTriedge(aB->border_a(), aB->border_b(), aB->border_c())
                          ) ;
    }
    else return EQUAL ;
  }

  Comparison_result CompareEventsDistanceToSeed ( Vertex_handle   aSeed
                                                , EventPtr const& aA
                                                , EventPtr const& aB
                                                ) const
  {
    if ( !AreTheSameTriple( aA->border_a(), aA->border_b(), aA->border_c()
                          , aB->border_a(), aB->border_b(), aB->border_c()
                          )
        )
    {
      if ( aSeed->is_skeleton() )
      {
        BorderTriple lTriple = GetSkeletonVertexDefiningBorders(aSeed);

        CGAL_SSBUILDER_TRACE(3
                            ,"Seed N" << aSeed->id() << " is a skeleton node,"
                            << " defined by: E" << lTriple.get<0>()->id()
                                       << ", E" << lTriple.get<1>()->id()
                                       << ", E" << lTriple.get<2>()->id()
                            );

        return Compare_sls_event_distance_to_seed_2<Traits>(mTraits)()( CreateTriedge(lTriple)
                                                                      , CreateTriedge(aA->border_a(), aA->border_b(), aA->border_c())
                                                                      , CreateTriedge(aB->border_a(), aB->border_b(), aB->border_c())
                                                                      ) ;
      }
      else
      {
        return Compare_sls_event_distance_to_seed_2<Traits>(mTraits)()( aSeed->point()
                                                                      , CreateTriedge(aA->border_a(), aA->border_b(), aA->border_c())
                                                                      , CreateTriedge(aB->border_a(), aB->border_b(), aB->border_c())
                                                                      ) ;
      }
    }
    else return EQUAL ;
  }
  
  bool AreEventsSimultaneous( EventPtr const& aX, EventPtr const& aY ) const
  {
    Halfedge_handle xa = aX->border_a() ;
    Halfedge_handle xb = aX->border_b() ;
    Halfedge_handle xc = aX->border_c() ;
    Halfedge_handle ya = aY->border_a() ;
    Halfedge_handle yb = aY->border_b() ;
    Halfedge_handle yc = aY->border_c() ;

    if ( HaveTwoInCommon(xa,xb,xc,ya,yb,yc) )
         return Are_sls_events_simultaneous_2<Traits>(mTraits)()( CreateTriedge(xa,xb,xc), CreateTriedge(ya,yb,yc)) ;
    else return false ;
  }

  bool IsNewEventInThePast( Halfedge_handle aBorderA
                          , Halfedge_handle aBorderB
                          , Halfedge_handle aBorderC
                          , Vertex_handle   aSeedNode
                          ) const
  {
    bool rResult = false ;

    Halfedge_handle lSeedBorderA, lSeedBorderB, lSeedBorderC ;

    boost::tie(lSeedBorderA,lSeedBorderB,lSeedBorderC) = GetSkeletonVertexDefiningBorders(aSeedNode) ;

    if ( !AreTheSameTriple(aBorderA,aBorderB,aBorderC,lSeedBorderA,lSeedBorderB,lSeedBorderC) )
    {
      if ( CompareEvents( CreateTriedge(aBorderA,aBorderB,aBorderC)
                        , CreateTriedge(lSeedBorderA,lSeedBorderB,lSeedBorderC)
                        ) == SMALLER
         )
      {
#ifdef CGAL_STRAIGHT_SKELETON_ENABLE_TRACE
        FT lTime1, lTime2 ;
        Point_2 lP1, lP2 ;
        boost::tie(lTime1,lP1) = ConstructEventTimeAndPoint(CreateTriedge(aBorderA,aBorderB,aBorderC));
        boost::tie(lTime2,lP2) = ConstructEventTimeAndPoint(CreateTriedge(lSeedBorderA,lSeedBorderB,lSeedBorderC));
        CGAL_SSBUILDER_TRACE(1,"New event for N" << aSeedNode->id() << ", with t=" << lTime1 << ", is in the past (current t="
                            << lTime2 << "). discarded." 
                            ) ;
#endif                            
        rResult = true ;
      }   
    }

    return rResult ;
  }

  boost::tuple<FT,Point_2> ConstructEventTimeAndPoint( iTriedge const& aTri ) const
  {
    return Construct_sls_event_time_and_point_2<Traits>(mTraits)()(aTri);
  }

  void SetEventTimeAndPoint( Event& aE )
  {
    FT lTime ; Point_2 lP ;
    boost::tie(lTime,lP) = ConstructEventTimeAndPoint( CreateTriedge(aE.border_a(),aE.border_b(),aE.border_c()) );
    aE.SetTimeAndPoint(lTime,lP);
  }

  void EraseBisector( Halfedge_handle aB )
  {
    mSSkel.SBase::edges_erase(aB);
  }

  BorderTriple GetDefiningBorders( Vertex_handle aA, Vertex_handle aB ) ;

  bool AreBisectorsCoincident ( Halfedge_const_handle aA, Halfedge_const_handle aB ) const ;

  void CollectSplitEvent( Vertex_handle    aNode
                        , Halfedge_handle  aReflexLBorder
                        , Halfedge_handle  aReflexRBorder
                        , Halfedge_handle  aOppositeBorder
                        ) ;

  void CollectSplitEvents( Vertex_handle aNode ) ;

  EventPtr FindEdgeEvent( Vertex_handle aLNode, Vertex_handle aRNode ) ;


  EventPtr FindVertexEvent( EventPtr aE0, Vertex_handle aOV ) ;
  EventPtr FindVertexEvent( EventPtr aE0 ) ;

  void HandleSimultaneousEdgeEvent( Vertex_handle aA, Vertex_handle aB ) ;

  void CollectNewEvents( Vertex_handle aNode ) ;
  void UpdatePQ( Vertex_handle aV ) ;
  void CreateInitialEvents();
  void CreateContourBisectors();
  void InitPhase();

  bool SetupVertexEventNode( Vertex_handle   aNode
                           , Halfedge_handle aDefiningBorderA
                           , Halfedge_handle aDefiningBorderB
                           );

  Vertex_handle LookupOnSLAV ( Halfedge_handle aOBorder, Event const& aEvent ) ;

  Vertex_handle_pair ConstructSplitEventNodes ( SplitEvent&  aEvent, Vertex_handle aOppR ) ;
  Vertex_handle      ConstructEdgeEventNode   ( EdgeEvent&   aEvent ) ;
  Vertex_handle_pair ConstructVertexEventNodes( VertexEvent& aEvent ) ;

  void HandleSplitEvent          ( EventPtr aEvent, Vertex_handle aOppR ) ;
  void HandleEdgeEvent           ( EventPtr aEvent ) ;
  void HandleVertexEvent         ( EventPtr aEvent ) ;
  void HandlePotentialSplitEvent ( EventPtr aEvent ) ;
  bool IsProcessed               ( EventPtr aEvent ) ;

  void Propagate();

  void MergeSplitNodes ( Vertex_handle_pair aSplitNodes ) ;

  void FinishUp();

  void Run();

private:

#ifdef CGAL_STRAIGHT_SKELETON_ENABLE_SHOW
  template<class Halfedge>
  void DrawBisector ( Halfedge aHalfedge )
  {
    SS_IO_AUX::ScopedSegmentDrawing draw_( aHalfedge->opposite()->vertex()->point()
                                         , aHalfedge->vertex()->point()
                                         , aHalfedge->is_inner_bisector() ? CGAL::BLUE  : CGAL::GREEN
                                         , aHalfedge->is_inner_bisector() ? "IBisector" : "CBisector"
                                         ) ;
    draw_.Release();
  }

#endif


private:

  //Input
  Traits mTraits ;

  typename Traits::Left_turn_2 Left_turn ;
  typename Traits::Collinear_2 Collinear ;

  std::vector<VertexWrapper> mWrappedVertices ;
  Halfedge_handle_vector     mDanglingBisectors ;
  Halfedge_handle_vector     mContourHalfedges ;

  std::list<Vertex_handle> mSLAV ;

  EventPtr_Vector  mSplitEvents ;
  SplitNodesVector mSplitNodes ;

  Event_compare mEventCompare ;

  int mVertexID ;
  int mEdgeID   ;
  int mEventID  ;
  int mStepID   ;

  PQ mPQ ;

  //Output
  SSkel mSSkel ;

public:

  template<class InputPointIterator>
  Straight_skeleton_builder_2& enter_contour ( InputPointIterator aBegin, InputPointIterator aEnd  )
  {
    CGAL_SSBUILDER_TRACE(0,"Inserting Connected Component of the Boundary....");

    if ( std::distance(aBegin,aEnd) >= 3 )
    {
      Halfedge_handle lFirstCCWBorder ;
      Halfedge_handle lPrevCCWBorder ;
      Halfedge_handle lNextCWBorder ;
      Vertex_handle   lFirstVertex ;
      Vertex_handle   lPrevVertex ;

      InputPointIterator lCurr = aBegin ;
      while ( lCurr != aEnd )
      {
        Halfedge_handle lCCWBorder = mSSkel.SBase::edges_push_back ( Halfedge(mEdgeID),Halfedge(mEdgeID+1)  );
        Halfedge_handle lCWBorder = lCCWBorder->opposite();
        mEdgeID += 2 ;

        mContourHalfedges.push_back(lCCWBorder);

        Vertex_handle lVertex = mSSkel.SBase::vertices_push_back( Vertex(mVertexID++,*lCurr) ) ;
        CGAL_SSBUILDER_TRACE(2,"Vertex: V" << lVertex->id() << " at " << lVertex->point() );
        mWrappedVertices.push_back( VertexWrapper(lVertex) ) ;

        Face_handle lFace = mSSkel.SBase::faces_push_back( Face() ) ;

        lCCWBorder->HBase::set_face    (lFace);
        lFace     ->FBase::set_halfedge(lCCWBorder);

        lVertex   ->VBase::set_halfedge(lCCWBorder);
        lCCWBorder->HBase::set_vertex  (lVertex);

        if ( lCurr == aBegin )
        {
          lFirstVertex    = lVertex ;
          lFirstCCWBorder = lCCWBorder ;
        }
        else
        {
          SetPrevInLAV(lVertex    ,lPrevVertex);
          SetNextInLAV(lPrevVertex,lVertex    );

          SetDefiningBorderA(lVertex    ,lCCWBorder);
          SetDefiningBorderB(lPrevVertex,lCCWBorder);

          lCWBorder->HBase::set_vertex(lPrevVertex);

          lCCWBorder    ->HBase::set_prev(lPrevCCWBorder);
          lPrevCCWBorder->HBase::set_next(lCCWBorder);

          lNextCWBorder->HBase::set_prev(lCWBorder);
          lCWBorder    ->HBase::set_next(lNextCWBorder);

          CGAL_SSBUILDER_TRACE(2,"CCW Border: E" << lCCWBorder->id() << ' ' << lPrevVertex->point() << " -> " << lVertex    ->point());
          CGAL_SSBUILDER_TRACE(2,"CW  Border: E" << lCWBorder ->id() << ' ' << lVertex    ->point() << " -> " << lPrevVertex->point() );

          CGAL_SSBUILDER_SHOW
          ( 
            SS_IO_AUX::ScopedSegmentDrawing draw_(lPrevVertex->point(),lVertex->point(), CGAL::RED, "Border" ) ;
            draw_.Release();
          )
        }

        ++ lCurr ;

        lPrevVertex    = lVertex ;
        lPrevCCWBorder = lCCWBorder ;
        lNextCWBorder  = lCWBorder ;
      }

      SetPrevInLAV(lFirstVertex,lPrevVertex );
      SetNextInLAV(lPrevVertex ,lFirstVertex);

      SetDefiningBorderA(lFirstVertex,lFirstCCWBorder);
      SetDefiningBorderB(lPrevVertex ,lFirstCCWBorder);

      lFirstCCWBorder->opposite()->HBase::set_vertex(lPrevVertex);

      CGAL_SSBUILDER_SHOW
      ( SS_IO_AUX::ScopedSegmentDrawing draw_(lPrevVertex->point(),lFirstVertex->point(), CGAL::RED, "Border" ) ;
        draw_.Release();
      )

      lFirstCCWBorder->HBase::set_prev(lPrevCCWBorder);
      lPrevCCWBorder ->HBase::set_next(lFirstCCWBorder);

      lPrevCCWBorder ->opposite()->HBase::set_prev(lFirstCCWBorder->opposite());
      lFirstCCWBorder->opposite()->HBase::set_next(lPrevCCWBorder ->opposite());

      CGAL_SSBUILDER_TRACE(2
                          , "CCW Border: E" << lFirstCCWBorder->id()
                          << ' ' << lPrevVertex ->point() << " -> " << lFirstVertex->point() << '\n'
                          << "CW  Border: E" << lFirstCCWBorder->opposite()->id()
                          << ' ' << lFirstVertex->point() << " -> " << lPrevVertex ->point()
                          );
    }

    for ( Vertex_iterator v = mSSkel.SBase::vertices_begin(); v != mSSkel.SBase::vertices_end(); ++ v )
    {
      mSLAV.push_back(static_cast<Vertex_handle>(v));
      Vertex_handle lPrev = GetPrevInLAV(v) ;
      Vertex_handle lNext = GetNextInLAV(v) ;

      bool lCollinear = Collinear( lPrev->point(),v->point(),lNext->point() ) ;
      if ( lCollinear || !Left_turn( lPrev->point(),v->point(),lNext->point() ) )
      {
        SetIsReflex(v);
        CGAL_SSBUILDER_TRACE(1,(lCollinear ? "COLLINEAR " : "Reflex ") << "vertex: N" << v->id() );
      }
    }

    return *this ;
  }

} ;

CGAL_END_NAMESPACE

#ifdef CGAL_CFG_NO_AUTOMATIC_TEMPLATE_INCLUSION
#  include <CGAL/Straight_skeleton_builder_2.C>
#endif


#endif // CGAL_STRAIGHT_SKELETON_BUILDER_2_H //
// EOF //
