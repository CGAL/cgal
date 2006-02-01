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
// file          : include/CGAL/Straight_skeleton_builder_2.h
// package       : Straight_skeleton_2 (1.1.0)
//
// author(s)     : Fernando Cacciola
// maintainer    : Fernando Cacciola <fernando_cacciola@hotmail>
// coordinator   : Fernando Cacciola <fernando_cacciola@hotmail>
//
// ============================================================================
#ifndef CGAL_STRAIGHT_SKELETON_BUILDER_2_H
#define CGAL_STRAIGHT_SKELETON_BUILDER_2_H 1

#include <list>

#include <boost/tuple/tuple.hpp>
#include <boost/intrusive_ptr.hpp>

#include <CGAL/Straight_skeleton_aux.h>
#include <CGAL/Straight_skeleton_2.h>
#include <CGAL/Straight_skeleton_builder_traits_2.h>
#include <CGAL/Straight_skeleton_builder_events_2.h>
#include <CGAL/enum.h>

CGAL_BEGIN_NAMESPACE

template<class Traits_, class Ssds_>
class Straight_skeleton_builder_2
{
public:

  typedef Traits_ Traits ;
  typedef Ssds_   Ssds ;

  typedef typename Ssds::Traits::Segment_2 Segment_2 ;

private :

  typedef typename Traits::FT      FT ;
  typedef typename Traits::Point_2 Point_2 ;

  typedef typename Ssds::Vertex   Vertex ;
  typedef typename Ssds::Halfedge Halfedge ;
  typedef typename Ssds::Face     Face ;

  typedef typename Ssds::Vertex_const_handle   Vertex_const_handle ;
  typedef typename Ssds::Halfedge_const_handle Halfedge_const_handle ;
  typedef typename Ssds::Face_const_handle     Face_const_handle ;

  typedef typename Ssds::Vertex_const_iterator   Vertex_const_iterator ;
  typedef typename Ssds::Halfedge_const_iterator Halfedge_const_iterator ;
  typedef typename Ssds::Face_const_iterator     Face_const_iterator ;

  typedef typename Ssds::Vertex_handle   Vertex_handle ;
  typedef typename Ssds::Halfedge_handle Halfedge_handle ;
  typedef typename Ssds::Face_handle     Face_handle ;

  typedef typename Ssds::Vertex_iterator   Vertex_iterator ;
  typedef typename Ssds::Halfedge_iterator Halfedge_iterator ;
  typedef typename Ssds::Face_iterator     Face_iterator ;

  typedef typename Ssds::size_type size_type ;

  typedef Straight_skeleton_builder_event_2<Ssds>        Event ;
  typedef Straight_skeleton_builder_edge_event_2<Ssds>   EdgeEvent ;
  typedef Straight_skeleton_builder_split_event_2<Ssds>  SplitEvent ;
  typedef Straight_skeleton_builder_vertex_event_2<Ssds> VertexEvent ;

  typedef boost::intrusive_ptr<Event> EventPtr ;

  typedef std::vector<EventPtr> EventPtr_Vector ;

  typedef typename EventPtr_Vector::iterator event_iterator ;

  typedef tuple<Halfedge_handle, Halfedge_handle, Halfedge_handle> BorderTriple ;

  typedef tuple<FT,FT,FT,FT> Edge ;

  typedef tuple<Edge,Edge,Edge> Edge_triple ;

  typedef typename Halfedge::SSBase SSBase;
  typedef typename Halfedge::Base   HBase;
  typedef typename Vertex::Base     VBase;
  typedef typename Face::Base       FBase;
  typedef typename Ssds::Base       SBase ;

  typedef Straight_skeleton_builder_2<Traits,Ssds> Self ;

public:

  Straight_skeleton_builder_2 ( Traits const& = Traits() ) ;

  template<class InputPointIterator>
  Straight_skeleton_builder_2& enter_contour ( InputPointIterator aBegin, InputPointIterator aEnd )
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

  Ssds construct_skeleton() ;

private :

  class Event_compare : public std::binary_function<bool,EventPtr,EventPtr>
  {
  public:

    Event_compare ( Self const& aBuilder ) : mBuilder(aBuilder) {}

    bool operator() ( EventPtr const& aA, EventPtr const& aB ) const
    {
      return mBuilder.CompareEvents(aA,aB) == SMALLER ;
    }

  private:

    Self const& mBuilder ;
  } ;

  typedef std::list<EventPtr> PQ ;

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
  } ;

private :

  inline Halfedge_handle GetDefiningBorderA ( Vertex_handle aV )
  {
    return mWrappedVertices[aV->id()].mDefiningBorderA ;
  }
  inline Halfedge_handle GetDefiningBorderB ( Vertex_handle aV )
  {
    return mWrappedVertices[aV->id()].mDefiningBorderB ;
  }
  inline void SetDefiningBorderA ( Vertex_handle aV, Halfedge_handle aH )
  {
    mWrappedVertices[aV->id()].mDefiningBorderA = aH ;
  }
  inline void SetDefiningBorderB ( Vertex_handle aV, Halfedge_handle aH )
  {
    mWrappedVertices[aV->id()].mDefiningBorderB = aH ;
  }

  static inline Edge GetEdge ( Halfedge_const_handle aH )
  {
    Point_2 s = aH->opposite()->vertex()->point() ;
    Point_2 t = aH->vertex()->point() ;
    return make_tuple(s.x(),s.y(),t.x(),t.y());
  }

  static inline Edge_triple GetEdgeTriple ( Halfedge_const_handle aE0
                                          , Halfedge_const_handle aE1
                                          , Halfedge_const_handle aE2
                                          )
  {
    return make_tuple(GetEdge(aE0),GetEdge(aE1),GetEdge(aE2));
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

  EventPtr PopEventFromPQ()
  {
    typename PQ::iterator f = std::min_element(mPQ.begin(), mPQ.end(), mEventCompare) ;
    EventPtr rR = *f ;
    mPQ.erase(f);
    return rR ;
  }

  bool ExistEvent ( Halfedge_const_handle aE0, Halfedge_const_handle aE1, Halfedge_const_handle aE2 ) const
  {
    return Exist_sls_event_2<Traits>(mTraits)()(GetEdgeTriple(aE0, aE1, aE2));
  }

  bool IsEventInsideOffsetZone( Halfedge_const_handle aReflexL
                              , Halfedge_const_handle aReflexR
                              , Halfedge_const_handle aOpposite
                              , Halfedge_const_handle aOppositePrev
                              , Halfedge_const_handle aOppositeNext
                              ) const
  {
    return Is_sls_event_inside_offset_zone_2<Traits>(mTraits)()( GetEdgeTriple(aReflexL     , aReflexR, aOpposite)
                                                               , GetEdgeTriple(aOppositePrev,aOpposite, aOppositeNext)
                                                               ) ;
  }

  Comparison_result CompareEvents ( EventPtr const& aA, EventPtr const& aB ) const
  {
    return Compare_sls_event_times_2<Traits>(mTraits)()( GetEdgeTriple(aA->border_a(), aA->border_b(), aA->border_c())
                                                       , GetEdgeTriple(aB->border_a(), aB->border_b(), aB->border_c())
                                                       ) ;
  }

  Comparison_result CompareEventsDistanceToSeed ( Vertex_handle   aSeed
                                                , EventPtr const& aA
                                                , EventPtr const& aB
                                                ) const
  {
    if ( aSeed->is_skeleton() )
    {
      Halfedge_handle lBorder0 = aSeed->halfedge()->face()->halfedge();
      Halfedge_handle lBorder1 = aSeed->halfedge()->opposite()->prev()->face()->halfedge();
      Halfedge_handle lBorder2 = aSeed->halfedge()->opposite()->prev()->opposite()->face()->halfedge();

      return Compare_sls_event_distance_to_seed_2<Traits>(mTraits)()( GetEdgeTriple(lBorder0,lBorder1,lBorder2)
                                                                    , GetEdgeTriple(aA->border_a(), aA->border_b(), aA->border_c())
                                                                    , GetEdgeTriple(aB->border_a(), aB->border_b(), aB->border_c())
                                                                    ) ;
    }
    else
    {
      return Compare_sls_event_distance_to_seed_2<Traits>(mTraits)()( aSeed->point()
                                                                    , GetEdgeTriple(aA->border_a(), aA->border_b(), aA->border_c())
                                                                    , GetEdgeTriple(aB->border_a(), aB->border_b(), aB->border_c())
                                                                    ) ;
    }
  }

  // Returns true if aE is in the set (aA,aB,aC)
  bool IsBorderInTriple( Halfedge_handle aE, Halfedge_handle aA, Halfedge_handle aB, Halfedge_handle aC ) const
  {
    return aE == aA || aE == aB || aE == aC ;
  }

  // Returns true if the intersection of the sets (aXA,aXB,aXC) and (aYA,aYB,aYC) has size exactly 2
  // (that is, both sets have 2 elements in common)
  bool HaveTwoInCommon( Halfedge_handle aXA, Halfedge_handle aXB, Halfedge_handle aXC
                      , Halfedge_handle aYA, Halfedge_handle aYB, Halfedge_handle aYC
                      ) const
  {
    int lC = IsBorderInTriple(aXA,aYA,aYB,aYC) ? 1 :0 ;
    lC    += IsBorderInTriple(aXB,aYA,aYB,aYC) ? 1 :0 ;
    lC    += IsBorderInTriple(aXC,aYA,aYB,aYC) ? 1 :0 ;
    return lC == 2 ;
  }

  // Returns the 0-base index of the one element from (aX[3]) NOT IN (aY[3])
  // NOTE: This function shall be called only when it is known that such an element exists
  // as 2 is returned by default without proper testing. That is, this function is for vertex-event analysis only.
  int GetUnique( Halfedge_handle aX[], Halfedge_handle aY[] ) const
  {
    return !IsBorderInTriple(aX[0],aY[0],aY[1],aY[2]) ? 0
             : !IsBorderInTriple(aX[1],aY[0],aY[1],aY[2]) ? 1
               : 2 ;
  }

  // Sorts the elements in the sets aX[2] and aY[3] returing (D0,D1,E0,E1)
  // where D0,D1 are unique elements in aX and aY respectively and E0,E1 are elements in common.
  // NOTE: This function shall only be called when it is known that thet sets aX and aY can indeed be sorted this way.
  // That is, this function is for vertex-event analysis only.
  tuple<Halfedge_handle,Halfedge_handle,Halfedge_handle,Halfedge_handle>
    SortTwoDistinctAndTwoEqual( Halfedge_handle aX[], Halfedge_handle aY[] ) const
  {
     int lUniqueX = GetUnique(aX,aY) ;
     int lUniqueY = GetUnique(aY,aX) ;
     int lCommon1 = ( lUniqueX + 1 ) % 3 ;
     int lCommon2 = ( lUniqueX + 2 ) % 3 ;
     return make_tuple(aX[lUniqueX],aY[lUniqueY],aX[lCommon1],aX[lCommon2]);
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
         return Are_sls_events_simultaneous_2<Traits>(mTraits)()( GetEdgeTriple(xa,xb,xc), GetEdgeTriple(ya,yb,yc)) ;
    else return false ;
  }

  void SetEventTimeAndPoint( Event& aE )
  {
    FT lTime ; Point_2 lP ;
    tie(lTime,lP) = Construct_sls_event_time_and_point_2<Traits>(mTraits)()( GetEdgeTriple(aE.border_a(), aE.border_b(), aE.border_c()) );
    aE.SetTimeAndPoint(lTime,lP);
  }

  void EraseBisector( Halfedge_handle aB )
  {
    mSS.SBase::edges_erase(aB);
  }

  BorderTriple GetDefiningBorders( Vertex_handle aA, Vertex_handle aB ) ;

  bool AreBisectorsCoincident ( Halfedge_const_handle aA, Halfedge_const_handle aB ) const ;

  void CollectSplitEvent( Vertex_handle    aNode
                        , Halfedge_handle  aReflexLBorder
                        , Halfedge_handle  aReflexRBorder
                        , Halfedge_handle  aOppositeBorder
                        , EventPtr_Vector& aCandidates
                        ) ;

  EventPtr FindEdgeEvent( Vertex_handle aLNode, Vertex_handle aRNode ) ;

  EventPtr FindSplitEvent( Vertex_handle aNode ) ;

  EventPtr ChooseBestNewEvent( Vertex_handle aNode, EventPtr_Vector& aEvents ) ;

  EventPtr FindBestNewEvent( Vertex_handle aNode ) ;

  void FindVertexEvents();

  void HandleSimultaneousEdgeEvent( Vertex_handle aA, Vertex_handle aB ) ;

  void CollectNewEvents( Vertex_handle aNode, EventPtr_Vector& aEvents ) ;
  void AddNewEvent( Vertex_handle aNode ) ;
  void UpdatePQ( Vertex_handle aV ) ;
  void CreateInitialEvents();
  void CreateContourBisectors();
  void InitPhase();

  bool SetupVertexEventNode( Vertex_handle   aNode
                           , Halfedge_handle aDefiningBorderA
                           , Halfedge_handle aDefiningBorderB
                           );

  Vertex_handle LookupOnSLAV ( Halfedge_handle aOBorder, Event const& aEvent ) ;

  Vertex_handle      ConstructEdgeEventNode   ( EdgeEvent&   aEvent ) ;
  Vertex_handle_pair ConstructSplitEventNodes ( SplitEvent&  aEvent ) ;
  Vertex_handle_pair ConstructVertexEventNodes( VertexEvent& aEvent ) ;

  void HandleEdgeEvent  ( EventPtr aEvent ) ;
  void HandleSplitEvent ( EventPtr aEvent ) ;
  void HandleVertexEvent( EventPtr aEvent ) ;

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


  //Internal

  typedef std::vector<Halfedge_handle> Halfedge_handle_vector ;

  typedef typename Halfedge_handle_vector::iterator Halfedge_handle_vector_iterator ;

  std::vector<VertexWrapper>   mWrappedVertices ;
  Halfedge_handle_vector       mDanglingBisectors ;
  Halfedge_handle_vector       mContourHalfedges ;

  EventPtr_Vector  mSplitEvents ;
  SplitNodesVector mSplitNodes ;

  Event_compare mEventCompare ;

  PQ mPQ ;

  int mVertexID ;
  int mEdgeID   ;
  int mEventID  ;
  int mStepID   ;

  //Output
  Ssds mSS ;
} ;

CGAL_END_NAMESPACE

#ifdef CGAL_CFG_NO_AUTOMATIC_TEMPLATE_INCLUSION
#  include <CGAL/Straight_skeleton_builder_2.C>
#endif


#endif // CGAL_STRAIGHT_SKELETON_BUILDER_2_H //
// EOF //
