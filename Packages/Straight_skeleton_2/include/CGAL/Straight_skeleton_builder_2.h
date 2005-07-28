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

  typedef Straight_skeleton_builder_event_2<Ssds>       Event ;
  typedef Straight_skeleton_builder_edge_event_2<Ssds>  EdgeEvent ;
  typedef Straight_skeleton_builder_split_event_2<Ssds> SplitEvent ;

  typedef boost::intrusive_ptr<Event> EventPtr ;

  typedef std::vector<EventPtr> EventPtr_Vector ;

  typedef boost::tuple<Halfedge_handle, Halfedge_handle, Halfedge_handle> BorderTriple ;

  typedef typename Halfedge::SSBase SSBase;
  typedef typename Halfedge::Base   HBase;
  typedef typename Vertex::Base     VBase;
  typedef typename Face::Base       FBase;
  typedef typename Ssds::Base       SBase ;
  
  typedef Straight_skeleton_builder_2<Traits,Ssds> Self ;
  
public:

  Straight_skeleton_builder_2 ( Traits const& = Traits() ) ;

  template<class InputPointIterator>
  Straight_skeleton_builder_2& insert_CCB ( InputPointIterator aBegin, InputPointIterator aEnd ) ;

  Ssds proceed() ;

private :

class Event_compare : public std::binary_function<bool,EventPtr,EventPtr>
  {
  public:

    Event_compare ( Self const& aBuilder ) : mBuilder(aBuilder) {}

    bool operator() ( EventPtr const& aA, EventPtr const& aB )  { return mBuilder.CompareEvents(aA,aB); }

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

    Vertex_handle mVertex ;
    bool          mIsReflex ;
    bool          mIsProcessed ;
    bool          mIsExcluded ;
    int           mPrev ;
    int           mNext ;
  } ;

  struct HalfedgeWrapper
  {
    HalfedgeWrapper( Halfedge_handle aHalfedge)
        : mHalfedge(aHalfedge), mPrevInCCB(-1), mNextInCCB(-1), mIsExcluded(false) {}

    Halfedge_handle mHalfedge ;
    int  mPrevInCCB ;
    int  mNextInCCB ;
    bool mIsExcluded ;
  } ;

private :

  static inline Halfedge_handle GetDefiningBorder0 ( Vertex_handle aV )
  {
    return aV->halfedge()->face()->halfedge();
  }
  static inline Halfedge_handle HasDefiningBorder1 ( Vertex_handle aV )
  {
    return handle_assigned(aV->halfedge()->opposite()->prev()->face()) ;
  }
  static inline Halfedge_handle GetDefiningBorder1 ( Vertex_handle aV )
  {
    return aV->halfedge()->opposite()->prev()->face()->halfedge();
  }
  static inline Halfedge_handle GetDefiningBorder2 ( Vertex_handle aV )
  {
    return aV->halfedge()->opposite()->prev()->opposite()->face()->halfedge();
  }
  static inline std::pair<Point_2,Point_2> GetSegment ( Halfedge_const_handle aH )
  {
    return std::make_pair(aH->opposite()->vertex()->point(),aH->vertex()->point());
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

  void Exclude ( Halfedge_handle aHalfedge )
  {
    mWrappedHalfedges[aHalfedge->id()].mIsExcluded = true ;
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

  void SetPrevInCCB ( Halfedge_handle aBorder, int aIdx )
  {
    mWrappedHalfedges[aBorder->id()].mPrevInCCB = aIdx ;
  }

  void SetNextInCCB ( Halfedge_handle aBorder, int aIdx )
  {
    mWrappedHalfedges[aBorder->id()].mNextInCCB = aIdx ;
  }

  Halfedge_handle GetPrevInCCB ( Halfedge_handle aBorder )
  {
    return mBorderHalfedges[ mWrappedHalfedges[aBorder->id()].mPrevInCCB ] ;
  }

  Halfedge_handle GetNextInCCB ( Halfedge_handle aBorder )
  {
    return mBorderHalfedges[ mWrappedHalfedges[aBorder->id()].mNextInCCB ] ;
  }

  EventPtr PopEventFromPQ()
  {
    typename PQ::iterator f = std::min_element(mPQ.begin(), mPQ.end(), mEventCompare) ;
    EventPtr rR = *f ;
    mPQ.erase(f);
    return rR ;
  }

  bool Exist_event ( Halfedge_const_handle aE0, Halfedge_const_handle aE1, Halfedge_const_handle aE2 )
  {
    return mTraits.exist_event_object()(GetSegment(aE0), GetSegment(aE1), GetSegment(aE2));
  }
  
  bool IsEventInsideOffsetZone( Halfedge_const_handle aReflexL     
                              , Halfedge_const_handle aReflexR
                              , Halfedge_const_handle aOpposite
                              , Halfedge_const_handle aOppositePrev
                              , Halfedge_const_handle aOppositeNext
  {
    return mTraits.is_event_inside_offset_zone_object()( GetSegment(aReflexL)
                                                       , GetSegment(aReflexR)
                                                       , GetSegment(aOpposite)
                                                       , GetSegment(aOppositePrev)
                                                       , GetSegment(aOppositeNext)
                                                       ) ;
  }
  
  Comparison_result CompareEvents ( EventPtr const& aA, EventPtr const& aB )
  {
    return mTraits.compare_event_times_object()( GetSegment(aA->border_a())
                                               , GetSegment(aA->border_b())
                                               , GetSegment(aA->border_c())
                                               , GetSegment(aB->border_a())
                                               , GetSegment(aB->border_b())
                                               , GetSegment(aB->border_c())
                                               ) ;
  }
  
  Comparison_result CompareEventsDistanceToSeed ( Vertex_handle   aSeed
                                                , EventPtr const& aA
                                                , EventPtr const& aB 
                                                )
  {
    if ( aSeed->is_inner() )
      return mTraits.compare_event_distance_to_seed_object()( GetSegment(GetDefiningBorder0(aSeed))
                                                            , GetSegment(GetDefiningBorder1(aSeed))
                                                            , GetSegment(GetDefiningBorder2(aSeed))
                                                            , GetSegment(aA->border_a())
                                                            , GetSegment(aA->border_b())
                                                            , GetSegment(aA->border_c())
                                                            , GetSegment(aB->border_a())
                                                            , GetSegment(aB->border_b())
                                                            , GetSegment(aB->border_c())
                                                            ) ;
    else
      return mTraits.compare_event_distance_to_seed_object()( aSeed->point()
                                                            , GetSegment(aA->border_a())
                                                            , GetSegment(aA->border_b())
                                                            , GetSegment(aA->border_c())
                                                            , GetSegment(aB->border_a())
                                                            , GetSegment(aB->border_b())
                                                            , GetSegment(aB->border_c())
                                                            ) ;
  }

  std::pair<Point_2,FT> ConstructEventPointAndTime( EventPtr const& aE )
  {
    return mTraits.construct_event_object()( GetSegment(aE->border_a())
                                           , GetSegment(aE->border_b())
                                           , GetSegment(aE->border_c())
                                           ); 
  }
     
  BorderTriple GetDefiningBorders( Vertex_handle aA, Vertex_handle aB ) ;

  bool AreBisectorsCoincident ( Halfedge_const_handle aA, Halfedge_const_handle aB ) const ;

  void CollectSplitEvent( Vertex_handle    aNode
                          ,Halfedge_handle  aReflexLBorder
                          ,Halfedge_handle  aReflexRBorder
                          ,Halfedge_handle  aOppositeBorder
                          ,EventPtr_Vector& aCandidates
                        ) ;

  EventPtr FindEdgeEvent( Vertex_handle aLNode, Vertex_handle aRNode ) ;

  EventPtr FindSplitEvent( Vertex_handle aNode ) ;

  EventPtr ChooseBestNewEvent( Vertex_handle aNode, EventPtr_Vector& aEvents ) ;

  EventPtr FindBestNewEvent( Vertex_handle aNode ) ;

  void HandleSimultaneousEdgeEvent( Vertex_handle aA, Vertex_handle aB ) ;

  void CollectNewEvents( Vertex_handle aNode, EventPtr_Vector& aEvents ) ;
  void AddNewEvent( Vertex_handle aNode ) ;
  void UpdatePQ( Vertex_handle aV ) ;
  void CreateInitialEvents();
  void CreateContourBisectors();
  void InitPhase();

  Vertex_handle LookupOnSLAV ( Halfedge_handle aOBorder, Event const& aEvent ) ;

  Vertex_handle ConstructEdgeEventNode( EdgeEvent& aEvent ) ;

  Vertex_handle_pair ConstructSplitEventNodes( SplitEvent& aEvent ) ;

  void HandleEdgeEvent       ( EdgeEvent&  aEvent ) ;
  void HandleSplitEvent      ( SplitEvent& aEvent ) ;
  void HandleMultiSplitEvent ( SplitEvent* aEvent ) ;

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
                                         , aHalfedge->is_contour_bisector() ? CGAL::GREEN : CGAL::BLUE
                                         , aHalfedge->is_contour_bisector() ? "CBisector" : "IBisector"
                                         ) ;
    draw_.Release();
  }

#endif

private:

  //Input
  Traits mTraits ;

  //Internal

  typedef std::vector<Halfedge_handle> HalfedgeCollection ;
  
  std::vector<Vertex_handle>   mReflexVertices ;
  std::vector<VertexWrapper>   mWrappedVertices ;
  std::vector<HalfedgeWrapper> mWrappedHalfedges ;
  HalfedgeCollection           mBorderHalfedges ;

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
