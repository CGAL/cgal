// Copyright (c) 2006 Fernando Luis Cacciola Carballal. All rights reserved.
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
#include <exception>
#include <string>
#include <map>

#include <boost/tuple/tuple.hpp>
#include <boost/intrusive_ptr.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/scoped_ptr.hpp>

#include <CGAL/algorithm.h>
#include <CGAL/Straight_skeleton_2/Straight_skeleton_aux.h>
#include <CGAL/Straight_skeleton_2/Straight_skeleton_builder_events_2.h>
#include <CGAL/Straight_skeleton_2.h>
#include <CGAL/Straight_skeleton_builder_traits_2.h>
#include <CGAL/HalfedgeDS_const_decorator.h>
#include <CGAL/enum.h>

CGAL_BEGIN_NAMESPACE

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

  void on_split_event_created( Vertex_const_handle const& node )  const {}

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
  
  typedef typename Traits::FT                  FT ;
  typedef typename Traits::Point_2             Point_2 ;
  typedef typename Traits::Segment_2           Segment_2 ;
  typedef typename Traits::Trisegment_2        Trisegment_2 ;
  typedef typename Traits::Seeded_trisegment_2 Seeded_trisegment_2 ;
  
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

  typedef std::vector<EventPtr>        EventPtr_Vector ;
  typedef std::vector<Halfedge_handle> Halfedge_handle_vector ;
  typedef std::vector<Vertex_handle>   Vertex_handle_vector ;

  typedef typename Halfedge_handle_vector::iterator Halfedge_handle_vector_iterator ;
  typedef typename Vertex_handle_vector  ::iterator Vertex_handle_vector_iterator ;
  typedef typename EventPtr_Vector       ::iterator event_iterator ;

  typedef Straight_skeleton_builder_2<Traits,SSkel,Visitor> Self ;
  
  typedef typename Halfedge::Base_base HBase_base ;
  typedef typename Halfedge::Base      HBase ;
  typedef typename Vertex::Base        VBase ;
  typedef typename Face::Base          FBase ;

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

  struct Halfedge_ID_compare : std::binary_function<bool,Halfedge_handle,Halfedge_handle>
  {
    bool operator() ( Halfedge_handle const& aA, Halfedge_handle const& aB ) const
    {
      return aA->id() < aB->id() ;
    }
  } ;
  
public:

  struct straight_skeleton_exception : std::runtime_error
  {
    straight_skeleton_exception( std::string what ) : std::runtime_error(what) {}  
  } ;
  
  Straight_skeleton_builder_2 ( Traits const& = Traits(), Visitor const& aVisitor = Visitor() ) ;

  // NOTE: The following public method is implemented here in this header file to support some broken compilers.
  // But it's right at the end of the class declaration becuause it needs all of the class.
  //
  // template<class InputPointIterator> Straight_skeleton_builder_2& enter_contour ( InputPointIterator aBegin, InputPointIterator aEnd ) ;

  
  SSkelPtr construct_skeleton() ;

private :

  void throw_error ( char const* what ) const ;

  
  class Event_compare : public std::binary_function<bool,EventPtr,EventPtr>
  {
  public:

    Event_compare ( Self const* aBuilder ) : mBuilder(aBuilder) {}

    bool operator() ( EventPtr const& aA, EventPtr const& aB ) const
    {
      return mBuilder->CompareEvents(aA,aB) == LARGER ;
    }

  private:

    Self const* mBuilder ;
  } ;

  typedef std::priority_queue<EventPtr,std::vector<EventPtr>,Event_compare> PQ ;

  typedef std::pair<Vertex_handle,Vertex_handle> Vertex_handle_pair ;

  typedef std::vector<Vertex_handle_pair> SplitNodesVector ;

  struct Vertex_data : public Ref_counted_base
  {
    Vertex_data ( Vertex_handle aVertex, Event_compare const& aComparer )
      :
        mVertex(aVertex)
      , mIsReflex(false)
      , mIsDegenerate(false)
      , mIsProcessed(false)
      , mIsExcluded(false)
      , mPrevInLAV(-1)
      , mNextInLAV(-1)
      , mNextSplitEventInMainPQ(false)
      , mSplitEvents(aComparer)
    {}

    Vertex_handle        mVertex ;
    bool                 mIsReflex ;
    bool                 mIsDegenerate ;
    bool                 mIsProcessed ;
    bool                 mIsExcluded ;
    int                  mPrevInLAV ;
    int                  mNextInLAV ;
    bool                 mNextSplitEventInMainPQ;
    PQ                   mSplitEvents ;
    Triedge              mTriedge ; 
    Seeded_trisegment_2  mSTrisegment ; // Skeleton nodes cache the seeded trisegment that defines the originating event
  } ;
  
  typedef boost::intrusive_ptr<Vertex_data> Vertex_data_ptr ;
  
private :

  inline Halfedge_handle validate( Halfedge_handle aH ) const
  {
    if ( !handle_assigned(aH) )
      throw_error("Unassigned halfedge handle") ;
    return aH ;
  }
  
  inline Vertex_handle validate( Vertex_handle aH ) const
  {
    if ( !handle_assigned(aH) )
      throw_error("Unassigned vertex handle") ;
    return aH ;
  }

  void InitVertexData( Vertex_handle aV )
  {
    mVertexData.push_back( Vertex_data_ptr( new Vertex_data(aV,mEventCompare) ) ) ;
  } 
  
  inline Vertex_data const& GetVertexData( Vertex_const_handle aV ) const { return *mVertexData[aV->id()]; }
  inline Vertex_data&       GetVertexData( Vertex_const_handle aV )       { return *mVertexData[aV->id()]; }
  
  Vertex_handle GetVertex ( int aIdx )
  {
    CGAL_precondition(aIdx>=0);
    return mVertexData[aIdx]->mVertex ;
  }

  inline Triedge const& GetTriedge ( Vertex_const_handle aV ) const
  {
    return GetVertexData(aV).mTriedge ;
  }
  inline void SetTriedge ( Vertex_handle aV, Triedge aTriedge )
  {
    GetVertexData(aV).mTriedge = aTriedge ;
  }

  // Contour nodes do not have a trisegment because they don't come from an event
  inline Seeded_trisegment_2 const& GetSeededTrisegment ( Vertex_handle aV ) const
  {
    return GetVertexData(aV).mSTrisegment ;
  }
  
  inline void SetSeededTrisegment ( Vertex_handle aV, Seeded_trisegment_2 const& aSTrisegment )
  {
    GetVertexData(aV).mSTrisegment = aSTrisegment ;
  }
  
  inline Segment_2 CreateSegment ( Halfedge_const_handle aH ) const
  {
    Point_2 s = aH->opposite()->vertex()->point() ;
    Point_2 t = aH->vertex()->point() ;
    return K().construct_segment_2_object()(s,t);
  }

  Trisegment_2 CreateTrisegment ( Triedge const& aTriedge ) const
  {
    boost::optional<Trisegment_2> r = Construct_ss_trisegment_2(mTraits)(CreateSegment(aTriedge.e0())
                                                                        ,CreateSegment(aTriedge.e1())
                                                                        ,CreateSegment(aTriedge.e2())
                                                                        );
                                                                        
    CGAL_STSKEL_BUILDER_TRACE(4,"Trisegment for " << aTriedge << "=" << ( !!r ? Trisegment_2() : *r ) ) ;
    
    if ( !r )
      throw_error("Unable to determine edges collinearity");
      
    return *r ;  
  }
  
  Seeded_trisegment_2 CreateSeededTrisegment ( Triedge const& aTriedge
                                             , Vertex_handle  aLSeed
                                             , Vertex_handle  aRSeed
                                             ) const
  {
    return Construct_ss_seeded_trisegment_2(mTraits)(CreateTrisegment(aTriedge)
                                                    ,handle_assigned(aLSeed) ? GetSeededTrisegment(aLSeed).event() : Trisegment_2::null()
                                                    ,handle_assigned(aRSeed) ? GetSeededTrisegment(aRSeed).event() : Trisegment_2::null()
                                                    );
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

  bool IsProcessed ( Vertex_handle aV )
  {
    return GetVertexData(aV).mIsProcessed ;
  }

  void AddSplitEvent ( Vertex_handle aV, EventPtr aEvent )
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
      PQ& lPQ = lData.mSplitEvents ;
      if ( !lPQ.empty() )
      {
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



  bool ExistEvent ( Seeded_trisegment_2 const& aS )
  {
    return Do_ss_event_exist_2(mTraits)(aS);
  }  
  
  bool IsOppositeEdgeFacingTheSplitSeed( Vertex_handle aSeed, Halfedge_handle aOpposite ) const
  {
    if ( aSeed->is_skeleton() )
         return Is_edge_facing_ss_node_2(mTraits)( GetSeededTrisegment(aSeed), CreateSegment(aOpposite) ) ;
    else return Is_edge_facing_ss_node_2(mTraits)( aSeed->point()            , CreateSegment(aOpposite) ) ;
  }
  
  bool IsSplitEventInsideOffsetZone( EventPtr const& aSplit
                                   , Triedge const&  aOppTriedge
                                   , Vertex_handle   aOppLSeed
                                   , Vertex_handle   aOppRSeed
                                   ) const
  {
    return Is_ss_event_inside_offset_zone_2(mTraits)( aSplit->strisegment()
                                                    , CreateSeededTrisegment(aOppTriedge, aOppLSeed, aOppRSeed )
                                                    ) ;
  }

  Comparison_result CompareEvents ( Seeded_trisegment_2 const& aA, Seeded_trisegment_2 const& aB ) const
  {
    return Compare_ss_event_times_2(mTraits)(aA,aB) ;
  }

  Comparison_result CompareEvents ( EventPtr const& aA, EventPtr const& aB ) const
  {
    return aA->triedge() != aB->triedge() ? CompareEvents( aA->strisegment(), aB->strisegment() ) : EQUAL ;
  }

  Comparison_result CompareEventsDistanceToSeed ( Vertex_handle   aSeed
                                                , EventPtr const& aA
                                                , EventPtr const& aB
                                                ) const
  {
    if ( aA->triedge() != aB->triedge() )
    {
      if ( aSeed->is_skeleton() )
           return Compare_ss_event_distance_to_seed_2(mTraits)( GetSeededTrisegment(aSeed), aA->strisegment(), aB->strisegment() ) ;
      else return Compare_ss_event_distance_to_seed_2(mTraits)( aSeed->point()            , aA->strisegment(), aB->strisegment() ) ;
    }
    else return EQUAL ;
  }
  
  bool AreEventsSimultaneous( Seeded_trisegment_2 const& x, Seeded_trisegment_2 const& y ) const
  {
    return Are_ss_events_simultaneous_2(mTraits)(x,y) ;
  }
  
  bool AreSkeletonNodesCoincident( Vertex_handle aX, Vertex_handle aY ) const
  {
    return AreEventsSimultaneous( GetSeededTrisegment(aX),  GetSeededTrisegment(aY) ) ;
  }
 
  bool IsNewEventInThePast( Seeded_trisegment_2 const& aTrisegment, Vertex_handle aSeedNode ) const
  {
    return aSeedNode->is_skeleton() ? CompareEvents( aTrisegment, GetSeededTrisegment(aSeedNode) ) == SMALLER 
                                    : false  ;
  }

  boost::tuple<FT,Point_2> ConstructEventTimeAndPoint( Seeded_trisegment_2 const& aS ) const
  {
    boost::optional< boost::tuple<FT,Point_2> > r = Construct_ss_event_time_and_point_2(mTraits)(aS);
    if ( !r )
      throw_error("Unable to compute skeleton node coordinates");
    return *r ;  
  }

  void SetEventTimeAndPoint( Event& aE )
  {
    FT      lTime ;
    Point_2 lP ;
    boost::tie(lTime,lP) = ConstructEventTimeAndPoint(aE.strisegment());
    
    aE.SetTimeAndPoint(lTime,lP);
  }

  void EraseBisector( Halfedge_handle aB )
  {
    CGAL_STSKEL_BUILDER_TRACE(1,"Dangling B" << aB->id() << " and B" << aB->opposite()->id() << " removed.");
    
    mSSkel->SSkel::Base::edges_erase(aB);
  }

  Triedge GetCommonTriedge( Vertex_handle aA, Vertex_handle aB ) ;

  bool AreBisectorsCoincident ( Halfedge_const_handle aA, Halfedge_const_handle aB ) const ;

  bool IsInverseSplitEventCoincident( Vertex_handle const&       aReflexOppN 
                                    , Triedge const&             aEventTriedge
                                    , Seeded_trisegment_2 const& aEventTrisegment
                                    ) ;
                                    
  EventPtr IsPseudoSplitEvent( EventPtr const& aEvent, Vertex_handle aOppN ) ;
  
  void CollectSplitEvent( Vertex_handle aNode, Triedge const& aTriedge ) ;

  void CollectSplitEvents( Vertex_handle aNode ) ;

  EventPtr FindEdgeEvent( Vertex_handle aLNode, Vertex_handle aRNode ) ;

  void HandleSimultaneousEdgeEvent( Vertex_handle aA, Vertex_handle aB ) ;

  void CollectNewEvents( Vertex_handle aNode ) ;
  void UpdatePQ( Vertex_handle aV ) ;
  void CreateInitialEvents();
  void CreateContourBisectors();
  void InitPhase();

  void SetupPseudoSplitEventNode( Vertex_handle   aNode
                                , Halfedge_handle aDefiningBorderA
                                , Halfedge_handle aDefiningBorderB
                                );

  Vertex_handle LookupOnSLAV ( Halfedge_handle aOBorder, EventPtr const& aEvent ) ;

  Vertex_handle      ConstructEdgeEventNode         ( EdgeEvent&   aEvent ) ;
  Vertex_handle_pair ConstructSplitEventNodes       ( SplitEvent&  aEvent, Vertex_handle aOppR ) ;
  Vertex_handle_pair ConstructPseudoSplitEventNodes ( PseudoSplitEvent& aEvent ) ;

  void HandleEdgeEvent               ( EventPtr aEvent ) ;
  void HandleSplitEvent              ( EventPtr aEvent, Vertex_handle aOppR ) ;
  void HandlePseudoSplitEvent        ( EventPtr aEvent ) ;
  void HandleSplitOrPseudoSplitEvent ( EventPtr aEvent ) ;
  
  void InsertNextSplitEventInPQ( Vertex_handle v ) ;
  void InsertNextSplitEventsInPQ() ;
  
  bool IsProcessed( EventPtr aEvent ) ;

  void Propagate();

  void MergeSplitNodes ( Vertex_handle_pair aSplitNodes ) ;
  
  void RelinkBisectorsAroundMultinode( Vertex_handle const& v0, Halfedge_handle_vector& aLinks ) ;
  
  void PreprocessMultinode( Multinode& aMN ) ;
                          
  void ProcessMultinode( Multinode&              aMN 
                       , Halfedge_handle_vector& rHalfedgesToRemove 
                       , Vertex_handle_vector&   rVerticesToRemove                                                             
                       ) ;

  MultinodePtr CreateMultinode( Halfedge_handle begin, Halfedge_handle end ) ;
  
  void MergeCoincidentNodes() ;
  
  bool FinishUp();

  bool Run();

private:

  //Input
  Traits  mTraits ;

  Visitor const& mVisitor ;

  std::vector<Vertex_data_ptr> mVertexData ;
  
  Vertex_handle_vector   mReflexVertices ;
  Halfedge_handle_vector mDanglingBisectors ;
  Halfedge_handle_vector mContourHalfedges ;

  std::list<Vertex_handle> mSLAV ;

  SplitNodesVector mSplitNodes ;

  Event_compare mEventCompare ;

  int mVertexID ;
  int mEdgeID   ;
  int mEventID  ;
  int mStepID   ;

  PQ mPQ ;

  //Output
  SSkelPtr mSSkel ;

private :  

  template<class InputPointIterator>
  void enter_valid_contour ( InputPointIterator aBegin, InputPointIterator aEnd )
  {
    CGAL_STSKEL_BUILDER_TRACE(0,"Inserting Connected Component of the Boundary....");

    Halfedge_handle lFirstCCWBorder ;
    Halfedge_handle lPrevCCWBorder ;
    Halfedge_handle lNextCWBorder ;
    Vertex_handle   lFirstVertex ;
    Vertex_handle   lPrevVertex ;

    InputPointIterator lCurr = aBegin ;
    InputPointIterator lPrev = aBegin ;
    
    int c = 0 ;
    
    while ( lCurr != aEnd )
    {
      Halfedge_handle lCCWBorder = mSSkel->SSkel::Base::edges_push_back ( Halfedge(mEdgeID),Halfedge(mEdgeID+1)  );
      Halfedge_handle lCWBorder = lCCWBorder->opposite();
      mEdgeID += 2 ;

      mContourHalfedges.push_back(lCCWBorder);

      Vertex_handle lVertex = mSSkel->SSkel::Base::vertices_push_back( Vertex(mVertexID++,*lCurr) ) ;
      CGAL_STSKEL_BUILDER_TRACE(1,"Vertex: V" << lVertex->id() << " at " << lVertex->point() );
      InitVertexData(lVertex);

      Face_handle lFace = mSSkel->SSkel::Base::faces_push_back( Face() ) ;

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

        SetTriedge( lPrevVertex, Triedge(lPrevCCWBorder,lCCWBorder) ) ;
        
        //SetDefiningBorderA(lVertex    ,lCCWBorder);
        //SetDefiningBorderB(lPrevVertex,lCCWBorder);

        lCWBorder->HBase_base::set_vertex(lPrevVertex);

        lCCWBorder    ->HBase_base::set_prev(lPrevCCWBorder);
        lPrevCCWBorder->HBase_base::set_next(lCCWBorder);

        lNextCWBorder->HBase_base::set_prev(lCWBorder);
        lCWBorder    ->HBase_base::set_next(lNextCWBorder);

        CGAL_STSKEL_BUILDER_TRACE(1,"CCW Border: E" << lCCWBorder->id() << ' ' << lPrevVertex->point() << " -> " << lVertex    ->point());
        CGAL_STSKEL_BUILDER_TRACE(1,"CW  Border: E" << lCWBorder ->id() << ' ' << lVertex    ->point() << " -> " << lPrevVertex->point() );
        
        mVisitor.on_contour_edge_entered(lCCWBorder);

      }

      lPrev = lCurr ;
      
      ++ lCurr ;

      lPrevVertex    = lVertex ;
      lPrevCCWBorder = lCCWBorder ;
      lNextCWBorder  = lCWBorder ;
    }

    SetPrevInLAV(lFirstVertex,lPrevVertex );
    SetNextInLAV(lPrevVertex ,lFirstVertex);

    //SetDefiningBorderA(lFirstVertex,lFirstCCWBorder);
    //SetDefiningBorderB(lPrevVertex ,lFirstCCWBorder);
    SetTriedge( lPrevVertex, Triedge(lPrevCCWBorder,lFirstCCWBorder) ) ;

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
  
public:
  
  // This compares INPUT vertices so the comparison must be exact and there is no need to filter it.
  struct AreVerticesEqual
  {
    bool operator() ( Point_2 const&x, Point_2 const& y ) const
    {
      return CGAL::compare_xy(x,y) == EQUAL ;
    }
  } ; 
  
  template<class InputPointIterator>
  Straight_skeleton_builder_2& enter_contour ( InputPointIterator aBegin, InputPointIterator aEnd, bool aCheckValidity = true )
  {
    if ( aCheckValidity )
    {
      typedef std::vector<Point_2> Point_vector ;
      typedef typename Point_vector::iterator Point_iterator ;
     
      // Remove coincident consecutive vertices
      Point_vector lList;
      std::unique_copy(aBegin,aEnd,std::back_inserter(lList),AreVerticesEqual());

      while ( lList.size() > 0 && CGAL::compare_xy(lList.front(),lList.back()) == EQUAL ) 
        lList.pop_back();
      
      if ( lList.size() >= 3 )
      {
        enter_valid_contour(lList.begin(),lList.end());
      }
      else
      {
        CGAL_STSKEL_BUILDER_TRACE(0,"Degenerate contour (less than 3 non-degenerate vertices).");
      }
    }
    else
    {
      enter_valid_contour(aBegin,aEnd);
    } 

    return *this ;
  }

} ;

CGAL_END_NAMESPACE

#include <CGAL/Straight_skeleton_2/Straight_skeleton_builder_2_impl.h>


#endif // CGAL_STRAIGHT_SKELETON_BUILDER_2_H //
// EOF //
