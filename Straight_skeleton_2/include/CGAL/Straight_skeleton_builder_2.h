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

#ifdef CGAL_STRAIGHT_SKELETON_STATS


#define CGAL_STSKEL_STATS_CODE(c) c
#else
#define CGAL_STSKEL_STATS_CODE(c)
#endif

CGAL_BEGIN_NAMESPACE

template<class SSkel_>
struct Dummy_straight_skeleton_builder_2_visitor
{
  typedef SSkel_ SSkel ;

  typedef typename SSkel::Halfedge_const_handle Halfedge_const_handle ;
  typedef typename SSkel::Vertex_const_handle   Vertex_const_handle ;

  void on_contour_edge_entered ( Halfedge_const_handle const& he ) const {}
                               
  void on_initialization_started( int size_of_vertices ) const {}
  
  void on_initial_events_collected( Vertex_const_handle const& v, bool is_reflex, bool is_degenerate )  const  {}
  
  void on_edge_event_created( Vertex_const_handle const& lnode
                            , Vertex_const_handle const& rnode
                            )  const {}

  void on_split_event_created( Vertex_const_handle const& node )  const {}

  void on_pseudo_split_event_created( Vertex_const_handle const& lnode
                                    , Vertex_const_handle const& rnode
                                    )  const {}
                                    
  void on_initialization_finished() const {}
  
  void on_propagation_started() const {}
  
  void on_anihiliation_event_processed ( Vertex_const_handle const& node0
                                       , Vertex_const_handle const& node1
                                       )  const  {}


  void on_edge_event_processed( Vertex_const_handle const& lseed
                              , Vertex_const_handle const& rseed
                              , Vertex_const_handle const& node
                              )  const {} 

  void on_split_event_processed( Vertex_const_handle const& seed
                               , Vertex_const_handle const& node0
                               , Vertex_const_handle const& node1
                               )  const {}

  void on_pseudo_split_event_processed( Vertex_const_handle const& lseed
                                      , Vertex_const_handle const& rseed
                                      , Vertex_const_handle const& node0
                                      , Vertex_const_handle const& node1
                                      )  const {}

  void on_vertex_processed( Vertex_const_handle const& node )  const {}
  
  void on_propagation_finished() const {}
  
  void on_cleanup_started( bool mergin_coincident_nodes ) const {}
  
  void on_cleanup_finished() const {}
  
  void on_algorithm_finished ( bool finished_ok ) const {}
  
  void on_error( char const* ) const {}
} ;


template<class Traits_, class SSkel_, class Visitor_ = Dummy_straight_skeleton_builder_2_visitor<SSkel_> >
class Straight_skeleton_builder_2
{
public:

  typedef Traits_ Traits ;
  typedef SSkel_  SSkel ;
  typedef Visitor_ Visitor ;

  typedef boost::shared_ptr<SSkel> SSkelPtr ;
  
private :

  typedef typename Traits::Kernel K ;
  
  typedef typename Traits::FT                FT ;
  typedef typename Traits::Point_2           Point_2 ;
  typedef typename Traits::Segment_2         Segment_2 ;
  typedef typename Traits::Triedge_2         Triedge_2 ;
  typedef typename Traits::Sorted_triedge_2  Sorted_triedge_2 ;
  
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

  typedef Straight_skeleton_builder_event_2             <SSkel,Traits> Event ;
  typedef Straight_skeleton_builder_edge_event_2        <SSkel,Traits> EdgeEvent ;
  typedef Straight_skeleton_builder_split_event_2       <SSkel,Traits> SplitEvent ;
  typedef Straight_skeleton_builder_pseudo_split_event_2<SSkel,Traits> PseudoSplitEvent ;

  typedef boost::intrusive_ptr<Event> EventPtr ;

  typedef std::vector<EventPtr>        EventPtr_Vector ;
  typedef std::vector<Halfedge_handle> Halfedge_handle_vector ;
  typedef std::vector<Vertex_handle>   Vertex_handle_vector ;

  typedef typename Halfedge_handle_vector::iterator Halfedge_handle_vector_iterator ;
  typedef typename Vertex_handle_vector  ::iterator Vertex_handle_vector_iterator ;
  typedef typename EventPtr_Vector       ::iterator event_iterator ;

  typedef boost::tuple<Halfedge_handle, Halfedge_handle, Halfedge_handle> BorderTriple ;

  typedef Straight_skeleton_builder_2<Traits,SSkel,Visitor> Self ;
  
  typedef typename Halfedge::Base_base HBase_base ;
  typedef typename Halfedge::Base      HBase ;
  typedef typename Vertex::Base        VBase ;
  typedef typename Face::Base          FBase ;
  
  struct Halfedge_ID_compare : std::binary_function<bool,Halfedge_handle,Halfedge_handle>
  {
    bool operator() ( Halfedge_handle const& aA, Halfedge_handle const& aB ) const
    {
      return aA->id() < aB->id() ;
    }
  } ;
  
  typedef std::map<Halfedge_handle,bool,Halfedge_ID_compare> Is_bond_map ;
    
  // Orders two halfedges pointing to a common vertex around it ccw
  struct Halfedge_compare_ccw : std::binary_function<bool,Halfedge_handle,Halfedge_handle>
  {
    bool operator() ( Halfedge_handle const& aA, Halfedge_handle const& aB ) const
    {
      Point_2 o = aA->vertex()->point();
      Point_2 a = aA->opposite()->vertex()->point();
      Point_2 b = aB->opposite()->vertex()->point();
      
      return K().compare_angle_with_x_axis_2_object()( K().construct_direction_2_object()( K().construct_vector_2_object()(a,o) )
                                                     , K().construct_direction_2_object()( K().construct_vector_2_object()(b,o) )
                                                     ) == SMALLER ;
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

  
  SSkelPtr construct_skeleton( bool aMergeCoincidentNodes = false ) ;

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

  struct VertexWrapper
  {
    VertexWrapper( Vertex_handle aVertex, Event_compare const& aComparer )
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

    Vertex_handle   mVertex ;
    bool            mIsReflex ;
    bool            mIsDegenerate ;
    bool            mIsProcessed ;
    bool            mIsExcluded ;
    int             mPrevInLAV ;
    int             mNextInLAV ;
    Halfedge_handle mDefiningBorderA ;
    Halfedge_handle mDefiningBorderB ;
    Halfedge_handle mDefiningBorderC ;
    bool            mNextSplitEventInMainPQ;
    PQ              mSplitEvents ;
  } ;
  
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

  inline Segment_2 CreateEdge ( Halfedge_const_handle aH ) const
  {
    Point_2 s = aH->opposite()->vertex()->point() ;
    Point_2 t = aH->vertex()->point() ;
    return K().construct_segment_2_object()(s,t);
  }

  Triedge_2 CreateTriedge ( Halfedge_const_handle aE0
                          , Halfedge_const_handle aE1
                          , Halfedge_const_handle aE2
                          ) const
  {
    return Construct_ss_triedge_2(mTraits)(CreateEdge(aE0),CreateEdge(aE1),CreateEdge(aE2));
  }
  
  Sorted_triedge_2 CreateSortedTriedge ( Triedge_2 const& aTriedge, Triedge_collinearity aC ) const
  {
    return Construct_ss_sorted_triedge_2(mTraits)(aTriedge,aC);
  }

  Triedge_collinearity GetCollinearity ( Triedge_2 const& aTriedge ) const
  {
    return Get_ss_triedge_collinearity_2(mTraits)(aTriedge);
  }  
  
  Triedge_2 CreateTriedge ( BorderTriple const& aTriple ) const
  {
    Halfedge_handle lE0, lE1, lE2 ;
    boost::tie(lE0,lE1,lE2) = aTriple ;
    return Construct_ss_triedge_2(mTraits)(CreateEdge(lE0),CreateEdge(lE1),CreateEdge(lE2));
  }
  
  
  Sorted_triedge_2 CreateSortedTriedge ( Halfedge_const_handle aE0
                                       , Halfedge_const_handle aE1
                                       , Halfedge_const_handle aE2
                                       ) const
  {
    Triedge_2 lTriedge = CreateTriedge(aE0,aE1,aE2);
    return CreateSortedTriedge(lTriedge, GetCollinearity(lTriedge) );
  }

  Sorted_triedge_2 CreateSortedTriedge( BorderTriple const& aTriple ) const
  {
    Triedge_2 lTriedge = CreateTriedge(aTriple);
    return CreateSortedTriedge(lTriedge, GetCollinearity(lTriedge) );
  }
    
  Vertex_handle GetVertex ( int aIdx )
  {
    CGAL_precondition(aIdx>=0);
    return mWrappedVertices[aIdx].mVertex ;
  }

  Vertex_handle GetPrevInLAV ( Vertex_handle aV )
  {
    return GetVertex ( mWrappedVertices[aV->id()].mPrevInLAV ) ;
  }

  Vertex_handle GetNextInLAV ( Vertex_handle aV )
  {
    return GetVertex ( mWrappedVertices[aV->id()].mNextInLAV ) ;
  }

  void SetPrevInLAV ( Vertex_handle aV, Vertex_handle aPrev )
  {
    mWrappedVertices[aV->id()].mPrevInLAV = aPrev->id();
  }

  void SetNextInLAV ( Vertex_handle aV, Vertex_handle aPrev )
  {
    mWrappedVertices[aV->id()].mNextInLAV = aPrev->id();
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
  bool IsExcluded ( Vertex_const_handle aVertex ) const
  {
    return mWrappedVertices[aVertex->id()].mIsExcluded ;
  }

  void SetIsReflex ( Vertex_handle aVertex )
  {
    mWrappedVertices[aVertex->id()].mIsReflex = true ;
  }

  bool IsReflex ( Vertex_handle aVertex )
  {
    return mWrappedVertices[aVertex->id()].mIsReflex ;
  }

  void SetIsDegenerate ( Vertex_handle aVertex )
  {
    mWrappedVertices[aVertex->id()].mIsDegenerate = true ;
  }

  bool IsDegenerate ( Vertex_handle aVertex )
  {
    return mWrappedVertices[aVertex->id()].mIsDegenerate ;
  }
  
  void SetIsProcessed ( Vertex_handle aVertex )
  {
    mWrappedVertices[aVertex->id()].mIsProcessed = true ;

    mVisitor.on_vertex_processed(aVertex);
  }

  bool IsProcessed ( Vertex_handle aVertex )
  {
    return mWrappedVertices[aVertex->id()].mIsProcessed ;
  }

  void AddSplitEvent ( Vertex_handle aVertex, EventPtr aEvent )
  {
    CGAL_STSKEL_BUILDER_TRACE(2, "V" << aVertex->id() << " PQ: " << *aEvent);
    mWrappedVertices[aVertex->id()].mSplitEvents.push(aEvent);
  }
  
  EventPtr PopNextSplitEvent ( Vertex_handle aVertex )
  {
    EventPtr rEvent ;
    VertexWrapper& lW = mWrappedVertices[aVertex->id()] ;
    if ( !lW.mNextSplitEventInMainPQ )
    {
      PQ& lPQ = lW.mSplitEvents ;
      if ( !lPQ.empty() )
      {
        rEvent = lPQ.top(); 
        lPQ.pop();
        lW.mNextSplitEventInMainPQ = true ;
      }
    }
    return rEvent ;
  }

  void AllowNextSplitEvent ( Vertex_handle aVertex )
  {
    mWrappedVertices[aVertex->id()].mNextSplitEventInMainPQ = false ;
  }  
  
  void InsertEventInPQ( EventPtr aEvent ) ;

  EventPtr PopEventFromPQ() ;

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
  
  // Returns true if the intersection of the sets (aXA,aXB,aXC) and (aYA,aYB,aYC) has size exactly 2
  // (that is, both sets have 2 elements in common)
  bool HaveTwoInCommon( BorderTriple aX, BorderTriple aY ) const
  {
    Halfedge_handle lXA, lXB, lXC, lYA, lYB, lYC ;
    boost::tie(lXA,lXB,lXC) = aX ;
    boost::tie(lYA,lYB,lYC) = aY ;
    return CountInCommon(lXA,lXB,lXC,lYA,lYB,lYC) == 2 ;
  }

  // Returns true if the sets of halfedges (aXA,aXB,aXC) and (aYA,aYB,aYC) are equivalent
  // (one is a permutation of the other)
  bool AreTheSameTriple( Halfedge_handle aXA, Halfedge_handle aXB, Halfedge_handle aXC
                       , Halfedge_handle aYA, Halfedge_handle aYB, Halfedge_handle aYC
                       ) const
  {
    return CountInCommon(aXA,aXB,aXC,aYA,aYB,aYC) == 3 ;
  }


  bool ExistEvent ( Sorted_triedge_2 aS )
  {
    return Do_ss_event_exist_2(mTraits)(aS);
  }  
  
  bool IsOppositeEdgeFacingTheSplitSeed( Vertex_const_handle aSeed, Halfedge_const_handle aOpposite ) const
  {
    if ( aSeed->is_skeleton() )
    {
      BorderTriple lTriple  = GetSkeletonVertexDefiningBorders(aSeed);
      Triedge_2    lTriedge = CreateTriedge(lTriple);

      return Is_edge_facing_ss_seed_2(mTraits)( CreateSortedTriedge(lTriedge,GetCollinearity(lTriedge))
                                              , CreateEdge(aOpposite)
                                              ) ;
    }
    else
    {
      return Is_edge_facing_ss_seed_2(mTraits)( aSeed->point(), CreateEdge(aOpposite) ) ;
    }
  }
  
  bool IsSplitEventInsideOffsetZone( EventPtr const& aSplit
                                   , Halfedge_const_handle aOppositePrev
                                   , Halfedge_const_handle aOpposite
                                   , Halfedge_const_handle aOppositeNext
                                   ) const
  {
    return Is_ss_event_inside_offset_zone_2(mTraits)( aSplit->sorted_triedge()
                                                    , CreateTriedge(aOppositePrev,aOpposite, aOppositeNext)
                                                    ) ;
  }

  Comparison_result CompareEvents ( Sorted_triedge_2 const& aA, Sorted_triedge_2 const& aB ) const
  {
  
    return Compare_ss_event_times_2(mTraits)(aA,aB) ;
  }

  Comparison_result CompareEvents ( EventPtr const& aA, EventPtr const& aB ) const
  {
    if ( !AreTheSameTriple( aA->border_a(), aA->border_b(), aA->border_c()
                          , aB->border_a(), aB->border_b(), aB->border_c()
                          )
        )
    {
      return CompareEvents( aA->sorted_triedge(), aB->sorted_triedge() ) ;
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
        BorderTriple lTriple  = GetSkeletonVertexDefiningBorders(aSeed);
        Triedge_2    lTriedge = CreateTriedge(lTriple);

        CGAL_STSKEL_BUILDER_TRACE(3
                                 ,"Seed N" << aSeed->id() << " is a skeleton node,"
                                 << " defined by: E" << lTriple.get<0>()->id()
                                            << ", E" << lTriple.get<1>()->id()
                                             << ", E" << lTriple.get<2>()->id()
                                 );

        return Compare_ss_event_distance_to_seed_2(mTraits)( CreateSortedTriedge(lTriedge,GetCollinearity(lTriedge))
                                                           , aA->sorted_triedge()
                                                           , aB->sorted_triedge()
                                                           ) ;
      }
      else
      {
        return Compare_ss_event_distance_to_seed_2(mTraits)( aSeed->point(), aA->sorted_triedge(), aB->sorted_triedge() ) ;
      }
    }
    else return EQUAL ;
  }
  
  bool AreEventsSimultaneous( Sorted_triedge_2 const& x, Sorted_triedge_2 const& y ) const
  {
    return Are_ss_events_simultaneous_2(mTraits)(x,y) ;
  }
  
  bool AreSkeletonNodesCoincident( Vertex_handle aX, Vertex_handle aY ) const
  {
    BorderTriple lBordersX = GetSkeletonVertexDefiningBorders(aX);
    BorderTriple lBordersY = GetSkeletonVertexDefiningBorders(aY);
    return Are_ss_events_simultaneous_2(mTraits)( CreateSortedTriedge(lBordersX), CreateSortedTriedge(lBordersY)) ;
  }
 
  bool IsNewEventInThePast( Halfedge_handle         aBorderA
                          , Halfedge_handle         aBorderB
                          , Halfedge_handle         aBorderC
                          , Sorted_triedge_2 const& aSortedTriedge1 // (aBorderA,aBorderB,aBorderC)
                          , Vertex_handle           aSeedNode
                          ) const
  {
    bool rResult = false ;

    Halfedge_handle lSeedBorderA, lSeedBorderB, lSeedBorderC ;

    boost::tie(lSeedBorderA,lSeedBorderB,lSeedBorderC) = GetSkeletonVertexDefiningBorders(aSeedNode) ;

    if ( !AreTheSameTriple(aBorderA,aBorderB,aBorderC,lSeedBorderA,lSeedBorderB,lSeedBorderC) )
    {
      Triedge_2        lTriedge2       = CreateTriedge(lSeedBorderA,lSeedBorderB,lSeedBorderC);
      Sorted_triedge_2 lSortedTriedge2 = CreateSortedTriedge(lTriedge2,GetCollinearity(lTriedge2));
      
      if ( CompareEvents( aSortedTriedge1, lSortedTriedge2 ) == SMALLER )
        rResult = true ;
    }

    return rResult ;
  }

  boost::tuple<FT,Point_2> ConstructEventTimeAndPoint( Sorted_triedge_2 const& aS ) const
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
    boost::tie(lTime,lP) = ConstructEventTimeAndPoint(aE.sorted_triedge());
    
    aE.SetTimeAndPoint(lTime,lP);
  }

  void EraseBisector( Halfedge_handle aB )
  {
    CGAL_STSKEL_BUILDER_TRACE(1,"Dangling B" << aB->id() << " and B" << aB->opposite()->id() << " removed.");
    
    mSSkel->SSkel::Base::edges_erase(aB);
  }

  BorderTriple GetDefiningBorders( Vertex_handle aA, Vertex_handle aB ) ;

  bool AreBisectorsCoincident ( Halfedge_const_handle aA, Halfedge_const_handle aB ) const ;

  bool IsInverseSplitEventCoincident( Vertex_handle    const& aReflexOppN 
                                    , Halfedge_handle  const& aReflexLBorder
                                    , Halfedge_handle  const& aReflexRBorder
                                    , Sorted_triedge_2 const& aEventSTriedge
                                    ) ;
                                    
  EventPtr IsPseudoSplitEvent( EventPtr const& aEvent, Vertex_handle aOppN ) ;
  
  void CollectSplitEvent( Vertex_handle    aNode
                        , Halfedge_handle  aReflexLBorder
                        , Halfedge_handle  aReflexRBorder
                        , Halfedge_handle  aOppositeBorder
                        ) ;

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

  void ClassifyBisectorsAroundMultinode( Vertex_handle const&         v0
                                       , Vertex_handle_vector const& aCluster
                                       , Is_bond_map&                rIsBond
                                       ) ;
                                       
  void ClassifyBisectorsAroundMultinode( Vertex_handle_vector const& aCluster
                                       , Is_bond_map&                rIsBond
                                       ) ;
  
  void RearrangeBisectorsAroundMultinode( Vertex_handle const& v0, Is_bond_map& rIsBond ) ;
  
  bool AreNodesConnected( Vertex_handle v, Vertex_handle u ) ;
  
  void MergeCoincidentNodes() ;
  
  bool FinishUp( bool aMergeCoincidentNodes );

  bool Run( bool aMergeCoincidentNodes );

private:

  //Input
  Traits  mTraits ;

  Visitor const& mVisitor ;

  std::vector<VertexWrapper> mWrappedVertices ;
  Vertex_handle_vector       mReflexVertices ;
  Halfedge_handle_vector     mDanglingBisectors ;
  Halfedge_handle_vector     mContourHalfedges ;

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
      mWrappedVertices.push_back( VertexWrapper(lVertex,mEventCompare) ) ;

      Face_handle lFace = mSSkel->SSkel::Base::faces_push_back( Face() ) ;

      ++ c ;

      lCCWBorder->HBase_base::set_face(lFace);
      lFace     ->FBase::set_halfedge(lCCWBorder);

      lVertex   ->VBase::set_halfedge(lCCWBorder);
      lCCWBorder->HBase_base::set_vertex(lVertex);

      if ( lCurr == aBegin )
      {
        lFirstVertex    = lVertex ;
        lFirstCCWBorder = lCCWBorder ;
      }
      else
      {
        SetPrevInLAV    (lVertex    ,lPrevVertex);
        SetNextInLAV    (lPrevVertex,lVertex    );

        SetDefiningBorderA(lVertex    ,lCCWBorder);
        SetDefiningBorderB(lPrevVertex,lCCWBorder);

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

    SetPrevInLAV    (lFirstVertex,lPrevVertex );
    SetNextInLAV    (lPrevVertex ,lFirstVertex);

    SetDefiningBorderA(lFirstVertex,lFirstCCWBorder);
    SetDefiningBorderB(lPrevVertex ,lFirstCCWBorder);

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
