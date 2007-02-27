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
#ifndef CGAL_POLYGON_OFFSET_BUILDER_2_H
#define CGAL_POLYGON_OFFSET_BUILDER_2_H 1

#include <vector>
#include <algorithm>

#include <boost/shared_ptr.hpp>
#include <boost/optional/optional.hpp>

#include <CGAL/Straight_skeleton_2/Straight_skeleton_aux.h>

#include <CGAL/Polygon_offset_builder_traits_2.h>

CGAL_BEGIN_NAMESPACE

template<class Ss_, class Traits_, class Container_>
class Polygon_offset_builder_2
{
public :

  typedef Ss_        Ss ;
  typedef Traits_    Traits ;
  typedef Container_ Container ;

  typedef typename Traits::Kernel K ;
  
  typedef typename Traits::FT      FT ;
  typedef typename Traits::Point_2 Point_2 ;

  typedef boost::shared_ptr<Container> ContainerPtr ;

  Polygon_offset_builder_2( Ss const& aSs, Traits const& aTraits = Traits() )  ;

  template<class OutputIterator>
  OutputIterator construct_offset_contours( FT aTime, OutputIterator aOut ) ;

private:

  typedef typename Ss::Halfedge_const_handle Halfedge_const_handle  ;
  typedef typename Ss::Vertex_const_handle   Vertex_const_handle  ;

  typedef std::vector<Halfedge_const_handle> Halfedge_vector ;

  typedef typename Traits::Segment_2           Segment_2 ;
  typedef typename Traits::Trisegment_2        Trisegment_2 ;
  typedef typename Traits::Seeded_trisegment_2 Seeded_trisegment_2 ;

  typedef CGAL_SS_i::Triedge<Halfedge_const_handle> Triedge ;
    
  bool handled_assigned( Halfedge_const_handle aH ) const
  {
    const Halfedge_const_handle cNull ;
    return aH != cNull ;
  }

  Halfedge_const_handle LocateHook( FT aTime, Halfedge_const_handle aBisector ) ;

  void AddOffsetVertex( FT aTime, Halfedge_const_handle aHook, ContainerPtr aPoly ) ;

  template<class OutputIterator>
  OutputIterator TraceOffsetPolygon( FT aTime, Halfedge_const_handle aHook, OutputIterator aOut ) ;

  Halfedge_const_handle LocateSeed( FT aTime ) ;

  bool IsVisited( Halfedge_const_handle aBisector ) { return mVisitedBisectors[aBisector->id()] != 0 ; }

  void Visit( Halfedge_const_handle aBisector ) { mVisitedBisectors[aBisector->id()] = 1 ; }

  inline Segment_2 CreateSegment ( Halfedge_const_handle aH ) const
  {
    Point_2 s = aH->opposite()->vertex()->point() ;
    Point_2 t = aH->vertex()->point() ;
    return K().construct_segment_2_object()(s,t);
  }

  Triedge GetTriedge( Vertex_const_handle aSeed ) const
  {
    typedef typename Ss::Vertex Vertex ;
    typedef typename Vertex::Defining_contour_halfedges_const_circulator circ ;
    
    circ cb = aSeed->defining_contour_halfedges_begin() ;
    circ c = cb ;

    Halfedge_const_handle lE[3] ;
    
    int d = 0 ;
    
    do
    {
      lE[d++] = *c++;
    }
    while( d < 3 && c != cb ) ;

    Triedge rTriedge(lE[0],lE[1],lE[2]);
    
    CGAL_POLYOFFSET_TRACE(3,"Triedge for " << v2str(*aSeed) << ": (E" << lE[0]->id() << ", E" << lE[1]->id() << ", E" << lE[2]->id() << ")" ) ;

    CGAL_postcondition(rTriedge.is_valid());
    
    return rTriedge ;
  }
  
  Trisegment_2 CreateTrisegment ( Halfedge_const_handle aE0
                                , Halfedge_const_handle aE1
                                , Halfedge_const_handle aE2
                                ) const
  {
    return *Construct_ss_trisegment_2(mTraits)(CreateSegment(aE0),CreateSegment(aE1),CreateSegment(aE2));
  }
  
  Trisegment_2 CreateTrisegment ( Triedge const& aTriedge ) const
  {
    return *Construct_ss_trisegment_2(mTraits)(CreateSegment(aTriedge.e0()),CreateSegment(aTriedge.e1()),CreateSegment(aTriedge.e2()));
  }
  
  Trisegment_2 CreateTrisegment ( Vertex_const_handle aSeed ) const
  {
    return aSeed->is_skeleton() ? CreateTrisegment(GetTriedge(aSeed)) : Trisegment_2::null() ;
  }
  
  Seeded_trisegment_2 CreateSeededTrisegment ( Halfedge_const_handle aE0
                                             , Halfedge_const_handle aE1
                                             , Halfedge_const_handle aE2
                                             , Vertex_const_handle   aLSeed
                                             , Vertex_const_handle   aRSeed
                                             ) const
  {
    return Construct_ss_seeded_trisegment_2(mTraits)(CreateTrisegment(aE0,aE1,aE2)
                                                    ,CreateTrisegment(aLSeed)
                                                    ,CreateTrisegment(aRSeed)
                                                    );
  }
  
  Seeded_trisegment_2 CreateSeededTrisegment ( Vertex_const_handle   aNode
                                             , Vertex_const_handle   aLSeed
                                             , Vertex_const_handle   aRSeed
                                             ) const
  {
    return Construct_ss_seeded_trisegment_2(mTraits)(CreateTrisegment(aNode)
                                                    ,CreateTrisegment(aLSeed)
                                                    ,CreateTrisegment(aRSeed)
                                                    );
  }
  
  Comparison_result Compare_offset_against_event_time( FT aT, Halfedge_const_handle aBisector, Halfedge_const_handle aNextBisector ) const
  {
    CGAL_assertion(aBisector->is_bisector());
    CGAL_assertion(handle_assigned(aBisector->opposite()));
    CGAL_assertion(aBisector->opposite()->is_bisector());
    CGAL_assertion(aNextBisector->is_bisector());
    CGAL_assertion(handle_assigned(aNextBisector->opposite()));
    CGAL_assertion(aNextBisector->opposite()->is_bisector());

    Halfedge_const_handle lBorderA = aBisector->defining_contour_edge();
    Halfedge_const_handle lBorderB = aBisector->opposite()->defining_contour_edge();
    Halfedge_const_handle lBorderC = aNextBisector->opposite()->defining_contour_edge();

    Vertex_const_handle lLSeed = aBisector->opposite()->vertex();
    Vertex_const_handle lRSeed = aBisector->vertex();
    
    CGAL_POLYOFFSET_TRACE(3,"Comparing offset " << aT << " against event [(E" << lBorderA->id() << ",E" << lBorderB->id() << ",E" << lBorderC->id() << ")(LS=" << v2str(*lLSeed) << ",RS=" << v2str(*lRSeed) << ")]");

    return Compare_offset_against_event_time_2(mTraits)(aT,CreateSeededTrisegment(lBorderA,lBorderB,lBorderC,lLSeed,lRSeed));
  }

  boost::optional<Point_2> Construct_offset_point( FT aT, Halfedge_const_handle aBisector ) const
  {
    CGAL_assertion(aBisector->is_bisector());
    CGAL_assertion(handle_assigned(aBisector->opposite()));
    CGAL_assertion(aBisector->opposite()->is_bisector());

    Halfedge_const_handle lBorderA = aBisector->defining_contour_edge();
    Halfedge_const_handle lBorderB = aBisector->opposite()->defining_contour_edge();
    
    // If aBisector is not a border bisector the offset point construction needs to get to the event
    // that generated the node closer to the border. That event is lSrcEvent 
    // (relevant only for inner bisectors)
    Seeded_trisegment_2 lSrcEvent ;
    if ( aBisector->is_inner_bisector() )
    {
      // aBisector might be pointing toward or away from the border, so we need to determine which
      // one of its endpoint nodes is closer to the border, that is, which one has the smaller event time.
      
      Seeded_trisegment_2 lSrcEventS, lSrcEventT ;
      
      Vertex_const_handle lNodeS = aBisector->opposite()->vertex();
      Vertex_const_handle lNodeT = aBisector->vertex();
      
      CGAL_assertion ( lNodeS->is_skeleton() ) ;
      CGAL_assertion ( lNodeT->is_skeleton() ) ;
        
      Vertex_const_handle lLSeedS = aBisector->prev()->opposite()->vertex();
      Vertex_const_handle lRSeedS = aBisector->opposite()->next()->vertex();
      lSrcEventS = CreateSeededTrisegment(lNodeS,lLSeedS,lRSeedS);
      
      Vertex_const_handle lLSeedT = aBisector->opposite()->prev()->opposite()->vertex();
      Vertex_const_handle lRSeedT = aBisector->next()->vertex();
      lSrcEventT = CreateSeededTrisegment(lNodeT,lLSeedT,lRSeedT);
      
      lSrcEvent = Compare_ss_event_times_2(mTraits)(lSrcEventS,lSrcEventT) != LARGER ? lSrcEventS : lSrcEventT ;
     
    }

    return Construct_offset_point_2(mTraits)(aT
                                            ,CreateSegment(lBorderA)
                                            ,CreateSegment(lBorderB)
                                            ,lSrcEvent
                                            );
  }

  void ResetVisitedBisectorsMap();

  Traits const&    mTraits ;
  std::vector<int> mVisitedBisectors;
  Halfedge_vector  mBorders ;
};

CGAL_END_NAMESPACE

#include <CGAL/Straight_skeleton_2/Polygon_offset_builder_2_impl.h>

#endif // CGAL_POLYGON_OFFSET_BUILDER_2_H //
// EOF //

