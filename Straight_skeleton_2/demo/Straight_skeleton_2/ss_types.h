// Copyright (c) 2005, 2006 Fernando Luis Cacciola Carballal. All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
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
#ifndef SSKEL_DEMO_SS_TYPES_H
#define SSKEL_DEMO_SS_TYPES_H

#include <CGAL/Straight_skeleton_builder_2.h>
#include <CGAL/Polygon_offset_builder_2.h>
#include <CGAL/compute_outer_frame_margin.h>

namespace demo
{

typedef CGAL::Straight_skeleton_2<K> SSkel;

typedef SSkel::Halfedge_iterator     Halfedge_iterator;
typedef SSkel::Vertex_handle         Vertex_handle;
typedef SSkel::Face_const_iterator   Face_const_iterator;
typedef SSkel::Halfedge_const_handle Halfedge_const_handle ;
typedef SSkel::Vertex_const_handle   Vertex_const_handle ;

typedef boost::shared_ptr<SSkel> SSkelPtr ;

extern void draw_point  ( Point const& v, CGAL::Color c ) ;
extern void draw_segment( Point const& s, Point const& t, CGAL::Color c ) ;
extern void wait_on_user() ;

struct Visitor
{
  void on_contour_edge_entered ( Halfedge_const_handle const& he ) const
  {
    draw_segment(he->opposite()->vertex()->point(),he->vertex()->point(),CGAL::RED);
  }

  void on_initialization_started( int size_of_vertices ) const
  {
    mTotalVertices = size_of_vertices ;
    mVertexCount0 = 0 ;
    mVertexCount1 = 0 ;
    mReflexVertexCount = 0 ;
    mDegenerateVertexCount = 0 ;
    mFoundEdgeEventCount = 0 ;
    mFoundSplitEventCount = 0 ;
    mProcessedEdgeEventCount = 0 ;
    mProcessedSplitEventCount = 0 ;
    mProcessedPseudoSplitEventCount = 0 ;
    mAnihiliationCount = 0 ;
    mStage = 0 ;
  }

  void on_initial_events_collected( Vertex_const_handle const& v, bool is_reflex, bool is_degenerate ) const
  {
    ++ mVertexCount0 ;
    if ( is_reflex )
      ++ mReflexVertexCount ;
    if ( is_degenerate )
      ++ mDegenerateVertexCount ;

    draw_point(v->point(),CGAL::BLACK );

    //printf("\rInitialization: %d/%d (%d%%)",mVertexCount0,mTotalVertices,(mVertexCount0*100/mTotalVertices));
  }

  void on_edge_event_created( Vertex_const_handle const& 
                            , Vertex_const_handle const& 
                            )  const
  {
    ++ mFoundEdgeEventCount ;
  }

  void on_split_event_created( Vertex_const_handle const&  ) const
  {
    ++ mFoundSplitEventCount ;
  }

  void on_pseudo_split_event_created( Vertex_const_handle const& 
                                    , Vertex_const_handle const& 
                                    ) const
  {
  }

  void on_initialization_finished() const { printf("\n"); ++ mStage ; }

  void on_propagation_started() const {}

  void on_anihiliation_event_processed ( Vertex_const_handle const& node0
                                       , Vertex_const_handle const& node1
                                       ) const
  {
    draw_segment(node0->point(),node1->point(),CGAL::BLACK);
    ++ mAnihiliationCount ;
  }

  void on_edge_event_processed( Vertex_const_handle const& lseed
                              , Vertex_const_handle const& rseed
                              , Vertex_const_handle const& node
                              ) const
  {
    draw_segment(lseed->point(),node->point(), CGAL::BLACK );
    draw_segment(rseed->point(),node->point(), CGAL::BLACK );

    ++ mProcessedEdgeEventCount ;
  }

  void on_split_event_processed( Vertex_const_handle const& seed
                               , Vertex_const_handle const& node0
                               , Vertex_const_handle const& 
                               ) const
  {
    draw_segment(seed->point(),node0->point(), CGAL::BLACK );
    ++ mProcessedSplitEventCount ;
  }

  void on_pseudo_split_event_processed( Vertex_const_handle const& lseed
                                      , Vertex_const_handle const& rseed
                                      , Vertex_const_handle const& node0
                                      , Vertex_const_handle const& node1
                                      ) const
  {
    draw_segment(lseed->point(),node0->point(), CGAL::BLACK );
    draw_segment(rseed->point(),node1->point(), CGAL::BLACK );
    ++ mProcessedPseudoSplitEventCount ;
  }

  void on_vertex_processed( Vertex_const_handle const& node ) const
  {
    if ( node->is_contour() )
    {
      ++ mVertexCount1 ;
      wait_on_user();
      //printf("\rPropagation: %d/%d (%d%%)",mVertexCount1,mTotalVertices,(mVertexCount1*100/mTotalVertices));
    }
  }

  void on_propagation_finished() const { printf("\n"); ++ mStage ; }

  void on_cleanup_started() const {}

  void on_cleanup_finished() const {}

  void on_error( char const* what ) const
  {
    std::cerr << what << std::endl ;
  }

  void on_algorithm_finished ( bool /* finished_ok */ ) const
  {
  }

  mutable int mTotalVertices ;
  mutable int mVertexCount0 ;
  mutable int mVertexCount1 ;
  mutable int mReflexVertexCount ;
  mutable int mDegenerateVertexCount ;
  mutable int mFoundEdgeEventCount ;
  mutable int mFoundSplitEventCount ;
  mutable int mProcessedEdgeEventCount ;
  mutable int mProcessedSplitEventCount ;
  mutable int mProcessedPseudoSplitEventCount ;
  mutable int mAnihiliationCount ;
  mutable int mStage ;

  void print_stats()
  {
    std::cout << "VertexCount                   =" << mTotalVertices << std::endl
              << "ReflexVertexCount             =" << mReflexVertexCount << std::endl
              << "DegenerateVertexCount         =" << mDegenerateVertexCount << std::endl
              << "FoundEdgeEventCount           =" << mFoundEdgeEventCount << std::endl
              << "FoundSplitEventCount          =" << mFoundSplitEventCount << std::endl
              << "ProcessedEdgeEventCount       =" << mProcessedEdgeEventCount << std::endl
              << "ProcessedSplitEventCount      =" << mProcessedSplitEventCount << std::endl
              << "ProcessedPseudoSplitEventCount=" << mProcessedPseudoSplitEventCount << std::endl
              << "AnihiliationCount             =" << mAnihiliationCount << std::endl ;
  }

} ;

typedef CGAL::Straight_skeleton_builder_traits_2<K>                         SSkelBuilderTraits;
typedef CGAL::Straight_skeleton_builder_2<SSkelBuilderTraits,SSkel,Visitor> SSkelBuilder;

typedef CGAL::Polygon_offset_builder_traits_2<K>                          OffsetBuilderTraits;
typedef CGAL::Polygon_offset_builder_2<SSkel,OffsetBuilderTraits,Polygon> OffsetBuilder;


}

#endif // SSKEL_DEMO_SS_TYPES_H
