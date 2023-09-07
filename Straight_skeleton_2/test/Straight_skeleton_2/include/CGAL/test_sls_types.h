// Copyright (c) 2005, 2006 Fernando Luis Cacciola Carballal. All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Fernando Cacciola <fernando_cacciola@ciudad.com.ar>
//
#ifndef CGAL_TEST_SLS_TYPES_H
#define CGAL_TEST_SLS_TYPES_H


#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
//#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
//#include <CGAL/Exact_predicates_exact_constructions_kernel_with_sqrt.h>

#include <CGAL/Straight_skeleton_builder_2.h>
#include <CGAL/Straight_skeleton_converter_2.h>
#include <CGAL/Polygon_offset_builder_2.h>
#include <CGAL/compute_outer_frame_margin.h>
#include <CGAL/HalfedgeDS_const_decorator.h>

#include <memory>
#include <vector>

#include <CGAL/Real_timer.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel   IK;
//typedef CGAL::Exact_predicates_exact_constructions_kernel     IK;
//typedef CGAL::Exact_predicates_exact_constructions_kernel_with_sqrt IK;


typedef CGAL::Exact_predicates_inexact_constructions_kernel   OK;
//typedef CGAL::Exact_predicates_exact_constructions_kernel     OK;


typedef IK::FT                         IFT ;
typedef OK::FT                         OFT ;

typedef IK::Point_2                    IPoint;
typedef OK::Point_2                    OPoint;

typedef CGAL::Aff_transformation_2<IK> ITransformation;

typedef std::vector<double>           Doubles ;

typedef CGAL::Segment_2<IK>           ISegment;
typedef std::vector<IPoint>           IPolygon;
typedef std::shared_ptr<IPolygon>   IPolygonPtr;
typedef std::vector<IPolygonPtr>      IRegion ;
typedef std::shared_ptr<IRegion>    IRegionPtr ;
typedef std::vector<IRegionPtr>       IRegions ;

typedef CGAL::Segment_2<OK>           OSegment;
typedef std::vector<OPoint>           OPolygon;
typedef std::shared_ptr<OPolygon>   OPolygonPtr;
typedef std::vector<OPolygonPtr>      ORegion ;
typedef std::shared_ptr<ORegion>    ORegionPtr ;
typedef std::vector<ORegionPtr>       ORegions ;

typedef CGAL::Straight_skeleton_2<IK>                             ISls;

typedef CGAL::Straight_skeleton_2<OK>                                     OSls;

typedef ISls::Halfedge_iterator       Halfedge_iterator;
typedef ISls::Vertex_handle           Vertex_handle;
typedef ISls::Face_const_iterator     Face_const_iterator;
typedef ISls::Halfedge_const_handle   Halfedge_const_handle ;
typedef ISls::Halfedge_const_iterator Halfedge_const_iterator ;
typedef ISls::Vertex_const_handle     Vertex_const_handle ;
typedef ISls::Vertex_const_iterator   Vertex_const_iterator ;

class VisitorBase
{

public:

  typedef void (*CheckTimeoutCallbackType) () ;

  VisitorBase( CheckTimeoutCallbackType aCheckTimeoutCallback ) : check_timeout(aCheckTimeoutCallback) {}

protected:

  CheckTimeoutCallbackType check_timeout ;
} ;

class ISlsBuilderVisitor : public VisitorBase
{
public:

  ISlsBuilderVisitor( CheckTimeoutCallbackType aCheckTimeoutCallback ) : VisitorBase(aCheckTimeoutCallback) {}

  void on_contour_edge_entered ( Halfedge_const_handle const& ) const {}

  void on_initialization_started( int /*size_of_vertices*/ ) const {}

  void on_initial_events_collected( Vertex_const_handle const& , bool , bool ) const { check_timeout(); }

  void on_edge_event_created( Vertex_const_handle const&, Vertex_const_handle const& ) const { check_timeout(); }

  void on_split_event_created( Vertex_const_handle const& ) const { check_timeout(); }

  void on_pseudo_split_event_created( Vertex_const_handle const&, Vertex_const_handle const& ) const { check_timeout(); }

  void on_initialization_finished() const { check_timeout(); }

  void on_propagation_started() const {}

  void on_anihiliation_event_processed ( Vertex_const_handle const&, Vertex_const_handle const& ) const { check_timeout(); }

  void on_edge_event_processed( Vertex_const_handle const&, Vertex_const_handle const&, Vertex_const_handle const& ) const { check_timeout(); }

  void on_split_event_processed( Vertex_const_handle const&, Vertex_const_handle const&, Vertex_const_handle const& ) const { check_timeout(); }

  void on_pseudo_split_event_processed( Vertex_const_handle const&, Vertex_const_handle const&, Vertex_const_handle const&, Vertex_const_handle const& ) const { check_timeout(); }

  void on_vertex_processed( Vertex_const_handle const& ) const { check_timeout(); }

  void on_propagation_finished() const { check_timeout(); }

  void on_cleanup_started() const {}

  void on_cleanup_finished() const { check_timeout(); }

  void on_error( char const* ) const { check_timeout(); }

  void on_algorithm_finished ( bool ) const {}

} ;


class IOffsetBuilderVisitor : public VisitorBase
{
public:

  IOffsetBuilderVisitor( CheckTimeoutCallbackType aCheckTimeoutCallback ) : VisitorBase(aCheckTimeoutCallback) {}

  void on_construction_started ( IFT ) const {}

  void on_offset_contour_started() const { check_timeout(); }

  void on_offset_point ( IPoint const&, Halfedge_const_handle ) const { check_timeout(); }

  IPoint on_offset_point_overflowed( Halfedge_const_handle ) const { return CGAL::ORIGIN ; }

  void on_offset_contour_finished ( bool ) const { check_timeout(); }

  void on_construction_finished () const {}

  void on_error( char const* ) const { check_timeout(); }

} ;

typedef CGAL::HalfedgeDS_const_decorator<ISls> Sls_const_decorator ;

typedef std::shared_ptr<ISls> ISlsPtr ;
typedef std::shared_ptr<OSls> OSlsPtr ;

typedef CGAL::Straight_skeleton_items_converter_2<ISls,OSls> SlsItemsConverter ;

typedef CGAL::Straight_skeleton_converter_2<ISls,OSls,SlsItemsConverter> SlsConverter ;

typedef CGAL::Straight_skeleton_builder_traits_2<IK>                                 ISlsBuilderTraits;
typedef CGAL::Straight_skeleton_builder_2<ISlsBuilderTraits,ISls,ISlsBuilderVisitor> ISlsBuilder;

typedef CGAL::Polygon_offset_builder_traits_2<OK>                                               OffsetBuilderTraits;
typedef CGAL::Polygon_offset_builder_2<OSls,OffsetBuilderTraits,OPolygon,IOffsetBuilderVisitor> OffsetBuilder;

#endif

