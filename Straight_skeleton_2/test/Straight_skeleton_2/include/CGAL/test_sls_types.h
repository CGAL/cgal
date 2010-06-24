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
// $URL: svn+ssh://fcacciola@scm.gforge.inria.fr/svn/cgal/trunk/Straight_skeleton_2/test/Straight_skeleton_2/include/CGAL/test_offset_builder_types.h $
// $Id: test_offset_builder_types.h 37121 2007-03-15 10:47:09Z hemmer $
//
// Author(s)     : Fernando Cacciola <fernando_cacciola@ciudad.com.ar>
//
#ifndef CGAL_TEST_SLS_TYPES_H
#define CGAL_TEST_SLS_TYPES_H

#include <CGAL/basic.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Straight_skeleton_builder_2.h>
#include <CGAL/Straight_skeleton_converter_2.h>
#include <CGAL/Polygon_offset_builder_2.h>
#include <CGAL/compute_outer_frame_margin.h>
#include <CGAL/HalfedgeDS_const_decorator.h>

#include <boost/shared_ptr.hpp>
#include <vector>

#include <CGAL/Real_timer.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel IK;
typedef CGAL::Exact_predicates_inexact_constructions_kernel OK;


typedef IK::FT                         IFT ;
typedef OK::FT                         OFT ;

typedef IK::Point_2                    IPoint;
typedef OK::Point_2                    OPoint;

typedef CGAL::Aff_transformation_2<IK> ITransformation;

typedef std::vector<double>           Doubles ;

typedef CGAL::Segment_2<IK>           ISegment;
typedef std::vector<IPoint>           IPolygon;
typedef std::vector<IFT>              IWeights;
typedef boost::shared_ptr<IPolygon>   IPolygonPtr;
typedef boost::shared_ptr<IWeights>   IWeightsPtr;
typedef std::vector<IWeightsPtr>      IWeightsPtrContainer ;

struct IWeightedPolygon 
{
  typedef IPolygonPtr polygon_ptr_type ;
  typedef IWeightsPtr weights_ptr_type ;
  
  IWeightedPolygon(IPolygonPtr aPolygonPtr, IWeightsPtr aWeightsPtr, bool aClosed = true) 
    : polygon(aPolygonPtr), weights(aWeightsPtr), closed(aClosed) {}
    
  IWeightedPolygon clone() const
  {
    return IWeightedPolygon( IPolygonPtr( new IPolygon(*polygon) ), IWeightsPtr( new IWeights(*weights) ), closed );
  }

  IPolygonPtr polygon; 
  IWeightsPtr weights;
  bool        closed;
};

IWeightedPolygon revert_weighted_polygon(IWeightedPolygon aWeightedPolygon)
{
  IPolygonPtr lPolygonPtr = IPolygonPtr( new IPolygon(aWeightedPolygon.polygon->rbegin(),aWeightedPolygon.polygon->rend()) ) ;  
  IWeightsPtr lWeightsPtr = IWeightsPtr( new IWeights(aWeightedPolygon.weights->rbegin(),aWeightedPolygon.weights->rend()) ) ;
  std::rotate(lWeightsPtr->begin(),lWeightsPtr->begin()+(lWeightsPtr->size()-1),lWeightsPtr->end()) ;
  return IWeightedPolygon(lPolygonPtr, lWeightsPtr, aWeightedPolygon.closed ) ;
}

IWeightedPolygon invert_weighted_polygon(IWeightedPolygon aWeightedPolygon)
{
  IWeightedPolygon rInverseWeightedPolygon = aWeightedPolygon.clone() ;
  for( IWeights::iterator wi = rInverseWeightedPolygon.weights->begin(), wi_end = rInverseWeightedPolygon.weights->end() ;
       wi != wi_end ; ++wi )
  {
    *wi *= -1.0;
  }
  return rInverseWeightedPolygon ;
}

typedef std::vector<IWeightedPolygon>          IWeightedBoundaries ;
typedef boost::shared_ptr<IWeightedBoundaries> IWeightedBoundariesPtr ;

typedef std::vector<IPolygonPtr >       IBoundaries ;
typedef boost::shared_ptr<IBoundaries>  IBoundariesPtr ;


IBoundariesPtr extract_polygons_view(IWeightedBoundaries const & aWeightedBoundaries)
{
  IBoundariesPtr rBoundaries = IBoundariesPtr(new IBoundaries());
  for( IWeightedBoundaries::const_iterator iter = aWeightedBoundaries.begin(), iter_end = aWeightedBoundaries.end(); iter != iter_end ; ++iter )
  {
    rBoundaries->push_back( iter->polygon );
  }
  return rBoundaries;
}


typedef CGAL::Segment_2<OK>            OSegment ;
typedef std::vector<OPoint>            OPolygon ;
typedef boost::shared_ptr<OPolygon>    OPolygonPtr ;
typedef std::vector<OPolygonPtr>       OBoundaries ;
typedef boost::shared_ptr<OBoundaries> OBoundariesPtr ;

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

  void on_initialization_started( std::size_t size_of_vertices ) const {}

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
  
  void on_offset_point ( IPoint const& ) const { check_timeout(); }

  IPoint on_offset_point_overflowed( Halfedge_const_handle ) const { return CGAL::ORIGIN ; }
  
  void on_offset_contour_finished ( bool ) const { check_timeout(); }
  
  void on_construction_finished () const {}
  
  void on_error( char const* ) const { check_timeout(); }

} ;

typedef CGAL::HalfedgeDS_const_decorator<ISls> Sls_const_decorator ;

typedef boost::shared_ptr<ISls> ISlsPtr ;
typedef boost::shared_ptr<OSls> OSlsPtr ;

typedef CGAL::Straight_skeleton_items_converter_2<ISls,OSls> SlsItemsConverter ;

typedef CGAL::Straight_skeleton_converter_2<ISls,OSls,SlsItemsConverter> SlsConverter ;

typedef CGAL::Straight_skeleton_builder_traits_2<IK>                                 ISlsBuilderTraits;
typedef CGAL::Straight_skeleton_builder_2<ISlsBuilderTraits,ISls,ISlsBuilderVisitor> ISlsBuilder;

typedef CGAL::Polygon_offset_builder_traits_2<OK>                                               OffsetBuilderTraits;
typedef CGAL::Polygon_offset_builder_2<OSls,OffsetBuilderTraits,OPolygon,IOffsetBuilderVisitor> OffsetBuilder;

#endif

