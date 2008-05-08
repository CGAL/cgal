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
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel_with_sqrt.h>

#include <CGAL/Straight_skeleton_builder_2.h>
#include <CGAL/Straight_skeleton_converter_2.h>
#include <CGAL/Polygon_offset_builder_2.h>
#include <CGAL/compute_outer_frame_margin.h>
#include <CGAL/HalfedgeDS_const_decorator.h>

#include <boost/shared_ptr.hpp>
#include <vector>

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
typedef boost::shared_ptr<IPolygon>   IPolygonPtr;
typedef std::vector<IPolygonPtr>      IRegion ;
typedef boost::shared_ptr<IRegion>    IRegionPtr ;
typedef std::vector<IRegionPtr>       IRegions ;

typedef CGAL::Segment_2<OK>           OSegment;
typedef std::vector<OPoint>           OPolygon;
typedef boost::shared_ptr<OPolygon>   OPolygonPtr;
typedef std::vector<OPolygonPtr>      ORegion ;
typedef boost::shared_ptr<ORegion>    ORegionPtr ;
typedef std::vector<ORegionPtr>       ORegions ;

typedef CGAL::Straight_skeleton_2<IK>                             ISls;
typedef CGAL::Straight_skeleton_builder_traits_2<IK>              ISlsBuilderTraits;
typedef CGAL::Straight_skeleton_builder_2<ISlsBuilderTraits,ISls> ISlsBuilder;

typedef CGAL::Straight_skeleton_2<OK>                                     OSls;
typedef CGAL::Polygon_offset_builder_traits_2<OK>                         OffsetBuilderTraits;
typedef CGAL::Polygon_offset_builder_2<OSls,OffsetBuilderTraits,OPolygon> OffsetBuilder;

typedef ISls::Halfedge_iterator       Halfedge_iterator;
typedef ISls::Vertex_handle           Vertex_handle;
typedef ISls::Face_const_iterator     Face_const_iterator;
typedef ISls::Halfedge_const_handle   Halfedge_const_handle ;
typedef ISls::Halfedge_const_iterator Halfedge_const_iterator ;
typedef ISls::Vertex_const_handle     Vertex_const_handle ;
typedef ISls::Vertex_const_iterator   Vertex_const_iterator ;

typedef CGAL::HalfedgeDS_const_decorator<ISls> Sls_const_decorator ;

typedef boost::shared_ptr<ISls> ISlsPtr ;
typedef boost::shared_ptr<OSls> OSlsPtr ;

typedef CGAL::Straight_skeleton_items_converter_2<ISls,OSls> SlsItemsConverter ;

typedef CGAL::Straight_skeleton_converter_2<ISls,OSls,SlsItemsConverter> SlsConverter ;


#endif

