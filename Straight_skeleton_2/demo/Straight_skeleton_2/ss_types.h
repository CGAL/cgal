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

#include <CGAL/Straight_skeleton_builder_2.h>
#include <CGAL/Polygon_offset_builder_2.h>
#include <CGAL/compute_outer_frame_margin.h>

namespace demo
{

typedef CGAL::Straight_skeleton_2<K>                                SSkel;
typedef CGAL::Straight_skeleton_builder_traits_2<K>                 SSkelBuilderTraits;
typedef CGAL::Straight_skeleton_builder_2<SSkelBuilderTraits,SSkel> SSkelBuilder;

typedef CGAL::Polygon_offset_builder_traits_2<K>                          OffsetBuilderTraits;
typedef CGAL::Polygon_offset_builder_2<SSkel,OffsetBuilderTraits,Polygon> OffsetBuilder;

typedef SSkel::Halfedge_iterator     Halfedge_iterator;
typedef SSkel::Vertex_handle         Vertex_handle;
typedef SSkel::Face_const_iterator   Face_const_iterator;
typedef SSkel::Halfedge_const_handle Halfedge_const_handle ;
typedef SSkel::Vertex_const_handle   Vertex_const_handle ;

typedef boost::shared_ptr<SSkel> SSkelPtr ;

}
