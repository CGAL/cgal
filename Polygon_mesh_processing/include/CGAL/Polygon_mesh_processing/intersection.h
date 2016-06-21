// Copyright (c) 2016 GeometryFactory (France).
// All rights reserved.
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
//
// Author(s)     : Sebastien Loriot

#ifndef CGAL_POLYGON_MESH_PROCESSING_INTERSECTION_H
#define CGAL_POLYGON_MESH_PROCESSING_INTERSECTION_H

#include <CGAL/Polygon_mesh_processing/internal/Corefinement/intersection_impl.h>

namespace CGAL {
namespace Polygon_mesh_processing{

///\todo use a functor similar to the one in split_graph_into_polylines
///      instead of the OutputIterator (+explicit construction from vector
///      that is document a default one)
template <class OutputIterator,
          class TriangleMesh
        //   class NamedParametersP,
        //   class NamedParametersQ,
        >
OutputIterator
surface_intersection(const TriangleMesh& tm1,
                     const TriangleMesh& tm2,
                     OutputIterator polyline_output,
                     bool throw_on_self_intersection=false)
{
  typedef typename boost::property_map<
    TriangleMesh,
    boost::vertex_point_t >::const_type VertexPointMap;
  Corefinement::Intersection_of_triangle_meshes<TriangleMesh,VertexPointMap>
    functor(tm1,
            tm2,
            get(boost::vertex_point, tm1),
            get(boost::vertex_point, tm2));
  functor(polyline_output, throw_on_self_intersection, true);
  return polyline_output;
}

} } //end of namespace CGAL::Polygon_mesh_processing

#endif // CGAL_POLYGON_MESH_PROCESSING_INTERSECTION_H
