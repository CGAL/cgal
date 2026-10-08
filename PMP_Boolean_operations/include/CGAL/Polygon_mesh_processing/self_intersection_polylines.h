// Copyright (c) 2024-2025 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Sébastien Loriot


#ifndef CGAL_POLYGON_MESH_PROCESSING_SELF_INTERSECTION_POLYLINES_H
#define CGAL_POLYGON_MESH_PROCESSING_SELF_INTERSECTION_POLYLINES_H

#include <CGAL/license/Polygon_mesh_processing/corefinement.h>

#include <CGAL/Polygon_mesh_processing/internal/Corefinement/intersection_impl.h>

#include <CGAL/boost/graph/named_params_helper.h>

namespace CGAL {
namespace Polygon_mesh_processing {
namespace experimental {

/**
 * \ingroup PMP_corefinement_grp
 *
 * computes the autointersection of triangles of `tm`. The output is a
 * set of polylines with all vertices but endpoints being of degree 2.
 *
 * @tparam TriangleMesh a model of `HalfedgeListGraph` and `FaceListGraph`
 * @tparam NamedParameters a sequence of \ref namedparameters
 * @tparam OutputIterator an output iterator in which `std::vector` of points
 *                        can be put. The point type is the one from the
 *                        vertex property map
 *
 * @param tm input triangulated surface mesh
 * @param polyline_output output iterator of polylines. Each polyline will be
 *        given as a vector of points
 * @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
 *
 * \cgalNamedParamsBegin
 *   \cgalParamNBegin{vertex_point_map}
 *     \cgalParamDescription{a property map associating points to the vertices of `tm`}
 *     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<TriangleMesh>::%vertex_descriptor`
 *                    as key type and `%Point_3` as value type}
 *     \cgalParamDefault{`boost::get(CGAL::vertex_point, tm)`}
 *   \cgalParamNEnd
 * \cgalNamedParamsEnd
 *
 */
template <class OutputIterator,
          class TriangleMesh,
          class NamedParameters = parameters::Default_named_parameters >
OutputIterator
self_intersection_polylines(const TriangleMesh& tm,
                            OutputIterator polyline_output,
                            const NamedParameters& np = parameters::default_values())
{
// Vertex point maps
  typedef typename GetVertexPointMap<TriangleMesh, NamedParameters>::const_type VPM;

  VPM vpm = parameters::choose_parameter(parameters::get_parameter(np, internal_np::vertex_point),
                                         get_const_property_map(CGAL::vertex_point, tm));

// surface intersection algorithm call
  typedef Corefinement::Default_surface_intersection_visitor<TriangleMesh, true> Visitor;
  Corefinement::Intersection_of_triangle_meshes<TriangleMesh, VPM, VPM, Visitor> functor(tm, vpm);

  polyline_output = functor(polyline_output, true);
  return polyline_output;
}

} // namespace experimental
} // namespace Polygon_mesh_processing
} // namespace CGAL

#endif // CGAL_POLYGON_MESH_PROCESSING_SELF_INTERSECTION_POLYLINES_H
