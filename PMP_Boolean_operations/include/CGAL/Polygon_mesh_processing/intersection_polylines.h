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


#ifndef CGAL_POLYGON_MESH_PROCESSING_INTERSECTION_POLYLINES_H
#define CGAL_POLYGON_MESH_PROCESSING_INTERSECTION_POLYLINES_H

#include <CGAL/license/Polygon_mesh_processing/corefinement.h>

#include <CGAL/Polygon_mesh_processing/internal/Corefinement/intersection_impl.h>

#include <CGAL/boost/graph/named_params_helper.h>

#include <type_traits>

namespace CGAL {
namespace Polygon_mesh_processing {

/**
 * \ingroup PMP_corefinement_grp
 *
 * computes the intersection of triangles of `tm1` and `tm2`. The output is a
 * set of polylines with all vertices but endpoints being of degree 2.
 *
 * \pre \link CGAL::Polygon_mesh_processing::does_self_intersect() `!CGAL::Polygon_mesh_processing::does_self_intersect(tm1)` \endlink
 * \pre \link CGAL::Polygon_mesh_processing::does_self_intersect() `!CGAL::Polygon_mesh_processing::does_self_intersect(tm2)` \endlink
 *
 * @tparam TriangleMesh a model of `MutableFaceGraph`, `HalfedgeListGraph` and `FaceListGraph`
 * @tparam NamedParameters1 a sequence of \ref bgl_namedparameters "Named Parameters"
 * @tparam NamedParameters2 a sequence of \ref bgl_namedparameters "Named Parameters"
 * @tparam OutputIterator an output iterator in which `std::vector` of points
 *                        can be put. The point type is the one from the
 *                        vertex property map
 *
 * @param tm1 first input triangulated surface mesh
 * @param tm2 second input triangulated surface mesh
 * @param polyline_output output iterator of polylines. Each polyline will be
 *        given as a vector of points
 * @param np1 an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
 * @param np2 an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
 *
 * \cgalNamedParamsBegin
 *   \cgalParamNBegin{vertex_point_map}
 *     \cgalParamDescription{a property map associating points to the vertices of `tm1` (`tm2`)}
 *     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<TriangleMesh>::%vertex_descriptor`
 *                    as key type and `%Point_3` as value type}
 *     \cgalParamDefault{`boost::get(CGAL::vertex_point, tm1 (tm2))`}
 *     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
 *                     should be available for the vertices of `tm1` (`tm2`)}
 *     \cgalParamExtra{Both vertex point maps must have the same value type}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{throw_on_self_intersection}
 *     \cgalParamDescription{If `true`, the set of triangles close to the intersection of `tm1` and `tm2` will be
 *                           checked for self-intersections and `Corefinement::Self_intersection_exception`
 *                           will be thrown if at least one self-intersection is found.}
 *     \cgalParamType{Boolean}
 *     \cgalParamDefault{`false`}
 *     \cgalParamExtra{`np1` only}
 *   \cgalParamNEnd
 * \cgalNamedParamsEnd
 *
 * \see `do_intersect()`
 */
template <class OutputIterator,
          class TriangleMesh,
          class NamedParameters1 = parameters::Default_named_parameters,
          class NamedParameters2 = parameters::Default_named_parameters >
OutputIterator
intersection_polylines(const TriangleMesh& tm1,
                       const TriangleMesh& tm2,
                       OutputIterator polyline_output,
                       const NamedParameters1& np1 = parameters::default_values(),
                       const NamedParameters2& np2 = parameters::default_values())
{
  const bool throw_on_self_intersection =
    parameters::choose_parameter(parameters::get_parameter(np1, internal_np::throw_on_self_intersection), false);

  typedef typename GetVertexPointMap<TriangleMesh, NamedParameters1>::const_type VPM1;
  typedef typename GetVertexPointMap<TriangleMesh, NamedParameters2>::const_type VPM2;

  static_assert(std::is_same<typename boost::property_traits<VPM1>::value_type,
                                      typename boost::property_traits<VPM2>::value_type>::value);

  VPM1 vpm1 = parameters::choose_parameter(parameters::get_parameter(np1, internal_np::vertex_point),
                                           get_const_property_map(CGAL::vertex_point, tm1));
  VPM2 vpm2 = parameters::choose_parameter(parameters::get_parameter(np2, internal_np::vertex_point),
                                           get_const_property_map(CGAL::vertex_point, tm2));

  Corefinement::Intersection_of_triangle_meshes<TriangleMesh, VPM1, VPM2>
    functor(tm1, tm2, vpm1, vpm2);

  // Fill non-manifold feature maps if provided
  functor.set_non_manifold_feature_map_1(parameters::get_parameter(np1, internal_np::non_manifold_feature_map));
  functor.set_non_manifold_feature_map_2(parameters::get_parameter(np2, internal_np::non_manifold_feature_map));

  return functor(polyline_output, throw_on_self_intersection, true);
}

#ifndef CGAL_NO_DEPRECATED_CODE

/// \ingroup PMP_corefinement_grp
/// \deprecated This function is deprecated since CGAL 6.2, use `CGAL::Polygon_mesh_processing::intersection_polylines()` instead.
template <class OutputIterator,
          class TriangleMesh,
          class NamedParameters1 = parameters::Default_named_parameters,
          class NamedParameters2 = parameters::Default_named_parameters >
CGAL_DEPRECATED
OutputIterator
surface_intersection(const TriangleMesh& tm1,
                     const TriangleMesh& tm2,
                     OutputIterator polyline_output,
                     const NamedParameters1& np1 = parameters::default_values(),
                     const NamedParameters2& np2 = parameters::default_values())
{
  return intersection_polylines(tm1, tm2, polyline_output, np1, np2);
}

#endif // CGAL_NO_DEPRECATED_CODE

} // namespace Polygon_mesh_processing
} // namespace CGAL

#endif // CGAL_POLYGON_MESH_PROCESSING_INTERSECTION_POLYLINES_H
