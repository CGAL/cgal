// Copyright (c) 2018  GeometryFactory Sarl (France).
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
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Liubomyr Piadyk

#ifndef CGAL_SURFACE_MESH_SEGMENTATION_APPROXIMATE_CONVEX_SEGMENTATION_H
#define CGAL_SURFACE_MESH_SEGMENTATION_APPROXIMATE_CONVEX_SEGMENTATION_H

#include <CGAL/license/Surface_mesh_segmentation.h>

/**
 * @file approximate_convex_segmentation.h
 * @brief The API which contains template functions for concavity value computation and approximate convex segmentation
 */
#include <CGAL/internal/Approximate_convex_segmentation/Approx_segmentation.h>
#include <CGAL/internal/Approximate_convex_segmentation/Concavity.h>
#include <CGAL/boost/graph/named_params_helper.h>

namespace CGAL
{

/*!
 * \ingroup PkgSurfaceSegmentation
 * @brief computes the concavity value of a set of faces that are assigned a given id in a triangle mesh.
 *
 * The concavity value of a set of faces is the largest distance of a vertex of the set of faces to its convex hull.
 * The distance is either computed using the closest point or the projection along the normal at the vertex,
 * depending on the value of the named parameter `use_closest_point`.
 *
 * @pre `is_triangle_mesh(mesh)`
 *
 * @tparam ConcurrencyTag enables sequential versus parallel algorithm (possible values are `Sequential_tag` and `Parallel_tag`)
 * @tparam TriangleMesh a model of `FaceListGraph`
 * @tparam FacePropertyMap a `ReadablePropertyMap` with `boost::graph_traits<TriangleMesh>::%face_descriptor` as key type and an integer as value type
 * @tparam DistancesPropertyMap a `WritablePropertyMap` with `boost::graph_traits<TriangleMesh>::%vertex_descriptor` as key type and a floating-point as value type
 *                              If `ConcurrencyTag` is `Parallel_tag`, then this type must support thread-safe calls to `put()` by concurrent threads.
 *                              An initial sequential pass will first set all values to `0` before values being updated by concurrent threads.
 * @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
 *
 * @param mesh segmented triangle mesh
 * @param face_ids property map with per face segment ids
 * @param segment_id id of the target segment on which concavity value is computed
 * @param[out] distances optional property map with per vertex squared distance to the convex hull
 * @param np optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
 *
 * \cgalNamedParamsBegin
 *    \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `mesh` \cgalParamEnd
 *    \cgalParamBegin{geom_traits} a geometric traits class instance \cgalParamEnd
 *    \cgalParamBegin{use_closest_point} if true, the concavity at each vertex is evaluated by using the distance to the closest point on the convex hull of the set of faces. If false, the distance to the first intersected point following the normal at each vertex is used. Default is false \cgalParamEnd
 *    \cgalParamBegin{face_index_map} a property map containing an index for each face initialized from 0 to `num_faces(graph)`
 *    \cgalParamBegin{vertex_index_map} a property map containing an index for each vertex initialized from 0 to `num_vertices(graph)`
 *    \cgalParamBegin{halfedge_index_map} a property map containing an index for each halfedge initialized from 0 to `num_halfedges(graph)`
 * \cgalNamedParamsEnd
 *
 * @return the concavity value of the set of faces
 */
template <class ConcurrencyTag, class TriangleMesh, class FacePropertyMap, class DistancesPropertyMap, class NamedParameters>
double
concavity_values(const TriangleMesh& mesh,
                 FacePropertyMap face_ids,
                 std::size_t segment_id,
                 DistancesPropertyMap distances,
                 const NamedParameters& np)
{
  using Vpm = typename CGAL::GetVertexPointMap<TriangleMesh, NamedParameters>::const_type;
  using Geom_traits = typename GetGeomTraits<TriangleMesh, NamedParameters>::type;

  using CGAL::parameters::get_parameter;
  using CGAL::parameters::choose_parameter;

  Geom_traits geom_traits = choose_parameter<Geom_traits>(get_parameter(np, internal_np::geom_traits));
  Vpm vpm = choose_parameter(get_parameter(np, CGAL::internal_np::vertex_point), get_const_property_map(CGAL::vertex_point, mesh));

  CGAL_precondition(is_triangle_mesh(mesh));

  bool use_closest_point = choose_parameter(get_parameter(np, internal_np::use_closest_point), false);
  //const double w2 = parameters::choose_parameter(parameters::get_parameter(np1, internal_np::point_to_plane_weight), 2);

  internal::Concavity<TriangleMesh, Vpm, Geom_traits, ConcurrencyTag> algorithm(mesh, vpm, geom_traits, use_closest_point);
  return algorithm.compute(face_ids, segment_id, distances);
}


#ifndef DOXYGEN_RUNNING
template <class ConcurrencyTag, class TriangleMesh, class FacePropertyMap, class DistancesPropertyMap>
double
concavity_values(const TriangleMesh& mesh,
                 FacePropertyMap face_ids,
                 std::size_t segment_id,
                 DistancesPropertyMap distances)
{
    return concavity_values<ConcurrencyTag>(mesh, face_ids, segment_id, distances, Polygon_mesh_processing::parameters::all_default());
}

template <class ConcurrencyTag, class TriangleMesh, class FacePropertyMap>
double
concavity_values(const TriangleMesh& mesh,
                 FacePropertyMap face_ids,
                 std::size_t segment_id)
{
    CGAL::Constant_property_map<typename boost::graph_traits<TriangleMesh>::vertex_descriptor, double > distances_pmap(0);

    return concavity_values<ConcurrencyTag>(mesh, face_ids, segment_id, distances_pmap, Polygon_mesh_processing::parameters::all_default());
}
#endif


/*!
 * \ingroup PkgSurfaceSegmentation
 * @brief computes the concavity value of a triangle mesh.
 *
 * The concavity value is the largest distance of a vertex to the convex hull of the mesh.
 * The distance is either computed using the closest point or the projection along the normal at the vertex,
 * depending on the value of the named parameter `use_closest_point`.
 * @pre `is_triangle_mesh(mesh)`
 *
 * @tparam ConcurrencyTag enables sequential versus parallel algorithm (possible values are `Sequential_tag` and `Parallel_tag`)
 * @tparam TriangleMesh a model of `FaceListGraph`
 * @tparam DistancesPropertyMap a `WritablePropertyMap` with the `boost::graph_traits<TriangleMesh>::%vertex_descriptor` key type and a floating-point as value type
 *                              If `ConcurrencyTag` is `Parallel_tag`, then this type must support thread-safe calls to `put()` by concurrent threads.
 *                              An initial sequential pass will first set all values to `0` before values being updated by concurrent threads.
 * @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
 *
 * @param mesh triangle mesh on which concavity value is computed
 * @param distances optional property map with per vertex squared distance to the convex hull
 * @param np optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
 *
 * \cgalNamedParamsBegin
 *    \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `mesh` \cgalParamEnd
 *    \cgalParamBegin{geom_traits} a geometric traits class instance \cgalParamEnd
 *    \cgalParamBegin{use_closest_point} if true, the concavity at each vertex is evaluated by using the distance to the closest point on the convex hull of the set of faces. If false, the distance to the first intersected point following the normal at each vertex is used. Default is false \cgalParamEnd
 * \cgalNamedParamsEnd
 *
 * @return the concavity value
 */
template <class ConcurrencyTag, class TriangleMesh, class DistancesPropertyMap, class NamedParameters>
double
concavity_values(const TriangleMesh& mesh,
                 DistancesPropertyMap distances,
                 const NamedParameters& np)
{
  using Vpm = typename CGAL::GetVertexPointMap<TriangleMesh, NamedParameters>::const_type;
  using Geom_traits = typename GetGeomTraits<TriangleMesh, NamedParameters>::type;


  using CGAL::parameters::get_parameter;
  using CGAL::parameters::choose_parameter;

  Geom_traits geom_traits = choose_parameter<Geom_traits>(get_parameter(np, internal_np::geom_traits));
  Vpm vpm = choose_parameter(get_parameter(np, CGAL::internal_np::vertex_point), get_const_property_map(CGAL::vertex_point, mesh));

  CGAL_precondition(is_triangle_mesh(mesh));

  bool use_closest_point = choose_parameter(get_parameter(np, internal_np::use_closest_point), false);

  internal::Concavity<TriangleMesh, Vpm, Geom_traits, ConcurrencyTag> algorithm(mesh, vpm, geom_traits, use_closest_point);
  return algorithm.compute(distances);
}


#ifndef DOXYGEN_RUNNING
template <class ConcurrencyTag, class TriangleMesh, class DistancesPropertyMap>
double
concavity_values(const TriangleMesh& mesh,
                 DistancesPropertyMap distances)
{
    return concavity_values<ConcurrencyTag>(mesh, distances, Polygon_mesh_processing::parameters::all_default());
}

template <class ConcurrencyTag, class TriangleMesh>
double
concavity_values(const TriangleMesh& mesh)
{
    CGAL::Constant_property_map<typename boost::graph_traits<TriangleMesh>::vertex_descriptor, double > distances_pmap(0);

    return concavity_values<ConcurrencyTag>(mesh, distances_pmap, Polygon_mesh_processing::parameters::all_default());
}
#endif


/*!
 * \ingroup PkgSurfaceSegmentation
 * @brief computes an approximate convex segmentation of a triangle mesh.
 *
 * This function fills a property map associating a segment-id to each face (in the range `[0, N-1])`.
 * `N` is the number of segments computed by the functions (greater are equal to `min_number_of_segments`).
 * The set of faces in the same segment defines an edge-connected patch with concavity value less are equal to `concavity_threshold`.
 *
 * @pre `is_triangle_mesh(mesh)`
 * @pre `faces(mesh).size() >= min_number_of_segments`
 *
 * @tparam ConcurrencyTag enables sequential versus parallel algorithm (possible values are `Sequential_tag` and `Parallel_tag`)
 * @tparam TriangleMesh a model of `FaceListGraph`
 * @tparam FacePropertyMap a `ReadWritePropertyMap` with the `boost::graph_traits<TriangleMesh>::%face_descriptor` key type and an integer as value type
 * @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
 *
 * @param mesh triangle mesh on which approximate convex segmentation is computed
 * @param face_ids property map with per face segment ids
 * @param concavity_threshold maximal concavity value of the set of faces in a segment
 * @param np optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
 *
 * \cgalNamedParamsBegin
 *    \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `mesh` \cgalParamEnd
 *    \cgalParamBegin{geom_traits} a geometric traits class instance model of Kernel \cgalParamEnd
 *    \cgalParamBegin{min_number_of_segments} minimal number of segments, default is 1 \cgalParamEnd
 *    \cgalParamBegin{convex_hulls_of_segments} a property map model of mutable `LvaluePropertyMap` with a key convertible to `std::size_t`
 *                                              and a model of `MutableFaceGraph` as value type. A good candidate is `boost::vector_property_map<TriangleMesh>`.
 *    \cgalParamEnd
 *    \cgalParamBegin{use_closest_point} if true, the concavity at each vertex is evaluated by using the distance to the closest point on the convex hull of the set of faces. If false, the distance to the first intersected point following the normal at each vertex is used. Default is false \cgalParamEnd
 *    \cgalParamBegin{segment_size_threshold} a bound on the length of the diagonal of the bounding box of the segments. It is a percentage of the length of the diagonal
 *                                            of the input mesh (valid values are within `[0,100`]). If different from `0`, any produced segment that has a bounding box
 *                                            with a diagonal shortest than the specified ratio will be merged with a neighbour segment regardless of the concavity constraint
 *                                            (but still minimizing the concavity while choosing). Default is 0%
 *    \cgalParamEnd
 * \cgalNamedParamsEnd
 *
 * @return number of segments computed
 */
template <class ConcurrencyTag, class TriangleMesh, class FacePropertyMap, class NamedParameters = parameters::Default_named_parameters>
std::size_t
approximate_convex_segmentation(const TriangleMesh& mesh,
                                FacePropertyMap face_ids,
                                double concavity_threshold,
                                const NamedParameters& np = parameters::default_values())
{
  using Vpm = typename CGAL::GetVertexPointMap<TriangleMesh, NamedParameters>::const_type;
  using Geom_traits = typename GetGeomTraits<TriangleMesh, NamedParameters>::type;

  using CGAL::parameters::get_parameter;
  using CGAL::parameters::choose_parameter;

  Geom_traits geom_traits = choose_parameter<Geom_traits>(get_parameter(np, internal_np::geom_traits));
  Vpm vpm = choose_parameter(get_parameter(np, CGAL::internal_np::vertex_point), get_const_property_map(CGAL::vertex_point, mesh));

  std::size_t min_number_of_segments = choose_parameter(get_parameter(np, internal_np::min_number_of_segments), 1);

  CGAL_precondition(is_triangle_mesh(mesh));
  CGAL_precondition(num_faces(mesh) >= min_number_of_segments);

  bool use_closest_point = choose_parameter(get_parameter(np, internal_np::use_closest_point), false);

  internal::Approx_segmentation<TriangleMesh, Vpm, Geom_traits, ConcurrencyTag> algorithm(mesh, vpm, geom_traits, use_closest_point);
  algorithm.segmentize(concavity_threshold, min_number_of_segments);

  double segment_size_threshold = choose_parameter(get_parameter(np, internal_np::segment_size_threshold), 0.);
  if (segment_size_threshold < 0) segment_size_threshold = 0;
  if (segment_size_threshold > 100) segment_size_threshold = 100;

  if (segment_size_threshold != 0)
  {
    algorithm.postprocess(min_number_of_segments, segment_size_threshold, concavity_threshold);
  }

  algorithm.fill_convex_hull_map(get_parameter(np, internal_np::convex_hulls_of_segments));

  return algorithm.result(face_ids);
}

} //namespace CGAL

#endif //CGAL_SURFACE_MESH_SEGMENTATION_APPROXIMATE_CONVEX_DECOMPOSITION_H
