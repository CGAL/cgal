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
// SPDX-License-Identifier: GPL-3.0+
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
#include <CGAL/boost/graph/named_function_params.h>
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
 * @tparam FacePropertyMap a `ReadWritePropertyMap` with the `boost::graph_traits<TriangleMesh>::%face_descriptor` key type and an integer value type
 * @tparam DistancesPropertyMap a `ReadWritePropertyMap` with the `boost::graph_traits<TriangleMesh>::%vertex_descriptor` key type and a floating-point value type
 * @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
 *
 * @param mesh segmented triangle mesh
 * @param face_ids property map with per face segment ids
 * @param segment_id id of the target segment on which concavity value is computed
 * @param[out] distances optional property map with per vertex distances to the convex hull
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
    typedef typename Polygon_mesh_processing::GetVertexPointMap<TriangleMesh, NamedParameters>::const_type Vpm;
    typedef typename Polygon_mesh_processing::GetGeomTraits<TriangleMesh, NamedParameters>::type Geom_traits;

    Geom_traits geom_traits = boost::choose_param(boost::get_param(np, internal_np::geom_traits), Geom_traits());
    Vpm vpm = boost::choose_param(boost::get_param(np, internal_np::vertex_point),
                                  get_const_property_map(boost::vertex_point, mesh));

    CGAL_precondition(is_triangle_mesh(mesh));
  
    bool use_closest_point = boost::choose_param(boost::get_param(np, internal_np::use_closest_point), false);

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
    CGAL::Static_property_map<typename boost::graph_traits<TriangleMesh>::vertex_descriptor, double > distances_pmap(0);
    
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
 * @tparam DistancesPropertyMap a `ReadWritePropertyMap` with the `boost::graph_traits<TriangleMesh>::%vertex_descriptor` key type and a floating-point as value type
 * @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
 *
 * @param mesh triangle mesh on which concavity value is computed
 * @param distances optional property map with per vertex distances to the convex hull
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
    typedef typename Polygon_mesh_processing::GetVertexPointMap<TriangleMesh, NamedParameters>::const_type Vpm;
    typedef typename Polygon_mesh_processing::GetGeomTraits<TriangleMesh, NamedParameters>::type Geom_traits;

    Geom_traits geom_traits = boost::choose_param(boost::get_param(np, internal_np::geom_traits), Geom_traits());
    Vpm vpm = boost::choose_param(boost::get_param(np, internal_np::vertex_point),
                                  get_const_property_map(boost::vertex_point, mesh));

    CGAL_precondition(is_triangle_mesh(mesh));
  
    bool use_closest_point = boost::choose_param(boost::get_param(np, internal_np::use_closest_point), false);

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
    CGAL::Static_property_map<typename boost::graph_traits<TriangleMesh>::vertex_descriptor, double > distances_pmap(0);
    
    return concavity_values<ConcurrencyTag>(mesh, distances_pmap, Polygon_mesh_processing::parameters::all_default());
}
#endif


/*!
 * \ingroup PkgSurfaceSegmentation
 * @brief computes an approximate convex segmentation of a triangle mesh.
 *
 * This function fills a property map associating a segment-id to each face (in the range [0, `N-1`]).
 * `N` is the number of segments computed by the functions (greater are equal to `min_number_of_segments`).
 * The set of faces in the same segment defines an edge-connected patch with concavity value less are equal to 'concavity_threshold'.
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
 *    \cgalParamBegin{segments_convex_hulls} an array filled up with the convex hulls of all segments \cgalParamEnd
 *    \cgalParamBegin{use_closest_point} if true, the concavity at each vertex is evaluated by using the distance to the closest point on the convex hull of the set of faces. If false, the distance to the first intersected point following the normal at each vertex is used. Default is false \cgalParamEnd
 *    \cgalParamBegin{postprocess_segments} if true, any produced segment that is smaller than `small_segment_threshold` will be merged with a neighbour segment regardless the concavity constraint). Default is false \cgalParamEnd
 *    \cgalParamBegin{small_segment_threshold} the minimal size of a segment postprocessing procedure is allowed to merge in percentage with regard to the length of diagonal of the input mesh. The value must be in the range [0, 100]. Default is 10.0% \cgalParamEnd
 * \cgalNamedParamsEnd
 *
 * @return number of segments computed
 */
template <class ConcurrencyTag, class TriangleMesh, class FacePropertyMap, class NamedParameters>
std::size_t
approximate_convex_segmentation(const TriangleMesh& mesh,
                                FacePropertyMap face_ids,
                                double concavity_threshold,
                                const NamedParameters& np)
{
  typedef typename Polygon_mesh_processing::GetVertexPointMap<TriangleMesh, NamedParameters>::const_type Vpm;
  typedef typename Polygon_mesh_processing::GetGeomTraits<TriangleMesh, NamedParameters>::type Geom_traits;

  Geom_traits geom_traits = boost::choose_param(boost::get_param(np, internal_np::geom_traits), Geom_traits());
  Vpm vpm = boost::choose_param(boost::get_param(np, internal_np::vertex_point),
                                get_const_property_map(boost::vertex_point, mesh));
  std::size_t min_number_of_segments = boost::choose_param(boost::get_param(np, internal_np::min_number_of_segments), 1);

  CGAL_precondition(is_triangle_mesh(mesh));
  CGAL_precondition(num_faces(mesh) >= min_number_of_segments);

  boost::vector_property_map<TriangleMesh> convex_hulls_pmap = boost::choose_param(boost::get_param(np, internal_np::segments_convex_hulls), boost::vector_property_map<TriangleMesh>());
  bool use_closest_point = boost::choose_param(boost::get_param(np, internal_np::use_closest_point), false);

  internal::Approx_segmentation<TriangleMesh, Vpm, Geom_traits, ConcurrencyTag> algorithm(mesh, vpm, geom_traits, use_closest_point);
  algorithm.segmentize(concavity_threshold, min_number_of_segments);

  bool postprocess_segments = boost::choose_param(boost::get_param(np, internal_np::postprocess_segments), false);
  double small_segment_threshold = boost::choose_param(boost::get_param(np, internal_np::small_segment_threshold), 10.);

  if (postprocess_segments)
  {
    algorithm.postprocess(min_number_of_segments, small_segment_threshold, concavity_threshold);
  }

  return algorithm.result(face_ids, convex_hulls_pmap);
}


#ifndef DOXYGEN_RUNNING
template <class ConcurrencyTag, class TriangleMesh, class FacePropertyMap>
std::size_t
approximate_convex_segmentation(const TriangleMesh& mesh,
                                FacePropertyMap face_ids,
                                double concavity_threshold)
{
  return approximate_convex_segmentation<ConcurrencyTag>(mesh, face_ids, concavity_threshold, Polygon_mesh_processing::parameters::all_default());
}
#endif


} //namespace CGAL

#endif //CGAL_SURFACE_MESH_SEGMENTATION_APPROXIMATE_CONVEX_DECOMPOSITION_H
