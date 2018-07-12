#ifndef CGAL_SURFACE_MESH_SEGMENTATION_APPROX_DECOMPOSITION_H

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

#define CGAL_SURFACE_MESH_SEGMENTATION_APPROX_DECOMPOSITION_H

#include <CGAL/license/Surface_mesh_segmentation.h>

/**
 * @file approx_decomposition.h
 * @brief The API which contains template functions for concavity value computation and approximate convex decomposition
 */
#include <CGAL/internal/Approximate_convex_decomposition/Approx_decomposition.h>
#include <CGAL/internal/Approximate_convex_decomposition/Concavity.h>
#include <CGAL/boost/graph/named_function_params.h>
#include <CGAL/boost/graph/named_params_helper.h>

namespace CGAL
{


/*!
 * \ingroup PkgSurfaceSegmentation
 * @brief Function computing the concavity value of a cluster described by an id in a triangle mesh.
 *
 * Cluster is a subset of connected faces in a triangle mesh.
 *
 * \note The function relies on the `Face_filtered_graph`, i.e. to compile this function with the `Polyhedron_3` mesh type it must provide indices for its components, for example `Polyhedron_items_with_id_3` can be used for such purpose.
 *
 * @pre `is_triangle_mesh(mesh)`
 *
 * @tparam ConcurrencyTag enables sequential versus parallel algorithm (possible values are `Sequential_tag` and `Parallel_tag`)
 * @tparam TriangleMesh a model of `FaceListGraph`
 * @tparam FacePropertyMap a `ReadWritePropertyMap` with the `boost::graph_traits<TriangleMesh>::%face_descriptor` key type and an integer value type
 * @tparam DistancesPropertyMap a `ReadWritePropertyMap` with the `boost::graph_traits<TriangleMesh>::%vertex_descriptor` key type and a floating-point value type
 * @tparam NamedParameters a sequence of \ref pmp_namedparameters "Named Parameters"
 *
 * @param mesh clustered triangle mesh
 * @param face_ids property map with per face cluster ids
 * @param cluster_id id of the target cluster on which concavity value is computed
 * @param[out] distances optional property map with per vertex distances to the convex hull
 * @param np optional sequence of \ref pmp_namedparameters "Named Parameters" among the ones listed below
 *
 * \cgalNamedParamsBegin
 *    \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `mesh` \cgalParamEnd
 *    \cgalParamBegin{geom_traits} a geometric traits class instance \cgalParamEnd
 * \cgalNamedParamsEnd
 *
 * @return concativity value, the largest distance from a vertex in a cluster to its projection on the convex hull of a triangle mesh
 */
template <class ConcurrencyTag, class TriangleMesh, class FacePropertyMap, class DistancesPropertyMap, class NamedParameters>
double
concavity_value(const TriangleMesh& mesh,
                FacePropertyMap face_ids,
                std::size_t cluster_id,
                DistancesPropertyMap distances,
                const NamedParameters& np)
{
    typedef typename Polygon_mesh_processing::GetVertexPointMap<TriangleMesh, NamedParameters>::const_type Vpm;
    typedef typename Polygon_mesh_processing::GetGeomTraits<TriangleMesh, NamedParameters>::type Geom_traits;

    Geom_traits geom_traits = boost::choose_param(boost::get_param(np, internal_np::geom_traits), Geom_traits());
    Vpm vpm = boost::choose_param(boost::get_param(np, internal_np::vertex_point),
                                  get_const_property_map(boost::vertex_point, mesh));

    CGAL_precondition(is_triangle_mesh(mesh));

    internal::Concavity<TriangleMesh, Vpm, Geom_traits, ConcurrencyTag> algorithm(mesh, vpm, geom_traits);
    return algorithm.compute(face_ids, cluster_id, distances);
}


#ifndef DOXYGEN_RUNNING
template <class ConcurrencyTag, class TriangleMesh, class FacePropertyMap, class DistancesPropertyMap>
double
concavity_value(const TriangleMesh& mesh,
                FacePropertyMap face_ids,
                std::size_t cluster_id,
                DistancesPropertyMap distances)
{
    return concavity_value<ConcurrencyTag>(mesh, face_ids, cluster_id, distances, Polygon_mesh_processing::parameters::all_default());
}

template <class ConcurrencyTag, class TriangleMesh, class FacePropertyMap>
double
concavity_value(const TriangleMesh& mesh,
                FacePropertyMap face_ids,
                std::size_t cluster_id)
{
    typedef boost::unordered_map<typename boost::graph_traits<TriangleMesh>::vertex_descriptor, double> Vertex_double_map;
    Vertex_double_map distances_map;
    boost::associative_property_map<Vertex_double_map> distances_pmap(distances_map);
    
    return concavity_value<ConcurrencyTag>(mesh, face_ids, cluster_id, distances_pmap, Polygon_mesh_processing::parameters::all_default());
}
#endif


/*!
 * \ingroup PkgSurfaceSegmentation
 * @brief Function computing the concavity value of a triangle mesh.
 *
 * @pre `is_triangle_mesh(mesh)`
 *
 * @tparam ConcurrencyTag enables sequential versus parallel algorithm (possible values are `Sequential_tag` and `Parallel_tag`)
 * @tparam TriangleMesh a model of `FaceListGraph`
 * @tparam DistancesPropertyMap a `ReadWritePropertyMap` with the `boost::graph_traits<TriangleMesh>::%vertex_descriptor` key type and a floating-point value type
 * @tparam NamedParameters a sequence of \ref pmp_namedparameters "Named Parameters"
 *
 * @param mesh triangle mesh on which concavity value is computed
 * @param[out] distances optional property map with per vertex distances to the convex hull
 * @param np optional sequence of \ref pmp_namedparameters "Named Parameters" among the ones listed below
 *
 * \cgalNamedParamsBegin
 *    \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `mesh` \cgalParamEnd
 *    \cgalParamBegin{geom_traits} a geometric traits class instance \cgalParamEnd
 * \cgalNamedParamsEnd
 *
 * @return concativity value, the largest distance from a vertex in a cluster to its projection onto the convex hull of a triangle mesh
 */
template <class ConcurrencyTag, class TriangleMesh, class DistancesPropertyMap, class NamedParameters>
double
concavity_value(const TriangleMesh& mesh,
                DistancesPropertyMap distances,
                const NamedParameters& np)
{
    typedef typename Polygon_mesh_processing::GetVertexPointMap<TriangleMesh, NamedParameters>::const_type Vpm;
    typedef typename Polygon_mesh_processing::GetGeomTraits<TriangleMesh, NamedParameters>::type Geom_traits;

    Geom_traits geom_traits = boost::choose_param(boost::get_param(np, internal_np::geom_traits), Geom_traits());
    Vpm vpm = boost::choose_param(boost::get_param(np, internal_np::vertex_point),
                                  get_const_property_map(boost::vertex_point, mesh));

    CGAL_precondition(is_triangle_mesh(mesh));

    internal::Concavity<TriangleMesh, Vpm, Geom_traits, ConcurrencyTag> algorithm(mesh, vpm, geom_traits);
    return algorithm.compute(distances);
}


#ifndef DOXYGEN_RUNNING
template <class ConcurrencyTag, class TriangleMesh, class DistancesPropertyMap>
double
concavity_value(const TriangleMesh& mesh,
                DistancesPropertyMap distances)
{
    return concavity_value<ConcurrencyTag>(mesh, distances, Polygon_mesh_processing::parameters::all_default());
}

template <class ConcurrencyTag, class TriangleMesh>
double
concavity_value(const TriangleMesh& mesh)
{
    typedef boost::unordered_map<typename boost::graph_traits<TriangleMesh>::vertex_descriptor, double> Vertex_double_map;
    Vertex_double_map distances_map;
    boost::associative_property_map<Vertex_double_map> distances_pmap(distances_map);
    
    return concavity_value<ConcurrencyTag>(mesh, distances_pmap, Polygon_mesh_processing::parameters::all_default());
}
#endif


/*!
 * \ingroup PkgSurfaceSegmentation
 * @brief Function computing approximate convex decomposition of a triangle mesh, which is the smallest set of clusters satisfying concavity constraint. 
 *
 * Clusters are subsets of connected faces in a triangle mesh which concavity values are less or equal to 'concavity_threshold'.
 *
 * This function fills a property map which associates a cluster-id (in the range [0, `number_of_clusters`-1]) to each face.
 *
 * \note For verbose output define CGAL_APPROX_DECOMPOSITION_VERBOSE before the include derective.
 *
 * @pre `is_triangle_mesh(mesh)`
 * @pre `num_face(mesh) >= min_number_of_clusters`
 *
 * @tparam ConcurrencyTag enables sequential versus parallel algorithm (possible values are `Sequential_tag` and `Parallel_tag`)
 * @tparam TriangleMesh a model of `FaceListGraph`
 * @tparam FacePropertyMap a `ReadWritePropertyMap` with the `boost::graph_traits<TriangleMesh>::%face_descriptor` key type and an integer value type
 * @tparam NamedParameters a sequence of \ref pmp_namedparameters "Named Parameters"
 *
 * @param mesh triangle mesh on which approximate convex decomposition is computed
 * @param[out] face_ids property map with per face cluster ids
 * @param concavity_threshold maximal concavity value of a cluster constraint
 * @param min_number_of_clusters optional minimal number of clusters constraint (default=1)
 * @param np optional sequence of \ref pmp_namedparameters "Named Parameters" among the ones listed below
 *
 * \cgalNamedParamsBegin
 *    \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `mesh` \cgalParamEnd
 *    \cgalParamBegin{geom_traits} a geometric traits class instance \cgalParamEnd
 * \cgalNamedParamsEnd
 *
 * @return number of clusters
 */
template <class ConcurrencyTag, class TriangleMesh, class FacePropertyMap, class NamedParameters>
std::size_t
approximate_convex_decomposition(const TriangleMesh& mesh,
                                 FacePropertyMap face_ids,
                                 double concavity_threshold,
                                 std::size_t min_number_of_clusters,
                                 const NamedParameters& np)
{
  typedef typename Polygon_mesh_processing::GetVertexPointMap<TriangleMesh, NamedParameters>::const_type Vpm;
  typedef typename Polygon_mesh_processing::GetGeomTraits<TriangleMesh, NamedParameters>::type Geom_traits;

  Geom_traits geom_traits = boost::choose_param(boost::get_param(np, internal_np::geom_traits), Geom_traits());
  Vpm vpm = boost::choose_param(boost::get_param(np, internal_np::vertex_point),
                                get_const_property_map(boost::vertex_point, mesh));

  CGAL_precondition(is_triangle_mesh(mesh));
  CGAL_precondition(num_faces(mesh) >= min_number_of_clusters);

  internal::Approx_decomposition<TriangleMesh, Vpm, Geom_traits, ConcurrencyTag> algorithm(mesh, vpm, geom_traits);
  return algorithm.decompose(face_ids, concavity_threshold, min_number_of_clusters);
}


#ifndef DOXYGEN_RUNNING
template <class ConcurrencyTag, class TriangleMesh, class FacePropertyMap>
std::size_t
approximate_convex_decomposition(const TriangleMesh& mesh,
                                 FacePropertyMap face_ids,
                                 double concavity_threshold,
                                 std::size_t min_number_of_clusters = 1)
{
  return approximate_convex_decomposition<ConcurrencyTag>(mesh, face_ids, concavity_threshold, min_number_of_clusters, Polygon_mesh_processing::parameters::all_default());
}
#endif


} //namespace CGAL

#endif //CGAL_SURFACE_MESH_SEGMENTATION_APPROX_DECOMPOSITION_H
