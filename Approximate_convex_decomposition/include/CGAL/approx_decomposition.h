#ifndef CGAL_APPROXIMATE_CONVEX_DECOMPOSITION_APPROX_DECOMPOSITION_H

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

#define CGAL_APPROXIMATE_CONVEX_DECOMPOSITION_APPROX_DECOMPOSITION_H

#include <CGAL/internal/Approximate_convex_decomposition/approx_decomposition.h>
#include <CGAL/internal/Approximate_convex_decomposition/concavity.h>
#include <CGAL/boost/graph/named_function_params.h>
#include <CGAL/boost/graph/named_params_helper.h>

namespace CGAL
{


/*
 * @brief Function computing concavity value of a cluster.
 *
 * Cluster is a subset of connected faces in a triangle mesh.
 * Concavity value is the largest distance from a vertex in a cluster to the projected point onto the convex hull of a mesh.
 *
 * @return concativity value.
 */
template <class ConcurrencyTag, class TriangleMesh, class FacePropertyMap, class NamedParameters>
double
concavity_value(const TriangleMesh& mesh,
                FacePropertyMap face_ids,
                std::size_t cluster_id,
                const NamedParameters& np)
{
    typedef typename Polygon_mesh_processing::GetVertexPointMap<TriangleMesh, NamedParameters>::const_type Vpm;
    typedef typename Polygon_mesh_processing::GetGeomTraits<TriangleMesh, NamedParameters>::type Geom_traits;

    Geom_traits geom_traits = boost::choose_param(boost::get_param(np, internal_np::geom_traits), Geom_traits());
    Vpm vpm = boost::choose_param(boost::get_param(np, internal_np::vertex_point),
                                  get_const_property_map(boost::vertex_point, mesh));


    CGAL_precondition(is_triangle_mesh(mesh));

    internal::Concavity<TriangleMesh, Vpm, Geom_traits, ConcurrencyTag> algorithm(mesh, vpm, geom_traits);
    return algorithm.compute(face_ids, cluster_id);
}


template <class ConcurrencyTag, class TriangleMesh, class FacePropertyMap>
double
concavity_value(const TriangleMesh& mesh,
                FacePropertyMap face_ids,
                std::size_t cluster_id)
{
    return concavity_value<ConcurrencyTag>(mesh, face_ids, cluster_id, Polygon_mesh_processing::parameters::all_default());
}


/*
 * @brief Function computing concavity value of a triangle mesh.
 *
 * Concavity value is the largest distance from a vertex in a mesh to the projected point onto the convex hull of a mesh.
 *
 * @return concativity value.
 */


template <class ConcurrencyTag, class TriangleMesh, class NamedParameters>
double
concavity_value(const TriangleMesh& mesh,
                const NamedParameters& np)
{
    typedef typename Polygon_mesh_processing::GetVertexPointMap<TriangleMesh, NamedParameters>::const_type Vpm;
    typedef typename Polygon_mesh_processing::GetGeomTraits<TriangleMesh, NamedParameters>::type Geom_traits;

    Geom_traits geom_traits = boost::choose_param(boost::get_param(np, internal_np::geom_traits), Geom_traits());
    Vpm vpm = boost::choose_param(boost::get_param(np, internal_np::vertex_point),
                                  get_const_property_map(boost::vertex_point, mesh));


    CGAL_precondition(is_triangle_mesh(mesh));

    internal::Concavity<TriangleMesh, Vpm, Geom_traits, ConcurrencyTag> algorithm(mesh, vpm, geom_traits);
    return algorithm.compute();
}


template <class ConcurrencyTag, class TriangleMesh>
double
concavity_value(const TriangleMesh& mesh)
{
    return concavity_value<ConcurrencyTag>(mesh, Polygon_mesh_processing::parameters::all_default());
}


/*
 * @brief Function computing the approximate convex decomposition of a triangle mesh.
 *
 * This function fills a property map which associates a cluster-id (in the range [0, 'number_of_clusters'-1]) to each face.
 * Clusters are subsets of connected faces in a triangle mesh which concavity values are less than 'concavity_threshold'.
 *
 * @return number of clusters.
 */
template <class ConcurrencyTag, class TriangleMesh, class FacePropertyMap, class NamedParameters>
std::size_t
convex_decomposition(const TriangleMesh& mesh,
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


template <class ConcurrencyTag, class TriangleMesh, class FacePropertyMap>
std::size_t
convex_decomposition(const TriangleMesh& mesh,
                     FacePropertyMap face_ids,
                     double concavity_threshold = 0.05,
                     std::size_t min_number_of_clusters = 1)
{
    return convex_decomposition<ConcurrencyTag>(mesh, face_ids, concavity_threshold, min_number_of_clusters, Polygon_mesh_processing::parameters::all_default());
}


} //namespace CGAL

#endif //CGAL_APPROXIMATE_CONVEX_DECOMPOSITION_APPROX_DECOMPOSITION_H
