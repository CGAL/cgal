// Copyright (c) 2023 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Hossam Saeed
//

// #ifndef CGAL_POLYGON_MESH_PROCESSING_<>
// #define CGAL_POLYGON_MESH_PROCESSING_<>

// #include <CGAL/license/Polygon_mesh_processing/<>>

#include <CGAL/assertions.h>
#include <CGAL/IO/Color.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Named_function_parameters.h>
#include <CGAL/property_map.h>
#include <CGAL/boost/graph/named_params_helper.h>
#include <Eigen/Eigenvalues>

#include <numeric>
#include <vector>
#include <queue>
#include <unordered_set>
#include <iostream>

namespace CGAL {

namespace Polygon_mesh_processing {

namespace internal {

template <typename GT>
struct IsotropicMetricCluster {
  typename GT::Vector_3 site_sum;
  typename GT::FT weight;
  typename GT::FT energy;

  IsotropicMetricCluster() : site_sum(0, 0, 0), weight(0), energy(0) {}

  void add_vertex(const typename GT::Vector_3& vertex_position, const typename GT::FT& weight = 1)
  {
    site_sum += vertex_position * weight;
    this->weight += weight;
  }

  void remove_vertex(const typename GT::Vector_3& vertex_position, const typename GT::FT& weight = 1)
  {
    site_sum -= vertex_position * weight;
    this->weight -= weight;
  }
};

template <typename PolygonMesh, /*ClusteringMetric,*/
            typename NamedParameters = parameters::Default_named_parameters>
typename boost::property_map<PolygonMesh, CGAL::dynamic_vertex_property_t<CGAL::IO::Color> >::type
    acvd_simplification(
    PolygonMesh& pmesh,
    const int& nb_vertices,
    const NamedParameters& np = parameters::default_values()
    // seed_randomization can be a named parameter
)
{
    typedef typename GetGeomTraits<PolygonMesh, NamedParameters>::type GT;
    typedef typename GetVertexPointMap<PolygonMesh, NamedParameters>::const_type Vertex_position_map;
    typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor Halfedge_descriptor;
    typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor Vertex_descriptor;
    typedef typename boost::property_map<PolygonMesh, CGAL::dynamic_vertex_property_t<CGAL::IO::Color> >::type VertexColorMap;
    typedef typename boost::property_map<PolygonMesh, CGAL::dynamic_vertex_property_t<int> >::type VertexClusterMap;


    using parameters::choose_parameter;
    using parameters::get_parameter;
    using parameters::is_default_parameter;

    Vertex_position_map vpm = choose_parameter(get_parameter(np, CGAL::vertex_point),
        get_property_map(CGAL::vertex_point, pmesh));

    // initial random clusters
    // property map from vertex_descriptor to cluster index
    boost::property_map<PolygonMesh, CGAL::dynamic_vertex_property_t<int>>::type
        vertex_clusters_pmap = get(CGAL::dynamic_vertex_property_t<int>(), pmesh);
    VertexClusterMap clusters = get(CGAL::dynamic_vertex_property_t<int>(), pmesh);
    std::vector<IsotropicMetricCluster<GT>> clusters_sites(nb_vertices);
    std::queue<Halfedge_descriptor> clusters_edges;

    srand(time(NULL));
    for(int ci = 0; ci < nb_vertices; ci++)
    {
        // random index
        int vi = rand() % num_vertices(pmesh);
        Vertex_descriptor vd = *(vertices(pmesh).begin() + vi);
        // TODO: check for cluster conflict at the same vertex
        put(vertex_clusters_pmap, vd, ci + 1); // TODO: should be ci but for now we start from 1 (can't set null value to -1)
        typename GT::Point_3 vp = get(vpm, vd);
        typename GT::Vector_3 vpv(vp.x(), vp.y(), vp.z());
        clusters_sites[ci].add_vertex(vpv);
        for (Halfedge_descriptor hd : halfedges_around_source(halfedge(vd, pmesh), pmesh))
            clusters_edges.push(hd);
    }

    // minimize the energy
    int nb_modifications = 0;
    int prev_nb_modifications = 0;

    while (nb_modifications != prev_nb_modifications) {
        Halfedge_descriptor hi = clusters_edges.front();
        clusters_edges.pop();

        Vertex_descriptor v1 = source(hi, pmesh);
        Vertex_descriptor v2 = target(hi, pmesh);

        int c1 = get(vertex_clusters_pmap, v1);
        int c2 = get(vertex_clusters_pmap, v2);

        prev_nb_modifications = nb_modifications;

        if (c1 == 0)
        {
            // expand cluster c2 (add v1 to c2)
            put(vertex_clusters_pmap, v1, c2);
            typename GT::Point_3 vp = get(vpm, v1);
            typename GT::Vector_3 vpv(vp.x(), vp.y(), vp.z());
            clusters_sites[c2].add_vertex(vpv);

            // add all halfedges around v1 except hi to the queue
            for (Halfedge_descriptor hd : halfedges_around_source(halfedge(v1, pmesh), pmesh))
                if (hd != hi)
                    clusters_edges.push(hd);
            nb_modifications++;

        }
        else if (c2 == 0)
        {
            // expand cluster c1 (add v2 to c1)
            put(vertex_clusters_pmap, v2, c1);
            typename GT::Point_3 vp = get(vpm, v2);
            typename GT::Vector_3 vpv(vp.x(), vp.y(), vp.z());
            clusters_sites[c1].add_vertex(vpv);
            // add all halfedges around v2 except hi to the queue
            for (Halfedge_descriptor hd : halfedges_around_source(halfedge(v2, pmesh), pmesh))
                if (hd != hi)
                    clusters_edges.push(hd);
            nb_modifications++;
        }
        else if (c1 == c2)
            continue; // no modification
        else
        {
            // compare the energy of the 3 cases
            continue;
        }
    }

    VertexColorMap vcm = get(CGAL::dynamic_vertex_property_t<CGAL::IO::Color>(), pmesh);

    for (Vertex_descriptor vd : vertices(pmesh))
    {
        int c = get(vertex_clusters_pmap, vd);
        if (!c) c = 0;
        std::cout << c << std::endl;
        CGAL::IO::Color color(255 - (c * 255 / nb_vertices), (c % 10) * 255 / 10, (c % 50) * 255 / 50);
        put(vcm, vd, color);
    }

    return vcm;
}


} // namespace internal

template <typename PolygonMesh,
          typename NamedParameters = parameters::Default_named_parameters>
  typename boost::property_map<PolygonMesh, CGAL::dynamic_vertex_property_t<CGAL::IO::Color> >::type
  acvd_isotropic_simplification(
    PolygonMesh& pmesh,
    const int& nb_vertices,
    const NamedParameters& np = parameters::default_values()
)
{
    return internal::acvd_simplification<PolygonMesh, /*IsotropicMetric,*/ NamedParameters>(
        pmesh,
        nb_vertices,
        np
    );
}

} // namespace Polygon_mesh_processing

} // namespace CGAL
