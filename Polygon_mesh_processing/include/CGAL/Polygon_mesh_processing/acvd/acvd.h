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
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>


#include <numeric>
#include <vector>
#include <queue>
#include <unordered_set>
#include <iostream>

namespace CGAL {

namespace Polygon_mesh_processing {

namespace internal {

template <typename GT>
struct ClusterData {
  typename GT::Vector_3 site_sum;
  typename GT::FT weight_sum;
  typename GT::FT energy;

  ClusterData() : site_sum(0, 0, 0), weight_sum(0), energy(0) {}

  void add_vertex(const typename GT::Vector_3 vertex_position, const typename GT::FT weight)
  {
    this->site_sum += vertex_position;
    this->weight_sum += weight;
  }

  void remove_vertex(const typename GT::Vector_3 vertex_position, const typename GT::FT weight)
  {
    this->site_sum -= vertex_position;
    this->weight_sum -= weight;
  }

  typename GT::FT compute_energy()
  {
    this->energy = - (this->site_sum).squared_length() / this->weight_sum;
    return this->energy;
  }
};

template <typename PolygonMesh, /*ClusteringMetric,*/
  typename NamedParameters = parameters::Default_named_parameters>
void acvd_simplification(
    PolygonMesh& pmesh,
    const int& nb_clusters,
    const NamedParameters& np = parameters::default_values()
    // seed_randomization can be a named parameter
  )
{
  typedef typename GetGeomTraits<PolygonMesh, NamedParameters>::type GT;
  typedef typename GetVertexPointMap<PolygonMesh, NamedParameters>::const_type Vertex_position_map;
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor Halfedge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor Vertex_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::face_descriptor Face_descriptor;
  typedef typename boost::property_map<PolygonMesh, CGAL::dynamic_vertex_property_t<CGAL::IO::Color> >::type VertexColorMap;
  typedef typename boost::property_map<PolygonMesh, CGAL::dynamic_vertex_property_t<int> >::type VertexClusterMap;
  typedef typename boost::property_map<PolygonMesh, CGAL::dynamic_vertex_property_t<typename GT::FT> >::type VertexWeightMap;

  using parameters::choose_parameter;
  using parameters::get_parameter;
  using parameters::is_default_parameter;

  Vertex_position_map vpm = choose_parameter(get_parameter(np, CGAL::vertex_point),
    get_property_map(CGAL::vertex_point, pmesh));

  // initial random clusters
  // property map from vertex_descriptor to cluster index
  VertexClusterMap vertex_clusters_pmap = get(CGAL::dynamic_vertex_property_t<int>(), pmesh);
  VertexWeightMap vertex_weight_pmap = get(CGAL::dynamic_vertex_property_t<typename GT::FT>(), pmesh);
  std::vector<ClusterData<GT>> clusters(nb_clusters + 1);
  std::queue<Halfedge_descriptor> clusters_edges_active;
  std::queue<Halfedge_descriptor> clusters_edges_new;

  int nb_vertices = num_vertices(pmesh);

  // compute vertex weights (dual area)
  for (Face_descriptor fd : faces(pmesh))
  {
    typename GT::FT weight = CGAL::Polygon_mesh_processing::face_area(fd, pmesh) / 3;
    for (Vertex_descriptor vd : vertices_around_face(halfedge(fd, pmesh), pmesh))
    {
      put(vertex_weight_pmap, vd, weight);
    }
  }

  srand(3);
  for (int ci = 0; ci < nb_clusters; ci++)
  {
    // random index
    int vi = rand() % num_vertices(pmesh);
    Vertex_descriptor vd = *(vertices(pmesh).begin() + vi);
    /*int vi;
    Vertex_descriptor vd;
    do {
      vi = ci * nb_vertices / nb_clusters;
      vd = *(vertices(pmesh).begin() + vi);
    } while (get(vertex_clusters_pmap, vd));*/

    // TODO: check for cluster conflict at the same vertex
    put(vertex_clusters_pmap, vd, ci + 1); // TODO: should be ci but for now we start from 1 (can't set null value to -1)
    typename GT::Point_3 vp = get(vpm, vd);
    typename GT::Vector_3 vpv(vp.x(), vp.y(), vp.z());
    clusters[ci].add_vertex(vpv, get(vertex_weight_pmap, vd));

    for (Halfedge_descriptor hd : halfedges_around_source(halfedge(vd, pmesh), pmesh))
      clusters_edges_active.push(hd);
  }

  int nb_modifications;

  do
  {
    nb_modifications = 0;

    while (clusters_edges_active.empty() == false) {
      Halfedge_descriptor hi = clusters_edges_active.front();
      clusters_edges_active.pop();

      Vertex_descriptor v1 = source(hi, pmesh);
      Vertex_descriptor v2 = target(hi, pmesh);

      int c1 = get(vertex_clusters_pmap, v1);
      int c2 = get(vertex_clusters_pmap, v2);

      if (c1 == 0)
      {
        // expand cluster c2 (add v1 to c2)
        put(vertex_clusters_pmap, v1, c2);
        typename GT::Point_3 vp1 = get(vpm, v1);
        typename GT::Vector_3 vpv(vp1.x(), vp1.y(), vp1.z());
        clusters[c2].add_vertex(vpv, get(vertex_weight_pmap, v1));
        clusters[c1].remove_vertex(vpv, get(vertex_weight_pmap, v1));


        // add all halfedges around v1 except hi to the queue
        for (Halfedge_descriptor hd : halfedges_around_source(halfedge(v1, pmesh), pmesh))
          if (hd != hi)
            clusters_edges_new.push(hd);
        nb_modifications++;

      }
      else if (c2 == 0)
      {
        // expand cluster c1 (add v2 to c1)
        put(vertex_clusters_pmap, v2, c1);
        typename GT::Point_3 vp2 = get(vpm, v2);
        typename GT::Vector_3 vpv(vp2.x(), vp2.y(), vp2.z());
        clusters[c1].add_vertex(vpv, get(vertex_weight_pmap, v2));
        clusters[c2].remove_vertex(vpv, get(vertex_weight_pmap, v2));


        // add all halfedges around v2 except hi to the queue
        for (Halfedge_descriptor hd : halfedges_around_source(halfedge(v2, pmesh), pmesh))
          if (hd != hi)
            clusters_edges_new.push(hd);
        nb_modifications++;
      }
      else if (c1 == c2)
      {
        clusters_edges_new.push(hi);
        continue; // no modification
      }
      else
      {
        // compare the energy of the 3 cases
        typename GT::Point_3 vp1 = get(vpm, v1);
        typename GT::Vector_3 vpv1(vp1.x(), vp1.y(), vp1.z());
        typename GT::Point_3 vp2 = get(vpm, v2);
        typename GT::Vector_3 vpv2(vp2.x(), vp2.y(), vp2.z());

        typename GT::FT e_no_change = clusters[c1].compute_energy() + clusters[c2].compute_energy();

        clusters[c1].remove_vertex(vpv1, get(vertex_weight_pmap, v1));
        clusters[c2].add_vertex(vpv1, get(vertex_weight_pmap, v1));

        typename GT::FT e_v1_to_c2 = clusters[c1].compute_energy() + clusters[c2].compute_energy();

        // reset to no change
        clusters[c1].add_vertex(vpv1, get(vertex_weight_pmap, v1));
        clusters[c2].remove_vertex(vpv1, get(vertex_weight_pmap, v1));

        // The effect of the following should always be reversed after the comparison
        clusters[c2].remove_vertex(vpv2, get(vertex_weight_pmap, v2));
        clusters[c1].add_vertex(vpv2, get(vertex_weight_pmap, v2));

        typename GT::FT e_v2_to_c1 = clusters[c1].compute_energy() + clusters[c2].compute_energy();

        if (e_v2_to_c1 < e_no_change && e_v2_to_c1 < e_v1_to_c2)
        {
          // move v2 to c1
          put(vertex_clusters_pmap, v2, c1);

          // cluster data is already updated

          // add all halfedges around v2 except hi to the queue
          for (Halfedge_descriptor hd : halfedges_around_source(halfedge(v2, pmesh), pmesh))
            if (hd != hi)
              clusters_edges_new.push(hd);
          nb_modifications++;
        }
        else if (e_v1_to_c2 < e_no_change)
        {
          // move v1 to c2
          put(vertex_clusters_pmap, v1, c2);

          // need to reset cluster data and then update
          clusters[c2].add_vertex(vpv2, get(vertex_weight_pmap, v2));
          clusters[c1].remove_vertex(vpv2, get(vertex_weight_pmap, v2));

          clusters[c1].remove_vertex(vpv1, get(vertex_weight_pmap, v1));
          clusters[c2].add_vertex(vpv1, get(vertex_weight_pmap, v1));

          // add all halfedges around v1 except hi to the queue
          for (Halfedge_descriptor hd : halfedges_around_source(halfedge(v1, pmesh), pmesh))
            if (hd != hi)
              clusters_edges_new.push(hd);
          nb_modifications++;
        }
        else
        {
            // no change but need to reset cluster data
            clusters[c2].add_vertex(vpv2, get(vertex_weight_pmap, v2));
            clusters[c1].remove_vertex(vpv2, get(vertex_weight_pmap, v2));

            clusters_edges_new.push(hi);
        }

        continue;
      }
    }
    clusters_edges_active.swap(clusters_edges_new);
  } while (nb_modifications > 0);


  VertexColorMap vcm = get(CGAL::dynamic_vertex_property_t<CGAL::IO::Color>(), pmesh);

  for (Vertex_descriptor vd : vertices(pmesh))
  {
    int c = get(vertex_clusters_pmap, vd);
    CGAL::IO::Color color(255 - (c * 255 / nb_clusters), (c * c % 7) * 255 / 7, (c * c * c % 31) * 255 / 31);
    std::cout << vd.idx() << " " << c << " " << color << std::endl;
    put(vcm, vd, color);
  }
  std::cout << "kak1" << std::endl;
  CGAL::IO::write_OFF("S52k_clustered_0.off", pmesh, CGAL::parameters::vertex_color_map(vcm));
  std::cout << "kak2" << std::endl;

}


} // namespace internal

template <typename PolygonMesh,
  typename NamedParameters = parameters::Default_named_parameters>
void  acvd_isotropic_simplification(
    PolygonMesh& pmesh,
    const int& nb_vertices,
    const NamedParameters& np = parameters::default_values()
  )
{
  internal::acvd_simplification<PolygonMesh, /*IsotropicMetric,*/ NamedParameters>(
    pmesh,
    nb_vertices,
    np
    );
}

} // namespace Polygon_mesh_processing

} // namespace CGAL
