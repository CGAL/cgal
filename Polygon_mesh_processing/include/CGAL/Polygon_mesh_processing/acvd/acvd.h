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
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/border.h>

#include <CGAL/Point_set_3/IO.h>
#include <CGAL/Point_set_3.h>

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
    this->site_sum += vertex_position * weight;
    this->weight_sum += weight;
  }

  void remove_vertex(const typename GT::Vector_3 vertex_position, const typename GT::FT weight)
  {
    this->site_sum -= vertex_position * weight;
    this->weight_sum -= weight;
  }

  typename GT::FT compute_energy()
  {
    this->energy = - (this->site_sum).squared_length() / this->weight_sum;
    return this->energy;
  }

  typename GT::Vector_3 compute_centroid()
  {
    if (this->weight_sum > 0)
      return (this->site_sum) / this->weight_sum;
    else
      return typename GT::Vector_3 (-1, -1, -1); // TODO: Change this
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
  typedef typename boost::property_map<PolygonMesh, CGAL::dynamic_halfedge_property_t<bool> >::type HalfedgeVisitedMap;

  using parameters::choose_parameter;
  using parameters::get_parameter;
  using parameters::is_default_parameter;

  Vertex_position_map vpm = choose_parameter(get_parameter(np, CGAL::vertex_point),
    get_property_map(CGAL::vertex_point, pmesh));

  // initial random clusters
  // property map from vertex_descriptor to cluster index
  VertexClusterMap vertex_cluster_pmap = get(CGAL::dynamic_vertex_property_t<int>(), pmesh);
  VertexWeightMap vertex_weight_pmap = get(CGAL::dynamic_vertex_property_t<typename GT::FT>(), pmesh);
  std::vector<ClusterData<GT>> clusters(nb_clusters);
  std::queue<Halfedge_descriptor> clusters_edges_active;
  std::queue<Halfedge_descriptor> clusters_edges_new;

  int nb_vertices = num_vertices(pmesh);

  // initialize vertex weights and clusters
  for (Vertex_descriptor vd : vertices(pmesh))
  {
    put(vertex_weight_pmap, vd, 0);
    put(vertex_cluster_pmap, vd, -1);
  }

  //typename GT::FT max_area = 0;
  //typename GT::FT min_area = std::numeric_limits<typename GT::FT>::max();

  // compute vertex weights (dual area)
  for (Face_descriptor fd : faces(pmesh))
  {
    typename GT::FT weight = abs(CGAL::Polygon_mesh_processing::face_area(fd, pmesh)) / 3;

    //max_area = std::max(max_area, weight * 8.0);
    //min_area = std::min(min_area, weight * 3.0);

    for (Vertex_descriptor vd : vertices_around_face(halfedge(fd, pmesh), pmesh))
    {
      typename GT::FT vertex_weight = get(vertex_weight_pmap, vd);
      vertex_weight += weight;
      put(vertex_weight_pmap, vd, vertex_weight);
    }
  }

  // srand(3);
  //srand(time(NULL));
  double avg_rand = 0;
  for (int ci = 0; ci < nb_clusters; ci++)
  {
    //// random index
    //int vi = rand() % num_vertices(pmesh);
    //Vertex_descriptor vd = *(vertices(pmesh).begin() + vi);
    int vi;
    Vertex_descriptor vd;
    do {
      vi = CGAL::get_default_random().get_int(0, num_vertices(pmesh));
      avg_rand += vi;
      vd = *(vertices(pmesh).begin() + vi);
    } while (get(vertex_cluster_pmap, vd) != -1);

    put(vertex_cluster_pmap, vd, ci);
    typename GT::Point_3 vp = get(vpm, vd);
    typename GT::Vector_3 vpv(vp.x(), vp.y(), vp.z());
    clusters[ci].add_vertex(vpv, get(vertex_weight_pmap, vd));

    for (Halfedge_descriptor hd : halfedges_around_source(vd, pmesh))
      clusters_edges_active.push(hd);

    /*if (ci % (nb_clusters / 5) == 0)
      std::cout << "rand ci" << ci << " " << vi << "\n";*/
  }
  avg_rand = avg_rand / nb_clusters;

  std::cout << "avg_rand: " << avg_rand << " nVertices: " << num_vertices(pmesh) << "\n";

  // frequency of each cluster
  std::vector<int> cluster_frequency (nb_clusters, 0);

  for (Vertex_descriptor vd : vertices(pmesh))
  {
    int c = get(vertex_cluster_pmap, vd);
    if (c != -1)  cluster_frequency[c]++;
  }

  int nb_empty = 0;
  for (int i = 0; i < nb_clusters; i++)
  {
    if (cluster_frequency[i] == 0)
    {
      nb_empty++;
    }
  }

  std::cout << "nb_empty before: " << nb_empty << std::endl;

  int nb_modifications = 0;

  do
  {
    nb_modifications = 0;

    while (clusters_edges_active.empty() == false) {
      Halfedge_descriptor hi = clusters_edges_active.front();
      clusters_edges_active.pop();

      Vertex_descriptor v1 = source(hi, pmesh);
      Vertex_descriptor v2 = target(hi, pmesh);

      int c1 = get(vertex_cluster_pmap, v1);
      int c2 = get(vertex_cluster_pmap, v2);

      if (c1 == -1)
      {
        // expand cluster c2 (add v1 to c2)
        put(vertex_cluster_pmap, v1, c2);
        typename GT::Point_3 vp1 = get(vpm, v1);
        typename GT::Vector_3 vpv(vp1.x(), vp1.y(), vp1.z());
        clusters[c2].add_vertex(vpv, get(vertex_weight_pmap, v1));

        // add all halfedges around v1 except hi to the queue
        for (Halfedge_descriptor hd : halfedges_around_source(v1, pmesh))
          //if (hd != hi && hd != opposite(hi, pmesh))
            clusters_edges_new.push(hd);
        nb_modifications++;
      }
      else if (c2 == -1)
      {
        // expand cluster c1 (add v2 to c1)
        put(vertex_cluster_pmap, v2, c1);
        typename GT::Point_3 vp2 = get(vpm, v2);
        typename GT::Vector_3 vpv(vp2.x(), vp2.y(), vp2.z());
        clusters[c1].add_vertex(vpv, get(vertex_weight_pmap, v2));

        // add all halfedges around v2 except hi to the queue
        for (Halfedge_descriptor hd : halfedges_around_source(v2, pmesh))
          //if (hd != hi && hd != opposite(hi, pmesh))
            clusters_edges_new.push(hd);
        nb_modifications++;
      }
      else if (c1 == c2)
      {
        clusters_edges_new.push(hi);
      }
      else
      {
        // compare the energy of the 3 cases
        typename GT::Point_3 vp1 = get(vpm, v1);
        typename GT::Vector_3 vpv1(vp1.x(), vp1.y(), vp1.z());
        typename GT::Point_3 vp2 = get(vpm, v2);
        typename GT::Vector_3 vpv2(vp2.x(), vp2.y(), vp2.z());
        typename GT::FT v1_weight = get(vertex_weight_pmap, v1);
        typename GT::FT v2_weight = get(vertex_weight_pmap, v2);

        typename GT::FT e_no_change = clusters[c1].compute_energy() + clusters[c2].compute_energy();

        clusters[c1].remove_vertex(vpv1, v1_weight);
        clusters[c2].add_vertex(vpv1, v1_weight);

        typename GT::FT e_v1_to_c2 = clusters[c1].compute_energy() + clusters[c2].compute_energy();

        typename GT::FT c1_weight_threshold = clusters[c1].weight_sum;

        // reset to no change
        clusters[c1].add_vertex(vpv1, v1_weight);
        clusters[c2].remove_vertex(vpv1, v1_weight);

        // The effect of the following should always be reversed after the comparison
        clusters[c2].remove_vertex(vpv2, v2_weight);
        clusters[c1].add_vertex(vpv2, v2_weight);

        typename GT::FT e_v2_to_c1 = clusters[c1].compute_energy() + clusters[c2].compute_energy();

        typename GT::FT c2_weight_threshold = clusters[c2].weight_sum;


        if (e_v2_to_c1 < e_no_change && e_v2_to_c1 < e_v1_to_c2 && c2_weight_threshold > 0)
        {
          // move v2 to c1
          put(vertex_cluster_pmap, v2, c1);

          // cluster data is already updated

          // add all halfedges around v2 except hi to the queue
          for (Halfedge_descriptor hd : halfedges_around_source(v2, pmesh))
            //if (hd != hi && hd != opposite(hi, pmesh))
              clusters_edges_new.push(hd);
          nb_modifications++;
        }
        else if (e_v1_to_c2 < e_no_change && c1_weight_threshold > 0)
        {
          // move v1 to c2
          put(vertex_cluster_pmap, v1, c2);

          // need to reset cluster data and then update
          clusters[c2].add_vertex(vpv2, v2_weight);
          clusters[c1].remove_vertex(vpv2, v2_weight);

          clusters[c1].remove_vertex(vpv1, v1_weight);
          clusters[c2].add_vertex(vpv1, v1_weight);

          // add all halfedges around v1 except hi to the queue
          for (Halfedge_descriptor hd : halfedges_around_source(halfedge(v1, pmesh), pmesh))
            //if (hd != hi && hd != opposite(hi, pmesh))
              clusters_edges_new.push(hd);
          nb_modifications++;
        }
        else
        {
            // no change but need to reset cluster data
            clusters[c2].add_vertex(vpv2, v2_weight);
            clusters[c1].remove_vertex(vpv2, v2_weight);

            clusters_edges_new.push(hi);
        }
      }
    }
    clusters_edges_active.swap(clusters_edges_new);
  } while (nb_modifications > 0);

  VertexColorMap vcm = get(CGAL::dynamic_vertex_property_t<CGAL::IO::Color>(), pmesh);

  // frequency of each cluster
  cluster_frequency = std::vector<int>(nb_clusters, 0);

  for (Vertex_descriptor vd : vertices(pmesh))
  {
    int c = get(vertex_cluster_pmap, vd);
    cluster_frequency[c]++;
    //CGAL::IO::Color color((c - min_area) * 255 / (max_area - min_area), 0, 0);
    CGAL::IO::Color color(255 - (c * 255 / nb_clusters), (c * c % 7) * 255 / 7, (c * c * c % 31) * 255 / 31);
    //std::cout << vd.idx() << " " << c << " " << color << std::endl;
    put(vcm, vd, color);
  }

  nb_empty = 0;
  for (int i = 0; i < nb_clusters; i++)
  {
    if (cluster_frequency[i] == 0)
    {
      nb_empty++;
    }
  }

  std::cout << "nb_empty: " << nb_empty << std::endl;

  std::cout << "kak1" << std::endl;
  std::string name = std::to_string(nb_clusters) + ".off";
  CGAL::IO::write_OFF(name, pmesh, CGAL::parameters::vertex_color_map(vcm));
  std::cout << "kak2" << std::endl;

  /// Construct new Mesh 
  std::vector<int> valid_cluster_map(nb_clusters, -1);
  std::vector<typename GT::Point_3> points;
  Point_set_3<typename GT::Point_3> point_set;

  std::vector<std::vector<int> > polygons;
  PolygonMesh simplified_mesh;

  for (int i = 0; i < nb_clusters; i++) //should i =1 ?
  {
    if (clusters[i].weight_sum > 0)
    {
      valid_cluster_map[i] = points.size();
      typename GT::Vector_3 center_v = clusters[i].compute_centroid();
      typename GT::Point_3 center_p(center_v.x(), center_v.y(), center_v.z());
      points.push_back(center_p);
      point_set.insert(center_p);
    }
  }
   
  // extract boundary cycles
  std::vector<Halfedge_descriptor> border_hedges;
  extract_boundary_cycles(pmesh, std::back_inserter(border_hedges));

  // loop over boundary loops
  for (Halfedge_descriptor hd : border_hedges)
  {
    Halfedge_descriptor hd1 = hd;

    int cb_first = -1;

    do
    {
      // 1- get the target and source vertices vt, vs
      // 2- if the target and source vertices are in different clusters, create a new vertex vb between them vb = (vt + vs) / 2
      // it is also added to the point_set
      // 3- make a new face with the new vertex vb and the centers of the clusters of vt and vs
      // 4- also make a new face with vb, the next vb, and the center of the cluster of vt

      Vertex_descriptor vt = target(hd1, pmesh);
      Vertex_descriptor vs = source(hd1, pmesh);

      int ct = get(vertex_cluster_pmap, vt);
      int cs = get(vertex_cluster_pmap, vs);

      if (ct != cs)
      {
        typename GT::Point_3 vt_p = get(vpm, vt);
        typename GT::Point_3 vs_p = get(vpm, vs);
        typename GT::Vector_3 vt_v(vt_p.x(), vt_p.y(), vt_p.z());
        typename GT::Vector_3 vs_v(vs_p.x(), vs_p.y(), vs_p.z());

        typename GT::Vector_3 vb_v = (vt_v + vs_v) / 2;
        typename GT::Point_3 vb_p(vb_v.x(), vb_v.y(), vb_v.z());

        points.push_back(vb_p);
        point_set.insert(vb_p);

        int cb = points.size() - 1;

        if (cb_first == -1)
          cb_first = cb;

        int ct_mapped = valid_cluster_map[ct], cs_mapped = valid_cluster_map[cs];

        if (ct_mapped != -1 && cs_mapped != -1)
        {
          std::vector<int>
          polygon = {ct_mapped, cb, cs_mapped};
          polygons.push_back(polygon);

          // after the loop, the last cb+1 should be modified to the first cb
          polygon = {cb, ct_mapped, cb + 1};
          polygons.push_back(polygon);
        }
      }
      hd1 = next(hd1, pmesh);
    } while (hd1 != hd);
    polygons[polygons.size() - 1][2] = cb_first;
  }

  name = std::to_string(nb_clusters) + "_points.off";
  CGAL::IO::write_point_set(name, point_set);

  for (Face_descriptor fd : faces(pmesh))
  {
    Halfedge_descriptor hd1 = halfedge(fd, pmesh);
    Vertex_descriptor v1 = source(hd1, pmesh);
    Halfedge_descriptor hd2 = next(hd1, pmesh);
    Vertex_descriptor v2 = source(hd2, pmesh);
    Halfedge_descriptor hd3 = next(hd2, pmesh);
    Vertex_descriptor v3 = source(hd3, pmesh);

    int c1 = get(vertex_cluster_pmap, v1);
    int c2 = get(vertex_cluster_pmap, v2);
    int c3 = get(vertex_cluster_pmap, v3);

    if (c1 != c2 && c1 != c3 && c2 != c3)
    {
      int c1_mapped = valid_cluster_map[c1], c2_mapped = valid_cluster_map[c2], c3_mapped = valid_cluster_map[c3];
      if (c1_mapped != -1 && c2_mapped != -1 && c3_mapped != -1)
      {
        std::vector<int> polygon = {c1_mapped, c2_mapped, c3_mapped};
        polygons.push_back(polygon);
      }
    }
  }

  std::cout << "are polygons a valid mesh ? : " << is_polygon_soup_a_polygon_mesh(polygons) << std::endl;
  polygon_soup_to_polygon_mesh(points, polygons, simplified_mesh);

  name = std::to_string(nb_clusters) + "_simped.off";
  CGAL::IO::write_OFF(name, simplified_mesh);
  std::cout << "kak3" << std::endl;

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
