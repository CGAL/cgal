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

/// TODO:
// #ifndef CGAL_<>
// #define CGAL_<>
// #include <CGAL/license/<>>

#include <CGAL/assertions.h>
#include <CGAL/IO/Color.h>
#include <CGAL/Named_function_parameters.h>
#include <CGAL/property_map.h>
#include <CGAL/boost/graph/named_params_helper.h>

#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/border.h>

#include <CGAL/Polygon_mesh_processing/interpolated_corrected_curvatures.h>
#include <CGAL/Polygon_mesh_processing/acvd/qem_metrics.h>

#include <CGAL/subdivision_method_3.h>

#include <CGAL/Point_set_3/IO.h>
#include <CGAL/Point_set_3.h>

#include <Eigen/Eigenvalues>

#include <numeric>
#include <vector>
#include <queue>
#include <unordered_set>
#include <iostream>

#define CGAL_CLUSTERS_TO_VERTICES_THRESHOLD 0.1
#define CGAL_WEIGHT_CLAMP_RATIO_THRESHOLD 10000

namespace CGAL {

namespace Polygon_mesh_processing {

namespace internal {

template <typename GT>
struct IsotropicClusterData {
  typename GT::Vector_3 site_sum;
  typename GT::FT weight_sum;
  typename GT::FT energy;
  size_t nb_vertices;

  IsotropicClusterData() : site_sum(0, 0, 0), weight_sum(0), energy(0), nb_vertices(0) {}

  void add_vertex(const typename GT::Vector_3 vertex_position, const typename GT::FT weight)
  {
    this->site_sum += vertex_position * weight;
    this->weight_sum += weight;
    this->nb_vertices++;
  }

  void remove_vertex(const typename GT::Vector_3 vertex_position, const typename GT::FT weight)
  {
    this->site_sum -= vertex_position * weight;
    this->weight_sum -= weight;
    this->nb_vertices--;
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

template <typename GT>
struct QEMClusterData {
  typename GT::Vector_3 site_sum;
  typename GT::Vector_3 q_centroid;
  typename GT::FT weight_sum;
  typename GT::FT energy;
  typename GT::FT qem[9];
  char rank_deficiency;

  QEMClusterData() : site_sum(0, 0, 0), weight_sum(0), energy(0), rank_deficiency(0), q_centroid(0, 0, 0) {
    for (int i = 0; i < 9; i++)
      qem[i] = 0;
  }

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

template <class PolygonMesh, class NamedParameters = parameters::Default_named_parameters>
void upsample_subdivision_property(PolygonMesh& pmesh, const NamedParameters& np = parameters::default_values()) {
  typedef typename GetGeomTraits<PolygonMesh, NamedParameters>::type GT;
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor Vertex_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor Halfedge_descriptor;
  typedef Constant_property_map<Vertex_descriptor, Principal_curvatures_and_directions<GT>> Default_principal_map;
  typedef typename internal_np::Lookup_named_param_def<internal_np::vertex_principal_curvatures_and_directions_map_t,
    NamedParameters,
    Default_principal_map>::type VPCDM;

  using parameters::choose_parameter;
  using parameters::get_parameter;
  using parameters::is_default_parameter;

  typedef typename CGAL::GetVertexPointMap<PolygonMesh, NamedParameters>::type VPM;
  VPM vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                         get_property_map(CGAL::vertex_point, pmesh));

  // get curvature related parameters
  const typename VPCDM vpcd_map =
    choose_parameter(get_parameter(np, internal_np::vertex_principal_curvatures_and_directions_map),
      Default_principal_map());

  // unordered_set of old vertices
  std::unordered_set<Vertex_descriptor> old_vertices;

  bool curvatures_available = !is_default_parameter<NamedParameters, internal_np::vertex_principal_curvatures_and_directions_map_t>::value;

  unsigned int step = choose_parameter(get_parameter(np, internal_np::number_of_iterations), 1);
  Upsample_mask_3<PolygonMesh,VPM> mask(&pmesh, vpm);

  for (unsigned int i = 0; i < step; i++){
    for (Vertex_descriptor vd : vertices(pmesh))
      old_vertices.insert(vd);

    Subdivision_method_3::internal::PTQ_1step(pmesh, vpm, mask);
    // interpolate curvature values
    if (curvatures_available)
    {
      for (Vertex_descriptor vd : vertices(pmesh))
      {
        if (old_vertices.find(vd) == old_vertices.end())
        {
          Principal_curvatures_and_directions<GT> pcd;
          pcd.min_curvature = 0;
          pcd.max_curvature = 0;
          pcd.min_direction = GT::Vector_3(0, 0, 0);
          pcd.max_direction = GT::Vector_3(0, 0, 0);
          for (Halfedge_descriptor hd : halfedges_around_target(vd, pmesh))
          {
            Vertex_descriptor v1 = source(hd, pmesh);
            if (old_vertices.find(v1) != old_vertices.end())
            {
              Principal_curvatures_and_directions<GT> pcd1 = get(vpcd_map, v1);
              pcd.min_curvature += pcd1.min_curvature;
              pcd.max_curvature += pcd1.max_curvature;
              pcd.min_direction += pcd1.min_direction;
              pcd.max_direction += pcd1.max_direction;
            }
          }
          pcd.min_curvature = pcd.min_curvature / 2;
          pcd.max_curvature = pcd.max_curvature / 2;
          pcd.min_direction = pcd.min_direction / sqrt(pcd.min_direction.squared_length());
          pcd.max_direction = pcd.max_direction / sqrt(pcd.max_direction.squared_length());
          put(vpcd_map, vd, pcd);
        }
      }
    }
  }
}

// provide a property map for principal curvatures as a named parameter for adaptive clustering
// provide a gradation factor as a named parameter for adaptive clustering
template <typename PolygonMesh,
          typename NamedParameters = parameters::Default_named_parameters>
std::pair<
  std::vector<typename GetGeomTraits<PolygonMesh, NamedParameters>::type::Point_3>,
  std::vector<std::vector<int>>
> acvd_isotropic(
    PolygonMesh& pmesh,
    const int nb_clusters,
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
  typedef typename boost::property_map<PolygonMesh, CGAL::dynamic_vertex_property_t<bool> >::type VertexVisitedMap;
  typedef typename boost::property_map<PolygonMesh, CGAL::dynamic_vertex_property_t<typename GT::FT> >::type VertexWeightMap;
  typedef Constant_property_map<Vertex_descriptor, Principal_curvatures_and_directions<GT>> Default_principal_map;
  typedef typename internal_np::Lookup_named_param_def<internal_np::vertex_principal_curvatures_and_directions_map_t,
    NamedParameters,
    Default_principal_map>::type Vertex_principal_curvatures_and_directions_map;

  using parameters::choose_parameter;
  using parameters::get_parameter;
  using parameters::is_default_parameter;

  Vertex_position_map vpm = choose_parameter(get_parameter(np, CGAL::vertex_point),
    get_property_map(CGAL::vertex_point, pmesh));

  // get curvature related parameters
  const typename GT::FT gradation_factor = choose_parameter(get_parameter(np, internal_np::gradation_factor), 0);
  const typename Vertex_principal_curvatures_and_directions_map vpcd_map =
    choose_parameter(get_parameter(np, internal_np::vertex_principal_curvatures_and_directions_map),
      Default_principal_map());

  if (gradation_factor > 0 &&
    is_default_parameter<NamedParameters, internal_np::vertex_principal_curvatures_and_directions_map_t>::value)
    interpolated_corrected_principal_curvatures_and_directions(pmesh, vpcd_map);


  // TODO: handle cases where the mesh is not a triangle mesh
  CGAL_precondition(CGAL::is_triangle_mesh(pmesh));

  // TODO: copy the mesh in order to not modify the original mesh
  //PolygonMesh pmesh = pmesh_org;
  int nb_vertices = num_vertices(pmesh);

  // To provide the functionality remeshing (not just simplification), we might need to
  // subdivide the mesh before clustering
  // in either case, nb_clusters <= nb_vertices * CGAL_CLUSTERS_TO_VERTICES_THRESHOLD

  // do the following while nb_clusters > nb_vertices * CGAL_CLUSTERS_TO_VERTICES_THRESHOLD
  // That is, because the subdivision steps heuristic is not 100% guaranteed to produce
  // the desired number of vertices.
  while (nb_clusters > nb_vertices * CGAL_CLUSTERS_TO_VERTICES_THRESHOLD)
  {
    double curr_factor = nb_clusters / (nb_vertices * CGAL_CLUSTERS_TO_VERTICES_THRESHOLD);
    int subdivide_steps = max((int)ceil(log(curr_factor) / log(4)), 0);

    //std::cout << "subdivide_steps: " << subdivide_steps << std::endl;

    if (subdivide_steps > 0)
    {
      if (gradation_factor == 0) // no adaptive clustering
        Subdivision_method_3::Upsample_subdivision(
          pmesh,
          CGAL::parameters::number_of_iterations(subdivide_steps)
        );
      else // adaptive clustering
        upsample_subdivision_property(
          pmesh,
          CGAL::parameters::number_of_iterations(subdivide_steps).vertex_principal_curvatures_and_directions_map(vpcd_map)
        );
      vpm = get_property_map(CGAL::vertex_point, pmesh);
      nb_vertices = num_vertices(pmesh);
    }
  }

  // initial random clusters
  // property map from vertex_descriptor to cluster index
  VertexClusterMap vertex_cluster_pmap = get(CGAL::dynamic_vertex_property_t<int>(), pmesh);
  VertexVisitedMap vertex_visited_pmap = get(CGAL::dynamic_vertex_property_t<bool>(), pmesh);
  VertexWeightMap vertex_weight_pmap = get(CGAL::dynamic_vertex_property_t<typename GT::FT>(), pmesh);
  std::vector<IsotropicClusterData<GT>> clusters(nb_clusters);
  std::queue<Halfedge_descriptor> clusters_edges_active;
  std::queue<Halfedge_descriptor> clusters_edges_new;

  // initialize vertex weights and clusters
  for (Vertex_descriptor vd : vertices(pmesh))
  {
    put(vertex_weight_pmap, vd, 0);
    put(vertex_visited_pmap, vd, false);
    put(vertex_cluster_pmap, vd, -1);
  }

  // compute vertex weights (dual area)
  typename GT::FT weight_avg = 0;
  for (Face_descriptor fd : faces(pmesh))
  {
    typename GT::FT weight = abs(CGAL::Polygon_mesh_processing::face_area(fd, pmesh)) / 3;

    for (Vertex_descriptor vd : vertices_around_face(halfedge(fd, pmesh), pmesh))
    {
      typename GT::FT vertex_weight = get(vertex_weight_pmap, vd);
      if (gradation_factor == 0) // no adaptive clustering
        vertex_weight += weight;
      else // adaptive clustering
      {
        typename GT::FT k1 = get(vpcd_map, vd).min_curvature;
        typename GT::FT k2 = get(vpcd_map, vd).max_curvature;
        typename GT::FT k_sq = (k1 * k1 + k2 * k2);
        vertex_weight += weight * (/*`eps * avg_curvature` instead of 1 +*/ pow(k_sq, gradation_factor / 2.0));  // /2.0 because k_sq is squared
      }
      weight_avg += vertex_weight;
      put(vertex_weight_pmap, vd, vertex_weight);
    }
  }
  weight_avg /= nb_vertices;

  // clamp the weights up and below by a ratio (like 10,000) * avg_weights
  for (Vertex_descriptor vd : vertices(pmesh))
  {
    typename GT::FT vertex_weight = get(vertex_weight_pmap, vd);
    if (vertex_weight > CGAL_WEIGHT_CLAMP_RATIO_THRESHOLD * weight_avg)
      put(vertex_weight_pmap, vd, CGAL_WEIGHT_CLAMP_RATIO_THRESHOLD * weight_avg);
    else if (vertex_weight < 1.0 / CGAL_WEIGHT_CLAMP_RATIO_THRESHOLD * weight_avg)
      put(vertex_weight_pmap, vd, 1.0 / CGAL_WEIGHT_CLAMP_RATIO_THRESHOLD * weight_avg);
  }

  for (int ci = 0; ci < nb_clusters; ci++)
  {
    int vi;
    Vertex_descriptor vd;
    do {
      vi = CGAL::get_default_random().get_int(0, num_vertices(pmesh));
      vd = *(vertices(pmesh).begin() + vi);
    } while (get(vertex_cluster_pmap, vd) != -1);

    put(vertex_cluster_pmap, vd, ci);
    typename GT::Point_3 vp = get(vpm, vd);
    typename GT::Vector_3 vpv(vp.x(), vp.y(), vp.z());
    clusters[ci].add_vertex(vpv, get(vertex_weight_pmap, vd));

    for (Halfedge_descriptor hd : halfedges_around_source(vd, pmesh))
      clusters_edges_active.push(hd);
  }

  int nb_modifications = 0;
  int nb_disconnected = 0;

  do
  {
    nb_disconnected = 0;
    do
    {
      nb_modifications = 0;

      while (!clusters_edges_active.empty()) {
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

          // reset to no change
          clusters[c1].add_vertex(vpv1, v1_weight);
          clusters[c2].remove_vertex(vpv1, v1_weight);

          // The effect of the following should always be reversed after the comparison
          clusters[c2].remove_vertex(vpv2, v2_weight);
          clusters[c1].add_vertex(vpv2, v2_weight);

          typename GT::FT e_v2_to_c1 = clusters[c1].compute_energy() + clusters[c2].compute_energy();

          if (e_v2_to_c1 < e_no_change && e_v2_to_c1 < e_v1_to_c2 && clusters[c2].nb_vertices > 0) // > 0 as 1 vertex was removed from c2
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
          else if (e_v1_to_c2 < e_no_change && clusters[c1].nb_vertices > 2) // > 2 as 1 vertex was added to c1
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
      // std::cout << "# Modifications: " << nb_modifications << "\n";

      clusters_edges_active.swap(clusters_edges_new);
    } while (nb_modifications > 0);

    // clean clusters here
    // the goal is to delete clusters with multiple connected components
    // for each cluster, do a BFS from a vertex in the cluster
    // we need to keep the largest connected component for each cluster
    // and set the other connected components to -1 (empty cluster), would also need to update clusters_edges_new

    std::vector<std::vector<std::vector<Vertex_descriptor>>> cluster_components(nb_clusters, std::vector<std::vector<Vertex_descriptor>>());

    std::queue<Vertex_descriptor> vertex_queue;

    // loop over vertices
    for (Vertex_descriptor vd : vertices(pmesh))
    {
      if (get(vertex_visited_pmap, vd)) continue;
      int c = get(vertex_cluster_pmap, vd);
      if (c != -1)
      {
        // first component of this cluster
        cluster_components[c].push_back(std::vector<Vertex_descriptor>());
        int component_i = cluster_components[c].size() - 1;

        // visited_clusters[c] = true;
        vertex_queue.push(vd);
        put(vertex_visited_pmap, vd, true);
        while (!vertex_queue.empty())
        {
          Vertex_descriptor v = vertex_queue.front();
          vertex_queue.pop();
          cluster_components[c][component_i].push_back(v);

          for (Halfedge_descriptor hd : halfedges_around_source(v, pmesh))
          {
            Vertex_descriptor v2 = target(hd, pmesh);
            int c2 = get(vertex_cluster_pmap, v2);
            if (c2 == c && !get(vertex_visited_pmap, v2))
            {
              vertex_queue.push(v2);
              put(vertex_visited_pmap, v2, true);
            }
          }
        }
      }
    }

    // loop over clusters
    for (int c = 0; c < nb_clusters; c++)
    {
      if (cluster_components[c].size() <= 1) continue; // only one component, no need to do anything
      nb_disconnected++;
      int max_component_size = 0;
      int max_component_index = -1;
      for (int component_i = 0; component_i < cluster_components[c].size(); component_i++)
      {
        if (cluster_components[c][component_i].size() > max_component_size)
        {
          max_component_size = cluster_components[c][component_i].size();
          max_component_index = component_i;
        }
      }
      // set cluster to -1 for all components except the largest one
      for (int component_i = 0; component_i < cluster_components[c].size(); component_i++)
      {
        if (component_i != max_component_index)
        {
          for (Vertex_descriptor vd : cluster_components[c][component_i])
          {
            put(vertex_cluster_pmap, vd, -1);
            // remove vd from cluster c
            typename GT::Point_3 vp = get(vpm, vd);
            typename GT::Vector_3 vpv(vp.x(), vp.y(), vp.z());
            clusters[c].remove_vertex(vpv, get(vertex_weight_pmap, vd));
            // add all halfedges around v except hi to the queue
            for (Halfedge_descriptor hd : halfedges_around_source(vd, pmesh))
            {
              // add hd to the queue if its target is not in the same cluster
              Vertex_descriptor v2 = target(hd, pmesh);
              int c2 = get(vertex_cluster_pmap, v2);
              if (c2 != c)
                clusters_edges_new.push(hd);
            }
          }
        }
      }
    }

  } while (nb_disconnected > 0);


  VertexColorMap vcm = get(CGAL::dynamic_vertex_property_t<CGAL::IO::Color>(), pmesh);

  // // frequency of each cluster
  // cluster_frequency = std::vector<int>(nb_clusters, 0);

  // for (Vertex_descriptor vd : vertices(pmesh))
  // {
  //   int c = get(vertex_cluster_pmap, vd);
  //   cluster_frequency[c]++;
  //   CGAL::IO::Color color(255 - (c * 255 / nb_clusters), (c * c % 7) * 255 / 7, (c * c * c % 31) * 255 / 31);
  //   put(vcm, vd, color);
  // }

  // std::string name = std::to_string(nb_clusters) + ".off";
  // CGAL::IO::write_OFF(name, pmesh, CGAL::parameters::vertex_color_map(vcm));

  /// Construct new Mesh
  std::vector<int> valid_cluster_map(nb_clusters, -1);
  std::vector<typename GT::Point_3> points;
  Point_set_3<typename GT::Point_3> point_set;

  std::vector<std::vector<int>> polygons;
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

  orient_polygon_soup(points, polygons);

  return std::make_pair(points, polygons);
}


} // namespace internal

/**
* \ingroup PMP_acvd_grp
*
* Performs uniform (isotropic) centroidal voronoi diagram simplification on a polygon mesh.
* This can also be used for remeshing by setting the number of clusters to the desired number of vertices.
* The number of clusters is the number of vertices in the output mesh.
*
* @tparam PolygonMesh a model of `FaceListGraph`
* @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters".
*
* @param pmesh input polygon mesh
* @param nb_vertices number of target vertices in the output mesh
* @param np optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
*        `GT` stands for the type of the object provided to the named parameter `geom_traits()`.
*
* \cgalNamedParamsBegin
*
*   \cgalParamNBegin{vertex_principal_curvatures_and_directions_map}
*     \cgalParamDescription{a property map associating principal curvatures and directions to the vertices of `pmesh`, used for adaptive clustering.}
*     \cgalParamType{a class model of `ReadWritePropertyMap` with
*                    `boost::graph_traits<PolygonMesh>::%vertex_descriptor`
*                    as key type and `Principal_curvatures_and_directions<GT>` as value type.}
*     \cgalParamExtra{If this parameter is omitted, but `gradation_factor` is not (and is > 0), an internal property map
*                     will be created and curvature values will be computed.}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{gradation_factor}
*     \cgalParamDescription{a factor used to gradate the weights of the vertices based on their curvature values.}
*     \cgalParamType{`GT::FT`}
*     \cgalParamDefault{0}
*     \cgalParamExtra{If this parameter is omitted, no adaptive clustering will be performed.}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{vertex_point_map}
*       \cgalParamDescription{a property map associating points to the vertices of `pmesh`.}
*       \cgalParamType{a class model of `ReadablePropertyMap` with
*                      `boost::graph_traits<PolygonMesh>::%vertex_descriptor`
*                      as key type and `GT::Point_3` as value type.}
*       \cgalParamDefault{`boost::get(CGAL::vertex_point, pmesh)`.}
*       \cgalParamExtra{If this parameter is omitted, an internal property map for
*                       `CGAL::vertex_point_t` must be available in `PolygonMesh`.}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{geom_traits}
*      \cgalParamDescription{an instance of a geometric traits class.}
*      \cgalParamType{a class model of `Kernel`}
*      \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`.}
*      \cgalParamExtra{The geometric traits class must be compatible with the vertex point type.}
*   \cgalParamNEnd
*
* \cgalNamedParamsEnd
*
* @pre only triangle meshes are supported for now
* @return a pair of vectors of points and polygons representing the simplified mesh as a polygon soup
*/

template <typename PolygonMesh,
  typename NamedParameters = parameters::Default_named_parameters>
std::pair<
  std::vector<typename GetGeomTraits<PolygonMesh, NamedParameters>::type::Point_3>,
  std::vector<std::vector<int>>
>  acvd_isotropic_simplification_polygon_soup(
    PolygonMesh& pmesh,
    const int& nb_vertices,
    const NamedParameters& np = parameters::default_values()
  )
{
  return internal::acvd_isotropic<PolygonMesh, NamedParameters>(
    pmesh,
    nb_vertices,
    np
  );
}

/**
* \ingroup PMP_acvd_grp
*
* Performs uniform (isotropic) centroidal voronoi diagram simplification on a polygon mesh.
* This can also be used for remeshing by setting the number of clusters to the desired number of vertices.
* The number of clusters is the number of vertices in the output mesh.
*
* @tparam PolygonMesh a model of `FaceListGraph`
* @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters".
*
* @param pmesh input polygon mesh
* @param nb_vertices number of target vertices in the output mesh
* @param np optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
*        `GT` stands for the type of the object provided to the named parameter `geom_traits()`.
*
* \cgalNamedParamsBegin
*
*   \cgalParamNBegin{vertex_principal_curvatures_and_directions_map}
*     \cgalParamDescription{a property map associating principal curvatures and directions to the vertices of `pmesh`, used for adaptive clustering.}
*     \cgalParamType{a class model of `ReadWritePropertyMap` with
*                    `boost::graph_traits<PolygonMesh>::%vertex_descriptor`
*                    as key type and `Principal_curvatures_and_directions<GT>` as value type.}
*     \cgalParamExtra{If this parameter is omitted, but `gradation_factor` is not (and is > 0), an internal property map
*                     will be created and curvature values will be computed.}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{gradation_factor}
*     \cgalParamDescription{a factor used to gradate the weights of the vertices based on their curvature values.}
*     \cgalParamType{`GT::FT`}
*     \cgalParamDefault{0}
*     \cgalParamExtra{If this parameter is omitted, no adaptive clustering will be performed.}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{vertex_point_map}
*       \cgalParamDescription{a property map associating points to the vertices of `pmesh`.}
*       \cgalParamType{a class model of `ReadablePropertyMap` with
*                      `boost::graph_traits<PolygonMesh>::%vertex_descriptor`
*                      as key type and `GT::Point_3` as value type.}
*       \cgalParamDefault{`boost::get(CGAL::vertex_point, pmesh)`.}
*       \cgalParamExtra{If this parameter is omitted, an internal property map for
*                       `CGAL::vertex_point_t` must be available in `PolygonMesh`.}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{geom_traits}
*      \cgalParamDescription{an instance of a geometric traits class.}
*      \cgalParamType{a class model of `Kernel`}
*      \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`.}
*      \cgalParamExtra{The geometric traits class must be compatible with the vertex point type.}
*   \cgalParamNEnd
*
* \cgalNamedParamsEnd
*
* @pre only triangle meshes are supported for now
* @return the simplified mesh as a PolygonMesh
*/

template <typename PolygonMesh,
  typename NamedParameters = parameters::Default_named_parameters>
PolygonMesh acvd_isotropic_simplification(
    PolygonMesh& pmesh,
    const int& nb_vertices,
    const NamedParameters& np = parameters::default_values()
  )
{
  auto ps = acvd_isotropic_simplification_polygon_soup(
    pmesh,
    nb_vertices,
    np
  );

  PolygonMesh simplified_mesh;
  polygon_soup_to_polygon_mesh(ps.first, ps.second, simplified_mesh);
  return simplified_mesh;
}

// template <typename PolygonMesh,
//   typename NamedParameters = parameters::Default_named_parameters>
// PolygonMesh acvd_qem_simplification(
//     PolygonMesh& pmesh,
//     const int& nb_vertices,
//     const NamedParameters& np = parameters::default_values()
//   )
// {
//   return internal::acve_qem<PolygonMesh, /*IsotropicMetric,*/ NamedParameters>(
//     pmesh,
//     nb_vertices,
//     np
//     );
// }

} // namespace Polygon_mesh_processing

} // namespace CGAL
