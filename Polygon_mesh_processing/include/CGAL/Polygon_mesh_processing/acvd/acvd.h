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
// Author(s)     : Hossam Saeed, Sébastien Valette, and Sébastien Loriot
//

#ifndef CGAL_PMP_ACVD_REMESHING_H
#define CGAL_PMP_ACVD_REMESHING_H

#include <CGAL/license/Polygon_mesh_processing/acvd.h>

#include <CGAL/assertions.h>
#include <CGAL/Named_function_parameters.h>
#include <CGAL/property_map.h>
#include <CGAL/boost/graph/named_params_helper.h>

#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/border.h>

#include <CGAL/Polygon_mesh_processing/interpolated_corrected_curvatures.h>

#include <CGAL/subdivision_method_3.h>
#include <CGAL/Random.h>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include <numeric>
#include <vector>
#include <queue>
#include <unordered_set>
#include <iostream>


#define CGAL_TO_QEM_MODIFICATION_THRESHOLD 1e-3
#define CGAL_WEIGHT_CLAMP_RATIO_THRESHOLD 10000

namespace CGAL {

namespace Polygon_mesh_processing {

namespace internal {

template <typename GT>
void compute_qem_face(const typename GT::Vector_3& p1, const typename GT::Vector_3& p2, const typename GT::Vector_3& p3, Eigen::Matrix<typename GT::FT, 4, 4>& quadric)
{
  typename GT::Vector_3 crossX1X2 = CGAL::cross_product(p1, p2);
  typename GT::Vector_3 crossX2X3 = CGAL::cross_product(p2, p3);
  typename GT::Vector_3 crossX3X1 = CGAL::cross_product(p3, p1);
  typename GT::FT determinantABC = CGAL::determinant(p1, p2, p3);

  typename GT::FT n[4] = {
    crossX1X2.x() + crossX2X3.x() + crossX3X1.x(),
    crossX1X2.y() + crossX2X3.y() + crossX3X1.y(),
    crossX1X2.z() + crossX2X3.z() + crossX3X1.z(),
    -determinantABC
  };

  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 4; ++j)
      quadric(i, j) = n[i] * n[j];
}

template <typename GT>
void compute_qem_vertex(std::vector<std::vector<typename GT::Vector_3>> cluster_tris, Eigen::Matrix<typename GT::FT, 4, 4>& quadric)
{
  quadric.setZero();

  for (int i = 0; i < cluster_tris.size(); ++i)
  {
    Eigen::Matrix<typename GT::FT, 4, 4> q;
    compute_qem_face<GT>(cluster_tris[i][0], cluster_tris[i][1], cluster_tris[i][2], q);
    quadric = quadric + q;
  }
}

template <typename GT>
typename GT::Vector_3 compute_displacement(const Eigen::Matrix<typename GT::FT, 4, 4> quadric, const typename GT::Point_3& p, int& rank_deficiency)
{
  typedef Eigen::Matrix<typename GT::FT, 3, 3> Matrix3d;

  int MaxNumberOfUsedSingularValues = 3;
  Matrix3d A;
  A(0, 0) = quadric(0, 0);
  A(0, 1) = A(1, 0) = quadric(0, 1);
  A(0, 2) = A(2, 0) = quadric(0, 2);
  A(1, 1) = quadric(1, 1);
  A(1, 2) = A(2, 1) = quadric(1, 2);
  A(2, 2) = quadric(2, 2);

  Matrix3d U, V;
  typename GT::Vector_3 w;
  Eigen::JacobiSVD<Matrix3d> svd(A, Eigen::ComputeFullU | Eigen::ComputeFullV);
  U = svd.matrixU();
  V = svd.matrixV();
  auto w_temp = svd.singularValues();
  w = typename GT::Vector_3(w_temp(0), w_temp(1), w_temp(2));

  // compute all eigen values absolute values
  typename GT::FT AbsolutesEigenValues[3];
  typename GT::FT maxW = -1.0;
  for (int j = 0; j < 3; ++j)
  {
    typename GT::FT AbsoluteEigenValue = fabs(w[j]);
    AbsolutesEigenValues[j] = AbsoluteEigenValue;
    if (AbsoluteEigenValue > maxW)
      maxW = AbsoluteEigenValue;
  }
  typename GT::FT invmaxW = 1.0 / maxW;

  Matrix3d tempMatrix, tempMatrix2;

  for (int i = 0; i < 3; ++i)
  {
    typename GT::FT LocalMaxW = -1;
    int IndexMax = -1;

    // find the remaining eigenvalue with highest absolute value
    for (int j = 0; j < 3; ++j)
    {
      if (LocalMaxW < AbsolutesEigenValues[j])
      {
        LocalMaxW = AbsolutesEigenValues[j];
        IndexMax = j;
      }
    }

    if ((AbsolutesEigenValues[IndexMax] * invmaxW > 1e-3)
      && (MaxNumberOfUsedSingularValues > 0))
    {
      // If this is true, then w[i] != 0, so this division is ok.
      double Inv = 1.0 / w[IndexMax];
      tempMatrix(IndexMax, 0) = U(0, IndexMax) * Inv;
      tempMatrix(IndexMax, 1) = U(1, IndexMax) * Inv;
      tempMatrix(IndexMax, 2) = U(2, IndexMax) * Inv;
    }
    else
    {
      tempMatrix(IndexMax, 0) = 0;
      tempMatrix(IndexMax, 1) = 0;
      tempMatrix(IndexMax, 2) = 0;
      ++rank_deficiency;
    }

    // set the eigenvalu to -2 to remove it from subsequent tests
    AbsolutesEigenValues[IndexMax] = -2;
    MaxNumberOfUsedSingularValues--;
  }

  tempMatrix2 = V * tempMatrix;
  Eigen::Matrix<typename GT::FT, 3, 1> p_temp, b;
  b << -quadric(0, 3), -quadric(1, 3), -quadric(2, 3);
  p_temp << p.x(), p.y(), p.z();
  Eigen::Matrix<typename GT::FT, 3, 1> displacement_temp = tempMatrix2 * (b - A * p_temp);
  typename GT::Vector_3 displacement = { displacement_temp(0), displacement_temp(1), displacement_temp(2) };
  return displacement;
}

template <typename GT>
void compute_representative_point(const Eigen::Matrix<typename GT::FT, 4, 4>& quadric, typename GT::Point_3& p, int& rank_deficiency)
{
  // average point
  typename GT::Vector_3 displacement = compute_displacement<GT>(quadric, p, rank_deficiency);
  p = p + displacement;
}


template <typename GT>
struct QEMClusterData {
  typename GT::Vector_3 site_sum;
  typename GT::FT weight_sum;
  Eigen::Matrix<typename GT::FT, 4, 4> quadric_sum;
  typename GT::Vector_3 representative_point_;
  typename GT::FT energy;
  bool modified = true;
  int last_modification_iteration = 0;
  size_t nb_vertices = 0;

  QEMClusterData() :
    site_sum(0, 0, 0),
    weight_sum(0),
    representative_point_(0, 0, 0),
    energy(0)
  {
    quadric_sum.setZero();
  }

  void add_vertex(const typename GT::Vector_3 vertex_position, const typename GT::FT weight, const Eigen::Matrix<typename GT::FT, 4, 4>& quadric)
  {
    this->site_sum += vertex_position * weight;
    this->weight_sum += weight;
    this->quadric_sum += quadric;
    ++this->nb_vertices;
    this->modified=true;
  }

  void add_vertex(const typename GT::Point_3 vertex_position, const typename GT::FT weight, const Eigen::Matrix<typename GT::FT, 4, 4>& quadric)
  {
    add_vertex( vertex_position - ORIGIN, weight, quadric );
  }

  void remove_vertex(const typename GT::Vector_3 vertex_position, const typename GT::FT weight, const Eigen::Matrix<typename GT::FT, 4, 4>& quadric)
  {
    this->site_sum -= vertex_position * weight;
    this->weight_sum -= weight;
    this->quadric_sum -= quadric;
    this->nb_vertices--;
    this->modified=true;
  }

  void compute_representative(bool qem)
  {
    if (this->weight_sum > 0)
    {
      if (qem)
      {
        int rank_deficiency = 0;
        typename GT::Point_3 point =  CGAL::ORIGIN + (this->site_sum / this->weight_sum);
        compute_representative_point<GT>(this->quadric_sum, point, rank_deficiency);
        this->representative_point_ = { point.x(), point.y(), point.z() };
      }
      else
        this->representative_point_ = (this->site_sum) / this->weight_sum;
    }
  }

  typename GT::FT compute_energy(bool qem)
  {
    if (modified)
    {
      compute_representative(qem);
      auto dot_product = GT().compute_scalar_product_3_object();

      this->energy = (this->representative_point_).squared_length() * this->weight_sum
                     - 2 * dot_product(this->representative_point_, this->site_sum);
    }
    return this->energy;
  }

  typename GT::Point_3 representative_point(bool qem)
  {
    if (modified) compute_representative(qem);
    return ORIGIN+representative_point_;
  }
};

template <class TriangleMesh, class NamedParameters = parameters::Default_named_parameters>
void upsample_subdivision_property(TriangleMesh& pmesh, const NamedParameters& np = parameters::default_values()) {
  typedef typename GetGeomTraits<TriangleMesh, NamedParameters>::type GT;
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor Vertex_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor Halfedge_descriptor;
  typedef Constant_property_map<Vertex_descriptor, Principal_curvatures_and_directions<GT>> Default_principal_map;
  typedef typename internal_np::Lookup_named_param_def<internal_np::vertex_principal_curvatures_and_directions_map_t,
    NamedParameters,
    Default_principal_map>::type VPCDM;

  using parameters::choose_parameter;
  using parameters::get_parameter;
  using parameters::is_default_parameter;

  typedef typename CGAL::GetVertexPointMap<TriangleMesh, NamedParameters>::type VPM;
  VPM vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                         get_property_map(CGAL::vertex_point, pmesh));

  // get curvature related parameters
  const VPCDM vpcd_map =
    choose_parameter(get_parameter(np, internal_np::vertex_principal_curvatures_and_directions_map),
      Default_principal_map());

  // unordered_set of old vertices
  std::unordered_set<Vertex_descriptor> old_vertices;

  bool curvatures_available = !is_default_parameter<NamedParameters, internal_np::vertex_principal_curvatures_and_directions_map_t>::value;

  unsigned int step = choose_parameter(get_parameter(np, internal_np::number_of_iterations), 1);
  Upsample_mask_3<TriangleMesh,VPM> mask(&pmesh, vpm);

  for (unsigned int i = 0; i < step; ++i){
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
          pcd.min_direction = typename GT::Vector_3(0, 0, 0);
          pcd.max_direction = typename GT::Vector_3(0, 0, 0);
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

template <typename TriangleMesh,
          typename NamedParameters = parameters::Default_named_parameters>
std::pair<
  std::vector<typename GetGeomTraits<TriangleMesh, NamedParameters>::type::Point_3>,
  std::vector<std::vector<int>>
> acvd_impl(
    TriangleMesh& pmesh,
    const int nb_clusters,
    const NamedParameters& np = parameters::default_values()
  )
{
  typedef typename GetGeomTraits<TriangleMesh, NamedParameters>::type GT;
  typedef typename Eigen::Matrix<typename GT::FT, 4, 4> Matrix4x4;
  typedef typename GetVertexPointMap<TriangleMesh, NamedParameters>::const_type Vertex_position_map;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor Halfedge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor Vertex_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor Face_descriptor;
  typedef typename boost::property_map<TriangleMesh, CGAL::dynamic_vertex_property_t<int> >::type VertexClusterMap;
  typedef typename boost::property_map<TriangleMesh, CGAL::dynamic_vertex_property_t<bool> >::type VertexVisitedMap;
  typedef typename boost::property_map<TriangleMesh, CGAL::dynamic_vertex_property_t<typename GT::FT> >::type VertexWeightMap;
  typedef typename boost::property_map<TriangleMesh, CGAL::dynamic_vertex_property_t<Matrix4x4> >::type VertexQuadricMap;
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
  const Vertex_principal_curvatures_and_directions_map vpcd_map =
    choose_parameter(get_parameter(np, internal_np::vertex_principal_curvatures_and_directions_map),
      Default_principal_map());

  // using QEM ?
  bool postprocessing_qem = choose_parameter(get_parameter(np, internal_np::postprocessing_qem), false);
  const bool use_qem_metric = choose_parameter(get_parameter(np, internal_np::use_qem_metric), false);
  if (use_qem_metric) postprocessing_qem = false;

  // if adaptive clustering
  if (gradation_factor > 0 &&
    is_default_parameter<NamedParameters, internal_np::vertex_principal_curvatures_and_directions_map_t>::value)
    interpolated_corrected_curvatures(pmesh, parameters::vertex_principal_curvatures_and_directions_map(vpcd_map));

  CGAL_precondition(CGAL::is_triangle_mesh(pmesh));

  const double clusters_to_vertices_threshold_ratio = choose_parameter(get_parameter(np, internal_np::clusters_to_vertices_threshold_ratio), 0.1);
  // TODO also call split_long_edges: max edge length shall not be larger than 3 * input mean edge length
  // TODO cheese crashes
  // TODO compute less qem things if not used
  // TODO use symmetric matrices?

  // TODO: copy the mesh in order to not modify the original mesh
  //TriangleMesh pmesh = pmesh_org;
  int nb_vertices = num_vertices(pmesh);


  // For remeshing, we might need to subdivide the mesh before clustering
  // This should always hold: nb_clusters <= nb_vertices * clusters_to_vertices_threshold_ratio
  // So do the following till the condition is met
  while (nb_clusters > nb_vertices * clusters_to_vertices_threshold_ratio)
  {
    typename GT::FT curr_factor = nb_clusters / (nb_vertices * clusters_to_vertices_threshold_ratio);
    int subdivide_steps = (std::max)((int)ceil(log(curr_factor) / log(4)), 0);

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

  // creating needed property maps
  VertexClusterMap vertex_cluster_pmap = get(CGAL::dynamic_vertex_property_t<int>(), pmesh, -1);
  VertexWeightMap vertex_weight_pmap = get(CGAL::dynamic_vertex_property_t<typename GT::FT>(), pmesh, typename GT::FT(0));
  Matrix4x4 zero_mat;  zero_mat.setZero();
  VertexQuadricMap vertex_quadric_pmap = get(CGAL::dynamic_vertex_property_t<Matrix4x4>(), pmesh, zero_mat);

  std::vector<QEMClusterData<GT>> clusters(nb_clusters);
  std::queue<Halfedge_descriptor> clusters_edges_active;
  std::queue<Halfedge_descriptor> clusters_edges_new;

  // compute vertex weights (dual area), and quadrics
  typename GT::FT weight_avg = 0;
  for (Face_descriptor fd : faces(pmesh))
  {
    typename GT::FT weight = abs(CGAL::Polygon_mesh_processing::face_area(fd, pmesh)) / 3;

    // get points of the face
    Halfedge_descriptor hd = halfedge(fd, pmesh);
    typename GT::Point_3 pi = get(vpm, source(hd, pmesh));
    typename GT::Vector_3 vp1(pi.x(), pi.y(), pi.z());
    hd = next(hd, pmesh);
    typename GT::Point_3 pj = get(vpm, source(hd, pmesh));
    typename GT::Vector_3 vp2(pj.x(), pj.y(), pj.z());
    hd = next(hd, pmesh);
    typename GT::Point_3 pk = get(vpm, source(hd, pmesh));
    typename GT::Vector_3 vp3(pk.x(), pk.y(), pk.z());

    // compute quadric for the face
    Matrix4x4 face_quadric;
    if (use_qem_metric)
      compute_qem_face<GT>(vp1, vp2, vp3, face_quadric);

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
        vertex_weight += weight * pow(k_sq, gradation_factor / 2.0);  // /2.0 because k_sq is squared
      }

      weight_avg += vertex_weight;
      put(vertex_weight_pmap, vd, vertex_weight);

      if (use_qem_metric)
      {
        Matrix4x4 vertex_quadric = get(vertex_quadric_pmap, vd);
        vertex_quadric += face_quadric;
        put(vertex_quadric_pmap, vd, vertex_quadric);
      }
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

  // randomly initialize clusters
  //TODO: std::lower_bound with vertex_weight_pmap
  for (int ci = 0; ci < nb_clusters; ++ci)
  {
    int vi;
    Vertex_descriptor vd;
    do {
      vi = CGAL::get_default_random().get_int(0, num_vertices(pmesh));
      vd = *(vertices(pmesh).begin() + vi); // TODO: bad with Polyhedron
    } while (get(vertex_cluster_pmap, vd) != -1);

    put(vertex_cluster_pmap, vd, ci);
    typename GT::Point_3 vp = get(vpm, vd);
    typename GT::Vector_3 vpv(vp.x(), vp.y(), vp.z());
    clusters[ci].add_vertex(vpv, get(vertex_weight_pmap, vd), get(vertex_quadric_pmap, vd));

    for (Halfedge_descriptor hd : halfedges_around_source(vd, pmesh))
      clusters_edges_active.push(hd);
  }

  // the energy minimization loop (clustering loop)
  int nb_modifications = 0;
  int nb_disconnected = 0;

  // Turned on once nb_modifications < nb_vertices * CGAL_TO_QEM_MODIFICATION_THRESHOLD
  bool qem_energy_minimization = false;
  int nb_loops = 0;

  QEMClusterData<GT> cluster1_v1_to_c2, cluster2_v1_to_c2, cluster1_v2_to_c1, cluster2_v2_to_c1;

  do
  {
    // reset cluster data
    for ( auto &cluster : clusters) cluster = QEMClusterData<GT>();
    for ( Vertex_descriptor v : vertices( pmesh ) ) {
      typename GT::FT v_weight = get(vertex_weight_pmap, v);
      Matrix4x4 v_qem = get (vertex_quadric_pmap, v);
      int cluster_id = get(vertex_cluster_pmap, v );
      if ( cluster_id != -1 )
        clusters[ cluster_id ].add_vertex( get(vpm, v), v_weight, v_qem);
    }

    int nb_iterations = -1;
    nb_disconnected = 0;
    do
    {
      ++nb_iterations;
      nb_modifications = 0;

      while (!clusters_edges_active.empty()) {
        Halfedge_descriptor hi = clusters_edges_active.front();
        clusters_edges_active.pop();

        Vertex_descriptor v1 = source(hi, pmesh);
        Vertex_descriptor v2 = target(hi, pmesh);

        int c1 = get(vertex_cluster_pmap, v1);
        int c2 = get(vertex_cluster_pmap, v2);
        if ( ( clusters[ c1 ].last_modification_iteration < nb_iterations - 1 ) &&
          ( clusters[ c2 ].last_modification_iteration < nb_iterations - 1 ) )
          {
            clusters_edges_new.push(hi);
            continue;

          }

        if (c1 == -1)
        {
          // expand cluster c2 (add v1 to c2)
          put(vertex_cluster_pmap, v1, c2);
          typename GT::Point_3 vp1 = get(vpm, v1);
          typename GT::Vector_3 vpv(vp1.x(), vp1.y(), vp1.z());
          clusters[c2].add_vertex(vpv, get(vertex_weight_pmap, v1), get(vertex_quadric_pmap, v1));
          clusters[c2].last_modification_iteration = nb_iterations;

          // add all halfedges around v1 except hi to the queue
          for (Halfedge_descriptor hd : halfedges_around_source(v1, pmesh))
            //TODO: if (hd != hi && hd != opposite(hi, pmesh))
            clusters_edges_new.push(hd);
          ++nb_modifications;
        }
        else if (c2 == -1)
        {
          // expand cluster c1 (add v2 to c1)
          put(vertex_cluster_pmap, v2, c1);
          typename GT::Point_3 vp2 = get(vpm, v2);
          typename GT::Vector_3 vpv(vp2.x(), vp2.y(), vp2.z());
          clusters[c1].add_vertex(vpv, get(vertex_weight_pmap, v2), get(vertex_quadric_pmap, v2));
          clusters[c1].last_modification_iteration = nb_iterations;

          // add all halfedges around v2 except hi to the queue
          for (Halfedge_descriptor hd : halfedges_around_source(v2, pmesh))
            //TODO: if (hd != hi && hd != opposite(hi, pmesh))
            clusters_edges_new.push(hd);
          ++nb_modifications;
        }
        else if ( c1 != c2 )
        {
          // topological test to avoid creating disconnected clusters
          auto is_topologically_valid_merge = [&](Halfedge_descriptor hv, int cluster_id)
          {
            return true;
            // TODO : result seems buggy with this ON
            std::stringstream ss;
            ss << target(hv, pmesh) << " in cluster " << cluster_id <<  "(" << get(vertex_cluster_pmap,source(hv, pmesh)) << ")\n";

            CGAL_assertion(get(vertex_cluster_pmap,target(hv, pmesh))==cluster_id);
            Halfedge_descriptor h=hv;
            bool in_cluster=false;
            int nb_cc_cluster=0;
            do{
              h=next(h, pmesh);

              int ci = get(vertex_cluster_pmap, target(h,pmesh));
              ss << "  " << target(h, pmesh) << " in cluster " << ci << "\n";
              if (in_cluster)
              {
                if (ci!=cluster_id) in_cluster=false;
              }
              else
              {
                if (ci==cluster_id)
                {
                  in_cluster=true;
                  if (++nb_cc_cluster>1) {
                    std::cout << "BOOM\n";
                    std::cout << ss.str();
                    return false;
                  }
                }
              }
              h=opposite(h, pmesh);
            }
            while(h!=hv);

            return true;
          };

          // compare the energy of the 3 cases
          typename GT::Point_3 vp1 = get(vpm, v1);
          typename GT::Vector_3 vpv1(vp1.x(), vp1.y(), vp1.z());
          typename GT::Point_3 vp2 = get(vpm, v2);
          typename GT::Vector_3 vpv2(vp2.x(), vp2.y(), vp2.z());
          typename GT::FT v1_weight = get(vertex_weight_pmap, v1);
          typename GT::FT v2_weight = get(vertex_weight_pmap, v2);
          Matrix4x4 v1_qem = get (vertex_quadric_pmap, v1);
          Matrix4x4 v2_qem = get (vertex_quadric_pmap, v2);

          cluster1_v2_to_c1 = clusters[c1];
          cluster2_v2_to_c1 = clusters[c2];
          cluster1_v1_to_c2 = clusters[c1];
          cluster2_v1_to_c2 = clusters[c2];
          cluster1_v2_to_c1.last_modification_iteration = nb_iterations;
          cluster2_v2_to_c1.last_modification_iteration = nb_iterations;
          cluster1_v1_to_c2.last_modification_iteration = nb_iterations;
          cluster2_v1_to_c2.last_modification_iteration = nb_iterations;

          typename GT::FT e_no_change = clusters[c1].compute_energy(qem_energy_minimization) +
                                        clusters[c2].compute_energy(qem_energy_minimization);
          typename GT::FT e_v1_to_c2 = (std::numeric_limits< double >::max)();
          typename GT::FT e_v2_to_c1 = (std::numeric_limits< double >::max)();

          if ( ( clusters[ c1 ].nb_vertices > 1 ) && ( !qem_energy_minimization || is_topologically_valid_merge(opposite(hi, pmesh), c1) ) )
          {
            cluster1_v1_to_c2.remove_vertex(vpv1, v1_weight, v1_qem);
            cluster2_v1_to_c2.add_vertex(vpv1, v1_weight, v1_qem);
            e_v1_to_c2 = cluster1_v1_to_c2.compute_energy(qem_energy_minimization) +
                         cluster2_v1_to_c2.compute_energy(qem_energy_minimization);
          }

          if ( ( clusters[ c2 ].nb_vertices > 1 ) && ( !qem_energy_minimization || is_topologically_valid_merge(hi, c2) ) )
          {
            cluster1_v2_to_c1.add_vertex(vpv2, v2_weight, v2_qem);
            cluster2_v2_to_c1.remove_vertex(vpv2, v2_weight, v2_qem);
            e_v2_to_c1 = cluster1_v2_to_c1.compute_energy(qem_energy_minimization) +
                         cluster2_v2_to_c1.compute_energy(qem_energy_minimization);
          }


          if (e_v2_to_c1 < e_no_change && e_v2_to_c1 < e_v1_to_c2)
          {
            // move v2 to c1
            put(vertex_cluster_pmap, v2, c1);
            clusters[c1] = cluster1_v2_to_c1;
            clusters[c2] = cluster2_v2_to_c1;
            // add all halfedges around v2 except hi to the queue
            for (Halfedge_descriptor hd : halfedges_around_source(v2, pmesh))
              //TODO: if (hd != hi && hd != opposite(hi, pmesh))
              clusters_edges_new.push(hd);
            ++nb_modifications;
          }
          else if (e_v1_to_c2 < e_no_change) // > 2 as 1 vertex was added to c1
          {
            // move v1 to c2
            put(vertex_cluster_pmap, v1, c2);
            clusters[c1] = cluster1_v1_to_c2;
            clusters[c2] = cluster2_v1_to_c2;
            // add all halfedges around v1 except hi to the queue
            for (Halfedge_descriptor hd : halfedges_around_source(halfedge(v1, pmesh), pmesh))
              //TODO: if (hd != hi && hd != opposite(hi, pmesh))
              clusters_edges_new.push(hd);
            ++nb_modifications;
          }
          else
          {
            // no change
            clusters_edges_new.push(hi);
          }
        }
      }
      std::cout << "# Modifications: " << nb_modifications << "\n";

      clusters_edges_active.swap(clusters_edges_new);
      if (use_qem_metric && !qem_energy_minimization && nb_modifications < nb_vertices * CGAL_TO_QEM_MODIFICATION_THRESHOLD)
      {
        qem_energy_minimization = true;
        break;
      }

    } while (nb_modifications > 0);

    // Disconnected clusters handling
    // the goal is to delete clusters with multiple connected components and only keep the largest connected component of each cluster
    // For each cluster, do a BFS from a vertex in the cluster

    std::vector<std::vector<std::vector<Vertex_descriptor>>> cluster_components(nb_clusters, std::vector<std::vector<Vertex_descriptor>>());
    std::queue<Vertex_descriptor> vertex_queue;

    // loop over vertices to compute cluster components
    VertexVisitedMap vertex_visited_pmap = get(CGAL::dynamic_vertex_property_t<bool>(), pmesh, false);
    for (Vertex_descriptor vd : vertices(pmesh))
    {
      if (get(vertex_visited_pmap, vd)) continue;
      int c = get(vertex_cluster_pmap, vd);
      if (c != -1)
      {
        cluster_components[c].push_back(std::vector<Vertex_descriptor>());
        int component_i = cluster_components[c].size() - 1;

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

    // loop over clusters to delete disconnected components except the largest one
    for (int c = 0; c < nb_clusters; ++c)
    {
      if (cluster_components[c].size() <= 1) continue; // only one component, no need to do anything
      ++nb_disconnected;
      std::size_t max_component_size = 0;
      std::size_t max_component_index = -1;
      for (std::size_t component_i = 0; component_i < cluster_components[c].size(); ++component_i)
      {
        if (cluster_components[c][component_i].size() > max_component_size)
        {
          max_component_size = cluster_components[c][component_i].size();
          max_component_index = component_i;
        }
      }
      // set cluster to -1 for all components except the largest one
      for (std::size_t component_i = 0; component_i < cluster_components[c].size(); ++component_i)
      {
        if (component_i != max_component_index)
        {
          for (Vertex_descriptor vd : cluster_components[c][component_i])
          {
            put(vertex_cluster_pmap, vd, -1);
            // remove vd from cluster c
            typename GT::Point_3 vp = get(vpm, vd);
            typename GT::Vector_3 vpv(vp.x(), vp.y(), vp.z());
            clusters[c].remove_vertex(vpv, get(vertex_weight_pmap, vd), get(vertex_quadric_pmap, vd));
            // add all halfedges around v except hi to the queue
            for (Halfedge_descriptor hd : halfedges_around_source(vd, pmesh))
            {
              // add hd to the queue if its target is not in the same cluster
              Vertex_descriptor v2 = target(hd, pmesh);
              int c2 = get(vertex_cluster_pmap, v2);
              if (c2 != c)
                clusters_edges_active.push(hd);
            }
          }
        }
      }
    }



    std::cout << "# nb_disconnected: " << nb_disconnected << "\n";
    ++nb_loops;

  } while (nb_disconnected > 0 || nb_loops < 2 );

  /// Construct new Mesh
  std::vector<int> valid_cluster_map(nb_clusters, -1);
  std::vector<typename GT::Point_3> points;

  std::vector<std::vector<int>> polygons;
  TriangleMesh simplified_mesh;


  for (int i = 0; i < nb_clusters; ++i)
  {
    if (clusters[i].weight_sum > 0)
    {
      valid_cluster_map[i] = points.size();
      points.push_back(clusters[i].representative_point(qem_energy_minimization));
    }
  }

  if (postprocessing_qem){
    // create a point for each cluster
    std::vector<Eigen::Matrix<typename GT::FT, 4, 4>> cluster_quadrics(clusters.size());

    // initialize quadrics
    for (int i = 0; i < nb_clusters; ++i)
      cluster_quadrics[i].setZero();

    // for storing the vertex_descriptor of each face
    std::vector<Vertex_descriptor> face_vertices;

    for (Face_descriptor fd : faces(pmesh)) {
      // get Vs for fd
      // compute qem from Vs->"vector_3"s
      // add to the 3 indices of the cluster
      Halfedge_descriptor hd = halfedge(fd, pmesh);
      do {
        Vertex_descriptor vd = target(hd, pmesh);
        face_vertices.push_back(vd);
        hd = next(hd, pmesh);
      } while (hd != halfedge(fd, pmesh));

      auto p_i = get(vpm, face_vertices[0]);
      typename GT::Vector_3 vec_1 = typename GT::Vector_3(p_i.x(), p_i.y(), p_i.z());
      p_i = get(vpm, face_vertices[1]);
      typename GT::Vector_3 vec_2 = typename GT::Vector_3(p_i.x(), p_i.y(), p_i.z());
      p_i = get(vpm, face_vertices[2]);
      typename GT::Vector_3 vec_3 = typename GT::Vector_3(p_i.x(), p_i.y(), p_i.z());

      int c_1 = get(vertex_cluster_pmap, face_vertices[0]);
      int c_2 = get(vertex_cluster_pmap, face_vertices[1]);
      int c_3 = get(vertex_cluster_pmap, face_vertices[2]);

      if (c_1 != -1 && c_2 != -1 && c_3 != -1)
      {
        Eigen::Matrix<typename GT::FT, 4, 4> q;
        compute_qem_face<GT>(vec_1, vec_2, vec_3, q);
        cluster_quadrics[c_1] += q;
        cluster_quadrics[c_2] += q;
        cluster_quadrics[c_3] += q;
      }

      face_vertices.clear();
    }

    int valid_index = 0;

    for (int i = 0; i < nb_clusters; ++i)
    {
      if (clusters[i].weight_sum > 0)
      {
        int rank_deficiency = 0;
        compute_representative_point<GT>(cluster_quadrics[i], points[valid_index], rank_deficiency);
        ++valid_index;
      }
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

  // create a triangle for each face with all vertices in 3 different clusters
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
* performs isotropic centroidal voronoi diagram remeshing on a triangle mesh. The remeshing is either uniform or adaptative.
*
* @tparam TriangleMesh a model of `FaceListGraph`
* @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters".
*
* @param tmesh input triangle mesh
* @param nb_vertices lower bound on the number of target vertices in the output mesh.
*                    In the case the mesh is not closed or if the number of points is too low
*                    and no manifold mesh could be produced with that budget of points, extra points
*                    are added to get a manifold output.
* @param np optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
*        `GT` stands for the type of the object provided to the named parameter `geom_traits()`.
*
* \cgalNamedParamsBegin
*
*   \cgalParamNBegin{vertex_principal_curvatures_and_directions_map}
*     \cgalParamDescription{a property map associating principal curvatures and directions to the vertices of `tmesh`, used for adaptive clustering.}
*     \cgalParamType{a class model of `ReadWritePropertyMap` with
*                    `boost::graph_traits<TriangleMesh>::%vertex_descriptor`
*                    as key type and `Principal_curvatures_and_directions<GT>` as value type.}
*     \cgalParamExtra{If this parameter is omitted, but `gradation_factor` is provided, an internal property map
*                     will be created and curvature values will be computed using the function `interpolated_corrected_curvatures()` will be called to initialize it.}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{gradation_factor}
*     \cgalParamDescription{a factor used to gradate the weights of the vertices based on their curvature values.}
*     \cgalParamType{`GT::FT`}
*     \cgalParamDefault{0}
*     \cgalParamExtra{If this parameter is omitted, no adaptive clustering will be performed.}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{postprocessing_qem}
*     \cgalParamDescription{indicates if a project step using QEM should be applied to cluster centers at the end of the minimization,
*                           in order for example to recover sharp features.
*                           This is a fast method but can result in visual issues or even self-intersections.}
*     \cgalParamType{bool}
*     \cgalParamDefault{false}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{use_qem_metric}
*     \cgalParamDescription{indicates if QEM should be applied during the minimization algorithm in order for example to recover sharp features.
*                           This is a slower than when using `postprocessing_qem(true)` but is more accurate.}
*     \cgalParamType{bool}
*     \cgalParamDefault{false}
*     \cgalParamExtra{If this parameter is `true` then postprocessing_qem will automatically set to `false`.}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{clusters_to_vertices_threshold_ratio}
*     \cgalParamDescription{a ratio used to subdivide the input mesh so that the mesh used as input for the clustering has enough vertices.
*                           More precisely, the number of vertices of the subdivided mesh should at least be the ratio times `nb_vertices`.}
*     \cgalParamType{`GT::FT`}
*     \cgalParamDefault{0.1}
*     \cgalParamExtra{A value between 0.1 and 0.01 is recommended, the smaller the better the approximation will be.}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{vertex_point_map}
*       \cgalParamDescription{a property map associating points to the vertices of `tmesh`.}
*       \cgalParamType{a class model of `ReadablePropertyMap` with
*                      `boost::graph_traits<TriangleMesh>::%vertex_descriptor`
*                      as key type and `GT::Point_3` as value type.}
*       \cgalParamDefault{`boost::get(CGAL::vertex_point, tmesh)`.}
*       \cgalParamExtra{If this parameter is omitted, an internal property map for
*                       `CGAL::vertex_point_t` must be available in `TriangleMesh`.}
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
* @return the simplified mesh as a TriangleMesh
*
* @todo implement manifold version
* @todo how to handle output vertex point map
*/
template <typename TriangleMesh,
          typename NamedParameters = parameters::Default_named_parameters>
TriangleMesh acvd_isotropic_remeshing(
    TriangleMesh& tmesh,
    const int& nb_vertices,
    const NamedParameters& np = parameters::default_values()
  )
{
  auto ps = internal::acvd_impl(tmesh, nb_vertices, np);

  CGAL_assertion(is_polygon_soup_a_polygon_mesh(ps.second));

  TriangleMesh simplified_mesh;
  polygon_soup_to_polygon_mesh(ps.first, ps.second, simplified_mesh);
  return simplified_mesh;
}

} // namespace Polygon_mesh_processing

} // namespace CGAL

#undef CGAL_WEIGHT_CLAMP_RATIO_THRESHOLD
#undef CGAL_TO_QEM_MODIFICATION_THRESHOLD

#endif // CGAL_PMP_ACVD_REMESHING_H
