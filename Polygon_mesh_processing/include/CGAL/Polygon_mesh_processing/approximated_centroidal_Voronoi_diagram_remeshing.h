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

#ifndef CGAL_ACVD_DOES_NOT_USE_INTERPOLATED_CORRECTED_CURVATURES
#include <CGAL/Polygon_mesh_processing/interpolated_corrected_curvatures.h>
#endif
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/border.h>
#ifdef CGAL_DEBUG_ACVD
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/IO/Color.h>
#endif

#include <CGAL/assertions.h>
#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/Named_function_parameters.h>
#include <CGAL/property_map.h>
#include <CGAL/Random.h>
#include <CGAL/determinant.h>
#include <CGAL/subdivision_method_3.h>
#include <CGAL/utility.h>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>

#include <array>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <numeric>
#include <queue>
#include <set>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

// Early convergence ratio which is used for the first two minimization steps.
// The final convergence is carried out entirely i.e. convergence is reached only when no more optimization can be performed.
#define CGAL_TO_QEM_MODIFICATION_THRESHOLD 1e-3

// Used to clamp curvature values:
//  - regions with zero curvature (although this is not really a degenerate case....)
//  - regions with too high computed curvature (sometimes due to numerical issues in curvature computation)
// Clamping happens this way ;
//  -The average weight weight_avg is computed
//      weights higher than CGAL_WEIGHT_CLAMP_RATIO_THRESHOLD * weight_avg are set to this value
//      weights lower than 1/( CGAL_WEIGHT_CLAMP_RATIO_THRESHOLD * weight_avg ) are set to this value
#define CGAL_WEIGHT_CLAMP_RATIO_THRESHOLD 10000

namespace CGAL {
namespace Polygon_mesh_processing {
namespace internal {

#ifdef CGAL_DEBUG_ACVD
template <class TriangleMesh, class ClusterMap>
void dump_mesh_with_cluster_colors(TriangleMesh tmesh, ClusterMap cluster_map, std::string fname)
{
  std::vector<CGAL::IO::Color> palette {{ CGAL::IO::red(),
                                          CGAL::IO::green(),
                                          CGAL::IO::blue(),
                                          CGAL::IO::purple(),
                                          CGAL::IO::orange(),
                                          CGAL::IO::deep_blue(),
                                          CGAL::IO::yellow(),
                                          CGAL::IO::violet(),
                                          CGAL::IO::gray(),
                                          CGAL::IO::white() }};

  auto vcm = tmesh.template add_property_map<typename TriangleMesh::Vertex_index, CGAL::IO::Color>("f:color").first;

  for (auto v : vertices(tmesh))
  {
    int cluster_id=get(cluster_map, v);
    if (cluster_id!=-1)
      put(vcm, v, palette[ cluster_id % palette.size() ]);
    else
      put(vcm, v, CGAL::IO::black());
  }

  std::ofstream out(fname);
  CGAL::IO::write_PLY(out, tmesh, CGAL::parameters::use_binary_mode(false)
                                                   .stream_precision(17)
                                                   .vertex_color_map(vcm));

}
#endif // CGAL_DEBUG_ACVD

template <typename GT>
void compute_qem_face(const typename GT::Vector_3& p1, const typename GT::Vector_3& p2, const typename GT::Vector_3& p3,
                      Eigen::Matrix<typename GT::FT, 4, 4>& quadric,
                      const GT& gt)
{
  auto cross_product = gt.construct_cross_product_vector_3_object();

  typename GT::Vector_3 crossX1X2 = cross_product(p1, p2);
  typename GT::Vector_3 crossX2X3 = cross_product(p2, p3);
  typename GT::Vector_3 crossX3X1 = cross_product(p3, p1);
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
typename GT::Vector_3 compute_displacement(const Eigen::Matrix<typename GT::FT, 4, 4>& quadric,
                                           const typename GT::Point_3& p,
                                           int& rank_deficiency)
{
  using Matrix3d = Eigen::Matrix<typename GT::FT, 3, 3>;

  int max_nb_of_singular_values_used = 3;
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

    if ((AbsolutesEigenValues[IndexMax] * invmaxW > 1e-3) &&
        (max_nb_of_singular_values_used > 0))
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
    --max_nb_of_singular_values_used;
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
void compute_representative_point(const Eigen::Matrix<typename GT::FT, 4, 4>& quadric,
                                  typename GT::Point_3& p,
                                  int& rank_deficiency)
{
  // average point
  typename GT::Vector_3 displacement = compute_displacement<GT>(quadric, p, rank_deficiency);
  p = p + displacement;
}

template <typename GT>
struct QEMClusterData
{
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

  void add_vertex(const typename GT::Vector_3& vertex_position,
                  const typename GT::FT weight,
                  const Eigen::Matrix<typename GT::FT, 4, 4>& quadric)
  {
    this->site_sum += vertex_position * weight;
    this->weight_sum += weight;
    this->quadric_sum += quadric;
    ++this->nb_vertices;
    this->modified = true;
  }

  void add_vertex(const typename GT::Point_3& vertex_position,
                  const typename GT::FT weight,
                  const Eigen::Matrix<typename GT::FT, 4, 4>& quadric)
  {
    add_vertex(vertex_position - ORIGIN, weight, quadric);
  }

  void remove_vertex(const typename GT::Vector_3& vertex_position,
                     const typename GT::FT weight,
                     const Eigen::Matrix<typename GT::FT, 4, 4>& quadric)
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
      {
        this->representative_point_ = (this->site_sum) / this->weight_sum;
      }
    }
  }

  typename GT::FT compute_energy(bool qem)
  {
    if (modified)
    {
      if (qem)
      {
        compute_representative(qem);
        auto dot_product = GT().compute_scalar_product_3_object();

        this->energy = (this->representative_point_).squared_length() * this->weight_sum
                       - 2 * dot_product(this->representative_point_, this->site_sum);
      }
      else
      {
        this->energy = - (this->site_sum).squared_length() / this->weight_sum;
      }
    }
    return this->energy;
  }

  typename GT::Point_3 representative_point(bool qem)
  {
    if (modified)
      compute_representative(qem);
    return ORIGIN + representative_point_;
  }
};

#ifndef CGAL_ACVD_DOES_NOT_USE_INTERPOLATED_CORRECTED_CURVATURES
template <class TriangleMesh, class VPCDM, class NamedParameters>
void upsample_subdivision_property(TriangleMesh& tmesh,
                                   VPCDM vpcd_map,
                                   const NamedParameters& np)
{
  using GT = typename GetGeomTraits<TriangleMesh, NamedParameters>::type;
  using vertex_descriptor = typename boost::graph_traits<TriangleMesh>::vertex_descriptor;
  using halfedge_descriptor = typename boost::graph_traits<TriangleMesh>::halfedge_descriptor;

  using parameters::choose_parameter;
  using parameters::get_parameter;
  using parameters::is_default_parameter;

  using VPM = typename CGAL::GetVertexPointMap<TriangleMesh, NamedParameters>::type;
  VPM vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                         get_property_map(CGAL::vertex_point, tmesh));

  // unordered_set of old vertices
  std::unordered_set<vertex_descriptor> old_vertices;

  unsigned int step = choose_parameter(get_parameter(np, internal_np::number_of_iterations), 1);
  Linear_mask_3<TriangleMesh,VPM> mask(&tmesh, vpm);

  for (unsigned int i = 0; i < step; ++i)
  {
    for (vertex_descriptor vd : vertices(tmesh))
      old_vertices.insert(vd);

    Subdivision_method_3::internal::PTQ_1step(tmesh, vpm, mask);
    // interpolate curvature values
    for (vertex_descriptor vd : vertices(tmesh))
    {
      if (old_vertices.find(vd) == old_vertices.end())
      {
        Principal_curvatures_and_directions<GT> pcd;
        pcd.min_curvature = 0;
        pcd.max_curvature = 0;
        pcd.min_direction = typename GT::Vector_3(0, 0, 0);
        pcd.max_direction = typename GT::Vector_3(0, 0, 0);
        for (halfedge_descriptor hd : halfedges_around_target(vd, tmesh))
        {
          vertex_descriptor v1 = source(hd, tmesh);
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
#endif

template <typename TriangleMesh,
          typename NamedParameters = parameters::Default_named_parameters>
std::pair<std::vector<typename GetGeomTraits<TriangleMesh, NamedParameters>::type::Point_3>,
          std::vector<std::array<std::size_t, 3> > >
acvd_impl(TriangleMesh& tmesh,
          int nb_clusters,
          const NamedParameters& np = parameters::default_values())
{
  using GT = typename GetGeomTraits<TriangleMesh, NamedParameters>::type;
  using FT = typename GT::FT;
  using Point_3 = typename GT::Point_3;
  using Vector_3 = typename GT::Vector_3;

  using Matrix4x4 = typename Eigen::Matrix<FT, 4, 4>;

  using Vertex_position_map = typename GetVertexPointMap<TriangleMesh, NamedParameters>::const_type;

  using halfedge_descriptor = typename boost::graph_traits<TriangleMesh>::halfedge_descriptor;
  using edge_descriptor = typename boost::graph_traits<TriangleMesh>::edge_descriptor;
  using vertex_descriptor = typename boost::graph_traits<TriangleMesh>::vertex_descriptor;
  using face_descriptor = typename boost::graph_traits<TriangleMesh>::face_descriptor;

  using VertexClusterMap = typename boost::property_map<TriangleMesh, CGAL::dynamic_vertex_property_t<int> >::type;
  using VertexVisitedMap = typename boost::property_map<TriangleMesh, CGAL::dynamic_vertex_property_t<bool> >::type;
  using VertexWeightMap = typename boost::property_map<TriangleMesh, CGAL::dynamic_vertex_property_t<FT> >::type;
  using VertexQuadricMap = typename boost::property_map<TriangleMesh, CGAL::dynamic_vertex_property_t<Matrix4x4> >::type;

#ifndef CGAL_ACVD_DOES_NOT_USE_INTERPOLATED_CORRECTED_CURVATURES
  using Default_principal_map = typename boost::property_map<TriangleMesh, CGAL::dynamic_vertex_property_t<Principal_curvatures_and_directions<GT>> >::type;
  using Vertex_principal_curvatures_and_directions_map =
    typename internal_np::Lookup_named_param_def<internal_np::vertex_principal_curvatures_and_directions_map_t,
                                                 NamedParameters, Default_principal_map>::type ;
#endif

  using parameters::choose_parameter;
  using parameters::get_parameter;
  using parameters::is_default_parameter;

  GT gt = choose_parameter<GT>(get_parameter(np, internal_np::geom_traits));
  Vertex_position_map vpm = choose_parameter(get_parameter(np, CGAL::vertex_point),
                                             get_property_map(CGAL::vertex_point, tmesh));

  // get curvature related parameters
#ifdef CGAL_ACVD_DOES_NOT_USE_INTERPOLATED_CORRECTED_CURVATURES
  static_assert(is_default_parameter<NamedParameters, internal_np::gradation_factor_t>::value, "gradation_factor option is disabled");
#endif
  const FT gradation_factor = choose_parameter(get_parameter(np, internal_np::gradation_factor), 0);

  // using QEM ?
  bool use_postprocessing_qem = choose_parameter(get_parameter(np, internal_np::use_postprocessing_qem), false);
  const bool use_qem_based_energy = choose_parameter(get_parameter(np, internal_np::use_qem_based_energy), false);
  if (use_qem_based_energy)
    use_postprocessing_qem = false;

#ifndef CGAL_ACVD_DOES_NOT_USE_INTERPOLATED_CORRECTED_CURVATURES
  Vertex_principal_curvatures_and_directions_map vpcd_map;

  // if adaptive clustering
  if (gradation_factor > 0)
  {
    if constexpr (is_default_parameter<NamedParameters, internal_np::vertex_principal_curvatures_and_directions_map_t>::value)
    {
      vpcd_map = get(CGAL::dynamic_vertex_property_t<Principal_curvatures_and_directions<GT>>(), tmesh);
      interpolated_corrected_curvatures(tmesh, parameters::vertex_principal_curvatures_and_directions_map(vpcd_map)
                                                          .vertex_point_map(vpm).geom_traits(gt));
    }
    else
      vpcd_map = get_parameter(np, internal_np::vertex_principal_curvatures_and_directions_map);
  }
#endif

  CGAL_precondition(CGAL::is_triangle_mesh(tmesh));

  const double vertex_count_ratio = choose_parameter(get_parameter(np, internal_np::vertex_count_ratio), 0.1);
  // TODO: (possible optimization to save memory) if qem is not used use std::conditionnal to store and compute less things in QEMClusterData
  // TODO: (possible optimization to save memort) use a std::array to store quadrics and assemble the Eigen matric only at the end

  std::size_t nb_vertices = vertices(tmesh).size();

  // For remeshing, we might need to subdivide the mesh before clustering
  // This should always hold: nb_clusters <= nb_vertices * vertex_count_ratio
  // So do the following until the condition is met
  if (gradation_factor == 0 && nb_clusters > nb_vertices * vertex_count_ratio)
  {
    std::vector<FT> lengths;
    lengths.reserve(num_edges(tmesh));
    FT cum = 0;
    int nbe = 0;
    for (edge_descriptor e : edges(tmesh))
    {
      lengths.push_back(edge_length(halfedge(e, tmesh), tmesh, parameters::vertex_point_map(vpm).geom_traits(gt)));
      cum += lengths.back();
      ++nbe;
    }
    FT threshold = 3. * cum / nbe;

    std::vector<std::pair<edge_descriptor, double>> to_split;
    int ei = -1;
    for (edge_descriptor e : edges(tmesh))
    {
      FT l = lengths[++ei];
      if (l >= threshold)
      {
        double nb_subsegments = std::ceil(l / threshold);
        to_split.emplace_back(e, nb_subsegments);
      }
    }

    while (!to_split.empty())
    {
      std::vector<std::pair<edge_descriptor, double>> to_split_new;
      for (auto [e, nb_subsegments] : to_split)
      {
        halfedge_descriptor h = halfedge(e, tmesh);
        Point_3 s = get(vpm, source(h,tmesh));
        Point_3 t = get(vpm, target(h,tmesh));

        for (double k=1; k<nb_subsegments; ++k)
        {
          // the new halfedge hnew pointing to the inserted vertex.
          // The new halfedge is followed by the old halfedge, i.e., next(hnew,g) == h.
          halfedge_descriptor hnew = Euler::split_edge(h, tmesh);
          put(vpm, target(hnew, tmesh), barycenter(t, k/nb_subsegments, s));
          for (halfedge_descriptor hhh : {hnew, opposite(h, tmesh)})
          {
            if (!is_border(hhh, tmesh))
            {
              halfedge_descriptor hf = Euler::split_face(hhh, next(next(hhh, tmesh), tmesh), tmesh);
              FT l = edge_length(hf, tmesh, parameters::vertex_point_map(vpm).geom_traits(gt));
              if (l >= threshold)
              {
                double nb_subsegments = std::ceil(l / threshold);
                to_split_new.emplace_back(edge(hf,tmesh), nb_subsegments);
              }
            }
          }
        }
      }
      to_split.swap(to_split_new);
    }
    nb_vertices = vertices(tmesh).size();
  }

  while (nb_clusters > nb_vertices * vertex_count_ratio)
  {
    FT curr_factor = nb_clusters / (nb_vertices * vertex_count_ratio);
    int subdivide_steps = (std::max)((int)ceil(log(curr_factor) / log(4)), 0);

    if (subdivide_steps > 0)
    {
      if (gradation_factor == 0) // no adaptive clustering
      {
        Subdivision_method_3::linear_subdivision(tmesh, CGAL::parameters::number_of_iterations(subdivide_steps)
                                                                         .vertex_point_map(vpm)
                                                                         .geom_traits(gt));
      }
#ifndef CGAL_ACVD_DOES_NOT_USE_INTERPOLATED_CORRECTED_CURVATURES
      else // adaptive clustering
      {
        upsample_subdivision_property(tmesh, vpcd_map,
                                      CGAL::parameters::number_of_iterations(subdivide_steps)
                                                       .vertex_principal_curvatures_and_directions_map(vpcd_map)
                                                       .vertex_point_map(vpm)
                                                       .geom_traits(gt));
      }
#endif
      vpm = get_property_map(CGAL::vertex_point, tmesh);
      nb_vertices = vertices(tmesh).size();
    }
  }

  // creating needed property maps
  VertexClusterMap vertex_cluster_pmap = get(CGAL::dynamic_vertex_property_t<int>(), tmesh, -1);
  VertexWeightMap vertex_weight_pmap = get(CGAL::dynamic_vertex_property_t<FT>(), tmesh, FT(0));
  Matrix4x4 zero_mat;
  zero_mat.setZero();
  VertexQuadricMap vertex_quadric_pmap = get(CGAL::dynamic_vertex_property_t<Matrix4x4>(), tmesh, zero_mat);

  std::vector<QEMClusterData<GT>> clusters(nb_clusters);
  std::queue<halfedge_descriptor> clusters_edges_active;
  std::queue<halfedge_descriptor> clusters_edges_new;

  // compute vertex weights (dual area), and quadrics
  FT weight_avg = 0;
  for (face_descriptor fd : faces(tmesh))
  {
    FT weight = abs(face_area(fd, tmesh, parameters::vertex_point_map(vpm)
                                                    .geom_traits(gt))) / 3;

    // get points of the face
    halfedge_descriptor hd = halfedge(fd, tmesh);
    Point_3 pi = get(vpm, source(hd, tmesh));
    Vector_3 vp1(pi.x(), pi.y(), pi.z());
    hd = next(hd, tmesh);
    Point_3 pj = get(vpm, source(hd, tmesh));
    Vector_3 vp2(pj.x(), pj.y(), pj.z());
    hd = next(hd, tmesh);
    Point_3 pk = get(vpm, source(hd, tmesh));
    Vector_3 vp3 = pk - ORIGIN;

    // compute quadric for the face
    Matrix4x4 face_quadric;
    if (use_qem_based_energy)
      compute_qem_face<GT>(vp1, vp2, vp3, face_quadric, gt);

    for (vertex_descriptor vd : vertices_around_face(halfedge(fd, tmesh), tmesh))
    {
      FT vertex_weight = get(vertex_weight_pmap, vd);

      if (gradation_factor == 0) // no adaptive clustering
        vertex_weight += weight;
#ifndef CGAL_ACVD_DOES_NOT_USE_INTERPOLATED_CORRECTED_CURVATURES
      else // adaptive clustering
      {
        FT k1 = get(vpcd_map, vd).min_curvature;
        FT k2 = get(vpcd_map, vd).max_curvature;
        FT k_sq = (k1 * k1 + k2 * k2);
        vertex_weight += weight * pow(k_sq, gradation_factor / 2.0);  // /2.0 because k_sq is squared
      }
#endif
      weight_avg += vertex_weight;
      put(vertex_weight_pmap, vd, vertex_weight);

      if (use_qem_based_energy)
      {
        Matrix4x4 vertex_quadric = get(vertex_quadric_pmap, vd);
        vertex_quadric += face_quadric;
        put(vertex_quadric_pmap, vd, vertex_quadric);
      }
    }
  }
  weight_avg /= nb_vertices;

  // clamp the weights up and below by a ratio (like 10,000) * avg_weights
  for (vertex_descriptor vd : vertices(tmesh))
  {
    FT vertex_weight = get(vertex_weight_pmap, vd);
    if (vertex_weight > CGAL_WEIGHT_CLAMP_RATIO_THRESHOLD * weight_avg)
      put(vertex_weight_pmap, vd, CGAL_WEIGHT_CLAMP_RATIO_THRESHOLD * weight_avg);
    else if (vertex_weight < 1.0 / CGAL_WEIGHT_CLAMP_RATIO_THRESHOLD * weight_avg)
      put(vertex_weight_pmap, vd, 1.0 / CGAL_WEIGHT_CLAMP_RATIO_THRESHOLD * weight_avg);
  }

  Random rnd = is_default_parameter<NamedParameters, internal_np::random_seed_t>::value
             ? get_default_random()
             : CGAL::Random(choose_parameter(get_parameter(np, internal_np::random_seed),0));

  // randomly initialize clusters
  //TODO: (possible optimization) use std::lower_bound with vertex_weight_pmap for better sampling
  //      moreover with non-random access iterators, the complexity will be bad
  for (int ci = 0; ci < nb_clusters; ++ci)
  {
    int vi;
    vertex_descriptor vd;
    do {
      vi = rnd.get_int(0, nb_vertices);
      vd = *(vertices(tmesh).begin() + vi);
    } while (get(vertex_cluster_pmap, vd) != -1);

    put(vertex_cluster_pmap, vd, ci);
    Point_3 vp = get(vpm, vd);
    Vector_3 vpv = vp - ORIGIN;
    clusters[ci].add_vertex(vpv, get(vertex_weight_pmap, vd), get(vertex_quadric_pmap, vd));

    for (halfedge_descriptor hd : halfedges_around_source(vd, tmesh))
      clusters_edges_active.push(hd);
  }

  // the energy minimization loop (clustering loop)
  int nb_modifications = 0;
  int nb_disconnected = 0;

  // Turned on once nb_modifications < nb_vertices * CGAL_TO_QEM_MODIFICATION_THRESHOLD
  bool qem_energy_minimization = false;
  int nb_loops = 0;

  QEMClusterData<GT> cluster1_v1_to_c2, cluster2_v1_to_c2, cluster1_v2_to_c1, cluster2_v2_to_c1;
  std::vector<bool> frozen_clusters(nb_clusters, false);
  do
  {
    do
    {
      // reset cluster data
      for (int ci=0; ci<nb_clusters; ++ci)
        if (!frozen_clusters[ci])
          clusters[ci] = QEMClusterData<GT>();

      for (vertex_descriptor v : vertices(tmesh))
      {
        FT v_weight = get(vertex_weight_pmap, v);
        Matrix4x4 v_qem = get(vertex_quadric_pmap, v);
        int cluster_id = get(vertex_cluster_pmap, v );
        if (cluster_id != -1 && !frozen_clusters[cluster_id])
          clusters[cluster_id].add_vertex(get(vpm, v), v_weight, v_qem);
      }

      CGAL_assertion_code(for (int ci=0; ci<nb_clusters; ++ci))
      CGAL_assertion(clusters[ci].nb_vertices != 0);

      int nb_iterations = -1;
      nb_disconnected = 0;

      auto edge_seen = get(CGAL::dynamic_edge_property_t<int>(), tmesh, -1);
      do
      {
        ++nb_iterations;
        nb_modifications = 0;

        while (!clusters_edges_active.empty())
        {
          halfedge_descriptor hi = clusters_edges_active.front();
          clusters_edges_active.pop();

          if (get(edge_seen, edge(hi, tmesh)) == nb_iterations)
            continue;
          put(edge_seen, edge(hi, tmesh), nb_iterations);

          vertex_descriptor v1 = source(hi, tmesh);
          vertex_descriptor v2 = target(hi, tmesh);

          // add all halfedges around v1 except hi to the queue
          auto push_vertex_edge_ring_to_queue = [&](vertex_descriptor v)
          {
            for (halfedge_descriptor hd : halfedges_around_source(v, tmesh))
              clusters_edges_new.push(hd);
          };

          int c1 = get(vertex_cluster_pmap, v1);
          int c2 = get(vertex_cluster_pmap, v2);

          if ((c1 != -1 && frozen_clusters[c1]) || (c2 != -1 && frozen_clusters[c2]))
          {
            clusters_edges_new.push(hi);
            continue;
          }

          if (c1 == c2)
            continue;

          if (c1 == -1)
          {
            // expand cluster c2 (add v1 to c2)
            put(vertex_cluster_pmap, v1, c2);
            Point_3 vp1 = get(vpm, v1);
            Vector_3 vpv(vp1.x(), vp1.y(), vp1.z());
            clusters[c2].add_vertex(vpv, get(vertex_weight_pmap, v1), get(vertex_quadric_pmap, v1));
            clusters[c2].last_modification_iteration = nb_iterations;
            push_vertex_edge_ring_to_queue(v1);
            ++nb_modifications;
          }
          else if (c2 == -1)
          {
            // expand cluster c1 (add v2 to c1)
            put(vertex_cluster_pmap, v2, c1);
            Point_3 vp2 = get(vpm, v2);
            Vector_3 vpv(vp2.x(), vp2.y(), vp2.z());
            clusters[c1].add_vertex(vpv, get(vertex_weight_pmap, v2), get(vertex_quadric_pmap, v2));
            clusters[c1].last_modification_iteration = nb_iterations;
            push_vertex_edge_ring_to_queue(v2);
            ++nb_modifications;
          }
          else
          {
            if ( ( clusters[ c1 ].last_modification_iteration < nb_iterations - 1 ) &&
                 ( clusters[ c2 ].last_modification_iteration < nb_iterations - 1 ) )
            {
              clusters_edges_new.push(hi);
              continue;
            }

            // topological test to avoid creating disconnected clusters
            auto is_topologically_valid_merge = [&](halfedge_descriptor hv, int cluster_id)
            {
              CGAL_assertion(get(vertex_cluster_pmap,target(hv, tmesh)) == cluster_id);

              halfedge_descriptor h = hv;
              bool in_cluster = false;
              int nb_cc_cluster = 0;
              do
              {
                h = next(h, tmesh);

                int ci = get(vertex_cluster_pmap, target(h,tmesh));
                if (in_cluster)
                {
                  if (ci != cluster_id)
                    in_cluster=false;
                }
                else
                {
                  if (ci == cluster_id)
                  {
                    in_cluster = true;
                    if (++nb_cc_cluster > 1)
                      return false;
                  }
                }
                h = opposite(h, tmesh);
              }
              while(h != hv);

              return true;
            };

            // compare the energy of the 3 cases
            Point_3 vp1 = get(vpm, v1);
            Vector_3 vpv1(vp1.x(), vp1.y(), vp1.z());
            Point_3 vp2 = get(vpm, v2);
            Vector_3 vpv2(vp2.x(), vp2.y(), vp2.z());
            FT v1_weight = get(vertex_weight_pmap, v1);
            FT v2_weight = get(vertex_weight_pmap, v2);
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

            FT e_no_change = clusters[c1].compute_energy(qem_energy_minimization) +
                             clusters[c2].compute_energy(qem_energy_minimization);
            FT e_v1_to_c2 = (std::numeric_limits< double >::max)();
            FT e_v2_to_c1 = (std::numeric_limits< double >::max)();

            if ( ( clusters[ c1 ].nb_vertices > 1 ) &&
                 ( !qem_energy_minimization || nb_loops < 2 ||
                   is_topologically_valid_merge(opposite(hi, tmesh), c1) ) )
            {
              cluster1_v1_to_c2.remove_vertex(vpv1, v1_weight, v1_qem);
              cluster2_v1_to_c2.add_vertex(vpv1, v1_weight, v1_qem);
              e_v1_to_c2 = cluster1_v1_to_c2.compute_energy(qem_energy_minimization) +
                           cluster2_v1_to_c2.compute_energy(qem_energy_minimization);
            }

            if ( ( clusters[ c2 ].nb_vertices > 1 ) &&
                 ( !qem_energy_minimization || nb_loops < 2 ||
                   is_topologically_valid_merge(hi, c2) ) )
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
              push_vertex_edge_ring_to_queue(v2);
              ++nb_modifications;
            }
            else if (e_v1_to_c2 < e_no_change) // > 2 as 1 vertex was added to c1
            {
              // move v1 to c2
              put(vertex_cluster_pmap, v1, c2);
              clusters[c1] = cluster1_v1_to_c2;
              clusters[c2] = cluster2_v1_to_c2;
              push_vertex_edge_ring_to_queue(v1);
              ++nb_modifications;
            }
            else
            {
              // no change
              clusters_edges_new.push(hi);
            }
          }
        }
#ifdef CGAL_DEBUG_ACVD
        std::cout << "# Modifications: " << nb_modifications << "\n";
#endif

        clusters_edges_active.swap(clusters_edges_new);
        if (use_qem_based_energy &&
            nb_modifications < nb_vertices * CGAL_TO_QEM_MODIFICATION_THRESHOLD &&
            nb_loops < 2)
        {
          qem_energy_minimization = true;
          break;
        }

      } while (nb_modifications > 0);

      // Disconnected clusters handling
      // the goal is to delete clusters with multiple connected components and only keep the largest connected component of each cluster
      // For each cluster, do a BFS from a vertex in the cluster

      std::vector<std::vector<std::vector<vertex_descriptor>>> cluster_components(nb_clusters);
      std::queue<vertex_descriptor> vertex_queue;

      // loop over vertices to compute cluster components
      VertexVisitedMap vertex_visited_pmap = get(CGAL::dynamic_vertex_property_t<bool>(), tmesh, false);
      for (vertex_descriptor vd : vertices(tmesh))
      {
        if (get(vertex_visited_pmap, vd))
          continue;

        int c = get(vertex_cluster_pmap, vd);
        if (c != -1)
        {
          cluster_components[c].emplace_back();

          vertex_queue.push(vd);
          put(vertex_visited_pmap, vd, true);
          while (!vertex_queue.empty())
          {
            vertex_descriptor v = vertex_queue.front();
            vertex_queue.pop();
            cluster_components[c].back().push_back(v);

            for (halfedge_descriptor hd : halfedges_around_source(v, tmesh))
            {
              vertex_descriptor v2 = target(hd, tmesh);
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
        if (cluster_components[c].size() <= 1)
          continue; // only one component, no need to do anything

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
            for (vertex_descriptor vd : cluster_components[c][component_i])
            {
              put(vertex_cluster_pmap, vd, -1);
              // remove vd from cluster c
              Point_3 vp = get(vpm, vd);
              Vector_3 vpv(vp.x(), vp.y(), vp.z());
              clusters[c].remove_vertex(vpv, get(vertex_weight_pmap, vd), get(vertex_quadric_pmap, vd));
              // add all halfedges around v except hi to the queue
              for (halfedge_descriptor hd : halfedges_around_source(vd, tmesh))
              {
                // add hd to the queue if its target is not in the same cluster
                vertex_descriptor v2 = target(hd, tmesh);
                int c2 = get(vertex_cluster_pmap, v2);
                if (c2 != c)
                  clusters_edges_active.push(hd);
              }
            }
          }
        }
      }

#ifdef CGAL_DEBUG_ACVD
      std::cout << "# nb_disconnected: " << nb_disconnected << "\n";
#endif
      ++nb_loops;

    } while (nb_disconnected > 0 || nb_loops < 3 );

    // Construct new Mesh
    std::vector<Point_3> points;
    points.reserve(clusters.size());
    std::vector<std::array<std::size_t, 3>> polygons;
    polygons.reserve(2*clusters.size()-2); // upper bound

    std::vector<std::size_t> valid_cluster_map(nb_clusters, -1);
    for (int i = 0; i < nb_clusters; ++i)
    {
      CGAL_assertion(clusters[i].weight_sum > 0);
      valid_cluster_map[i] = points.size();
      points.push_back(clusters[i].representative_point(qem_energy_minimization));
    }

    if (use_postprocessing_qem)
    {
      // create a point for each cluster
      std::vector<Eigen::Matrix<FT, 4, 4>> cluster_quadrics(clusters.size());

      // initialize quadrics
      for (int i = 0; i < nb_clusters; ++i)
        cluster_quadrics[i].setZero();

      for (face_descriptor fd : faces(tmesh))
      {
        std::array<Vector_3,3> vecs;
        std::array<int,3> cids;
        int i=0;
        // get Vs for fd
        // compute qem from Vs->"vector_3"s
        // add to the 3 indices of the cluster
        halfedge_descriptor hd = halfedge(fd, tmesh);
        do {
          vertex_descriptor vd = target(hd, tmesh);
          cids[i]=get(vertex_cluster_pmap, vd);
          if (cids[i]==-1)
            break;
          const Point_3& p_i=get(vpm, vd);
          vecs[i++]=p_i - ORIGIN;
          hd = next(hd, tmesh);
        } while (hd != halfedge(fd, tmesh));

        if (i!=3)
          continue;

        Eigen::Matrix<FT, 4, 4> q;
        compute_qem_face(vecs[0], vecs[1], vecs[2], q, gt);
        cluster_quadrics[cids[0]] += q;
        cluster_quadrics[cids[1]] += q;
        cluster_quadrics[cids[2]] += q;
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
    std::vector<halfedge_descriptor> border_hedges;
    extract_boundary_cycles(tmesh, std::back_inserter(border_hedges));

    // loop over boundary loops
    for (halfedge_descriptor hd : border_hedges)
    {
      halfedge_descriptor hd1 = hd;

      int cb_first = -1;

      do
      {
        // 1- get the target and source vertices vt, vs
        // 2- if the target and source vertices are in different clusters, create a new vertex vb between them vb = (vt + vs) / 2
        // 3- make a new face with the new vertex vb and the centers of the clusters of vt and vs
        // 4- also make a new face with vb, the next vb, and the center of the cluster of vt

        vertex_descriptor vt = target(hd1, tmesh);
        vertex_descriptor vs = source(hd1, tmesh);

        int ct = get(vertex_cluster_pmap, vt);
        int cs = get(vertex_cluster_pmap, vs);

        if (ct != cs)
        {
          Point_3 vt_p = get(vpm, vt);
          Point_3 vs_p = get(vpm, vs);
          Vector_3 vt_v(vt_p.x(), vt_p.y(), vt_p.z());
          Vector_3 vs_v(vs_p.x(), vs_p.y(), vs_p.z());

          Vector_3 vb_v = (vt_v + vs_v) / 2;
          Point_3 vb_p(vb_v.x(), vb_v.y(), vb_v.z());

          points.push_back(vb_p);

          std::size_t cb = points.size() - 1;

          if (cb_first == -1)
            cb_first = cb;

          std::size_t ct_mapped = valid_cluster_map[ct], cs_mapped = valid_cluster_map[cs];

          if (ct_mapped != std::size_t(-1) && cs_mapped != std::size_t(-1))
          {
            std::array<std::size_t, 3> polygon = {ct_mapped, cb, cs_mapped};
            polygons.push_back(polygon);

            // after the loop, the last cb+1 should be modified to the first cb
            polygon = {cb, ct_mapped, cb + 1};
            polygons.push_back(polygon);
          }
        }
        hd1 = next(hd1, tmesh);
      } while (hd1 != hd);
      polygons[polygons.size() - 1][2] = cb_first;
    }

    // create a triangle for each face with all vertices in 3 different clusters
    for (face_descriptor fd : faces(tmesh))
    {
      halfedge_descriptor hd1 = halfedge(fd, tmesh);
      vertex_descriptor v1 = source(hd1, tmesh);
      halfedge_descriptor hd2 = next(hd1, tmesh);
      vertex_descriptor v2 = source(hd2, tmesh);
      halfedge_descriptor hd3 = next(hd2, tmesh);
      vertex_descriptor v3 = source(hd3, tmesh);

      int c1 = get(vertex_cluster_pmap, v1);
      int c2 = get(vertex_cluster_pmap, v2);
      int c3 = get(vertex_cluster_pmap, v3);

      if (c1 != c2 && c1 != c3 && c2 != c3)
      {
        std::size_t c1_mapped = valid_cluster_map[c1], c2_mapped = valid_cluster_map[c2], c3_mapped = valid_cluster_map[c3];
        if (c1_mapped !=  std::size_t(-1) && c2_mapped !=  std::size_t(-1) && c3_mapped !=  std::size_t(-1))
        {
          std::array<std::size_t,3> polygon = {c1_mapped, c2_mapped, c3_mapped};
          polygons.push_back(polygon);
        }
      }
    }

#ifdef CGAL_DEBUG_ACVD
    static int kkk = -1;
    CGAL::IO::write_polygon_soup("/tmp/soup_" + std::to_string(++kkk) + ".off", points, polygons);
    dump_mesh_with_cluster_colors(tmesh, vertex_cluster_pmap, "/tmp/cluster_" + std::to_string(kkk) + ".ply");
#endif

    // detect non-manifold edges in the output
    std::vector<std::unordered_map<std::size_t, std::size_t> > edge_map(points.size());
    for (const std::array<std::size_t, 3> & p : polygons)
    {
      for (int i=0; i<3; ++i)
      {
        std::pair<std::size_t, std::size_t> edge = make_sorted_pair(p[i], p[(i+1)%3]);
        edge_map[edge.first].emplace(edge.second,0).first->second+=1;
      }
    }

    std::vector<std::size_t> nm_clusters;
    for (std::size_t i=0; i<points.size(); ++i)
    {
      for ( auto [j, size] : edge_map[i] )
      {
        if (size > 2)
        {
#ifdef CGAL_DEBUG_ACVD
          std::cout << "non-manifold edge : " << i << " " << j << "\n";
#endif
          nm_clusters.push_back(i);
          nm_clusters.push_back(j);
        }
      }
    }


    // detect isolated graph edges
    for (edge_descriptor ed : edges(tmesh))
    {
      int c1 = get(vertex_cluster_pmap, source(ed, tmesh));
      int c2 = get(vertex_cluster_pmap, target(ed, tmesh));
#ifdef CGAL_DEBUG_ACVD
      if (c1 == -1 || c2 == -1)
        throw std::runtime_error("non assigned vertex");
#endif

      if (c1 == c2)
        continue;
      if (c2 < c1)
        std::swap(c1, c2);

      if (edge_map[c1].emplace(c2,0).second)
      {
#ifdef CGAL_DEBUG_ACVD
        std::cout << "isolated edge " << c1 << " " << c2 << "\n";
        std::cout << "   " << kkk << " " <<  points[c1] << " " << points[c2] << "\n";
#endif
        nm_clusters.push_back(c1);
        nm_clusters.push_back(c2);
      }
    }

    std::sort(nm_clusters.begin(), nm_clusters.end());
    nm_clusters.erase(std::unique(nm_clusters.begin(), nm_clusters.end()), nm_clusters.end());

    // detect non-manifold vertices
    std::vector< std::vector< std::pair<std::size_t, std::size_t> > > link_edges(points.size());
    for (const std::array<std::size_t, 3> & p : polygons)
    {
      for (int i=0; i<3; ++i)
        link_edges[p[i]].emplace_back(p[(i+1)%3], p[(i+2)%3]);
    }

    using Graph = boost::adjacency_list <boost::vecS, boost::vecS, boost::undirectedS>;
    for (std::size_t i=0; i< points.size(); ++i)
    {
      if (std::binary_search(nm_clusters.begin(), nm_clusters.end(), i)) continue;
      std::vector<Graph::vertex_descriptor> descriptors(points.size(), Graph::null_vertex());

      Graph graph;
      for (const auto& p : link_edges[i])
      {
        if (descriptors[p.first] == Graph::null_vertex())
          descriptors[p.first] = add_vertex(graph);
        if (descriptors[p.second] == Graph::null_vertex())
          descriptors[p.second] = add_vertex(graph);
        add_edge(descriptors[p.first], descriptors[p.second], graph);
      }

      // TODO: (possible optimisation) add an id to vertices in Graph and use a random acccess pmap
      std::map<typename Graph::vertex_descriptor, int> the_map;

      if (boost::connected_components(graph, boost::make_assoc_property_map(the_map)) > 1)
      {
#ifdef CGAL_DEBUG_ACVD
        std::cout << "non-manifold vertex " << i << "\n";
#endif
        nm_clusters.push_back(i);
      }
    }

    if (nm_clusters.empty())
    {
      // dump_mesh_with_cluster_colors(tmesh, vertex_cluster_pmap, "/tmp/debug.ply");
      // CGAL::IO::write_polygon_soup("/tmp/soup.off", points, polygons);

      return std::make_pair(points, polygons);
    }


    std::vector<std::size_t> one_ring;
    for (std::size_t nmi : nm_clusters)
    {
      for ( auto [n, s] : edge_map[nmi] )
        one_ring.push_back(n);
    }

    std::set<std::size_t> nm_clusters_picked; // TODO: (possible optimization) use cluster data instead of the set

    for (vertex_descriptor v : vertices(tmesh))
    {
      int c = get(vertex_cluster_pmap, v);
      if (std::binary_search(nm_clusters.begin(), nm_clusters.end(), c))
      {
        if (!nm_clusters_picked.insert(c).second)
          continue;

        if (clusters[c].nb_vertices == 1)
          continue;

        put(vertex_cluster_pmap, v, clusters.size());
        CGAL_assertion(get(vertex_cluster_pmap, v) == (int) clusters.size());
        clusters.emplace_back();
      }

      if (nm_clusters_picked.size() == nm_clusters.size())
        break;
    }

    frozen_clusters = std::vector<bool>(nb_clusters, true);
    for (std::size_t nmi : nm_clusters)
      frozen_clusters[nmi] = false;
    for (std::size_t nmi : one_ring)
      frozen_clusters[nmi] = false;

    nb_clusters = clusters.size();
    frozen_clusters.resize(nb_clusters, false);
    nb_loops = 0;
    qem_energy_minimization = false;

  } while(true);
}


} // namespace internal

/**
* \ingroup PkgPolygonMeshProcessingRef
*
* performs Approximated Centroidal Voronoi Diagram (ACVD) remeshing on a triangle mesh. The remeshing is either uniform or adaptative.
* Note that there is no guarantee that the output mesh will have the same topology as the input.
* See \ref acvdrem for more details on the algorithm.
*
* @tparam TriangleMesh a model of `FaceListGraph` and `MutableFaceGraph`
* @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters".
*
* @param tmesh triangle mesh to be remeshed
* @param nb_vertices lower bound on the target number of vertices in the output mesh.
*                    The requested number of vertices in the output will be respected except if the input mesh is not closed
*                    (extra vertices will be used on the boundary), or if the number of points is too low
*                    and no manifold mesh can be produced with that budget of points (extra points are added to get a manifold output).
* @param np optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below.
*        `GT` stands for the type of the traits object provided to the named parameter `geom_traits()`.
*
* \cgalNamedParamsBegin
*
*   \cgalParamNBegin{vertex_principal_curvatures_and_directions_map}
*     \cgalParamDescription{a property map associating principal curvatures and directions to the vertices of `tmesh`, used for adaptive clustering.}
*     \cgalParamType{a class model of `ReadWritePropertyMap` with
*                    `boost::graph_traits<TriangleMesh>::%vertex_descriptor`
*                    as key type and `Principal_curvatures_and_directions<GT>` as value type.}
*     \cgalParamDefault{If this parameter is omitted, but `gradation_factor` is provided, an internal property map
*                       will be created and curvature values will be computed using the function `interpolated_corrected_curvatures()`.}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{gradation_factor}
*     \cgalParamDescription{a factor used to gradate the weights of the vertices based on their curvature values.
*                           The larger the value is, the more the curvature will impact the distribution of output vertices
*                           on the surface. The original paper recommends the value `1.5`.}
*     \cgalParamType{`GT::FT`}
*     \cgalParamDefault{0}
*     \cgalParamExtra{If this parameter is omitted, no adaptive clustering will be performed.}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{use_postprocessing_qem}
*     \cgalParamDescription{indicates if a projection step using quadric error metric should be applied to cluster centers at the end of the minimization,
*                           in order for example to recover sharp features.
*                           This is a fast method but can result in badly shaped triangles and even self-intersections.}
*     \cgalParamType{Boolean}
*     \cgalParamDefault{false}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{use_qem_based_energy}
*     \cgalParamDescription{indicates if quadric error metric should be applied during the minimization algorithm in order for example to recover sharp features.
*                           This is slower than using `use_postprocessing_qem(true)`, but it is more accurate.}
*     \cgalParamType{Boolean}
*     \cgalParamDefault{false}
*     \cgalParamExtra{If this parameter is `true` then `use_postprocessing_qem` will be automatically set to `false`.}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{vertex_count_ratio}
*     \cgalParamDescription{a ratio used to control the subdivision of the input mesh in case it does not have enough vertices compared to `nb_vertices`.
*                           More precisely, the number of vertices of the input mesh should at least be the ratio times `nb_vertices`.
*                           If not, the mesh will first be subdivided until the aforementioned criterium is met.}
*     \cgalParamType{`GT::FT`}
*     \cgalParamDefault{0.1}
*     \cgalParamExtra{A value between 0.1 and 0.01 is recommended, the smaller the better the approximation will be, but it will increase the runtime.}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{vertex_point_map}
*       \cgalParamDescription{a property map associating points to the vertices of `tmesh`.}
*       \cgalParamType{a class model of `ReadWritePropertyMap` with
*                      `boost::graph_traits<TriangleMesh>::%vertex_descriptor`
*                      as key type and `GT::Point_3` as value type.}
*       \cgalParamDefault{`boost::get(CGAL::vertex_point, tmesh)`.}
*       \cgalParamExtra{If this parameter is omitted, an internal property map for
*                       `CGAL::vertex_point_t` must be available in `TriangleMesh`.}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{geom_traits}
*      \cgalParamDescription{an instance of a geometric traits class.}
*      \cgalParamType{a class model of `Kernel`, with `GT::FT` being either `float` or `double`.}
*      \cgalParamDefault{a \cgal Kernel deduced from the vertex point type, using `CGAL::Kernel_traits`.}
*      \cgalParamExtra{The geometric traits class must be compatible with the vertex point type.}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{random_seed}
*     \cgalParamDescription{a value to seed the random number generator}
*     \cgalParamType{unsigned int}
*     \cgalParamDefault{a value generated with `std::time()`}
*   \cgalParamNEnd
*
* \cgalNamedParamsEnd
*
* @return `true` if `nb_vertices` was sufficiently large for remeshing the input, and `false` if more points were used
*
*/
template <typename TriangleMesh,
          typename NamedParameters = parameters::Default_named_parameters>
bool approximated_centroidal_Voronoi_diagram_remeshing(TriangleMesh& tmesh,
                                                       std::size_t nb_vertices,
                                                       const NamedParameters& np = parameters::default_values())
{
  auto ps = internal::acvd_impl(tmesh, nb_vertices, np);
  CGAL_assertion(is_polygon_soup_a_polygon_mesh(ps.second));

  auto vpm = parameters::choose_parameter(
                parameters::get_parameter(np, CGAL::vertex_point),
                get_property_map(CGAL::vertex_point, tmesh));

  remove_all_elements(tmesh);
  polygon_soup_to_polygon_mesh(ps.first, ps.second, tmesh, parameters::vertex_point_map(vpm));
  return ps.first.size() == nb_vertices;
}

} // namespace Polygon_mesh_processing
} // namespace CGAL

#undef CGAL_WEIGHT_CLAMP_RATIO_THRESHOLD
#undef CGAL_TO_QEM_MODIFICATION_THRESHOLD

#endif // CGAL_PMP_ACVD_REMESHING_H
