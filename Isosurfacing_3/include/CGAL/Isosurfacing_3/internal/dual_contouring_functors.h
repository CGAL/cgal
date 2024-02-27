// Copyright (c) 2022-2024 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Daniel Zint
//                 Julian Stahl
//                 Mael Rouxel-Labbé

#ifndef CGAL_ISOSURFACING_3_INTERNAL_DUAL_CONTOURING_FUNCTORS_H
#define CGAL_ISOSURFACING_3_INTERNAL_DUAL_CONTOURING_FUNCTORS_H

#include <CGAL/license/Isosurfacing_3.h>

#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/Named_function_parameters.h>
#include <CGAL/Container_helper.h>
#include <CGAL/Origin.h>

#ifdef CGAL_EIGEN3_ENABLED
#include <CGAL/Eigen_matrix.h>
#include <CGAL/Eigen_vector.h>

#include <Eigen/SVD>
#endif

#ifdef CGAL_LINKED_WITH_TBB
#if TBB_INTERFACE_VERSION < 12010 && !defined(TBB_PREVIEW_CONCURRENT_ORDERED_CONTAINERS)
#define CGAL_HAS_DEFINED_TBB_PREVIEW_CONCURRENT_ORDERED_CONTAINERS
#define TBB_PREVIEW_CONCURRENT_ORDERED_CONTAINERS 1
#endif
#include <tbb/concurrent_map.h>
#endif
#include <iostream>
#include <map>
#include <mutex>
#include <vector>

namespace CGAL {
namespace Isosurfacing {
namespace internal {

template <typename Domain,
          typename EdgeToPointIDMap,
          typename PointRange,
          typename GradientRange>
bool cell_position_QEM(const typename Domain::Cell_descriptor& c,
                       const Domain& domain,
                       const bool constrain_to_cell,
                       const EdgeToPointIDMap& edge_to_point_id,
                       const PointRange& edge_points,
                       const GradientRange& edge_gradients,
                       typename Domain::Geom_traits::Point_3& p)
{
  using Geom_traits = typename Domain::Geom_traits;
  using FT = typename Geom_traits::FT;
  using Point_3 = typename Geom_traits::Point_3;
  using Vector_3 = typename Geom_traits::Vector_3;

  using Eigen_vector_3 = Eigen_vector<FT, 3>;
  using Eigen_matrix_3 = Eigen_matrix<FT, 3, 3>;
  using Eigen_vector_x = Eigen_vector<FT>;
  using Eigen_matrix_x = Eigen_matrix<FT>;

  typename Geom_traits::Compute_x_3 x_coord = domain.geom_traits().compute_x_3_object();
  typename Geom_traits::Compute_y_3 y_coord = domain.geom_traits().compute_y_3_object();
  typename Geom_traits::Compute_z_3 z_coord = domain.geom_traits().compute_z_3_object();
  typename Geom_traits::Construct_point_3 point = domain.geom_traits().construct_point_3_object();

  // compute edge intersections
  std::vector<Point_3> cell_edge_intersections;
  std::vector<Vector_3> cell_edge_intersection_normals;

  for(const auto& edge : domain.cell_edges(c))
  {
    const auto it = edge_to_point_id.find(edge);
    if(it == edge_to_point_id.end())
      continue;

    cell_edge_intersections.push_back(edge_points[it->second]);
    cell_edge_intersection_normals.push_back(edge_gradients[it->second]);
  }

  const std::size_t en = cell_edge_intersections.size();
  if(en == 0)
    return false;

#ifdef CGAL_ISOSURFACING_3_DC_FUNCTORS_DEBUG
  std::cout << "Points and normals: " << std::endl;
  for(std::size_t i=0; i<cell_edge_intersections.size(); ++i)
  {
    std::cout << x_coord(cell_edge_intersections[i]) << " "
              << y_coord(cell_edge_intersections[i]) << " "
              << z_coord(cell_edge_intersections[i]) << " "
              << x_coord(cell_edge_intersections[i]) + x_coord(cell_edge_intersection_normals[i]) << " "
              << y_coord(cell_edge_intersections[i]) + y_coord(cell_edge_intersection_normals[i]) << " "
              << z_coord(cell_edge_intersections[i]) + z_coord(cell_edge_intersection_normals[i]) << std::endl;
  }
#endif

  FT x_min, y_min, z_min, x_max, y_max, z_max;
  x_min = y_min = z_min = std::numeric_limits<FT>::max(); // @todo domain.span()
  x_max = y_max = z_max = std::numeric_limits<FT>::lowest();
  FT x(0), y(0), z(0);

  typename Domain::Cell_vertices vertices = domain.cell_vertices(c);
  for(const auto& v : vertices)
  {
    const Point_3& cp = domain.point(v);
    if(constrain_to_cell)
    {
      x_min = (std::min<FT>)(x_min, x_coord(cp));
      y_min = (std::min<FT>)(y_min, y_coord(cp));
      z_min = (std::min<FT>)(z_min, z_coord(cp));

      x_max = (std::max<FT>)(x_max, x_coord(cp));
      y_max = (std::max<FT>)(y_max, y_coord(cp));
      z_max = (std::max<FT>)(z_max, z_coord(cp));
    }
  }

  for(const auto& ep : cell_edge_intersections)
  {
    x += x_coord(ep);
    y += y_coord(ep);
    z += z_coord(ep);
  }

  Point_3 com = point(x / en, y / en, z / en);

#ifdef CGAL_ISOSURFACING_3_DC_FUNCTORS_DEBUG
  std::cout << "cell: " << x_min << " " << y_min << " " << z_min << " " << x_max << " " << y_max << " " << z_max << std::endl;
  std::cout << "COM " << com << std::endl;
#endif

  // SVD QEM
  Eigen_matrix_3 A;
  A.setZero();
  Eigen_vector_3 rhs;
  rhs.setZero();
  for(std::size_t i=0; i<cell_edge_intersections.size(); ++i)
  {
    Eigen_vector_3 n_k;
    n_k.set(0, x_coord(cell_edge_intersection_normals[i]));
    n_k.set(1, y_coord(cell_edge_intersection_normals[i]));
    n_k.set(2, z_coord(cell_edge_intersection_normals[i]));

    Eigen_vector_3 p_k;
    p_k.set(0, x_coord(cell_edge_intersections[i]));
    p_k.set(1, y_coord(cell_edge_intersections[i]));
    p_k.set(2, z_coord(cell_edge_intersections[i]));

    const FT d_k = n_k.transpose() * p_k;

    // have to cast to cast to the underlying Eigen type because the type
    // of 'n_k * n_k.transpose()' is Eigen::Product< ... > etc.,
    // so the double conversion is not implicit
    Eigen_matrix_3 A_k = typename Eigen_matrix_3::EigenType(n_k * n_k.transpose());
    Eigen_vector_3 b_k;
    b_k = d_k * n_k;
    A += A_k;
    rhs += b_k;
  }

  Eigen::JacobiSVD<typename Eigen_matrix_x::EigenType> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);

  // Lindstrom's paper, "Out-of-Core Simplification of Large Polygonal Models": 1e-3
  // Ju's paper, "Dual Contouring of Hermite Data": 1e-1
  svd.setThreshold(1e-1);

  Eigen_vector_3 x_hat;
  x_hat << x_coord(com), y_coord(com), z_coord(com);

  // Lindstrom formula for QEM new position for singular matrices
  Eigen_vector_x v_svd;
  v_svd = x_hat + svd.solve(rhs - A * x_hat);


  // bounding box
  if(constrain_to_cell)
  {
    v_svd[0] = std::clamp<FT>(v_svd[0], x_min, x_max);
    v_svd[1] = std::clamp<FT>(v_svd[1], y_min, y_max);
    v_svd[2] = std::clamp<FT>(v_svd[2], z_min, z_max);
  }

  p = point(v_svd[0], v_svd[1], v_svd[2]);

#ifdef CGAL_ISOSURFACING_3_DC_FUNCTORS_DEBUG
  std::cout << "CGAL QEM POINT: " << v_svd[0] << " " << v_svd[1] << " " << v_svd[2] << std::endl;
  std::cout << "CGAL clamped QEM POINT: " << p[0] << " " << p[1] << " " << p[2] << std::endl;
  std::cout << "--" << std::endl;
#endif

  return true;
}

template <typename Domain,
          typename EdgeToPointIDMap,
          typename CellToPointIDMap,
          typename PolygonRange>
void generate_face(const typename Domain::Edge_descriptor& e,
                   const Domain& domain,
                   const typename Domain::Geom_traits::FT isovalue,
                   const bool do_not_triangulate_faces,
                   const EdgeToPointIDMap& edge_to_point_id,
                   const CellToPointIDMap& cell_to_point_id,
                   std::mutex& mutex,
                   PolygonRange& polygons)
{
  using FT = typename Domain::Geom_traits::FT;

  using Cell_descriptor = typename Domain::Cell_descriptor;

  const auto& vertices = domain.incident_vertices(e);

  // @todo this check could be avoided for QEM: active edges are in `edge_to_point_id`
  const FT val_0 = domain.value(vertices[0]);
  const FT val_1 = domain.value(vertices[1]);
  if((val_0 <= isovalue) == (val_1 <= isovalue))
    return;

  std::vector<std::size_t> vertex_ids;

  const auto& icells = domain.incident_cells(e);
  for(const Cell_descriptor& c : icells)
  {
    auto it = cell_to_point_id.find(c);
    if(it == cell_to_point_id.end())
      continue;

    vertex_ids.push_back(it->second);
  }

  if(vertex_ids.size() < 3)
    return;

  if(val_0 > val_1)
    std::reverse(vertex_ids.begin(), vertex_ids.end());

  // @todo? filter degenerate faces?

  if(do_not_triangulate_faces)
  {
    std::lock_guard<std::mutex> lock(mutex);
    polygons.emplace_back();
    CGAL::internal::resize(polygons.back(), vertex_ids.size());
    std::copy(vertex_ids.begin(), vertex_ids.end(), std::begin(polygons.back()));
  }
  else
  {
    auto it = edge_to_point_id.find(e);
    if(it == edge_to_point_id.end())
    {
      CGAL_assertion(false);
      return;
    }
    const std::size_t ei = it->second;

    std::lock_guard<std::mutex> lock(mutex);
    for(std::size_t i=0; i<vertex_ids.size(); ++i)
      polygons.push_back({ei, vertex_ids[i], vertex_ids[(i + 1) % vertex_ids.size()]});
  }
}

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

namespace DC_Strategy {

struct QEM {};
struct Centroid_of_edge_intersections {};
struct Cell_center {};

} // namespace DC_Strategy

template <typename ConcurrencyTag,
          typename Domain,
          typename Dual_contouring_strategy>
class Dual_contourer;

template <typename ConcurrencyTag,
          typename Domain>
class Dual_contourer<ConcurrencyTag, Domain, DC_Strategy::QEM>
{
  using Geom_traits = typename Domain::Geom_traits;
  using FT = typename Geom_traits::FT;
  using Point_3 = typename Geom_traits::Point_3;
  using Vector_3 = typename Geom_traits::Vector_3;

  using Vertex_descriptor = typename Domain::Vertex_descriptor;
  using Edge_descriptor = typename Domain::Edge_descriptor;
  using Cell_descriptor = typename Domain::Cell_descriptor;

  std::mutex m_mutex;

public:
  template<typename PointRange, typename PolygonRange, typename NamedParameters>
  void operator()(const Domain& domain,
                  const FT isovalue,
                  PointRange& points,
                  PolygonRange& polygons,
                  const NamedParameters& np = parameters::default_values())
  {
    using parameters::choose_parameter;
    using parameters::get_parameter;

    // Otherwise the `edge_to_point_id` map might be messed up
    CGAL_precondition(points.empty());
    CGAL_precondition(polygons.empty());

    const bool constrain_to_cell = choose_parameter(get_parameter(np, internal_np::constrain_to_cell), false);

    const bool do_not_triangulate_faces =
      choose_parameter(get_parameter(np, internal_np::do_not_triangulate_faces), false);

    using Edge_to_point_ID_map = std::unordered_map<Edge_descriptor, std::size_t>;
    using Cell_to_point_ID_map = std::unordered_map<Cell_descriptor, std::size_t>;

    Edge_to_point_ID_map edge_to_point_id;
    Cell_to_point_ID_map cell_to_point_id;

    std::vector<Point_3> edge_points;
    std::vector<Vector_3> edge_gradients;

    // ---------------------------------------------------------------------------------------------
    // construct the intersection of the surface at active edges
    auto edge_positioner = [&](const Edge_descriptor& e)
    {
      const auto& evs = domain.incident_vertices(e);
      const Vertex_descriptor& v0 = evs[0];
      const Vertex_descriptor& v1 = evs[1];
      const Point_3& p0 = domain.point(v0);
      const Point_3& p1 = domain.point(v1);
      const FT val0 = domain.value(v0);
      const FT val1 = domain.value(v1);

      Point_3 p;
      bool res = domain.construct_intersection(p0, p1, val0, val1, isovalue, p);
      if(!res)
        return;

      const Vector_3 g = domain.gradient(p);

      std::lock_guard<std::mutex> lock(m_mutex);
      edge_to_point_id[e] = edge_points.size();
      edge_points.push_back(p);
      edge_gradients.push_back(g);
    };
    domain.template for_each_edge<ConcurrencyTag>(edge_positioner);

#ifdef CGAL_ISOSURFACING_3_DC_FUNCTORS_DEBUG
    typename Geom_traits::Compute_x_3 x_coord = domain.geom_traits().compute_x_3_object();
    typename Geom_traits::Compute_y_3 y_coord = domain.geom_traits().compute_y_3_object();
    typename Geom_traits::Compute_z_3 z_coord = domain.geom_traits().compute_z_3_object();

    std::ofstream out_active_edges("active_edges.polylines");
    for(const auto& ei : edge_to_point_id)
    {
      const Edge_descriptor& e = ei.first;
      const auto& evs = domain.incident_vertices(e);
      const Vertex_descriptor& v0 = evs[0];
      const Vertex_descriptor& v1 = evs[1];
      const Point_3& p0 = domain.point(v0);
      const Point_3& p1 = domain.point(v1);

      out_active_edges << "2 " << x_coord(p0) << " " << y_coord(p0) << " " << z_coord(p0) << " " << x_coord(p1) << " " << y_coord(p1) << " " << z_coord(p1) << std::endl;
    }

    std::ofstream out_edge_intersections("edge_intersections.polylines");
    for(const auto& ei : edge_to_point_id)
    {
      const Point_3& p = edge_points.at(ei.second);
      const Vector_3& g = edge_gradients.at(ei.second);

      out_edge_intersections << "2 " << x_coord(p) << " " << y_coord(p) << " " << z_coord(p) << " " << x_coord(p) + x_coord(g) << " " << y_coord(p) + y_coord(g) << " " << z_coord(p) + z_coord(g) << std::endl;
    }
#endif

    if(!do_not_triangulate_faces)
      points.insert(points.end(), edge_points.begin(), edge_points.end());

    // ---------------------------------------------------------------------------------------------
    // create a vertex for each cell that has at least one active edge
    auto cell_positioner = [&](const Cell_descriptor& c)
    {
      Point_3 p;
      if(cell_position_QEM(c, domain, constrain_to_cell, edge_to_point_id,
                           edge_points, edge_gradients, p))
      {
        std::lock_guard<std::mutex> lock(m_mutex); // @todo useless if sequential
        cell_to_point_id[c] = points.size();
        points.push_back(p);
      }
    };
    domain.template for_each_cell<ConcurrencyTag>(cell_positioner);

    // ---------------------------------------------------------------------------------------------
    // connect vertices around edges to form faces
    auto face_generator = [&](const Edge_descriptor& e)
    {
      generate_face(e, domain, isovalue, do_not_triangulate_faces,
                    edge_to_point_id, cell_to_point_id, m_mutex, polygons);
    };
    domain.template for_each_edge<ConcurrencyTag>(face_generator);

#ifdef CGAL_ISOSURFACING_3_DC_FUNCTORS_DEBUG
    std::cout << points.size() << " points" << std::endl;
    std::cout << polygons.size() << " polygons" << std::endl;
#endif
  }
};

template <typename ConcurrencyTag,
          typename Domain>
class Dual_contourer<ConcurrencyTag, Domain, DC_Strategy::Centroid_of_edge_intersections>
{
  using Geom_traits = typename Domain::Geom_traits;
  using FT = typename Geom_traits::FT;
  using Point_3 = typename Geom_traits::Point_3;
  using Vector_3 = typename Geom_traits::Vector_3;

  using Vertex_descriptor = typename Domain::Vertex_descriptor;
  using Edge_descriptor = typename Domain::Edge_descriptor;
  using Cell_descriptor = typename Domain::Cell_descriptor;

  std::mutex m_mutex;

public:
  template<typename PointRange, typename PolygonRange, typename NamedParameters>
  void operator()(const Domain& domain,
                  const typename Geom_traits::FT isovalue,
                  PointRange& points,
                  PolygonRange& polygons,
                  const NamedParameters& np = parameters::default_values())
  {
    using parameters::choose_parameter;
    using parameters::get_parameter;

    // Otherwise the `edge_to_point_id` map might be messed up
    CGAL_precondition(points.empty());
    CGAL_precondition(polygons.empty());

    bool do_not_triangulate_faces =
      choose_parameter(get_parameter(np, internal_np::do_not_triangulate_faces), false);

    using Edge_to_point_ID_map = std::unordered_map<Edge_descriptor, std::size_t>;
    using Cell_to_point_ID_map = std::unordered_map<Cell_descriptor, std::size_t>;

    Edge_to_point_ID_map edge_to_point_id;
    Cell_to_point_ID_map cell_to_point_id;

    std::vector<Point_3> edge_points;

    // ---------------------------------------------------------------------------------------------
    auto edge_positioner = [&](const Edge_descriptor& e)
    {
      const auto& evs = domain.incident_vertices(e);
      const Vertex_descriptor& v0 = evs[0];
      const Vertex_descriptor& v1 = evs[1];
      const Point_3& p0 = domain.point(v0);
      const Point_3& p1 = domain.point(v1);
      const FT val0 = domain.value(v0);
      const FT val1 = domain.value(v1);

      Point_3 p;
      bool res = domain.construct_intersection(p0, p1, val0, val1, isovalue, p);
      if(!res)
        return;

      std::lock_guard<std::mutex> lock(m_mutex);
      edge_to_point_id[e] = edge_points.size();
      edge_points.push_back(p);
    };
    domain.template for_each_edge<ConcurrencyTag>(edge_positioner);

    if(!do_not_triangulate_faces)
      points.insert(points.end(), edge_points.begin(), edge_points.end());

    // ---------------------------------------------------------------------------------------------
    auto cell_positioner = [&](const Cell_descriptor& c)
    {
      typename Geom_traits::Compute_x_3 x_coord = domain.geom_traits().compute_x_3_object();
      typename Geom_traits::Compute_y_3 y_coord = domain.geom_traits().compute_y_3_object();
      typename Geom_traits::Compute_z_3 z_coord = domain.geom_traits().compute_z_3_object();
      typename Geom_traits::Construct_point_3 point = domain.geom_traits().construct_point_3_object();

      // compute edge intersections
      std::vector<Point_3> edge_intersections;
      for(const Edge_descriptor& e : domain.cell_edges(c))
      {
        const auto it = edge_to_point_id.find(e);
        if(it == edge_to_point_id.end())
          continue;

        edge_intersections.push_back(edge_points[it->second]); // @todo could avoid copying
      }

      const std::size_t en = edge_intersections.size();
      if(en == 0)
        return;

      FT x = 0, y = 0, z = 0;
      for(const Point_3& p : edge_intersections)
      {
        x += x_coord(p);
        y += y_coord(p);
        z += z_coord(p);
      }

      const Point_3 p = point(x / en, y / en, z / en);

      std::lock_guard<std::mutex> lock(m_mutex);
      cell_to_point_id[c] = points.size();
      points.push_back(p);
    };
    domain.template for_each_cell<ConcurrencyTag>(cell_positioner);

    // ---------------------------------------------------------------------------------------------
    auto face_generator = [&](const Edge_descriptor& e)
    {
      generate_face(e, domain, isovalue, do_not_triangulate_faces,
                    edge_to_point_id, cell_to_point_id, m_mutex, polygons);
    };
    domain.template for_each_edge<ConcurrencyTag>(face_generator);
  }
};

template <typename ConcurrencyTag,
          typename Domain>
class Dual_contourer<ConcurrencyTag, Domain, DC_Strategy::Cell_center>
{
  using Geom_traits = typename Domain::Geom_traits;
  using FT = typename Geom_traits::FT;
  using Point_3 = typename Geom_traits::Point_3;
  using Vector_3 = typename Geom_traits::Vector_3;

  using Vertex_descriptor = typename Domain::Vertex_descriptor;
  using Edge_descriptor = typename Domain::Edge_descriptor;
  using Cell_descriptor = typename Domain::Cell_descriptor;

  std::mutex m_mutex;

public:
  template<typename PointRange, typename PolygonRange, typename NamedParameters>
  void operator()(const Domain& domain,
                  const FT isovalue,
                  PointRange& points,
                  PolygonRange& polygons,
                  const NamedParameters& np = parameters::default_values())
  {
    using parameters::choose_parameter;
    using parameters::get_parameter;

    // Otherwise the `edge_to_point_id` map might be messed up
    CGAL_precondition(points.empty());
    CGAL_precondition(polygons.empty());

    bool do_not_triangulate_faces =
      choose_parameter(get_parameter(np, internal_np::do_not_triangulate_faces), false);

    using Edge_to_point_ID_map = std::unordered_map<Edge_descriptor, std::size_t>;
    using Cell_to_point_ID_map = std::unordered_map<Cell_descriptor, std::size_t>;

    Edge_to_point_ID_map edge_to_point_id;
    Cell_to_point_ID_map cell_to_point_id;

    std::vector<Point_3> edge_points;

    // ---------------------------------------------------------------------------------------------
    auto edge_positioner = [&](const Edge_descriptor& e)
    {
      typename Geom_traits::Compute_x_3 x_coord = domain.geom_traits().compute_x_3_object();
      typename Geom_traits::Compute_y_3 y_coord = domain.geom_traits().compute_y_3_object();
      typename Geom_traits::Compute_z_3 z_coord = domain.geom_traits().compute_z_3_object();
      typename Geom_traits::Construct_point_3 point = domain.geom_traits().construct_point_3_object();

      const auto& evs = domain.incident_vertices(e);
      const Vertex_descriptor& v0 = evs[0];
      const Vertex_descriptor& v1 = evs[1];

      const FT val_0 = domain.value(v0);
      const FT val_1 = domain.value(v1);
      if((val_0 <= isovalue) == (val_1 <= isovalue))
        return;

      Point_3 p = point((x_coord(domain.point(v0)) + x_coord(domain.point(v1))) / FT(2),
                        (y_coord(domain.point(v0)) + y_coord(domain.point(v1))) / FT(2),
                        (z_coord(domain.point(v0)) + z_coord(domain.point(v1))) / FT(2));

      std::lock_guard<std::mutex> lock(m_mutex);
      edge_to_point_id[e] = edge_points.size();
      edge_points.push_back(p);
    };
    domain.template for_each_edge<ConcurrencyTag>(edge_positioner);

    if(!do_not_triangulate_faces)
      points.insert(points.end(), edge_points.begin(), edge_points.end());

    // ---------------------------------------------------------------------------------------------
    auto cell_positioner = [&](const Cell_descriptor& c)
    {
      typename Geom_traits::Compute_x_3 x_coord = domain.geom_traits().compute_x_3_object();
      typename Geom_traits::Compute_y_3 y_coord = domain.geom_traits().compute_y_3_object();
      typename Geom_traits::Compute_z_3 z_coord = domain.geom_traits().compute_z_3_object();
      typename Geom_traits::Construct_point_3 point = domain.geom_traits().construct_point_3_object();

      typename Domain::Cell_vertices vertices = domain.cell_vertices(c);
      const std::size_t cn = vertices.size();

      bool all_smaller = true;
      bool all_greater = true;
      for(const auto& v : vertices)
      {
        const bool b = (domain.value(v) <= isovalue);
        all_smaller = all_smaller && b;
        all_greater = all_greater && !b;
      }

      if(all_smaller || all_greater)
        return;

      FT x(0), y(0), z(0);
      for(const auto& v : vertices)
      {
        const Point_3& cp = domain.point(v);
        x += x_coord(cp);
        y += y_coord(cp);
        z += z_coord(cp);
      }

      // set point to cell center
      Point_3 p = point(x / cn, y / cn, z / cn);

      std::lock_guard<std::mutex> lock(m_mutex);
      cell_to_point_id[c] = points.size();
      points.push_back(p);
    };
    domain.template for_each_cell<ConcurrencyTag>(cell_positioner);

    // ---------------------------------------------------------------------------------------------
    auto face_generator = [&](const Edge_descriptor& e)
    {
      generate_face(e, domain, isovalue, do_not_triangulate_faces,
                    edge_to_point_id, cell_to_point_id, m_mutex, polygons);
    };
    domain.template for_each_edge<ConcurrencyTag>(face_generator);
  }
};

} // namespace internal
} // namespace Isosurfacing
} // namespace CGAL

#endif // CGAL_ISOSURFACING_3_INTERNAL_DUAL_CONTOURING_FUNCTORS_H
