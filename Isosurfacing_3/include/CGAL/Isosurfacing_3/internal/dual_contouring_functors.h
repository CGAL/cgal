// Copyright (c) 2022-2023 INRIA Sophia-Antipolis (France).
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

#ifndef CGAL_ISOSURFACING_3_INTERNAL_DUAL_CONTOURING_FUNCTORS_H
#define CGAL_ISOSURFACING_3_INTERNAL_DUAL_CONTOURING_FUNCTORS_H

#include <CGAL/license/Isosurfacing_3.h>

#include <CGAL/centroid.h>
#include <CGAL/Origin.h>

#ifdef CGAL_EIGEN3_ENABLED
#include <CGAL/Eigen_matrix.h>
#include <CGAL/Eigen_vector.h>

#include <Eigen/SVD>
#endif

#include <map>
#include <mutex>
#include <vector>

namespace CGAL {
namespace Isosurfacing {
namespace internal {
namespace Positioning {

/*
 * \ingroup PkgIsosurfacing3Ref
 *
 * \brief computes the vertex position for a point in Dual Contouring
 *        using Quadric Error Metrics and the SVD pseudo inverse.
 *
 * \tparam constrain_to_cell clamp vertex position to the bounding box of the cell
 */
template <bool constrain_to_cell = false>
class QEM_SVD
{
public:
  /*
   * \brief computes the vertex position for a point in Dual Contouring.
   *
   * \tparam Domain must be a model of `IsosurfacingDomainWithGradient_3`.
   *
   * \param domain the domain providing input data and its topology
   * \param isovalue value of the isosurface
   * \param cell the cell within the domain for which the vertex position is computed
   * \param point the point position of the vertex that belongs to that cell
   *
   * \return `true` if the voxel intersects the isosurface, `false` otherwise
   */
  template <typename Domain>
  bool position(const Domain& domain,
                const typename Domain::Geom_traits::FT isovalue,
                const typename Domain::Cell_descriptor& cell,
                typename Domain::Geom_traits::Point_3& p) const
  {
    using Geom_traits = typename Domain::Geom_traits;
    using FT = typename Geom_traits::FT;
    using Point_3 = typename Geom_traits::Point_3;
    using Vector_3 = typename Geom_traits::Vector_3;

    using Eigen_vector_3 = Eigen_vector<FT, 3>;
    using Eigen_vector_x = Eigen_vector<FT>;
    using Eigen_matrix_3 = Eigen_matrix<FT, 3, 3>;
    using Eigen_matrix_x = Eigen_matrix<FT>;

    using Vertex_descriptor = typename Domain::Vertex_descriptor;

    auto x_coord = domain.geom_traits().compute_x_3_object();
    auto y_coord = domain.geom_traits().compute_y_3_object();
    auto z_coord = domain.geom_traits().compute_z_3_object();
    auto point = domain.geom_traits().construct_point_3_object();

    typename Domain::Cell_vertices vertices = domain.cell_vertices(cell);
    const std::size_t cn = vertices.size();

    // compute edge intersections
    std::vector<Point_3> edge_intersections;
    std::vector<Vector_3> edge_intersection_normals;

    for(const auto& edge : domain.cell_edges(cell))
    {
      const auto& edge_vertices = domain.incident_vertices(edge);
      const Vertex_descriptor& v0 = edge_vertices[0];
      const Vertex_descriptor& v1 = edge_vertices[1];

      const FT& val0 = domain.value(v0);
      const FT& val1 = domain.value(v1);

      const Point_3& p0 = domain.point(v0);
      const Point_3& p1 = domain.point(v1);

      if((val0 <= isovalue) != (val1 <= isovalue))
      {
        // current edge is intersected by the isosurface
        const FT u = (val0 - isovalue) / (val0 - val1);
        const Point_3 p_lerp = point((FT(1) - u) * x_coord(p0) + u * x_coord(p1),
                                     (FT(1) - u) * y_coord(p0) + u * y_coord(p1),
                                     (FT(1) - u) * z_coord(p0) + u * z_coord(p1));
        edge_intersections.push_back(p_lerp);
        edge_intersection_normals.push_back(domain.gradient(p_lerp));
      }
    }

    if(edge_intersections.empty())
      return false;

    FT x_min, y_min, z_min, x_max, y_max, z_max;
    FT x(0), y(0), z(0);

    for(const auto& v : vertices)
    {
      const Point_3 cp = domain.point(v);
      x += x_coord(cp);
      y += y_coord(cp);
      z += z_coord(cp);

      if constexpr(constrain_to_cell)
      {
        x_min = (std::min<FT>)(x_min, x_coord(cp));
        y_min = (std::min<FT>)(y_min, y_coord(cp));
        z_min = (std::min<FT>)(z_min, z_coord(cp));

        x_max = (std::max<FT>)(x_max, x_coord(cp));
        y_max = (std::max<FT>)(y_max, y_coord(cp));
        z_max = (std::max<FT>)(z_max, z_coord(cp));
      }
    }

    p = point(x / cn, y / cn, z / cn);

    // SVD QEM
    Eigen_matrix_3 A;
    A.setZero();
    Eigen_vector_3 rhs;
    rhs.setZero();
    for(std::size_t i=0; i<edge_intersections.size(); ++i)
    {
      Eigen_vector_3 n_k;
      n_k.set(0, x_coord(edge_intersection_normals[i]));
      n_k.set(1, y_coord(edge_intersection_normals[i]));
      n_k.set(2, z_coord(edge_intersection_normals[i]));

      Eigen_vector_3 p_k;
      p_k.set(0, x_coord(edge_intersections[i]));
      p_k.set(1, y_coord(edge_intersections[i]));
      p_k.set(2, z_coord(edge_intersections[i]));

      const FT d_k = n_k.transpose() * p_k;

      // have to cast to cast to the underlying type because the type of 'n_k * n_k.transpose()' is Eigen::Product< ... > etc., so the double conversion is not implicit
      Eigen_matrix_3 A_k = typename Eigen_matrix_3::EigenType(n_k * n_k.transpose());
      Eigen_vector_3 b_k;
      b_k = d_k * n_k;
      A += A_k;
      rhs += b_k;
    }

    Eigen::JacobiSVD<typename Eigen_matrix_x::EigenType> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);

    // set threshold as in Peter Lindstrom's paper, "Out-of-Core Simplification of Large Polygonal Models"
    svd.setThreshold(1e-3);

    // Init x hat
    Eigen_vector_3 x_hat;
    x_hat << x_coord(p), y_coord(p), z_coord(p);

    // Lindstrom formula for QEM new position for singular matrices
    Eigen_vector_x v_svd;
    v_svd = x_hat + svd.solve(rhs - A * x_hat);

    // bounding box
    if constexpr(constrain_to_cell)
    {
      v_svd[0] = std::clamp<FT>(v_svd[0], x_min, x_max);
      v_svd[1] = std::clamp<FT>(v_svd[1], y_min, y_max);
      v_svd[2] = std::clamp<FT>(v_svd[2], z_min, z_max);
    }

    p = point(v_svd[0], v_svd[1], v_svd[2]);

    return true;
  }
};

/**
 * \brief returns the cell's center.
 */
class Cell_center
{
public:
  /*
   * \brief computes the vertex position for a point in Dual Contouring.
   *
   * \tparam Domain must be a model of `IsosurfacingDomainWithGradient_3`.
   *
   * \param domain the domain providing input data and its topology
   * \param isovalue value of the isosurface
   * \param cell the cell within the domain for which the vertex position is computed
   * \param point the point position of the vertex that belongs to that cell
   *
   * \return `true` if the voxel intersects the isosurface, `false` otherwise
   */
  template <typename Domain>
  bool position(const Domain& domain,
                const typename Domain::Geom_traits::FT isovalue,
                const typename Domain::Cell_descriptor& vh,
                typename Domain::Geom_traits::Point_3& point) const
  {
    using FT = typename Domain::Geom_traits::FT;
    using Point_3 = typename Domain::Geom_traits::Vector_3;

    using Vertex_descriptor = typename Domain::Vertex_descriptor;

    typename Domain::Cell_vertices vertices = domain.cell_vertices(vh);

    std::vector<Point_3> pos(vertices.size());
    std::transform(vertices.begin(), vertices.end(), pos.begin(),
                   [&](const Vertex_descriptor& v) { return domain.point(v); });

    // set point to cell center
    point = CGAL::centroid(pos.begin(), pos.end(), CGAL::Dimension_tag<0>());

    bool all_smaller = true;
    bool all_greater = true;
    for(const auto& v : vertices)
    {
      const bool b = (domain.value(v) <= isovalue);
      all_smaller = all_smaller && b;
      all_greater = all_greater && !b;
    }

    if(all_smaller || all_greater)
      return false;

    return true;
  }
};

/*
 * \brief computes the centroid of all cell edge intersections with the isosurface.
 */
class Centroid_of_edge_intersections
{
public:
  /*
   * \brief computes the vertex position for a point in Dual Contouring.
   *
   * \tparam Domain must be a model of `IsosurfacingDomainWithGradient_3`.
   *
   * \param domain the domain providing input data and its topology
   * \param isovalue value of the isosurface
   * \param cell the cell within the domain for which the vertex position is computed
   * \param point the point position of the vertex that belongs to that cell
   *
   * \return `true` if the voxel intersects the isosurface, `false` otherwise
   */
  template <typename Domain>
  bool position(const Domain& domain,
                const typename Domain::Geom_traits::FT isovalue,
                const typename Domain::Cell_descriptor& cell,
                typename Domain::Geom_traits::Point_3& point) const
  {
    using FT = typename Domain::Geom_traits::FT;
    using Point_3 = typename Domain::Geom_traits::Point_3;
    using Vector_3 = typename Domain::Geom_traits::Vector_3;

    using Vertex_descriptor = typename Domain::Vertex_descriptor;
    using Edge_descriptor = typename Domain::Edge_descriptor;

    typename Domain::Cell_vertices vertices = domain.cell_vertices(cell);

    // compute edge intersections
    std::vector<Point_3> edge_intersections;

    for(const Edge_descriptor& edge : domain.cell_edges(cell))
    {
      const auto& edge_vertices = domain.incident_vertices(edge);
      const Vertex_descriptor& v0 = edge_vertices[0];
      const Vertex_descriptor& v1 = edge_vertices[1];

      const FT val0 = domain.value(v0);
      const FT val1 = domain.value(v1);

      const Point_3& p0 = domain.point(v0);
      const Point_3& p1 = domain.point(v1);

      if((val0 <= isovalue) != (val1 <= isovalue))
      {
        // current edge is intersected by the isosurface
        const FT u = (val0 - isovalue) / (val0 - val1);
        const Point_3 p_lerp = CGAL::ORIGIN + ((FT(1) - u) * (p0 - CGAL::ORIGIN) + u * (p1 - CGAL::ORIGIN));
        edge_intersections.push_back(p_lerp);
      }
    }

    if(edge_intersections.empty())
      return false;

    point = CGAL::centroid(edge_intersections.begin(), edge_intersections.end(),
                            CGAL::Dimension_tag<0>());  // set point to center of edge intersections

    return true;
  }
};

} // namespace Positioning

template <typename Domain,
          typename Positioning>
class Dual_contouring_vertex_positioning
{
private:
  using FT = typename Domain::Geom_traits::FT;
  using Point_3 = typename Domain::Geom_traits::Point_3;

  using Cell_descriptor = typename Domain::Cell_descriptor;

public:
  Dual_contouring_vertex_positioning(const Domain& domain,
                                     const FT isovalue,
                                     const Positioning& positioning)
    : domain(domain),
      isovalue(isovalue),
      positioning(positioning),
      points_counter(0)
  { }

  void operator()(const Cell_descriptor& v)
  {
    // compute dc-vertices
    Point_3 p; // fixme: initialize?
    if(positioning.position(domain, isovalue, v, p))
    {
      std::lock_guard<std::mutex> lock(mutex);
      map_voxel_to_point[v] = p;
      map_voxel_to_point_id[v] = points_counter++;
    }
  }

  // private: // @todo
  const Domain& domain;
  const FT isovalue;
  const Positioning& positioning;

  std::map<Cell_descriptor, std::size_t> map_voxel_to_point_id;
  std::map<Cell_descriptor, Point_3> map_voxel_to_point;
  std::size_t points_counter;

  std::mutex mutex;
};

template <typename Domain>
class Dual_contouring_face_generation
{
private:
  using FT = typename Domain::Geom_traits::FT;

  using Edge_descriptor = typename Domain::Edge_descriptor;
  using Cell_descriptor = typename Domain::Cell_descriptor;

public:
  Dual_contouring_face_generation(const Domain& domain,
                                  const FT isovalue)
    : domain(domain),
      isovalue(isovalue)
  { }

  void operator()(const Edge_descriptor& e)
  {
    // save all faces
    const auto& vertices = domain.incident_vertices(e);
    const FT s0 = domain.value(vertices[0]);
    const FT s1 = domain.value(vertices[1]);

    if(s0 <= isovalue && s1 > isovalue)
    {
      const auto& voxels = domain.incident_cells(e);

      std::lock_guard<std::mutex> lock(mutex);
      faces[e].insert(faces[e].begin(), voxels.begin(), voxels.end());
    }
    else if(s1 <= isovalue && s0 > isovalue)
    {
      const auto& voxels = domain.incident_cells(e);

      std::lock_guard<std::mutex> lock(mutex);
      faces[e].insert(faces[e].begin(), voxels.rbegin(), voxels.rend());
    }
  }

  // private: // @todo
  std::map<Edge_descriptor, std::vector<Cell_descriptor> > faces;

  const Domain& domain;
  const FT isovalue;

  std::mutex mutex;
};

} // namespace internal
} // namespace Isosurfacing
} // namespace CGAL

#endif // CGAL_ISOSURFACING_3_INTERNAL_DUAL_CONTOURING_FUNCTORS_H
