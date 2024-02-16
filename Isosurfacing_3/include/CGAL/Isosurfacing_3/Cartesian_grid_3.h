// Copyright (c) 2022-2024 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Julian Stahl
//                 Mael Rouxel-Labbé

#ifndef CGAL_ISOSURFACING_3_CARTESIAN_GRID_3_H
#define CGAL_ISOSURFACING_3_CARTESIAN_GRID_3_H

#include <CGAL/license/Isosurfacing_3.h>

#include <CGAL/Isosurfacing_3/internal/partition_traits_Cartesian_grid.h>
#include <CGAL/Isosurfacing_3/IO/OBJ.h>

#include <CGAL/assertions.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/Named_function_parameters.h>

#include <array>
#include <type_traits>

namespace CGAL {
namespace Isosurfacing {

/**
 * \ingroup IS_Partitions_helpers_grp
 *
 * A policy to choose whether grid vertex positions should be cached, or recomputed at each access.
 *
 * \tparam Tag a tag that is either `Tag_true` (positions are cached) or `Tag_false` (positions are not cached).
 */
template <typename Tag>
struct Grid_vertex_memory_policy : public Tag { };

/**
 * \ingroup IS_Partitions_helpers_grp
 *
 * A convenience alias for the policy that caches grid vertex positions.
 */
using Cache_positions = Grid_vertex_memory_policy<Tag_true>;

/**
 * \ingroup IS_Partitions_helpers_grp
 *
 * A convenience alias for the policy that does not cache grid vertex positions.
 */
using Do_not_cache_positions = Grid_vertex_memory_policy<Tag_false>;

/**
 * \ingroup IS_Partitions_grp
 *
 * \cgalModels{Partition_3}
 *
 * \brief The class `Cartesian_grid_3` represents a 3D %Cartesian grid, that is the partition of
 * an iso-cuboid into identical iso-cuboidal cells.
 *
 * The class `Cartesian_grid_3` is one of the possible space partitioning data structures
 * that can be used along with value and gradient fields to make up a domain.
 *
 * \tparam GeomTraits must be a model of `IsosurfacingTraits_3`.
 * \tparam MemoryPolicy whether the geometric positions of the grid vertices are stored or not.
 *                      Possible values are `CGAL::Isosurfacing::Cache_positions` and `CGAL::Isosurfacing::Do_not_cache_positions`.
 *
 * \sa `CGAL::Isosurfacing::Marching_cubes_domain_3()`
 * \sa `CGAL::Isosurfacing::Dual_contouring_domain_3()`
 */
template <typename GeomTraits,
          typename MemoryPolicy = Do_not_cache_positions> // @todo actually implement it
class Cartesian_grid_3
{
public:
  using Geom_traits = GeomTraits;
  using FT = typename Geom_traits::FT;
  using Point_3 = typename Geom_traits::Point_3;
  using Vector_3 = typename Geom_traits::Vector_3;
  using Iso_cuboid_3 = typename Geom_traits::Iso_cuboid_3;

private:
  Iso_cuboid_3 m_bbox;
  std::array<std::size_t, 3> m_sizes;
  std::array<FT, 3> m_spacing;

  Geom_traits m_gt;

private:
  void compute_spacing()
  {
    auto x_coord = m_gt.compute_x_3_object();
    auto y_coord = m_gt.compute_y_3_object();
    auto z_coord = m_gt.compute_z_3_object();
    auto vertex = m_gt.construct_vertex_3_object();

    // calculate grid spacing
    const Point_3& min_p = vertex(m_bbox, 0);
    const Point_3& max_p = vertex(m_bbox, 7);
    const FT x_span = x_coord(max_p) - x_coord(min_p);
    const FT y_span = y_coord(max_p) - y_coord(min_p);
    const FT z_span = z_coord(max_p) - z_coord(min_p);

    const FT d_x = x_span / (m_sizes[0] - 1);
    const FT d_y = y_span / (m_sizes[1] - 1);
    const FT d_z = z_span / (m_sizes[2] - 1);

    m_spacing = make_array(d_x, d_y, d_z);
  }

public:
  /*!
   * \brief Default constructor
   */
  Cartesian_grid_3()
    : m_bbox{Point_3{0, 0, 0}, Point_3{0, 0, 0}},
      m_sizes{2, 2, 2},
      m_spacing{0, 0, 0},
      m_gt{Geom_traits()}
  { }

  /**
   * \brief creates a %Cartesian grid with `xdim * ydim * zdim` grid vertices.
   *
   * The grid covers the space described by a bounding box.
   *
   * \param bbox the bounding box of the grid
   * \param xdim the number of grid vertices in the `x` direction
   * \param ydim the number of grid vertices in the `y` direction
   * \param zdim the number of grid vertices in the `z` direction
   * \param gt the geometric traits
   *
   * \pre `xdim`, `ydim`, and `zdim` are (strictly) positive.
   */
  Cartesian_grid_3(const Iso_cuboid_3& bbox,
                   const std::size_t xdim,
                   const std::size_t ydim,
                   const std::size_t zdim,
                   const Geom_traits& gt = Geom_traits())
    : m_bbox{bbox},
      m_sizes{xdim, ydim, zdim},
      m_gt{gt}
  {
    compute_spacing();
  }

  /**
   * \brief creates a %Cartesian grid with `xdim * ydim * zdim` grid vertices.
   *
   * The grid covers the space described by a bounding box, itself described through two diagonal corners.
   *
   * \param p the lowest corner of the bounding box of the grid
   * \param q the upper corner of the bounding box of the grid
   * \param xdim the number of grid vertices in the `x` direction
   * \param ydim the number of grid vertices in the `y` direction
   * \param zdim the number of grid vertices in the `z` direction
   * \param gt the geometric traits
   *
   * \pre `p` is lexicographically (strictly) smaller than `q`
   * \pre `xdim`, `ydim`, and `zdim` are (strictly) positive.
   */
  Cartesian_grid_3(const Point_3& p, const Point_3& q,
                   const std::size_t xdim,
                   const std::size_t ydim,
                   const std::size_t zdim,
                   const Geom_traits& gt = Geom_traits())
    : Cartesian_grid_3{Iso_cuboid_3{p, q}, xdim, ydim, zdim, gt}
  { }

  /**
   * \brief creates a %Cartesian grid using a prescribed grid step.
   *
   * The grid covers the space described by a bounding box.
   *
   * \param bbox the bounding box of the grid
   * \param spacing the dimension of the paving cell, in the `x`, `y`, and `z` directions, respectively.
   * \param gt the geometric traits
   *
   * \pre the diagonal of `bbox` has length a multiple of `spacing`
   */
  Cartesian_grid_3(const Iso_cuboid_3& bbox,
                   const std::array<FT, 3>& spacing,
                   const Geom_traits& gt = Geom_traits())
    : m_bbox{bbox},
      m_spacing{spacing},
      m_gt{gt}
  {
    auto x_coord = gt.compute_x_3_object();
    auto y_coord = gt.compute_y_3_object();
    auto z_coord = gt.compute_z_3_object();
    auto vertex = gt.construct_vertex_3_object();
    auto vector = gt.construct_vector_3_object();

    const Point_3& min_p = vertex(bbox, 0);
    const Point_3& max_p = vertex(bbox, 7);
    const FT x_span = x_coord(max_p) - x_coord(min_p);
    const FT y_span = y_coord(max_p) - y_coord(min_p);
    const FT z_span = z_coord(max_p) - z_coord(min_p);

    m_sizes[0] = static_cast<std::size_t>(std::ceil(x_span / spacing[0])) + 1;
    m_sizes[1] = static_cast<std::size_t>(std::ceil(y_span / spacing[1])) + 1;
    m_sizes[2] = static_cast<std::size_t>(std::ceil(z_span / spacing[2])) + 1;
  }

  /**
   * \brief creates a %Cartesian grid using a prescribed grid step.
   *
   * The grid covers the space described by a bounding box, itself described through two diagonal corners.
   *
   * \param p the lowest corner of the bounding box of the grid
   * \param q the upper corner of the bounding box of the grid
   * \param spacing the dimension of the paving cell, in the `x`, `y`, and `z` directions, respectively.
   * \param gt the geometric traits
   *
   * \pre `p` is lexicographically (strictly) smaller than `q`
   * \pre the diagonal of the bounding box has length a multiple of `spacing`
   */
  Cartesian_grid_3(const Point_3& p, const Point_3& q,
                   const std::array<FT, 3>& spacing,
                   const Geom_traits& gt = Geom_traits())
    : Cartesian_grid_3{Iso_cuboid_3{p, q}, spacing, gt}
  { }

public:
  /**
   * \return the geometric traits class
   */
  const Geom_traits& geom_traits() const
  {
    return m_gt;
  }

  /**
   * \return the bounding box of the %Cartesian grid
   */
  const Iso_cuboid_3& bbox() const { return m_bbox; }

  /**
   * sets the bounding box of the %Cartesian grid and recomputes the spacing.
   */
  void set_bbox(const Iso_cuboid_3& bbox)
  {
    m_bbox = bbox;
    compute_spacing();
  }

  /**
   * \return the spacing of the %Cartesian grid, that is a vector whose coordinates are
   *         the grid steps in the `x`, `y`, and `z` directions, respectively
   */
  const std::array<FT, 3>& spacing() const { return m_spacing; }

  /**
   * \return the number of grid vertices in the `x` direction
   */
  std::size_t xdim() const { return m_sizes[0]; }

  /**
   * \return the number of grid vertices in the `y` direction
   */
  std::size_t ydim() const { return m_sizes[1]; }

  /**
   * \return the number of grid vertices in the `z` direction
   */
  std::size_t zdim() const { return m_sizes[2]; }

  /**
   * sets the number of grid vertices in the `x`, `y`, and `z` directions, respectively,
   * and recomputes the spacing.
   */
  void set_sizes(const std::size_t xdim, const std::size_t ydim, const std::size_t zdim)
  {
    m_sizes = {xdim, ydim, zdim};
    compute_spacing();
  }

public:
  /**
   * \brief gets the canonical index of a grid cell given its indices.
  */
  std::size_t linear_index(const std::size_t i,
                           const std::size_t j,
                           const std::size_t k) const
  {
    CGAL_precondition(i < m_sizes[0] && j < m_sizes[1] && k < m_sizes[2]);
    return (k * m_sizes[1] + j) * m_sizes[0] + i;
  }

public:
  /**
   * \brief gets the index of the grid cell that contains a given point.
   *
   * \param p the point to be located
   *
   * \return the index of the grid cell that contains `p`
   *
   * \pre `p` is inside the bounding box of the grid.
   */
  std::array<std::size_t, 3> index(const Point_3& p) const
  {
    auto x_coord = m_gt.compute_x_3_object();
    auto y_coord = m_gt.compute_y_3_object();
    auto z_coord = m_gt.compute_z_3_object();
    auto vertex = m_gt.construct_vertex_3_object();

    const Point_3& min_p = vertex(m_bbox, 0);
    std::size_t i = (x_coord(p) - x_coord(min_p)) / m_spacing[0];
    std::size_t j = (y_coord(p) - y_coord(min_p)) / m_spacing[1];
    std::size_t k = (z_coord(p) - z_coord(min_p)) / m_spacing[2];

    // @todo check this
    if(i == xdim() - 1)
      --i;
    if(j == ydim() - 1)
      --j;
    if(k == zdim() - 1)
      --k;

    return {i, j, k};
  }

  // Geometry
public:
  /**
   * \brief returns the geometric position of the grid vertex described by a set of indices.
   *
   * Depending on the value of the template parameter `cache_points`, positions might not be stored
   * but calculated using the lowest corner of the bounding box and grid spacing.
   *
   * \param i the index in the `x` direction
   * \param j the index in the `y` direction
   * \param k the index in the `z` direction
   *
   * \pre `i < xdim()` and `j < ydim()` and `k < zdim()`
   */
  Point_3 point(const std::size_t i,
                const std::size_t j,
                const std::size_t k) const
  {
    typename Geom_traits::Compute_x_3 x_coord = m_gt.compute_x_3_object();
    typename Geom_traits::Compute_y_3 y_coord = m_gt.compute_y_3_object();
    typename Geom_traits::Compute_z_3 z_coord = m_gt.compute_z_3_object();
    typename Geom_traits::Construct_point_3 point = m_gt.construct_point_3_object();
    typename Geom_traits::Construct_vertex_3 vertex = m_gt.construct_vertex_3_object();

    const Point_3& min_p = vertex(m_bbox, 0);
    return point(x_coord(min_p) + i * m_spacing[0],
                 y_coord(min_p) + j * m_spacing[1],
                 z_coord(min_p) + k * m_spacing[2]);
  }
};

} // namespace Isosurfacing
} // namespace CGAL

#endif // CGAL_ISOSURFACING_3_CARTESIAN_GRID_3_H
