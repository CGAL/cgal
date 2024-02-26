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

namespace internal {

template <typename GeomTraits,
          typename MemoryPolicy = Do_not_cache_positions>
struct Cartesian_grid_position
{
  using Point_3 = typename GeomTraits::Point_3;
  using Vector_3 = typename GeomTraits::Vector_3;
  using Iso_cuboid_3 = typename GeomTraits::Iso_cuboid_3;

  Cartesian_grid_position() { } // just for compilation

  Cartesian_grid_position(const Iso_cuboid_3& /*span*/,
                          const std::array<std::size_t, 3>& /*dims*/,
                          const Vector_3& /*spacing*/)
  { }

  template <typename Grid>
  Point_3 operator()(const std::size_t i,
                     const std::size_t j,
                     const std::size_t k,
                     const Grid& g) const
  {
    typename GeomTraits::Compute_x_3 x_coord = g.geom_traits().compute_x_3_object();
    typename GeomTraits::Compute_y_3 y_coord = g.geom_traits().compute_y_3_object();
    typename GeomTraits::Compute_z_3 z_coord = g.geom_traits().compute_z_3_object();
    typename GeomTraits::Construct_point_3 point = g.geom_traits().construct_point_3_object();
    typename GeomTraits::Construct_vertex_3 vertex = g.geom_traits().construct_vertex_3_object();

    const Point_3& min_p = vertex(g.span(), 0);
    return point(x_coord(min_p) + i * g.spacing()[0],
                 y_coord(min_p) + j * g.spacing()[1],
                 z_coord(min_p) + k * g.spacing()[2]);
  }
};

template <typename GeomTraits>
struct Cartesian_grid_position<GeomTraits, Cache_positions>
{
  using Point_3 = typename GeomTraits::Point_3;
  using Vector_3 = typename GeomTraits::Vector_3;
  using Iso_cuboid_3 = typename GeomTraits::Iso_cuboid_3;

  std::vector<Point_3> m_points;

  Cartesian_grid_position() { } // just for compilation

  Cartesian_grid_position(const Iso_cuboid_3& span,
                          const std::array<std::size_t, 3>& dims,
                          const Vector_3& spacing)
  {
    m_points.reserve(dims[0] * dims[1] * dims[2]);
    for(std::size_t k=0; k<dims[2]; ++k)
      for(std::size_t j=0; j<dims[1]; ++j)
        for(std::size_t i=0; i<dims[0]; ++i)
          m_points.emplace_back(span.min() + Vector_3(i * spacing[0], j * spacing[1], k * spacing[2]));
  }

  template <typename Grid>
  const Point_3& operator()(const std::size_t i,
                            const std::size_t j,
                            const std::size_t k,
                            const Grid& g) const
  {
    const std::size_t linear_index = g.linear_index(i, j, k);
    CGAL_precondition(linear_index < m_points.size());
    return m_points[linear_index];
  }
};

} // namespace internal

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
          typename MemoryPolicy = Do_not_cache_positions>
class Cartesian_grid_3
{
public:
  using Geom_traits = GeomTraits;
  using FT = typename Geom_traits::FT;
  using Point_3 = typename Geom_traits::Point_3;
  using Vector_3 = typename Geom_traits::Vector_3;
  using Iso_cuboid_3 = typename Geom_traits::Iso_cuboid_3;

  using Positioner = internal::Cartesian_grid_position<GeomTraits, MemoryPolicy>;

private:
  Iso_cuboid_3 m_span;
  std::array<std::size_t, 3> m_dims;
  Vector_3 m_spacing;

  Positioner m_positioner;
  Geom_traits m_gt;

private:
  void initialize_spacing()
  {
    typename Geom_traits::Compute_x_3 x_coord = m_gt.compute_x_3_object();
    typename Geom_traits::Compute_y_3 y_coord = m_gt.compute_y_3_object();
    typename Geom_traits::Compute_z_3 z_coord = m_gt.compute_z_3_object();
    typename Geom_traits::Construct_vector_3 vector = m_gt.construct_vector_3_object();
    typename Geom_traits::Construct_vertex_3 vertex = m_gt.construct_vertex_3_object();

    // calculate grid spacing
    const Point_3& min_p = vertex(m_span, 0);
    const Point_3& max_p = vertex(m_span, 7);
    const FT x_span = x_coord(max_p) - x_coord(min_p);
    const FT y_span = y_coord(max_p) - y_coord(min_p);
    const FT z_span = z_coord(max_p) - z_coord(min_p);

    const FT d_x = x_span / (m_dims[0] - 1);
    const FT d_y = y_span / (m_dims[1] - 1);
    const FT d_z = z_span / (m_dims[2] - 1);

    m_spacing = vector(d_x, d_y, d_z);

    m_positioner = Positioner { m_span, m_dims, m_spacing };
  }

  void initialize_dimensions()
  {
    typename Geom_traits::Compute_x_3 x_coord = m_gt.compute_x_3_object();
    typename Geom_traits::Compute_y_3 y_coord = m_gt.compute_y_3_object();
    typename Geom_traits::Compute_z_3 z_coord = m_gt.compute_z_3_object();
    typename Geom_traits::Construct_vertex_3 vertex = m_gt.construct_vertex_3_object();

    const Point_3& min_p = vertex(m_span, 0);
    const Point_3& max_p = vertex(m_span, 7);
    const FT x_span = x_coord(max_p) - x_coord(min_p);
    const FT y_span = y_coord(max_p) - y_coord(min_p);
    const FT z_span = z_coord(max_p) - z_coord(min_p);

    m_dims[0] = static_cast<std::size_t>(std::ceil(x_span / m_spacing[0])) + 1;
    m_dims[1] = static_cast<std::size_t>(std::ceil(y_span / m_spacing[1])) + 1;
    m_dims[2] = static_cast<std::size_t>(std::ceil(z_span / m_spacing[2])) + 1;

    m_positioner = Positioner { m_span, m_dims, m_spacing };
  }

public:
  /*!
   * \brief Default constructor
   */
  Cartesian_grid_3()
    : m_span{Point_3{0, 0, 0}, Point_3{0, 0, 0}},
      m_dims{2, 2, 2},
      m_spacing{0, 0, 0},
      m_gt{Geom_traits()}
  { }

  /**
   * \brief creates a %Cartesian grid with `xdim * ydim * zdim` grid vertices.
   *
   * The grid covers the space described by an iso-cuboid.
   *
   * \param span the geometric span of the grid
   * \param dimensions the number of grid vertices in the `x`, `y`, and `z` directions
   * \param gt the geometric traits
   *
   * \pre all dimensions are (strictly) positive.
   */
  Cartesian_grid_3(const Iso_cuboid_3& span,
                   const std::array<std::size_t, 3>& dimensions,
                   const Geom_traits& gt = Geom_traits())
    : m_span{span},
      m_dims{dimensions},
      m_gt{gt}
  {
    initialize_spacing();
  }

  /**
   * \brief creates a %Cartesian grid with `xdim * ydim * zdim` grid vertices.
   *
   * The grid covers the space described by an iso-cuboid,
   * itself described through two diagonal corners.
   *
   * \param p the lexicographically smallest corner of the iso-cuboid
   * \param q the lexicographically largest corner of the iso-cuboid
   * \param dimensions the number of grid vertices in the `x`, `y`, and `z` directions
   * \param gt the geometric traits
   *
   * \pre `p` is lexicographically (strictly) smaller than `q`
   * \pre all dimensions are (strictly) positive.
   */
  Cartesian_grid_3(const Point_3& p, const Point_3& q,
                   const std::array<std::size_t, 3>& dimensions,
                   const Geom_traits& gt = Geom_traits())
    : Cartesian_grid_3{Iso_cuboid_3{p, q}, dimensions, gt}
  { }

  /**
   * \brief creates a %Cartesian grid using a prescribed grid step.
   *
   * The grid covers the space described by an iso-cuboid.
   *
   * \param span the geometric span of the grid
   * \param spacing the dimension of the paving cell, in the `x`, `y`, and `z` directions
   * \param gt the geometric traits
   *
   * \pre the diagonal of `span` has length a multiple of `spacing`
   */
  Cartesian_grid_3(const Iso_cuboid_3& span,
                   const Vector_3& spacing,
                   const Geom_traits& gt = Geom_traits())
    : m_span{span},
      m_spacing{spacing},
      m_gt{gt}
  {
    initialize_dimensions();
  }

  /**
   * \brief creates a %Cartesian grid using a prescribed grid step.
   *
   * The grid covers the space described by an iso-cuboid, itself described through two diagonal corners.
   *
   * \param p the lexicographically smallest corner of the iso-cuboid
   * \param q the lexicographically largest corner of the iso-cuboid
   * \param spacing the dimension of the paving cell, in the `x`, `y`, and `z` directions, respectively.
   * \param gt the geometric traits
   *
   * \pre `p` is lexicographically (strictly) smaller than `q`
   * \pre the diagonal of the iso-cuboid has length a multiple of `spacing`
   */
  Cartesian_grid_3(const Point_3& p, const Point_3& q,
                   const Vector_3& spacing,
                   const Geom_traits& gt = Geom_traits())
    : Cartesian_grid_3{Iso_cuboid_3{p, q}, spacing, gt}
  { }

public:
  /**
   * returns the geometric traits class
   */
  const Geom_traits& geom_traits() const
  {
    return m_gt;
  }

  /**
   * returns an iso-cuboid representing the geometric span of the %Cartesian grid
   */
  const Iso_cuboid_3& span() const { return m_span; }

  /**
   * returns the number of grid vertices in the `x` direction
   */
  std::size_t xdim() const { return m_dims[0]; }

  /**
   * returns the number of grid vertices in the `y` direction
   */
  std::size_t ydim() const { return m_dims[1]; }

  /**
   * returns the number of grid vertices in the `z` direction
   */
  std::size_t zdim() const { return m_dims[2]; }

  /**
   * returns the spacing of the %Cartesian grid, that is a vector whose coordinates are
   * the grid steps in the `x`, `y`, and `z` directions, respectively
   */
  const Vector_3& spacing() const { return m_spacing; }

public:
  /**
   * \brief returns the index of a grid cell given its indices.
  */
  std::size_t linear_index(const std::size_t i,
                           const std::size_t j,
                           const std::size_t k) const
  {
    CGAL_precondition(i < m_dims[0] && j < m_dims[1] && k < m_dims[2]);
    return (k * m_dims[1] + j) * m_dims[0] + i;
  }

public:
  /**
   * \brief returns the index of the grid cell that contains a given point.
   *
   * \param p the point to be located
   *
   * \pre `p` is inside the grid.
   */
  std::array<std::size_t, 3> index(const Point_3& p) const
  {
    typename Geom_traits::Compute_x_3 x_coord = m_gt.compute_x_3_object();
    typename Geom_traits::Compute_y_3 y_coord = m_gt.compute_y_3_object();
    typename Geom_traits::Compute_z_3 z_coord = m_gt.compute_z_3_object();
    typename Geom_traits::Construct_vertex_3 vertex = m_gt.construct_vertex_3_object();

    const Point_3& min_p = vertex(m_span, 0);
    std::size_t i = std::size_t((x_coord(p) - x_coord(min_p)) / x_coord(m_spacing));
    std::size_t j = std::size_t((y_coord(p) - y_coord(min_p)) / y_coord(m_spacing));
    std::size_t k = std::size_t((z_coord(p) - z_coord(min_p)) / z_coord(m_spacing));

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
   * but calculated on-the-fly.
   *
   * \param i the index in the `x` direction
   * \param j the index in the `y` direction
   * \param k the index in the `z` direction
   *
   * \pre `i < xdim()` and `j < ydim()` and `k < zdim()`
   */
  decltype(auto) /*Point_3*/ point(const std::size_t i,
                                   const std::size_t j,
                                   const std::size_t k) const
  {
    return m_positioner(i, j, k, *this);
  }
};

} // namespace Isosurfacing
} // namespace CGAL

#endif // CGAL_ISOSURFACING_3_CARTESIAN_GRID_3_H
