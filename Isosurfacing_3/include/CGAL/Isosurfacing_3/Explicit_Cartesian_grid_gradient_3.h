// Copyright (c) 2022-2023 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Julian Stahl

#ifndef CGAL_ISOSURFACING_3_EXPLICIT_CARTESIAN_GRID_GRADIENT_3_H
#define CGAL_ISOSURFACING_3_EXPLICIT_CARTESIAN_GRID_GRADIENT_3_H

#include <CGAL/license/Isosurfacing_3.h>

#include <CGAL/Bbox_3.h>

namespace CGAL {
namespace Isosurfacing {

/**
 * \ingroup IS_Domains_grp
 *
 * \brief Class template for a gradient that is stored in a %Cartesian grid.
 *
 * \details The gradient at any point is calculated using trilinear interpolation.
 *
 * \tparam GeomTraits must be a model of `IsosurfacingTraits_3`.
 */
#ifdef DOXYGEN_RUNNING // Allow more than just Cartesian_grid_3
template <template <typename GeomTraits> class Cartesian_grid_3>
class Explicit_Cartesian_grid_gradient_3
#else
template <typename Grid>
class Explicit_Cartesian_grid_gradient_3
#endif
{
public:
  using Geom_traits = typename Grid::Geom_traits;
  using FT = typename Geom_traits::FT;
  using Point_3 = typename Geom_traits::Point_3;
  using Vector_3 = typename Geom_traits::Vector_3;

private:
  const Grid& m_grid;

public:
  /**
   * \brief creates a new instance of this gradient.
   *
   * \param grid the %Cartesian grid that stores the gradient
   */
  Explicit_Cartesian_grid_gradient_3(const Grid& grid)
    : m_grid(grid)
  { }

  /**
   * \brief evaluates the gradient at a point in space.
   *
   * \param p the position at which the gradient is computed
   */
  Vector_3 operator()(const Point_3& p) const
  {
    typename Geom_traits::Compute_x_3 x_coord = m_grid.geom_traits().compute_x_3_object();
    typename Geom_traits::Compute_y_3 y_coord = m_grid.geom_traits().compute_y_3_object();
    typename Geom_traits::Compute_z_3 z_coord = m_grid.geom_traits().compute_z_3_object();
    typename Geom_traits::Construct_scaled_vector_3 scale = m_grid.geom_traits().construct_scaled_vector_3_object();
    typename Geom_traits::Construct_sum_of_vectors_3 sum = m_grid.geom_traits().construct_sum_of_vectors_3_object();

    // trilinear interpolation of stored gradients
    const Bbox_3& bbox = m_grid.bbox();
    const Vector_3& spacing = m_grid.spacing();

    // calculate min index including border case
    std::size_t min_i = (x_coord(p) - bbox.xmin()) / x_coord(spacing);
    std::size_t min_j = (y_coord(p) - bbox.ymin()) / y_coord(spacing);
    std::size_t min_k = (z_coord(p) - bbox.zmin()) / z_coord(spacing);

    if(min_i == m_grid.xdim() - 1)
      --min_i;

    if(min_j == m_grid.ydim() - 1)
      --min_j;

    if(min_k == m_grid.zdim() - 1)
      --min_k;

    // calculate coordinates of min index
    const FT min_x = min_i * x_coord(spacing) + bbox.xmin();
    const FT min_y = min_j * y_coord(spacing) + bbox.ymin();
    const FT min_z = min_k * z_coord(spacing) + bbox.zmin();

    // interpolation factors between 0 and 1
    const FT f_i = (x_coord(p) - min_x) / x_coord(spacing);
    const FT f_j = (y_coord(p) - min_y) / y_coord(spacing);
    const FT f_k = (z_coord(p) - min_z) / z_coord(spacing);

    // read the gradient at all 8 corner points
    const Vector_3& g000 = m_grid.gradient(min_i + 0, min_j + 0, min_k + 0);
    const Vector_3& g001 = m_grid.gradient(min_i + 0, min_j + 0, min_k + 1);
    const Vector_3& g010 = m_grid.gradient(min_i + 0, min_j + 1, min_k + 0);
    const Vector_3& g011 = m_grid.gradient(min_i + 0, min_j + 1, min_k + 1);
    const Vector_3& g100 = m_grid.gradient(min_i + 1, min_j + 0, min_k + 0);
    const Vector_3& g101 = m_grid.gradient(min_i + 1, min_j + 0, min_k + 1);
    const Vector_3& g110 = m_grid.gradient(min_i + 1, min_j + 1, min_k + 0);
    const Vector_3& g111 = m_grid.gradient(min_i + 1, min_j + 1, min_k + 1);

    // interpolate along all axes by weighting the corner points
    const Vector_3 g0 = scale(g000, (1 - f_i) * (1 - f_j) * (1 - f_k));
    const Vector_3 g1 = scale(g001, (1 - f_i) * (1 - f_j) * f_k);
    const Vector_3 g2 = scale(g010, (1 - f_i) * f_j * (1 - f_k));
    const Vector_3 g3 = scale(g011, (1 - f_i) * f_j * f_k);
    const Vector_3 g4 = scale(g100, f_i * (1 - f_j) * (1 - f_k));
    const Vector_3 g5 = scale(g101, f_i * (1 - f_j) * f_k);
    const Vector_3 g6 = scale(g110, f_i * f_j * (1 - f_k));
    const Vector_3 g7 = scale(g111, f_i * f_j * f_k);

    // add weighted corners
    return sum(g0, sum(g1, sum(g2, sum(g3, sum(g4, sum(g5, sum(g6, g7)))))));
  }
};

} // namespace Isosurfacing
} // namespace CGAL

#endif // CGAL_ISOSURFACING_3_EXPLICIT_CARTESIAN_GRID_GRADIENT_3_H
