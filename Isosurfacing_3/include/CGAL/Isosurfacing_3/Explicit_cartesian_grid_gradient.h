// Copyright (c) 2022 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Julian Stahl

#ifndef CGAL_ISOSURFACING_3_EXPLICIT_CARTESIAN_GRID_GRADIENT_H
#define CGAL_ISOSURFACING_3_EXPLICIT_CARTESIAN_GRID_GRADIENT_H

#include <CGAL/license/Isosurfacing_3.h>

#include <CGAL/Isosurfacing_3/Cartesian_grid_3.h>

#include <CGAL/Bbox_3.h>

#include <memory>

namespace CGAL {
namespace Isosurfacing {

/**
 * \ingroup PkgIsosurfacing3Ref
 *
 * \brief Class template for a gradient that is stored in a `Cartesian_grid_3`.
 *
 * \details The gradient at any point is calculated using trilinear interpolation.
 *
 * \tparam GeomTraits must be a model of ``.
 */
template <typename GeomTraits>
class Explicit_cartesian_grid_gradient
{
public:
  using Geom_traits = GeomTraits;
  using FT = typename Geom_traits::FT;
  using Point = typename Geom_traits::Point_3;
  using Vector = typename Geom_traits::Vector_3;

  using Grid = Cartesian_grid_3<Geom_traits>;

public:
  /**
   * \ingroup PkgIsosurfacing3Ref
   *
   * \brief creates a new instance of this gradient.
   *
   * \param grid the Cartesian grid that stores the gradient
   */
  Explicit_cartesian_grid_gradient(const Grid& grid)
    : grid(grid)
  { }

  /**
   * \ingroup PkgIsosurfacing3Ref
   *
   * \brief evaluates the gradient at a point in space.
   *
   * \param point the position at which the gradient is computed
   */
  Vector operator()(const Point& point) const
  {
    // trilinear interpolation of stored gradients
    const Bbox_3& bbox = grid.get_bbox();
    const Vector& spacing = grid.get_spacing();

    // calculate min index including border case
    std::size_t min_i = (point.x() - bbox.xmin()) / spacing.x();
    std::size_t min_j = (point.y() - bbox.ymin()) / spacing.y();
    std::size_t min_k = (point.z() - bbox.zmin()) / spacing.z();
    if(min_i == grid.xdim() - 1)
      --min_i;

    if(min_j == grid.ydim() - 1)
      --min_j;

    if(min_k == grid.zdim() - 1)
      --min_k;

    // calculate coordinates of min index
    const FT min_x = min_i * spacing.x() + bbox.xmin();
    const FT min_y = min_j * spacing.y() + bbox.ymin();
    const FT min_z = min_k * spacing.z() + bbox.zmin();

    // interpolation factors between 0 and 1
    const FT f_i = (point.x() - min_x) / spacing.x();
    const FT f_j = (point.y() - min_y) / spacing.y();
    const FT f_k = (point.z() - min_z) / spacing.z();

    // read the gradient at all 8 corner points
    const Vector g000 = grid.gradient(min_i + 0, min_j + 0, min_k + 0);
    const Vector g001 = grid.gradient(min_i + 0, min_j + 0, min_k + 1);
    const Vector g010 = grid.gradient(min_i + 0, min_j + 1, min_k + 0);
    const Vector g011 = grid.gradient(min_i + 0, min_j + 1, min_k + 1);
    const Vector g100 = grid.gradient(min_i + 1, min_j + 0, min_k + 0);
    const Vector g101 = grid.gradient(min_i + 1, min_j + 0, min_k + 1);
    const Vector g110 = grid.gradient(min_i + 1, min_j + 1, min_k + 0);
    const Vector g111 = grid.gradient(min_i + 1, min_j + 1, min_k + 1);

    // interpolate along all axes by weighting the corner points
    const Vector g0 = g000 * (1 - f_i) * (1 - f_j) * (1 - f_k);
    const Vector g1 = g001 * (1 - f_i) * (1 - f_j) * f_k;
    const Vector g2 = g010 * (1 - f_i) * f_j * (1 - f_k);
    const Vector g3 = g011 * (1 - f_i) * f_j * f_k;
    const Vector g4 = g100 * f_i * (1 - f_j) * (1 - f_k);
    const Vector g5 = g101 * f_i * (1 - f_j) * f_k;
    const Vector g6 = g110 * f_i * f_j * (1 - f_k);
    const Vector g7 = g111 * f_i * f_j * f_k;

    // add weighted corners
    return g0 + g1 + g2 + g3 + g4 + g5 + g6 + g7;
  }

private:
  const Grid& grid;
};

} // namespace Isosurfacing
} // namespace CGAL

#endif // CGAL_ISOSURFACING_3_EXPLICIT_CARTESIAN_GRID_GRADIENT_H
