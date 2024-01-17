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
 * \details The gradient at a query point is calculated by trilinear interpolation.
 *
 * \warning This class keeps a pointer to the `grid` object, hence users must ensure that
 *          the lifetime of the `grid` object exceeds that of this gradient object.
 *
 * \tparam Grid must be a `CGAL::Isosurfacing::Cartesian_grid_3` whose `GeomTraits` template parameter
 *              is a model of `IsosurfacingTraits_3`.
*/
template <typename Grid> // allow more than just Cartesian_grid_3
class Explicit_Cartesian_grid_gradient_3
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
    typename Geom_traits::Construct_vector_3 vector = m_grid.geom_traits().construct_vector_3_object();

    // trilinear interpolation of stored gradients
    const Bbox_3& bbox = m_grid.bbox();
    const std::array<FT, 3>& spacing = m_grid.spacing();

    // calculate min index including border case
    std::size_t min_i = (x_coord(p) - bbox.xmin()) / spacing[0];
    std::size_t min_j = (y_coord(p) - bbox.ymin()) / spacing[1];
    std::size_t min_k = (z_coord(p) - bbox.zmin()) / spacing[2];

    if(min_i == m_grid.xdim() - 1)
      --min_i;

    if(min_j == m_grid.ydim() - 1)
      --min_j;

    if(min_k == m_grid.zdim() - 1)
      --min_k;

    // calculate coordinates of min index
    const FT min_x = min_i * spacing[0] + bbox.xmin();
    const FT min_y = min_j * spacing[1] + bbox.ymin();
    const FT min_z = min_k * spacing[2] + bbox.zmin();

    // interpolation factors between 0 and 1
    const FT f_i = (x_coord(p) - min_x) / spacing[0];
    const FT f_j = (y_coord(p) - min_y) / spacing[1];
    const FT f_k = (z_coord(p) - min_z) / spacing[2];

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
    const FT lambda000 = (FT(1) - f_i) * (FT(1) - f_j) * (FT(1) - f_k);
    const FT lambda001 = (FT(1) - f_i) * (FT(1) - f_j) * f_k;
    const FT lambda010 = (FT(1) - f_i) * f_j * (FT(1) - f_k);
    const FT lambda011 = (FT(1) - f_i) * f_j * f_k;
    const FT lambda100 = f_i * (FT(1) - f_j) * (FT(1) - f_k);
    const FT lambda101 = f_i * (FT(1) - f_j) * f_k;
    const FT lambda110 = f_i * f_j * (FT(1) - f_k);
    const FT lambda111 = f_i * f_j * f_k;

    // add weighted corners
    return vector( x_coord(g000) * lambda000 + x_coord(g001) * lambda001 +
                   x_coord(g010) * lambda010 + x_coord(g011) * lambda011 +
                   x_coord(g100) * lambda100 + x_coord(g101) * lambda101 +
                   x_coord(g110) * lambda110 + x_coord(g111) * lambda111,
                   y_coord(g000) * lambda000 + y_coord(g001) * lambda001 +
                   y_coord(g010) * lambda010 + y_coord(g011) * lambda011 +
                   y_coord(g100) * lambda100 + y_coord(g101) * lambda101 +
                   y_coord(g110) * lambda110 + y_coord(g111) * lambda111,
                   z_coord(g000) * lambda000 + z_coord(g001) * lambda001 +
                   z_coord(g010) * lambda010 + z_coord(g011) * lambda011 +
                   z_coord(g100) * lambda100 + z_coord(g101) * lambda101 +
                   z_coord(g110) * lambda110 + z_coord(g111) * lambda111 );
  }
};

} // namespace Isosurfacing
} // namespace CGAL

#endif // CGAL_ISOSURFACING_3_EXPLICIT_CARTESIAN_GRID_GRADIENT_3_H
