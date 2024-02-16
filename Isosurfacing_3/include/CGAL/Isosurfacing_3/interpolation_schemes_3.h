// Copyright (c) 2022-2024 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé

#ifndef CGAL_ISOSURFACING_3_INTERNAL_INTERPOLATION_SCHEMES_3_H
#define CGAL_ISOSURFACING_3_INTERNAL_INTERPOLATION_SCHEMES_3_H

#include <CGAL/license/Isosurfacing_3.h>

#include <CGAL/enum.h>

#include <array>
#include <vector>

namespace CGAL {
namespace Isosurfacing {

/*!
 * \ingroup IS_Fields_helpers_grp
 *
 * \cgalModels{InterpolationScheme_3}
 *
 * The class `Trilinear_interpolation` is the standard interpolation scheme to extrapolate
 * data defined only at vertices of a %Cartesian grid.
 *
 * \tparam Grid must be `CGAL::Isosurfacing::Cartesian_grid_3<GeomTraits>`, with `GeomTraits`
 *              a model of `IsosurfacingTraits_3`
 */
template <typename Grid>
class Trilinear_interpolation
{
public:
  using Geom_traits = typename Grid::Geom_traits;
  using FT = typename Geom_traits::FT;
  using Point_3 = typename Geom_traits::Point_3;
  using Vector_3 = typename Geom_traits::Vector_3;
  using Iso_cuboid_3 = typename Geom_traits::Iso_cuboid_3;

public:
  /*!
   * \brief interpolates the values at a given point using trilinear interpolation.
   *
   * \param p the point at which to interpolate the values
   * \param g the grid
   * \param values the continuous field of scalar values, defined over the bounding box of `g`
   *
   * \return the interpolated value at point `p`
   */
  FT interpolate_values(const Point_3& p,
                        const Grid& g,
                        const std::vector<FT>& values) const
{
    typename Geom_traits::Compute_x_3 x_coord = g.geom_traits().compute_x_3_object();
    typename Geom_traits::Compute_y_3 y_coord = g.geom_traits().compute_y_3_object();
    typename Geom_traits::Compute_z_3 z_coord = g.geom_traits().compute_z_3_object();
    typename Geom_traits::Construct_vertex_3 vertex = g.geom_traits().construct_vertex_3_object();

    // trilinear interpolation of stored values
    const Iso_cuboid_3& bbox = g.bbox();
    const std::array<FT, 3>& spacing = g.spacing();

    // calculate min index including border case
    const Point_3& min_p = vertex(bbox, 0);
    std::size_t i = (x_coord(p) - x_coord(min_p)) / spacing[0];
    std::size_t j = (y_coord(p) - y_coord(min_p)) / spacing[1];
    std::size_t k = (z_coord(p) - z_coord(min_p)) / spacing[2];

    // @todo check this
    if(i == g.xdim() - 1)
      --i;
    if(j == g.ydim() - 1)
      --j;
    if(k == g.zdim() - 1)
      --k;

    // calculate coordinates of min index
    const FT min_x = x_coord(min_p) + i * spacing[0];
    const FT min_y = y_coord(min_p) + j * spacing[1];
    const FT min_z = z_coord(min_p) + k * spacing[2];

    // interpolation factors between 0 and 1
    const FT f_i = std::clamp<FT>((x_coord(p) - min_x) / spacing[0], 0, 1);
    const FT f_j = std::clamp<FT>((y_coord(p) - min_y) / spacing[1], 0, 1);
    const FT f_k = std::clamp<FT>((z_coord(p) - min_z) / spacing[2], 0, 1);

    // read the value at all 8 corner points
    const FT g000 = values[g.linear_index(i + 0, j + 0, k + 0)];
    const FT g001 = values[g.linear_index(i + 0, j + 0, k + 1)];
    const FT g010 = values[g.linear_index(i + 0, j + 1, k + 0)];
    const FT g011 = values[g.linear_index(i + 0, j + 1, k + 1)];
    const FT g100 = values[g.linear_index(i + 1, j + 0, k + 0)];
    const FT g101 = values[g.linear_index(i + 1, j + 0, k + 1)];
    const FT g110 = values[g.linear_index(i + 1, j + 1, k + 0)];
    const FT g111 = values[g.linear_index(i + 1, j + 1, k + 1)];

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
    return g000 * lambda000 + g001 * lambda001 +
           g010 * lambda010 + g011 * lambda011 +
           g100 * lambda100 + g101 * lambda101 +
           g110 * lambda110 + g111 * lambda111;
  }

  /*!
   * \brief interpolates the gradients at a given point using trilinear interpolation.
   *
   * \param p the point at which to interpolate the gradients
   * \param g the grid
   * \param gradients the continuous field of vector values, defined over the bounding box of `g`
   *
   * \return the interpolated value at point `p`
   */
  Vector_3 interpolate_gradients(const Point_3& p,
                                 const Grid& g,
                                 const std::vector<Vector_3>& gradients) const
  {
    typename Geom_traits::Compute_x_3 x_coord = g.geom_traits().compute_x_3_object();
    typename Geom_traits::Compute_y_3 y_coord = g.geom_traits().compute_y_3_object();
    typename Geom_traits::Compute_z_3 z_coord = g.geom_traits().compute_z_3_object();
    typename Geom_traits::Construct_vector_3 vector = g.geom_traits().construct_vector_3_object();
    typename Geom_traits::Construct_vertex_3 vertex = g.geom_traits().construct_vertex_3_object();

    // trilinear interpolation of stored gradients
    const Iso_cuboid_3& bbox = g.bbox();
    const std::array<FT, 3>& spacing = g.spacing();

    // calculate min index including border case
    const Point_3& min_p = vertex(bbox, 0);
    std::size_t i = static_cast<std::size_t>((x_coord(p) - x_coord(min_p)) / spacing[0]);
    std::size_t j = static_cast<std::size_t>((y_coord(p) - y_coord(min_p)) / spacing[1]);
    std::size_t k = static_cast<std::size_t>((z_coord(p) - z_coord(min_p)) / spacing[2]);

    if(i == g.xdim() - 1) // dim is the point number
      --i;
    if(j == g.ydim() - 1)
      --j;
    if(k == g.zdim() - 1)
      --k;

    // calculate coordinates of min index
    const FT min_x = i * spacing[0] + x_coord(min_p);
    const FT min_y = j * spacing[1] + y_coord(min_p);
    const FT min_z = k * spacing[2] + z_coord(min_p);

    // interpolation factors between 0 and 1
    const FT f_i = std::clamp<FT>((x_coord(p) - min_x) / spacing[0], 0, 1);
    const FT f_j = std::clamp<FT>((y_coord(p) - min_y) / spacing[1], 0, 1);
    const FT f_k = std::clamp<FT>((z_coord(p) - min_z) / spacing[2], 0, 1);

    // read the value at all 8 corner points
    const Vector_3& g000 = gradients[g.linear_index(i + 0, j + 0, k + 0)];
    const Vector_3& g001 = gradients[g.linear_index(i + 0, j + 0, k + 1)];
    const Vector_3& g010 = gradients[g.linear_index(i + 0, j + 1, k + 0)];
    const Vector_3& g011 = gradients[g.linear_index(i + 0, j + 1, k + 1)];
    const Vector_3& g100 = gradients[g.linear_index(i + 1, j + 0, k + 0)];
    const Vector_3& g101 = gradients[g.linear_index(i + 1, j + 0, k + 1)];
    const Vector_3& g110 = gradients[g.linear_index(i + 1, j + 1, k + 0)];
    const Vector_3& g111 = gradients[g.linear_index(i + 1, j + 1, k + 1)];

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
    return vector(x_coord(g000) * lambda000 + x_coord(g001) * lambda001 +
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
                  z_coord(g110) * lambda110 + z_coord(g111) * lambda111);
  }
};

// This can be used for example when we have implicit functions for data (values & gradients),
// but use an interpolated values field as to store data.
template <typename Grid>
class Function_evaluation
{
  using Geom_traits = typename Grid::Geom_traits;
  using FT = typename Geom_traits::FT;
  using Point_3 = typename Geom_traits::Point_3;
  using Vector_3 = typename Geom_traits::Vector_3;

  std::function<FT(const Point_3&)> m_value_fn;
  std::function<Vector_3(const Point_3&)> m_gradient_fn;

public:
  template <typename ValueFunction, typename GradientFunction>
  Function_evaluation(const ValueFunction& value_fn,
                      const GradientFunction& gradient_fn = [](const Point_3&) -> Vector_3 { return CGAL::NULL_VECTOR; })
    : m_value_fn{value_fn},
      m_gradient_fn{gradient_fn}
  { }

public:
  FT interpolate_values(const Point_3& p, const Grid&, const std::vector<FT>&) const
  {
    return m_value_fn(p);
  }

  Vector_3 interpolate_gradients(const Point_3& p, const Grid&, const std::vector<Vector_3>&) const
  {
    return m_gradient_fn(p);
  }
};

} // namespace Isosurfacing
} // namespace CGAL

#endif // CGAL_ISOSURFACING_3_INTERNAL_INTERPOLATION_SCHEMES_3_H
