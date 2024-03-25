// Copyright (c) 2022-2024 INRIA Sophia-Antipolis (France), GeometryFactory (France).
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

#ifndef CGAL_ISOSURFACING_3_INTERPOLATED_DISCRETE_GRADIENTS_3_H
#define CGAL_ISOSURFACING_3_INTERPOLATED_DISCRETE_GRADIENTS_3_H

#include <CGAL/license/Isosurfacing_3.h>

#include <CGAL/Isosurfacing_3/internal/partition_traits.h>
#include <CGAL/Isosurfacing_3/interpolation_schemes_3.h>
#include <CGAL/Isosurfacing_3/Finite_difference_gradient_3.h>

#include <vector>

namespace CGAL {
namespace Isosurfacing {

/**
 * \ingroup IS_Domain_helpers_grp
 *
 * \cgalModels{IsosurfacingValueField_3}
 *
 * \brief Class template for a gradient that is calculated using discrete values and trilinear interpolation.
 *
 * \tparam Grid must be `CGAL::Isosurfacing::Cartesian_grid_3<GeomTraits>`, with `GeomTraits` a model of `IsosurfacingTraits_3`
 * \tparam InterpolationScheme must be a model of `IsosurfacingInterpolationScheme_3`
 */
template <typename Grid,
          typename InterpolationScheme = Trilinear_interpolation<Grid> >
class Interpolated_discrete_gradients_3
{
  using Geom_traits = typename Grid::Geom_traits;
  using FT = typename Geom_traits::FT;
  using Point_3 = typename Geom_traits::Point_3;
  using Vector_3 = typename Geom_traits::Vector_3;

  using vertex_descriptor = typename partition_traits<Grid>::vertex_descriptor;

private:
  const Grid& m_grid;
  const InterpolationScheme m_interpolation;

  std::vector<Vector_3> m_gradients;

public:
  Interpolated_discrete_gradients_3(const Grid& grid,
                                    const InterpolationScheme& interpolation = InterpolationScheme())
    : m_grid{grid},
      m_interpolation{interpolation}
  {
    // pre-allocate memory
    const std::size_t nv = grid.xdim() * grid.ydim() * grid.zdim();
    m_gradients.resize(nv);
  }

  /// computes (using finite difference) and stores gradients at the vertices of the grid.
  /// \tparam ValueField must be a model of `IsosurfacingValueField_3`
  /// \param values a field of values whose gradient are being computed
  template <typename ValueField>
  void compute_discrete_gradients(const ValueField& values)
  {
    const FT step = CGAL::approximate_sqrt(m_grid.spacing().squared_length()) * 0.01; // finite difference step
    Finite_difference_gradient_3<Geom_traits> g(values, step);

   for(std::size_t i=0; i<m_grid.xdim(); ++i)
      for(std::size_t j=0; j<m_grid.ydim(); ++j)
        for(std::size_t k=0; k<m_grid.zdim(); ++k)
          m_gradients[m_grid.linear_index(i, j, k)] = g(m_grid.point(i,j,k));
  }

public:
  /**
   * \brief returns the gradient stored at the grid vertex described by a set of indices.
   *
   * \note This function can be used to set the gradient at a grid vertex.
   *
   * \param i the index in the `x` direction
   * \param j the index in the `y` direction
   * \param k the index in the `z` direction
   *
   * \pre `i < xdim()` and `j < ydim()` and `k < zdim()`
   */
  Vector_3& operator()(const std::size_t i,
                       const std::size_t j,
                       const std::size_t k)
  {
    return m_gradients[m_grid.linear_index(i, j, k)];
  }

  /**
   * \brief returns the gradient stored at the grid vertex described by a set of indices.
   *
   * \param i the index in the `x` direction
   * \param j the index in the `y` direction
   * \param k the index in the `z` direction
   *
   * \pre `i < xdim()` and `j < ydim()` and `k < zdim()`
   */
  const Vector_3& operator()(const std::size_t i,
                             const std::size_t j,
                             const std::size_t k) const
  {
    return m_gradients[m_grid.linear_index(i, j, k)];
  }

  /*!
   * returns the gradient at a given point `p`.
   */
  Vector_3 operator()(const Point_3& p) const
  {
    return m_interpolation.interpolate_gradients(p, m_grid, m_gradients);
  }
};

} // namespace Isosurfacing
} // namespace CGAL

#endif // CGAL_ISOSURFACING_3_INTERPOLATED_DISCRETE_GRADIENTS_3_H
