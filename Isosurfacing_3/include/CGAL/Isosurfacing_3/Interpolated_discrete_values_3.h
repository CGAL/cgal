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

#ifndef CGAL_ISOSURFACING_3_INTERPOLATED_DISCRETE_VALUES_3_H
#define CGAL_ISOSURFACING_3_INTERPOLATED_DISCRETE_VALUES_3_H

#include <CGAL/license/Isosurfacing_3.h>

#include <CGAL/Isosurfacing_3/internal/partition_traits.h>

#include <CGAL/Isosurfacing_3/interpolation_schemes_3.h>

#include <array>
#include <vector>

namespace CGAL {
namespace Isosurfacing {

/**
 * \ingroup IS_Domain_helpers_grp
 *
 * \cgalModels{IsosurfacingValueField_3}
 *
 * \brief Class template for a field of values that are calculated using discrete values and interpolation.
 *
 * \tparam Grid must be `CGAL::Isosurfacing::Cartesian_grid_3<GeomTraits>`, with `GeomTraits` a model of `IsosurfacingTraits_3`
 * \tparam InterpolationScheme must be a model of `IsosurfacingInterpolationScheme_3`
 */
template <typename Grid,
          typename InterpolationScheme = Trilinear_interpolation<Grid> >
class Interpolated_discrete_values_3
{
  using Geom_traits = typename Grid::Geom_traits;
  using FT = typename Geom_traits::FT;
  using Point_3 = typename Geom_traits::Point_3;

  using Vertex_descriptor = typename partition_traits<Grid>::Vertex_descriptor;

private:
  const Grid& m_grid;
  const InterpolationScheme m_interpolation;

  std::vector<FT> m_values;

public:
  Interpolated_discrete_values_3(const Grid& grid,
                                 const InterpolationScheme& interpolation = InterpolationScheme())
    : m_grid{grid},
      m_interpolation{interpolation}
  {
    // pre-allocate memory
    const std::size_t nv = grid.xdim() * grid.ydim() * grid.zdim();
    m_values.resize(nv);
  }

public:
  /**
   * \brief returns the scalar value stored at the grid vertex described by a set of indices.
   *
   * \note This function can be used to set the value at a grid vertex.
   *
   * \param i the index in the `x` direction
   * \param j the index in the `y` direction
   * \param k the index in the `z` direction
   */
  FT& operator()(const std::size_t i,
                 const std::size_t j,
                 const std::size_t k)
  {
    const std::size_t id = m_grid.linear_index(i, j, k);
    if(id >= m_values.size())
      m_values.resize(id + 1);
    return m_values[id];
  }

  /**
   * \brief returns the scalar value stored at the grid vertex described by a set of indices.
   *
   * \param i the index in the `x` direction
   * \param j the index in the `y` direction
   * \param k the index in the `z` direction
   *
   * \pre `i < xdim()` and `j < ydim()` and `k < zdim()`
   */
  FT operator()(const std::size_t i,
                const std::size_t j,
                const std::size_t k) const
  {
    CGAL_precondition(i < m_grid.xdim() && j < m_grid.ydim() && k < m_grid.zdim());
    return m_values[m_grid.linear_index(i, j, k)];
  }

  /*!
   * returns the value at vertex `v`.
   */
  FT operator()(const Vertex_descriptor& v) const
  {
    return this->operator()(v[0], v[1], v[2]);
  }

  /*!
   * returns the value at point `p`.
   */
  FT operator()(const Point_3& p) const
  {
    return m_interpolation.interpolate_values(p, m_grid, m_values);
  }
};

} // namespace Isosurfacing
} // namespace CGAL

#endif // CGAL_ISOSURFACING_3_INTERPOLATED_DISCRETE_VALUES_3_H
