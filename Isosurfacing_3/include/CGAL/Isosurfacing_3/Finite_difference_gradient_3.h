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

#ifndef CGAL_ISOSURFACING_3_FINITE_DIFFERENCE_GRADIENT_3_H
#define CGAL_ISOSURFACING_3_FINITE_DIFFERENCE_GRADIENT_3_H

#include <CGAL/license/Isosurfacing_3.h>

#include <CGAL/number_utils.h>

#include <functional>

namespace CGAL {
namespace Isosurfacing {

/**
 * \ingroup IS_Fields_grp
 *
 * \cgalModels{IsosurfacingGradientField_3}
 *
 * \brief Class template for a gradient that is calculated using finite differences.
 *
 * \details This gradient function evaluates a value function at six points that are
 *          a given distance `delta` away from the queried point along the %Cartesian axes.
 *
 * \tparam GeomTraits must be a model of `IsosurfacingTraits_3`.
 */
template <typename GeomTraits>
class Finite_difference_gradient_3
{
public:
  using Geom_traits = GeomTraits;
  using FT = typename Geom_traits::FT;
  using Point_3 = typename Geom_traits::Point_3;
  using Vector_3 = typename Geom_traits::Vector_3;

private:
  const std::function<FT(const Point_3&)> m_function;
  const FT m_delta, m_half_step_inv;

  GeomTraits m_gt;

public:
  /**
   * \brief creates a new instance of this gradient class.
   *
   * \tparam ValueFunction must be a model of `ValueFunction_3`.
   *
   * \param function the function giving the scalar value at each point
   * \param delta the distance between samples for calculating the finite differences
   * \param gt the geometric traits class
   */
  template <typename ValueFunction>
  Finite_difference_gradient_3(const ValueFunction& function,
                               const FT delta,
                               const Geom_traits& gt = Geom_traits())
    : m_function{function},
      m_delta{delta},
      m_half_step_inv{FT{1} / (FT{2} * m_delta)},
      m_gt{gt}
  { }

  /**
   * \brief evaluates the gradient at a point in 3D space.
   *
   * \param p the position at which the gradient is computed
   */
  Vector_3 operator()(const Point_3& p) const
  {
    typename Geom_traits::Compute_x_3 x_coord = m_gt.compute_x_3_object();
    typename Geom_traits::Compute_y_3 y_coord = m_gt.compute_y_3_object();
    typename Geom_traits::Compute_z_3 z_coord = m_gt.compute_z_3_object();
    typename Geom_traits::Construct_point_3 point = m_gt.construct_point_3_object();
    typename Geom_traits::Construct_vector_3 vector = m_gt.construct_vector_3_object();

    // compute the gradient by sampling the function with finite differences
    // at six points with distance delta around the query point
    const FT x = x_coord(p), y = y_coord(p), z = z_coord(p);

    const Point_3 p0 = point(x + m_delta, y, z);
    const Point_3 p1 = point(x - m_delta, y, z);
    const Point_3 p2 = point(x, y + m_delta, z);
    const Point_3 p3 = point(x, y - m_delta, z);
    const Point_3 p4 = point(x, y, z + m_delta);
    const Point_3 p5 = point(x, y, z - m_delta);

    const FT gx = (m_function(p0) - m_function(p1)) * m_half_step_inv;
    const FT gy = (m_function(p2) - m_function(p3)) * m_half_step_inv;
    const FT gz = (m_function(p4) - m_function(p5)) * m_half_step_inv;

    const FT n = CGAL::approximate_sqrt(CGAL::square(gx) + CGAL::square(gy) + CGAL::square(gz));

    return vector(gx / n, gy / n, gz / n);
  }
};

} // namespace Isosurfacing
} // namespace CGAL

#endif // CGAL_ISOSURFACING_3_FINITE_DIFFERENCE_GRADIENT_3_H
