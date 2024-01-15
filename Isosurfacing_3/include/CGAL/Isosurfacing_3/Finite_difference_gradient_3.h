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

#ifndef CGAL_ISOSURFACING_3_FINITE_DIFFERENCE_GRADIENT_3_H
#define CGAL_ISOSURFACING_3_FINITE_DIFFERENCE_GRADIENT_3_H

#include <CGAL/license/Isosurfacing_3.h>

namespace CGAL {
namespace Isosurfacing {

/**
 * \ingroup IS_Domain_helpers_grp
 *
 * \brief Class template for a gradient that is calculated using finite differences.
 *
 * \details This gradient function evaluates an implicit value function at six points
 *          that are a given distance `delta` away from the queried point along the Cartesian axes.
 *
 * \tparam GeomTraits must be a model of `Kernel`.
 * \tparam PointFunction the type of the implicit function. It must be a model of `CopyConstructible` and implement
 *                       `GeomTraits::FT operator()(const GeomTraits::Point_3& point) const`.
 */
template <typename GeomTraits,
          typename PointFunction>
class Finite_difference_gradient_3
{
public:
  using Geom_traits = GeomTraits;
  using FT = typename Geom_traits::FT;
  using Point_3 = typename Geom_traits::Point_3;
  using Vector_3 = typename Geom_traits::Vector_3;

  using Point_function = PointFunction;

private:
  const Point_function m_func;
  const FT m_delta, m_den;

  const GeomTraits m_gt;

public:
  /**
   * \brief creates a new instance of this gradient.
   *
   * \param func the implicit function giving the value of the implicit function at each discretization point
   * \param delta the distance for calculating the finite differences
   * \param gt the geometric traits class
   */
  Finite_difference_gradient_3(const Point_function& func,
                               const FT delta = 0.001,
                               const Geom_traits& gt = Geom_traits())
    : m_func{func},
      m_delta{delta},
      m_den{FT{1} / (FT{2} * m_delta)},
      m_gt{gt}
  { }

  /**
   * \brief evaluates the gradient at a point in space.
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

    const FT gx = (m_func(p0) - m_func(p1)) * m_den;
    const FT gy = (m_func(p2) - m_func(p3)) * m_den;
    const FT gz = (m_func(p4) - m_func(p5)) * m_den;

    return vector(gx, gy, gz);
  }
};

} // namespace Isosurfacing
} // namespace CGAL

#endif // CGAL_ISOSURFACING_3_FINITE_DIFFERENCE_GRADIENT_3_H
