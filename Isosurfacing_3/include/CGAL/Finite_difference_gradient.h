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

#ifndef CGAL_FINITE_DIFFERENCE_GRADIENT_H
#define CGAL_FINITE_DIFFERENCE_GRADIENT_H

#include <CGAL/license/Isosurfacing_3.h>

namespace CGAL {
namespace Isosurfacing {

/**
 * \ingroup PkgIsosurfacing3Ref
 *
 * \brief Class template for a gradient that is calculated using finite differences.
 *
 * \details This gradient function evaluates an implicit value function at six points
 *          that are a given distance `delta` away from the queried point.
 *
 * \tparam GeomTraits the traits for this gradient.
 *
 * \tparam PointFunction the type of the implicit function. It must implement `GeomTraits::FT operator()(const
 * GeomTraits::Point& point) const`.
 */
template <typename GeomTraits,
          typename PointFunction>
class Finite_difference_gradient
{
public:
  using Geom_traits = GeomTraits;
  using FT = typename Geom_traits::FT;
  using Point = typename Geom_traits::Point_3;
  using Vector = typename Geom_traits::Vector_3;

public:
  /**
   * \ingroup PkgIsosurfacing3Ref
   *
   * \brief Create a new instance of this gradient.
   *
   * \param point_function the function with a point as argument
   * \param delta the distance for calculating the finite differences
   */
  Finite_difference_gradient(const PointFunction& point_function,
                             const FT delta = 0.001)
    : func(point_function),
      delta(delta)
  { }

  /**
   * \ingroup PkgIsosurfacing3Ref
   *
   * \brief Evaluate the gradient at a point in space.
   *
   * \param point the point at which the gradient is computed
   */
  Vector operator()(const Point& point) const
  {
    // compute the gradient by sampling the function with finite differences
    // at six points with distance delta around the query point
    const Point p0 = point + Vector(delta, 0, 0);
    const Point p1 = point - Vector(delta, 0, 0);
    const Point p2 = point + Vector(0, delta, 0);
    const Point p3 = point - Vector(0, delta, 0);
    const Point p4 = point + Vector(0, 0, delta);
    const Point p5 = point - Vector(0, 0, delta);

    const FT gx = (func(p0) - func(p1)) / (2 * delta);
    const FT gy = (func(p2) - func(p3)) / (2 * delta);
    const FT gz = (func(p4) - func(p5)) / (2 * delta);

    return Vector(gx, gy, gz);
  }

private:
  const PointFunction func;
  FT delta;
};

} // namespace Isosurfacing
} // namespace CGAL

#endif // CGAL_FINITE_DIFFERENCE_GRADIENT_H
