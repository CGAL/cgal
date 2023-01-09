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

#ifndef CGAL_ISOSURFACING_3_IMPLICIT_CARTESIAN_GRID_DOMAIN_H
#define CGAL_ISOSURFACING_3_IMPLICIT_CARTESIAN_GRID_DOMAIN_H

#include <CGAL/license/Isosurfacing_3.h>

#include <CGAL/Isosurfacing_3/internal/Base_domain.h>
#include <CGAL/Isosurfacing_3/internal/Cartesian_grid_geometry.h>
#include <CGAL/Isosurfacing_3/internal/Grid_topology.h>
#include <CGAL/Isosurfacing_3/internal/Implicit_function_with_geometry.h>
#include <CGAL/Isosurfacing_3/Zero_gradient.h>

#include <CGAL/Bbox_3.h>

namespace CGAL {
namespace Isosurfacing {

/**
 * \ingroup PkgIsosurfacing3Ref
 *
 * \cgalModels `IsosurfacingDomainWithGradient`
 *
 * \brief A domain that represents a Cartesian grid that discretizes an implicit function.
 *
 * \tparam GeomTraits must be a model of ``.
 * \tparam PointFunction the type of the implicit function. It must implement
 *                       `GeomTraits::FT operator()(const GeomTraits::Point& point) const`.
 * \tparam Gradient_ the type of the gradient functor. It must implement
 *                   `GeomTraits::Vector operator()(const GeomTraits::Point& point) const`.
 */
template <typename GeomTraits,
          typename PointFunction,
          typename Gradient_>
using Implicit_cartesian_grid_domain =
  internal::Base_domain<GeomTraits,
                        internal::Grid_topology,
                        internal::Cartesian_grid_geometry<GeomTraits>,
                        internal::Implicit_function_with_geometry<GeomTraits,
                                                                  internal::Cartesian_grid_geometry<GeomTraits>,
                                                                  PointFunction>,
                        Gradient_>;

/**
 * \ingroup PkgIsosurfacing3Ref
 *
 * \brief creates a domain from an implicit function that can be used as input for isosurfacing algorithms.
 *
 * \details The implicit function will be evaluated on the grid vertices of the virtual grid
 * defined by the bounding box and spacing. By not storing any function values implicitly,
 * fewer memory accesses are required in comparison to an `Explicit_cartesian_grid_domain`.
 *
 * \tparam GeomTraits must be a model of ``.
 * \tparam PointFunction the type of the implicit function. It must implement
 *                       `GeomTraits::FT operator()(const GeomTraits::Point& point) const`.
 * \tparam Gradient_ the type of the gradient functor. It must implement
 *                   `GeomTraits::Vector operator()(const GeomTraits::Point& point) const`.
 *
 * \param bbox a bounding box that specifies the size of the functions domain
 * \param spacing the distance between discretized points on the function
 * \param point_function the function with a point as argument
 * \param gradient a function that describes the gradient of the data
 *
 * \return a new `Implicit_cartesian_grid_domain`
 */
template <typename GeomTraits,
          typename PointFunction,
          typename Gradient_ = Zero_gradient>
Implicit_cartesian_grid_domain<GeomTraits, PointFunction, Gradient_>
create_implicit_cartesian_grid_domain(const Bbox_3& bbox,
                                      const typename GeomTraits::Vector_3& spacing,
                                      const PointFunction& point_function,
                                      const Gradient_& gradient = Gradient_())
{
  using Domain = Implicit_cartesian_grid_domain<GeomTraits, PointFunction, Gradient_>;

  using Topology = typename Domain::Topology;
  using Geometry = typename Domain::Geometry;
  using Function = typename Domain::Function;
  using Gradient = typename Domain::Gradient;
  using Point_function = PointFunction;

  // calculate grid dimensions
  const std::size_t size_i = std::ceil(bbox.x_span() / spacing.x()) + 1;
  const std::size_t size_j = std::ceil(bbox.y_span() / spacing.y()) + 1;
  const std::size_t size_k = std::ceil(bbox.z_span() / spacing.z()) + 1;

  const typename GeomTraits::Vector_3 offset(bbox.xmin(), bbox.ymin(), bbox.zmin());

  const Topology topo { size_i, size_j, size_k };
  const Geometry geom { offset, spacing };
  const Point_function point_func { point_function };
  const Function func { geom, point_func };
  const Gradient grad { gradient };

  return { topo, geom, func, grad };
}

} // namespace Isosurfacing
} // namespace CGAL

#endif // CGAL_ISOSURFACING_3_IMPLICIT_CARTESIAN_GRID_DOMAIN_H
