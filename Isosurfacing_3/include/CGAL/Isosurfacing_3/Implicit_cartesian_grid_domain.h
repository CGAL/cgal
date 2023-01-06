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
 * \brief A domain that respesents a cartesian grid that discretizes an implicit function. It is a model of the concept
 * `IsosurfacingDomainWithGradient`.
 *
 * \tparam GeomTraits the traits type
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
 * \brief Creates a domain from an implicit function that can be used as input for isosurfacing algorithms.
 *
 * \details The implicit function will be evaluated on the grid points of the virtual grid
 * defined by the bounding box and spacing. By not storing any function values implicitly
 * less memory accesses are required in comparison to an `Explicit_cartesian_grid_domain`.
 *
 * \tparam GeomTraits the traits type
 * \tparam PointFunction the type of the implicit function. It must implement `GeomTraits::FT operator()(const
 * GeomTraits::Point& point) const`.
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
          typename Gradient_ = Zero_gradient<GeomTraits> >
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
  using Point_function = typename Function::element_type::Point_function;

  // calculate grid dimensions
  const std::size_t size_i = std::ceil(bbox.x_span() / spacing.x()) + 1;
  const std::size_t size_j = std::ceil(bbox.y_span() / spacing.y()) + 1;
  const std::size_t size_k = std::ceil(bbox.z_span() / spacing.z()) + 1;

  const typename GeomTraits::Vector_3 offset(bbox.xmin(), bbox.ymin(), bbox.zmin());

  // create copies as shared_ptr for safe memory management
  const Topology topo = std::make_shared<Topology::element_type>(size_i, size_j, size_k);
  const Geometry geom = std::make_shared<Geometry::element_type>(offset, spacing);
  const Point_function point_func = std::make_shared<Point_function::element_type>(point_function);
  const Function func = std::make_shared<Function::element_type>(geom, point_func);
  const Gradient grad = std::make_shared<Gradient::element_type>(gradient);

  return Domain(topo, geom, func, grad);
}

} // namespace Isosurfacing
} // namespace CGAL

#endif // CGAL_ISOSURFACING_3_IMPLICIT_CARTESIAN_GRID_DOMAIN_H
