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

#ifndef CGAL_ISOSURFACING_3_IMPLICIT_CARTESIAN_GRID_DOMAIN_3_H
#define CGAL_ISOSURFACING_3_IMPLICIT_CARTESIAN_GRID_DOMAIN_3_H

#include <CGAL/license/Isosurfacing_3.h>

#include <CGAL/Isosurfacing_3/internal/Isosurfacing_domain_3.h>
#include <CGAL/Isosurfacing_3/internal/Grid_topology_3.h>
#include <CGAL/Isosurfacing_3/internal/Implicit_Cartesian_grid_geometry_3.h>
#include <CGAL/Isosurfacing_3/internal/Implicit_function_with_geometry.h>
#include <CGAL/Isosurfacing_3/Zero_gradient.h>

#include <CGAL/assertions.h>
#include <CGAL/Bbox_3.h>

namespace CGAL {
namespace Isosurfacing {

/**
 * \ingroup IS_Domains_grp
 *
 * \cgalModels `IsosurfacingDomainWithGradient_3`
 *
 * \brief A domain that represents a %Cartesian grid that discretizes an implicit function.
 *
 * \tparam GeomTraits must be a model of `IsosurfacingTraits_3`.
 * \tparam ImplicitFunction the type of the implicit function. It must be a model of `CopyConstructible` and implement
 *                          `GeomTraits::FT operator()(const GeomTraits::Point_3& point) const`.
 * \tparam Gradient the type of the gradient functor. It must be a model of `CopyConstructible` and implement
 *                  `GeomTraits::Vector_3 operator()(const GeomTraits::Point_3& point) const`.
 *
 * \sa `CGAL::Isosurfacing::create_implicit_Cartesian_grid_domain()`
 */
#ifdef DOXYGEN_RUNNING // Otherwise it shows what is behind "using" in the doc...
template <typename GeomTraits,
          typename ImplicitFunction,
          typename Gradient = Zero_gradient>
using Implicit_Cartesian_grid_domain_3 = unspecified_type;
#else
template <typename GeomTraits,
          typename ImplicitFunction,
          typename Gradient = Zero_gradient,
          typename Topology = internal::Grid_topology_3,
          typename Geometry = internal::Implicit_Cartesian_grid_geometry_3<GeomTraits>,
          typename Function = internal::Implicit_function_with_geometry<Geometry, ImplicitFunction> >
using Implicit_Cartesian_grid_domain_3 =
  internal::Isosurfacing_domain_3<GeomTraits,
                                  Topology,
                                  Geometry,
                                  Function,
                                  Gradient>;
#endif

/**
 * \ingroup IS_Domains_grp
 *
 * \brief creates a domain from an implicit function that can be used as input for isosurfacing algorithms.
 *
 * \details The implicit function is evaluated at the vertices of the virtual grid
 * defined by the bounding box and the spacing value. By not storing any function values explicitely,
 * less overall memory is required in comparison to an `Explicit_Cartesian_grid_domain_3`.
 *
 * \tparam GeomTraits must be a model of `IsosurfacingTraits_3`.
 * \tparam ImplicitFunction the type of the implicit function. It must be a model of `CopyConstructible` and implement
 *                          `GeomTraits::FT operator()(const GeomTraits::Point_3& point) const`.
 * \tparam Gradient the type of the gradient functor. It must be a model of `CopyConstructible` and implement
 *                  `GeomTraits::Vector_3 operator()(const GeomTraits::Point_3& point) const`.
 *
 * \param bbox a bounding box that specifies the dimensions of the implicit function's domain
 * \param spacing the distance between discretization points
 * \param point_function the implicit function giving the value of the implicit function at each discretization point
 * \param grad a function giving the value of the gradient of the implicit function at each discretization point
 * \param gt an instance of geometric traits
 *
 * \return a new object of type `CGAL::Implicit_Cartesian_grid_domain_3`
 *
 * \pre `spacing != CGAL::NULL_VECTOR`
 */
template <typename GeomTraits,
          typename ImplicitFunction,
          typename Gradient = Zero_gradient>
Implicit_Cartesian_grid_domain_3<GeomTraits, ImplicitFunction, Gradient>
create_implicit_Cartesian_grid_domain(const Bbox_3& bbox,
                                      const Vector_3<GeomTraits>& spacing,
                                      const ImplicitFunction& point_function,
                                      const Gradient& grad = Gradient(),
                                      const GeomTraits& gt = GeomTraits())
{
  using Domain = Implicit_Cartesian_grid_domain_3<GeomTraits, ImplicitFunction, Gradient>;

  using Topology = typename Domain::Topology;
  using Geometry = typename Domain::Geometry;
  using Function = typename Domain::Function;

  typename GeomTraits::Compute_x_3 x_coord = gt.compute_x_3_object();
  typename GeomTraits::Compute_y_3 y_coord = gt.compute_y_3_object();
  typename GeomTraits::Compute_z_3 z_coord = gt.compute_z_3_object();

  // calculate grid dimensions
  const std::size_t size_i = std::ceil(bbox.x_span() / x_coord(spacing)) + 1;
  const std::size_t size_j = std::ceil(bbox.y_span() / y_coord(spacing)) + 1;
  const std::size_t size_k = std::ceil(bbox.z_span() / z_coord(spacing)) + 1;

  CGAL_precondition(size_i > 0 && size_j > 0 && size_k > 0);

  // @fixme recompute the spacing?

  const typename GeomTraits::Vector_3 offset{bbox.xmin(), bbox.ymin(), bbox.zmin()};

  const Topology topo { size_i, size_j, size_k };
  const Geometry geom { offset, spacing };
  const Function func { geom, point_function };

  return Domain{ topo, geom, func, grad, gt };
}

} // namespace Isosurfacing
} // namespace CGAL

#endif // CGAL_ISOSURFACING_3_IMPLICIT_CARTESIAN_GRID_DOMAIN_3_H
