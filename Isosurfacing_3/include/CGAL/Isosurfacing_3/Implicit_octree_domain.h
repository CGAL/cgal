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

#ifndef CGAL_ISOSURFACING_3_IMPLICIT_OCTREE_DOMAIN_H
#define CGAL_ISOSURFACING_3_IMPLICIT_OCTREE_DOMAIN_H

#include <CGAL/license/Isosurfacing_3.h>

#include <CGAL/Isosurfacing_3/internal/Base_domain.h>
#include <CGAL/Isosurfacing_3/internal/Implicit_function_with_geometry.h>
#include <CGAL/Isosurfacing_3/internal/Octree_geometry.h>
#include <CGAL/Isosurfacing_3/internal/Octree_topology.h>
#include <CGAL/Isosurfacing_3/internal/Octree_wrapper.h>
#include <CGAL/Isosurfacing_3/Zero_gradient.h>

namespace CGAL {
namespace Isosurfacing {

/*
 * \ingroup PkgIsosurfacing3Ref
 *
 * \cgalModels `IsosurfacingDomainWithGradient`
 *
 * \brief A domain that respesents an octree that discretizes an implicit function.
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
using Implicit_octree_domain =
  internal::Base_domain<GeomTraits,
                        internal::Octree_topology<GeomTraits>,
                        internal::Octree_geometry<GeomTraits>,
                        internal::Implicit_function_with_geometry<GeomTraits,
                                                                  internal::Octree_geometry<GeomTraits>,
                                                                  PointFunction>,
                        Gradient_>;

/*
 * \ingroup PkgIsosurfacing3Ref
 *
 * \brief creates a domain from an octree that can be used as input for isosurfacing algorithms.
 *
 * \details The implicit function will be evaluated on the octree vertices.
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
Implicit_octree_domain<GeomTraits, PointFunction, Gradient_>
create_implicit_octree_domain(const internal::Octree_wrapper<GeomTraits>& octree,
                              const PointFunction& point_function,
                              const Gradient_& gradient = Gradient_())
{
  using Domain = Implicit_octree_domain<GeomTraits, PointFunction, Gradient_>;

  using Topology = typename Domain::Topology;
  using Geometry = typename Domain::Geometry;
  using Function = typename Domain::Function;
  using Gradient = typename Domain::Gradient;
  using Point_function = PointFunction;
  using Octree = internal::Octree_wrapper<GeomTraits>;

  const Topology topo { octree };
  const Geometry geom { octree };
  const Point_function point_func { point_function };
  const Function func { geom, point_func };
  const Gradient grad { gradient };

  return { topo, geom, func, grad };
}

} // namespace Isosurfacing
} // namespace CGAL

#endif // CGAL_ISOSURFACING_3_IMPLICIT_OCTREE_DOMAIN_H
