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

#ifndef CGAL_ISOSURFACING_3_EXPLICIT_CARTESIAN_GRID_DOMAIN_H
#define CGAL_ISOSURFACING_3_EXPLICIT_CARTESIAN_GRID_DOMAIN_H

#include <CGAL/license/Isosurfacing_3.h>

#include <CGAL/Isosurfacing_3/internal/Base_domain.h>
#include <CGAL/Isosurfacing_3/internal/Cartesian_grid_geometry.h>
#include <CGAL/Isosurfacing_3/internal/Grid_topology.h>

#include <CGAL/Isosurfacing_3/Cartesian_grid_3.h>
#include <CGAL/Isosurfacing_3/Zero_gradient.h>

namespace CGAL {
namespace Isosurfacing {

/**
 * \ingroup PkgIsosurfacing3Ref
 *
 * \cgalModels `IsosurfacingDomainWithGradient`
 *
 * \brief A domain that represents an explicitly stored Cartesian grid.
 *
 * \tparam GeomTraits must be a model of ``.
 * \tparam Gradient_ the type of the gradient functor. It must implement
 *                   `GeomTraits::Vector operator()(const GeomTraits::Point& point) const`.
 */
template <typename GeomTraits,
          typename Gradient_>
using Explicit_cartesian_grid_domain =
  internal::Base_domain<GeomTraits,
                        internal::Grid_topology,
                        internal::Cartesian_grid_geometry<GeomTraits>,
                        Cartesian_grid_3<GeomTraits>,
                        Gradient_>;

/**
 * \ingroup PkgIsosurfacing3Ref
 *
 * @todo maybe the grid should be a template and a concept
 *
 * \brief Creates a domain from a `Cartesian_grid_3` that can be used as input for isosurfacing algorithms.
 *
 * \tparam GeomTraits must be a model of ``.
 * \tparam Gradient_ the type of the gradient functor. It must implement
 *                   `GeomTraits::Vector operator()(const GeomTraits::Point& point) const`.
 *
 * \param grid the Cartesian grid containing input data
 * \param gradient a function that describes the gradient of the data
 *
 * \return a new `Explicit_cartesian_grid_domain`
 */
template <typename GeomTraits,
          typename Gradient_ = Zero_gradient>
Explicit_cartesian_grid_domain<GeomTraits, Gradient_>
create_explicit_cartesian_grid_domain(const Cartesian_grid_3<GeomTraits>& grid,
                                      const Gradient_& gradient = Gradient_())
{
  using Domain = Explicit_cartesian_grid_domain<GeomTraits, Gradient_>;

  using Topology = typename Domain::Topology;
  using Geometry = typename Domain::Geometry;
  using Function = typename Domain::Function;
  using Gradient = typename Domain::Gradient;

  const std::size_t size_i = grid.xdim();
  const std::size_t size_j = grid.ydim();
  const std::size_t size_k = grid.zdim();

  const Bbox_3& bbox = grid.get_bbox();
  const typename GeomTraits::Vector_3 offset(bbox.xmin(), bbox.ymin(), bbox.zmin());
  const typename GeomTraits::Vector_3 spacing = grid.get_spacing();

  const Topology topo { size_i, size_j, size_k };
  const Geometry geom { offset, spacing };
  const Function func { grid };
  const Gradient grad { gradient };

  return { topo, geom, func, grad };
}

} // namespace Isosurfacing
} // namespace CGAL

#endif // CGAL_ISOSURFACING_3_EXPLICIT_CARTESIAN_GRID_DOMAIN_H
