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

#ifndef CGAL_EXPLICIT_CARTESIAN_GRID_DOMAIN_H
#define CGAL_EXPLICIT_CARTESIAN_GRID_DOMAIN_H

#include <CGAL/license/Isosurfacing_3.h>

#include <CGAL/Cartesian_grid_3.h>
#include <CGAL/Isosurfacing_3/internal/Base_domain.h>
#include <CGAL/Isosurfacing_3/internal/Cartesian_grid_geometry.h>
#include <CGAL/Isosurfacing_3/internal/Grid_topology.h>
#include <CGAL/Zero_gradient.h>

namespace CGAL {
namespace Isosurfacing {

/**
 * \ingroup PkgIsosurfacing3Ref
 *
 * \brief A domain that respesents an explicitly stored cartesian grid. It is a model of the concept
 * `IsosurfacingDomainWithGradient`.
 *
 * \tparam GeomTraits the traits type
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
 * \brief Creates a domain from a `Cartesian_grid_3` that can be used as input for isosurfacing algorithms.
 *
 * \details
 *
 * \tparam GeomTraits the traits type
 *
 * \param grid the %Cartesian grid containing input data
 * \param gradient a function that describes the gradient of the data
 *
 * \return a new `Explicit_cartesian_grid_domain`
 */
template <typename GeomTraits,
          typename Gradient_ = Zero_gradient<GeomTraits> >
Explicit_cartesian_grid_domain<GeomTraits, Gradient_>
create_explicit_cartesian_grid_domain(const std::shared_ptr<Cartesian_grid_3<GeomTraits> > grid,
                                      const Gradient_& gradient = Gradient_())
{
  using Domain = Explicit_cartesian_grid_domain<GeomTraits, Gradient_>;

  using Topology = typename Domain::Topology ;
  using Geometry = typename Domain::Geometry;
  using Function = typename Domain::Function;
  using Gradient = typename Domain::Gradient;

  const std::size_t size_i = grid->xdim();
  const std::size_t size_j = grid->ydim();
  const std::size_t size_k = grid->zdim();

  const Bbox_3& bbox = grid->get_bbox();
  const typename GeomTraits::Vector_3 offset(bbox.xmin(), bbox.ymin(), bbox.zmin());
  const typename GeomTraits::Vector_3 spacing = grid->get_spacing();

  // create copies as shared_ptr for safe memory management
  const Topology topo = std::make_shared<Topology::element_type>(size_i, size_j, size_k);
  const Geometry geom = std::make_shared<Geometry::element_type>(offset, spacing);
  const Function func = grid;
  const Gradient grad = std::make_shared<Gradient::element_type>(gradient);

  return Domain(topo, geom, func, grad);
}

} // namespace Isosurfacing
} // namespace CGAL

#endif // CGAL_EXPLICIT_CARTESIAN_GRID_DOMAIN_H
