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

#ifndef CGAL_ISOSURFACING_3_EXPLICIT_CARTESIAN_GRID_DOMAIN_3_H
#define CGAL_ISOSURFACING_3_EXPLICIT_CARTESIAN_GRID_DOMAIN_3_H

#include <CGAL/license/Isosurfacing_3.h>

#include <CGAL/Isosurfacing_3/internal/Isosurfacing_domain_3.h>
#include <CGAL/Isosurfacing_3/internal/Explicit_Cartesian_grid_function.h>
#include <CGAL/Isosurfacing_3/internal/Explicit_Cartesian_grid_geometry_3.h>
#include <CGAL/Isosurfacing_3/internal/Grid_topology_3.h>
#include <CGAL/Isosurfacing_3/Zero_gradient.h>

namespace CGAL {
namespace Isosurfacing {

/**
 * \ingroup IS_Domains_grp
 *
 * \cgalModels{IsosurfacingDomainWithGradient_3}
 *
 * \brief A domain that represents an explicitly stored %Cartesian grid.
 *
 * \tparam GeomTraits must be a model of `IsosurfacingTraits_3`.
 * \tparam Gradient the type of the gradient functor. It must be a model of `CopyConstructible` and implement
 *                  `GeomTraits::Vector_3 operator()(const GeomTraits::Point_3& point) const`.
 *
 * \sa `CGAL::Isosurfacing::create_explicit_Cartesian_grid_domain()`
 */
#ifdef DOXYGEN_RUNNING // Allow more than a Cartesian_grid_3
template <template <typename GeomTraits> class CGAL::Isosurfacing::Cartesian_grid_3,
          typename Gradient = Zero_gradient>
using Explicit_Cartesian_grid_domain_3 = unspecified_type;
#else
template <typename Grid,
          typename Gradient = Zero_gradient,
          typename Topology = internal::Grid_topology_3,
          typename Geometry = internal::Explicit_Cartesian_grid_geometry_3<Grid>,
          typename Function = internal::Explicit_Cartesian_grid_function<Grid> >
using Explicit_Cartesian_grid_domain_3 =
  internal::Isosurfacing_domain_3<typename Grid::Geom_traits,
                                  Topology,
                                  Geometry,
                                  Function,
                                  Gradient>;
#endif

/**
 * \ingroup IS_Domains_grp
 *
 * \brief Creates a domain that can be used as input for isosurfacing algorithms.
 *
 * \warning As the domain will keep a pointer to the `grid` object, users must ensure that
 *          the lifetime of the `grid` object exceeds that of the object returned by this function.
 *
 * \tparam GeomTraits must be a model of `IsosurfacingTraits_3`.
 * \tparam Gradient the type of the gradient functor. It must be a model of `CopyConstructible` and implement
 *                  `GeomTraits::Vector_3 operator()(const GeomTraits::Point_3& point) const`.
 *
 * \param grid the %Cartesian grid containing input data
 * \param grad a function giving the value of the gradient of the implicit function at each discretization point
 *
 * \return a new `CGAL::Explicit_Cartesian_grid_domain_3`
 */
#ifdef DOXYGEN_RUNNING // Allow more than Cartesian_grid_3
template <typename GeomTraits,
          typename Gradient = Zero_gradient>
Explicit_Cartesian_grid_domain_3<Grid, Gradient>
create_explicit_Cartesian_grid_domain(const CGAL::Isosurfacing::Cartesian_grid_3<GeomTraits>& grid,
                                      const Gradient& grad = Gradient())
#else
// Actual code enables passing more than just a Cartesian_grid_3
template <typename Grid,
          typename Gradient = Zero_gradient>
Explicit_Cartesian_grid_domain_3<Grid, Gradient>
create_explicit_Cartesian_grid_domain(const Grid& grid,
                                      const Gradient& grad = Gradient())
#endif
{
  using Domain = Explicit_Cartesian_grid_domain_3<Grid, Gradient>;

  using Topology = typename Domain::Topology;
  using Geometry = typename Domain::Geometry;
  using Function = typename Domain::Function;

  const std::size_t size_i = grid.xdim();
  const std::size_t size_j = grid.ydim();
  const std::size_t size_k = grid.zdim();

  const Topology topo { size_i, size_j, size_k };
  const Geometry geom { grid };
  const Function func { grid };

  return Domain{ topo, geom, func, grad, grid.geom_traits() };
}

} // namespace Isosurfacing
} // namespace CGAL

#endif // CGAL_ISOSURFACING_3_EXPLICIT_CARTESIAN_GRID_DOMAIN_3_H
