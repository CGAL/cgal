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
 * \cgalModels{IsosurfacingDomain_3,IsosurfacingDomainWithGradient_3}
 *
 * \brief A domain that represents an explicitly stored %Cartesian grid.
 *
 * \warning The domain keeps a pointer to the `grid` object, hence users must ensure that
 *          the lifetime of the `grid` object exceeds that of the object returned by this function.
 *
 * \tparam Grid must be a `CGAL::Isosurfacing::Cartesian_grid_3` whose `GeomTraits` template parameter
 *              is a model of `IsosurfacingTraits_3`.
 * \tparam Gradient the type of the gradient functor. It must be a model of `CopyConstructible`
 *                  and implement `%Grid::GeomTraits::Vector_3 operator()(const %Grid::GeomTraits::Point_3& point) const`.
 *
 * \sa `CGAL::Isosurfacing::create_explicit_Cartesian_grid_domain()`
 */
template <typename Grid, // to allow more than a Cartesian_grid_3
          typename Gradient = Zero_gradient
#ifndef DOXYGEN_RUNNING // Do not document Topology, Geometry, Function
          , typename Topology = internal::Grid_topology_3
          , typename Geometry = internal::Explicit_Cartesian_grid_geometry_3<Grid>
          , typename Function = internal::Explicit_Cartesian_grid_function<Grid>
#endif
          >
class Explicit_Cartesian_grid_domain_3
#ifndef DOXYGEN_RUNNING
  : public internal::Isosurfacing_domain_3<typename Grid::Geom_traits,
                                           Topology, Geometry, Function, Gradient>
#endif
{
private:
  using Base = internal::Isosurfacing_domain_3<typename Grid::Geom_traits,
                                               Topology, Geometry, Function, Gradient>;

public:
  /**
   * \brief creates a domain that can be used as input for isosurfacing algorithms.
   *
   * \param grid the %Cartesian grid containing input data
   * \param gradient a function giving the value of the gradient at each discretization point
   */
  Explicit_Cartesian_grid_domain_3(const Grid& grid,
                                   const Gradient& gradient = Gradient())
    : Base(Topology { grid.xdim(), grid.ydim(), grid.zdim() },
           Geometry { grid },
           Function { grid },
           gradient,
           grid.geom_traits())
  {
  }
};

/**
 * \ingroup IS_Domains_grp
 *
 * \brief creates a domain that can be used as input for isosurfacing algorithms.
 *
 * \warning The domain keeps a pointer to the `grid` object, hence users must ensure that
 *          the lifetime of the `grid` object exceeds that of the object returned by this function.
 *
 * \tparam Grid must be a `CGAL::Isosurfacing::Cartesian_grid_3` whose `GeomTraits` template parameter
 *              is a model of `IsosurfacingTraits_3`.
 * \tparam Gradient the type of the gradient functor. It must be a model of `CopyConstructible`
 *                  and implement `%Grid::Geom_traits::Vector_3 operator()(const GeomTraits::Point_3& point) const
 *
 * \param grid the %Cartesian grid containing input data
 * \param gradient a function giving the value of the gradient of the implicit function at each discretization point
 */
template <typename Grid, // allow passing more than just a Cartesian_grid_3
          typename Gradient = Zero_gradient>
Explicit_Cartesian_grid_domain_3<Grid, Gradient>
create_explicit_Cartesian_grid_domain(const Grid& grid,
                                      const Gradient& gradient = Gradient())
{
  return { grid, gradient };
}

} // namespace Isosurfacing
} // namespace CGAL

#endif // CGAL_ISOSURFACING_3_EXPLICIT_CARTESIAN_GRID_DOMAIN_3_H
