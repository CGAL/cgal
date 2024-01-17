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

#ifndef CGAL_ISOSURFACING_3_IMPLICIT_OCTREE_DOMAIN_H
#define CGAL_ISOSURFACING_3_IMPLICIT_OCTREE_DOMAIN_H

#include <CGAL/license/Isosurfacing_3.h>

#include <CGAL/Isosurfacing_3/internal/Isosurfacing_domain_3.h>
#include <CGAL/Isosurfacing_3/internal/Implicit_function_with_geometry.h>
#include <CGAL/Isosurfacing_3/internal/Octree_geometry.h>
#include <CGAL/Isosurfacing_3/internal/Octree_topology.h>
#include <CGAL/Isosurfacing_3/internal/Octree_wrapper.h>
#include <CGAL/Isosurfacing_3/Zero_gradient.h>

namespace CGAL {
namespace Isosurfacing {

/*
 * \ingroup IS_Domains_grp
 *
 * \cgalModels{IsosurfacingDomainWithGradient_3}
 *
 * \brief A domain that represents an octree that discretizes an implicit function.
 *
 * \tparam Octree must be a model of `...`.
 * \tparam ImplicitFunction the type of the implicit function. It must be a model of CopyConstructible implement
 *                          `GeomTraits::FT operator()(const GeomTraits::Point_3& point) const`.
 * \tparam Gradient the type of the gradient functor. It must be a model of `CopyConstructible` and implement
 *                  `GeomTraits::Vector_3 operator()(const GeomTraits::Point_3& point) const`.
 */
template <typename Octree,
          typename ImplicitFunction,
          typename Gradient = Zero_gradient,
          typename Topology = internal::Octree_topology<Octree>,
          typename Geometry = internal::Octree_geometry<Octree> >
using Implicit_octree_domain =
  internal::Isosurfacing_domain_3<typename Octree::Geom_traits,
                                  Topology,
                                  Geometry,
                                  internal::Implicit_function_with_geometry<Geometry, ImplicitFunction>,
                                  Gradient>;

/*
 * \ingroup IS_Domains_grp
 *
 * \brief creates a domain from an octree that can be used as input for isosurfacing algorithms.
 *
 * \details The implicit function will be evaluated on the octree vertices.
 *
 * \tparam GeomTraits must be a model of `Kernel`.
 * \tparam ImplicitFunction the type of the implicit function. It must be a model of `CopyConstructible` and implement
 *                          `GeomTraits::FT operator()(const GeomTraits::Point_3& point) const`.
 * \tparam Gradient the type of the gradient functor. It must be a model of `CopyConstructible` and implement
 *                  `GeomTraits::Vector_3 operator()(const GeomTraits::Point_3& point) const`.
 *
 * \param octree an octree
 * \param point_function the function with a point as argument
 * \param gradient a function that describes the gradient of the data
 *
 * \return a new `CGAL::Implicit_octree_domain`
 */
#ifdef DOXYGEN_RUNNING // allow more than the Octree_wrapper
template <typename GeomTraits,
          typename ImplicitFunction,
          typename Gradient = Zero_gradient>
auto create_implicit_octree_domain(const internal::Octree_wrapper<Octree<GeomTraits> >& octree,
                                   const ImplicitFunction& point_func,
                                   const Gradient& grad = Gradient())
{
  using Octree = internal::Octree_wrapper<Octree<GeomTraits> >;
#else
template <typename Octree,
          typename ImplicitFunction,
          typename Gradient = Zero_gradient>
auto create_implicit_octree_domain(const Octree& octree,
                                   const ImplicitFunction& point_func,
                                   const Gradient& gradient = Gradient())
{
#endif
  using Domain = Implicit_octree_domain<Octree, ImplicitFunction, Gradient>;

  using Topology = typename Domain::Topology;
  using Geometry = typename Domain::Geometry;
  using Function = typename Domain::Function;

  const Topology topo { octree };
  const Geometry geom { octree };
  const Function func { geom, point_func };

  // @fixme Octree_wrapper's geom_traits() isn't octree's geom_traits()...
  return Domain{ topo, geom, func, gradient, octree.geom_traits() };
}

} // namespace Isosurfacing
} // namespace CGAL

#endif // CGAL_ISOSURFACING_3_IMPLICIT_OCTREE_DOMAIN_H
