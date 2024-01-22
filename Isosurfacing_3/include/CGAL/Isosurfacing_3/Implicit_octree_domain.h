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
 * \warning The domain keeps a pointer to the `octree` object, hence users must ensure that
 *          the lifetime of the `octree` object exceeds that of this object.
 *
 * \tparam Octree must be a `CGAL::Octree<GeomTraits>`.
 * \tparam ImplicitFunction the type of the implicit function. It must be a model of CopyConstructible implement
 *                          `GeomTraits::FT operator()(const GeomTraits::Point_3& point) const`.
 * \tparam Gradient the type of the gradient functor. It must be a model of `CopyConstructible` and implement
 *                  `GeomTraits::Vector_3 operator()(const GeomTraits::Point_3& point) const`.
 */
template <typename Octree,
          typename ImplicitFunction,
          typename Gradient = Zero_gradient
#ifndef DOXYGEN_RUNNING // Do not document Topology, Geometry, Function
          , typename Topology = internal::Octree_topology<Octree>
          , typename Geometry = internal::Octree_geometry<Octree>
          , typename Function = internal::Implicit_function_with_geometry<Geometry, ImplicitFunction>
#endif
          >
class Implicit_octree_domain
#ifndef DOXYGEN_RUNNING
  : public internal::Isosurfacing_domain_3<typename Octree::Geom_traits,
                                           Topology, Geometry, Function, Gradient>
#endif
{
  using Base = internal::Isosurfacing_domain_3<typename Octree::Geom_traits,
                                               Topology, Geometry, Function, Gradient>;

public:
  Implicit_octree_domain(const Octree& octree,
                         const ImplicitFunction& point_function,
                         const Gradient& gradient = Gradient())
    : Base(Topology { octree },
           Geometry { octree },
           Function { Geometry { octree }, point_function },
           gradient,
           octree.geom_traits())
  {
  }
};

/*
 * \ingroup IS_Domains_grp
 *
 * \brief creates a domain from an octree that can be used as input for isosurfacing algorithms.
 *
 * \warning The domain keeps a pointer to the `octree` object, hence users must ensure that
 *          the lifetime of the `octree` object exceeds that of the object returned by this function.
 *
 * \tparam Octree must be a `CGAL::Octree<GeomTraits>`.
 * \tparam ImplicitFunction the type of the implicit function. It must be a model of `CopyConstructible`
 *                          and implement `GeomTraits::FT operator()(const GeomTraits::Point_3& point) const`.
 * \tparam Gradient the type of the gradient functor. It must be a model of `CopyConstructible` and implement
 *                  `Octree::GeomTraits::Vector_3 operator()(const GeomTraits::Point_3& point) const`.
 *
 * \param octree an octree
 * \param point_function the function with a point as argument
 * \param gradient a function that describes the gradient of the data
 *
 * \return a new `CGAL::Implicit_octree_domain`
 */
template <typename Octree,
          typename ImplicitFunction,
          typename Gradient = Zero_gradient>
Implicit_octree_domain<Octree, ImplicitFunction, Gradient>
create_implicit_octree_domain(const Octree& octree,
                              const ImplicitFunction& point_function,
                              const Gradient& gradient = Gradient())
{
  return { octree, point_function, gradient };
}

} // namespace Isosurfacing
} // namespace CGAL

#endif // CGAL_ISOSURFACING_3_IMPLICIT_OCTREE_DOMAIN_H
