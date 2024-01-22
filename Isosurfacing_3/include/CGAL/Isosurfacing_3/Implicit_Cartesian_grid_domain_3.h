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
 * \cgalModels{IsosurfacingDomain_3,IsosurfacingDomainWithGradient_3}
 *
 * \brief A domain that represents a %Cartesian grid that discretizes an implicit function.
 *
 * \tparam GeomTraits must be a model of `IsosurfacingTraits_3`.
 * \tparam ImplicitFunction the type of the implicit function. It must be a model of `CopyConstructible`
 *                          and implement `GeomTraits::FT operator()(const GeomTraits::Point_3& point) const`.
 * \tparam Gradient the type of the gradient functor. It must be a model of `CopyConstructible` and implement
 *                  `GeomTraits::Vector_3 operator()(const GeomTraits::Point_3& point) const`.
 *
 * \sa `CGAL::Isosurfacing::create_implicit_Cartesian_grid_domain()`
 */
template <typename GeomTraits,
          typename ImplicitFunction,
          typename Gradient = Zero_gradient
#ifndef DOXYGEN_RUNNING // Do not document Topology, Geometry, Function
          , typename Topology = internal::Grid_topology_3
          , typename Geometry = internal::Implicit_Cartesian_grid_geometry_3<GeomTraits>
          , typename Function = internal::Implicit_function_with_geometry<Geometry, ImplicitFunction>
#endif
          >
class Implicit_Cartesian_grid_domain_3
#ifndef DOXYGEN_RUNNING
  : public internal::Isosurfacing_domain_3<GeomTraits, Topology, Geometry, Function, Gradient>
#endif
{
private:
  using Base = internal::Isosurfacing_domain_3<GeomTraits, Topology, Geometry, Function, Gradient>;

public:
  /**
   * \brief creates a domain from an implicit function.
   *
   * \details The implicit function is evaluated at the vertices of the virtual grid
   * defined by the bounding box and the spacing value. By not storing any function values explicitely,
   * less overall memory is required in comparison to an `Explicit_Cartesian_grid_domain_3`.
   *
   * \tparam GeomTraits must be a model of `IsosurfacingTraits_3`.
   * \tparam ImplicitFunction the type of the implicit function. It must be a model of `CopyConstructible`
   *                          and implement `GeomTraits::FT operator()(const GeomTraits::Point_3& point) const`.
   * \tparam Gradient the type of the gradient functor. It must be a model of `CopyConstructible` and implement
   *                  `GeomTraits::Vector_3 operator()(const GeomTraits::Point_3& point) const`.
   *
   * \param bbox a bounding box that specifies the dimensions of the implicit function's domain
   * \param spacing the distance between discretization points
   * \param point_function the implicit function giving the value of the implicit function at each discretization point
   * \param gradient a function giving the value of the gradient of the implicit function at each discretization point
   * \param gt an instance of geometric traits
   *
   * \pre `spacing != CGAL::NULL_VECTOR`
   */
  Implicit_Cartesian_grid_domain_3(const CGAL::Bbox_3& bbox,
                                   const typename GeomTraits::Vector_3& spacing,
                                   const ImplicitFunction& point_function,
                                   const Gradient& gradient = Gradient(),
                                   const GeomTraits& gt = GeomTraits())
    : Base(Topology { std::ceil(bbox.x_span() / gt.compute_x_3_object()(spacing)) + 1,
                      std::ceil(bbox.y_span() / gt.compute_y_3_object()(spacing)) + 1,
                      std::ceil(bbox.z_span() / gt.compute_z_3_object()(spacing)) + 1 },
           Geometry { typename GeomTraits::Vector_3 { bbox.xmin(), bbox.ymin(), bbox.zmin() }, spacing },
           Function { Geometry { typename GeomTraits::Vector_3{bbox.xmin(), bbox.ymin(), bbox.zmin()}, spacing },
                      point_function },
           gradient,
           gt)
  {
  }
};

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
 * \tparam ImplicitFunction the type of the implicit function. It must be a model of `CopyConstructible`
 *                          and implement `GeomTraits::FT operator()(const GeomTraits::Point_3& point) const`.
 * \tparam Gradient the type of the gradient functor. It must be a model of `CopyConstructible` and implement
 *                  `GeomTraits::Vector_3 operator()(const GeomTraits::Point_3& point) const`.
 *
 * \param bbox a bounding box that specifies the dimensions of the implicit function's domain
 * \param spacing the distance between discretization points
 * \param point_function the implicit function giving the value of the implicit function at each discretization point
 * \param gradient a function giving the value of the gradient of the implicit function at each discretization point
 * \param gt an instance of geometric traits
 *
 * \return a new instance of `CGAL::Isosurfacing::Implicit_Cartesian_grid_domain_3`
 *
 * \pre `spacing != CGAL::NULL_VECTOR`
 */
template <typename GeomTraits,
          typename ImplicitFunction,
          typename Gradient = Zero_gradient>
Implicit_Cartesian_grid_domain_3<GeomTraits, ImplicitFunction, Gradient>
create_implicit_Cartesian_grid_domain(const CGAL::Bbox_3& bbox,
                                      const typename GeomTraits::Vector_3& spacing,
                                      const ImplicitFunction& point_function,
                                      const Gradient& gradient = Gradient(),
                                      const GeomTraits& gt = GeomTraits())
{
  return { bbox, spacing, point_function, gradient, gt };
}

// @todo add an undocumented convenience overload with Vector_3<GeomTraits> to match CGAL kernels
// without having to provide the kernel in the call like f<kernel>(...)

} // namespace Isosurfacing
} // namespace CGAL

#endif // CGAL_ISOSURFACING_3_IMPLICIT_CARTESIAN_GRID_DOMAIN_3_H
