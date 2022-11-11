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

#ifndef CGAL_IMPLICIT_CARTESIAN_GRID_DOMAIN_H
#define CGAL_IMPLICIT_CARTESIAN_GRID_DOMAIN_H

#include <CGAL/Default_gradients.h>
#include <CGAL/Isosurfacing_3/internal/Base_domain.h>
#include <CGAL/Isosurfacing_3/internal/Cartesian_grid_geometry.h>
#include <CGAL/Isosurfacing_3/internal/Grid_topology.h>
#include <CGAL/Isosurfacing_3/internal/Implicit_function_with_geometry.h>
#include <CGAL/license/Isosurfacing_3.h>

namespace CGAL {
namespace Isosurfacing {

template <class GeomTraits, typename PointFunction, typename Gradient_>
using Implicit_cartesian_grid_domain =
    Base_domain<GeomTraits, Grid_topology, Cartesian_grid_geometry<GeomTraits>,
                Implicit_function_with_geometry<GeomTraits, Cartesian_grid_geometry<GeomTraits>, PointFunction>,
                Gradient_>;

/**
 * \ingroup PkgIsosurfacing3Ref
 *
 * \brief Creates a domain from a Cartesian_grid_3 that can be used as input for isosurfacing algorithms.
 *
 * \details
 *
 * \tparam GeomTraits the traits type
 * \tparam PointFunction the type of the implicit function
 *
 * \param bbox a bounding box that specifies the size of the functions domain
 * \param spacing the distance between discretized points on the function
 * \param point_function the function with a point as argument
 * \param gradient a function that describes the gradient of the data
 */
template <class GeomTraits, typename PointFunction, typename Gradient_ = Zero_gradient<GeomTraits>>
Implicit_cartesian_grid_domain<GeomTraits, PointFunction, Gradient_> create_implicit_cartesian_grid_domain(
    const Bbox_3& bbox, const typename GeomTraits::Vector_3& spacing, const PointFunction& point_function,
    const Gradient_& gradient = Gradient_()) {

    typedef Implicit_cartesian_grid_domain<GeomTraits, PointFunction, Gradient_> Domain;
    typedef typename Domain::Topology Topology;
    typedef typename Domain::Geometry Geometry;
    typedef typename Domain::Function Function;
    typedef typename Domain::Gradient Gradient;
    typedef typename Function::element_type::Point_function Point_function;

    const std::size_t size_i = std::ceil(bbox.x_span() / spacing.x()) + 1;
    const std::size_t size_j = std::ceil(bbox.y_span() / spacing.y()) + 1;
    const std::size_t size_k = std::ceil(bbox.z_span() / spacing.z()) + 1;

    const typename GeomTraits::Vector_3 offset(bbox.xmin(), bbox.ymin(), bbox.zmin());

    const Topology topo = std::make_shared<Topology::element_type>(size_i, size_j, size_k);
    const Geometry geom = std::make_shared<Geometry::element_type>(offset, spacing);
    const Point_function point_func = std::make_shared<Point_function::element_type>(point_function);
    const Function func = std::make_shared<Function::element_type>(geom, point_func);
    const Gradient grad = std::make_shared<Gradient::element_type>(gradient);

    return Domain(topo, geom, func, grad);
}

}  // namespace Isosurfacing
}  // namespace CGAL

#endif  // CGAL_IMPLICIT_CARTESIAN_GRID_DOMAIN_H
