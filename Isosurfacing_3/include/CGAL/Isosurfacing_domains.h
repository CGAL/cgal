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

#ifndef CGAL_ISOSURFACING_DOMAINS_H
#define CGAL_ISOSURFACING_DOMAINS_H

#include <CGAL/Bbox_3.h>
#include <CGAL/Cartesian_grid_3.h>
#include <CGAL/Isosurfacing_3/internal/Explicit_cartesian_grid_domain.h>
#include <CGAL/Isosurfacing_3/internal/Implicit_cartesian_grid_domain.h>
#include <CGAL/Isosurfacing_3/internal/Implicit_octree_domain.h>
#include <CGAL/Octree_wrapper.h>
#include <CGAL/license/Isosurfacing_3.h>

namespace CGAL {
namespace Isosurfacing {

template <class GeomTraits, typename Function, typename Gradient>
Implicit_cartesian_grid_domain_with_gradient<GeomTraits, Function, Gradient> create_implicit_cartesian_grid_domain(
    const Bbox_3& bbox, const typename GeomTraits::Vector_3& spacing, const Function& func, const Gradient& grad) {

    const std::size_t size_i = std::ceil(bbox.x_span() / spacing.x()) + 1;
    const std::size_t size_j = std::ceil(bbox.y_span() / spacing.y()) + 1;
    const std::size_t size_k = std::ceil(bbox.z_span() / spacing.z()) + 1;

    const typename GeomTraits::Vector_3 offset(bbox.xmin(), bbox.ymin(), bbox.zmin());

    return {size_i, size_j, size_k, offset, spacing, func, grad};
}

template <class GeomTraits, typename Function>
Implicit_cartesian_grid_domain<GeomTraits, Function> create_implicit_cartesian_grid_domain(
    const Bbox_3& bbox, const typename GeomTraits::Vector_3& spacing, const Function& func) {

    const std::size_t size_i = std::ceil(bbox.x_span() / spacing.x()) + 1;
    const std::size_t size_j = std::ceil(bbox.y_span() / spacing.y()) + 1;
    const std::size_t size_k = std::ceil(bbox.z_span() / spacing.z()) + 1;

    const typename GeomTraits::Vector_3 offset(bbox.xmin(), bbox.ymin(), bbox.zmin());

    return {size_i, size_j, size_k, offset, spacing, func};
}

template <class GeomTraits, typename Gradient>
Explicit_cartesian_grid_domain_with_gradient<GeomTraits, Gradient> create_explicit_cartesian_grid_domain(
    const Cartesian_grid_3<GeomTraits>& grid, const Gradient& grad) {

    const std::size_t size_i = grid.xdim();
    const std::size_t size_j = grid.ydim();
    const std::size_t size_k = grid.zdim();

    const Bbox_3& bbox = grid.get_bbox();
    const typename GeomTraits::Vector_3 offset(bbox.xmin(), bbox.ymin(), bbox.zmin());
    const typename GeomTraits::Vector_3 spacing = grid.get_spacing();

    return {size_i, size_j, size_k, offset, spacing, grid, grad};
}

template <class GeomTraits>
Explicit_cartesian_grid_domain<GeomTraits> create_explicit_cartesian_grid_domain(
    const Cartesian_grid_3<GeomTraits>& grid) {

    const std::size_t size_i = grid.xdim();
    const std::size_t size_j = grid.ydim();
    const std::size_t size_k = grid.zdim();

    const Bbox_3& bbox = grid.get_bbox();
    const typename GeomTraits::Vector_3 offset(bbox.xmin(), bbox.ymin(), bbox.zmin());
    const typename GeomTraits::Vector_3 spacing = grid.get_spacing();

    return {size_i, size_j, size_k, offset, spacing, grid};
}

template <class GeomTraits, typename Function, typename Gradient>
Implicit_octree_domain_with_gradient<GeomTraits, Function, Gradient> create_implicit_octree_domain(
    const Octree_wrapper<GeomTraits>& octree, const Function& func, const Gradient& grad) {

    return {octree, func, grad};
}

template <class GeomTraits, typename Function>
Implicit_octree_domain<GeomTraits, Function> create_implicit_octree_domain(const Octree_wrapper<GeomTraits>& octree,
                                                                           const Function& func) {

    return {octree, func};
}

}  // namespace Isosurfacing
}  // namespace CGAL

#endif  // CGAL_ISOSURFACING_DOMAINS_H
