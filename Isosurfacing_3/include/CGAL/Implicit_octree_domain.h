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

#ifndef CGAL_IMPLICIT_OCTREE_DOMAIN_H
#define CGAL_IMPLICIT_OCTREE_DOMAIN_H

#include <CGAL/Default_gradients.h>
#include <CGAL/Isosurfacing_3/internal/Base_domain.h>
#include <CGAL/Isosurfacing_3/internal/Implicit_function_with_geometry.h>
#include <CGAL/Isosurfacing_3/internal/Octree_geometry.h>
#include <CGAL/Isosurfacing_3/internal/Octree_topology.h>
#include <CGAL/Octree_wrapper.h>
#include <CGAL/license/Isosurfacing_3.h>

namespace CGAL {
namespace Isosurfacing {


template <class GeomTraits, typename PointFunction, typename Gradient_>
using Implicit_octree_domain =
    Base_domain<GeomTraits, Octree_topology<GeomTraits>, Octree_geometry<GeomTraits>,
                Implicit_function_with_geometry<GeomTraits, Octree_geometry<GeomTraits>, PointFunction>, Gradient_>;

template <class GeomTraits, typename PointFunction, typename Gradient_ = Zero_gradient<GeomTraits>>
Implicit_octree_domain<GeomTraits, PointFunction, Gradient_> create_implicit_octree_domain(
    const Octree_wrapper<GeomTraits>& octree, const PointFunction& point_function,
    const Gradient_& gradient = Gradient_()) {

    typedef Implicit_octree_domain<GeomTraits, PointFunction, Gradient_> Domain;
    typedef typename Domain::Topology Topology;
    typedef typename Domain::Geometry Geometry;
    typedef typename Domain::Function Function;
    typedef typename Domain::Gradient Gradient;
    typedef typename Function::element_type::Point_function Point_function;
    typedef typename Topology::element_type::Octree Octree;

    const Octree oct = std::make_shared<Octree::element_type>(octree);
    const Topology topo = std::make_shared<Topology::element_type>(oct);
    const Geometry geom = std::make_shared<Geometry::element_type>(oct);
    const Point_function point_func = std::make_shared<Point_function::element_type>(point_function);
    const Function func = std::make_shared<Function::element_type>(geom, point_func);
    const Gradient grad = std::make_shared<Gradient::element_type>(gradient);

    return Domain(topo, geom, func, grad);
}

}  // namespace Isosurfacing
}  // namespace CGAL

#endif  // CGAL_IMPLICIT_OCTREE_DOMAIN_H
