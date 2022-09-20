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

template <class GeomTraits, typename Function_, typename Gradient_>
class Implicit_octree_domain_with_gradient
    : public Base_domain<GeomTraits, Octree_topology<GeomTraits>, Octree_geometry<GeomTraits>,
                         Implicit_function_with_geometry<GeomTraits, Octree_geometry<GeomTraits>, Function_>,
                         Gradient_> {
public:
    typedef GeomTraits Geom_traits;

    typedef Octree_wrapper<Geom_traits> Octree;

    typedef Octree_topology<Geom_traits> Topology;
    typedef Octree_geometry<Geom_traits> Geometry;
    typedef Function_ Function_with_point;
    typedef Implicit_function_with_geometry<Geom_traits, Geometry, Function_with_point> Function;
    typedef Gradient_ Gradient;

public:
    Implicit_octree_domain_with_gradient(const Octree& octree, const Function_with_point& func_with_point,
                                         const Gradient& grad)
        : topo(octree), geom(octree), func(geom, func_with_point), Base_domain(topo, geom, func, grad) {}

private:
    Topology topo;
    Geometry geom;
    Function func;
};

template <class GeomTraits, typename Function_>
class Implicit_octree_domain
    : public Implicit_octree_domain_with_gradient<GeomTraits, Function_, Zero_gradient<GeomTraits>> {
public:
    Implicit_octree_domain(const Octree& octree, const Function_with_point& func)
        : Implicit_cartesian_grid_domain_with_gradient(octree, func, grad) {}

private:
    Gradient grad;
};

}  // namespace Isosurfacing
}  // namespace CGAL

#endif  // CGAL_IMPLICIT_OCTREE_DOMAIN_H
