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

template <class GeomTraits, typename Function_, typename Gradient_>
class Implicit_cartesian_grid_domain_with_gradient
    : public Base_domain<GeomTraits, Grid_topology, Cartesian_grid_geometry<GeomTraits>,
                         Implicit_function_with_geometry<GeomTraits, Cartesian_grid_geometry<GeomTraits>, Function_>,
                         Gradient_> {
public:
    typedef GeomTraits Geom_traits;
    typedef typename Geom_traits::Vector_3 Vector;

    typedef Grid_topology Topology;
    typedef Cartesian_grid_geometry<Geom_traits> Geometry;
    typedef Function_ Function_with_point;
    typedef Implicit_function_with_geometry<Geom_traits, Geometry, Function_with_point> Function;
    typedef Gradient_ Gradient;

    typedef Base_domain<Geom_traits, Topology, Geometry, Function, Gradient> Base;

public:
    Implicit_cartesian_grid_domain_with_gradient(const std::size_t size_i, const std::size_t size_j,
                                                 const std::size_t size_k, const Vector& offset, const Vector& spacing,
                                                 const Function_with_point& func_with_point, const Gradient& grad)
        : topo(size_i, size_j, size_k),
          geom(offset, spacing),
          func(geom, func_with_point),
          Base(topo, geom, func, grad) {}

private:
    Topology topo;
    Geometry geom;
    Function func;
};

template <class GeomTraits, typename Function_>
class Implicit_cartesian_grid_domain
    : public Implicit_cartesian_grid_domain_with_gradient<GeomTraits, Function_, Zero_gradient<GeomTraits>> {
public:
    typedef GeomTraits Geom_traits;
    typedef typename Geom_traits::Vector_3 Vector;

    typedef Function_ Function_with_point;
    typedef Zero_gradient<Geom_traits> Gradient;

    typedef Implicit_cartesian_grid_domain_with_gradient<Geom_traits, Function_, Gradient> Base;

public:
    Implicit_cartesian_grid_domain(const std::size_t size_i, const std::size_t size_j, const std::size_t size_k,
                                   const Vector& offset, const Vector& spacing, const Function_with_point& func)
        : Base(size_i, size_j, size_k, offset, spacing, func, grad) {}

private:
    Gradient grad;
};

}  // namespace Isosurfacing
}  // namespace CGAL

#endif  // CGAL_IMPLICIT_CARTESIAN_GRID_DOMAIN_H
