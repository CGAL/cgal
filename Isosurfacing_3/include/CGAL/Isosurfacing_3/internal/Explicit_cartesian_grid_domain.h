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

#ifndef CGAL_EXPLICIT_CARTESIAN_GRID_DOMAIN_H
#define CGAL_EXPLICIT_CARTESIAN_GRID_DOMAIN_H

#include <CGAL/Cartesian_grid_3.h>
#include <CGAL/Default_gradients.h>
#include <CGAL/Isosurfacing_3/internal/Base_domain.h>
#include <CGAL/Isosurfacing_3/internal/Cartesian_grid_geometry.h>
#include <CGAL/Isosurfacing_3/internal/Grid_topology.h>
#include <CGAL/license/Isosurfacing_3.h>

namespace CGAL {
namespace Isosurfacing {

template <class GeomTraits, typename Gradient_>
class Explicit_cartesian_grid_domain_with_gradient
    : public Base_domain<GeomTraits, Grid_topology, Cartesian_grid_geometry<GeomTraits>, Cartesian_grid_3<GeomTraits>,
                         Gradient_> {
public:
    typedef GeomTraits Geom_traits;
    typedef typename Geom_traits::Vector_3 Vector;

    typedef Grid_topology Topology;
    typedef Cartesian_grid_geometry<Geom_traits> Geometry;
    typedef Cartesian_grid_3<Geom_traits> Function;
    typedef Gradient_ Gradient;

public:
    Explicit_cartesian_grid_domain_with_gradient(const std::size_t size_i, const std::size_t size_j,
                                                 const std::size_t size_k, const Vector& offset, const Vector& spacing,
                                                 const Function& grid, const Gradient& grad)
        : topo(size_i, size_j, size_k), geom(offset, spacing), Base_domain(topo, geom, grid, grad) {}

private:
    Topology topo;
    Geometry geom;
};

template <class GeomTraits>
class Explicit_cartesian_grid_domain
    : public Explicit_cartesian_grid_domain_with_gradient<GeomTraits, Zero_gradient<GeomTraits>> {
public:
    typedef GeomTraits Geom_traits;
    typedef typename Geom_traits::Vector_3 Vector;

    typedef Cartesian_grid_3<Geom_traits> Function;
    typedef Zero_gradient<Geom_traits> Gradient;

public:
    Explicit_cartesian_grid_domain(const std::size_t size_i, const std::size_t size_j, const std::size_t size_k,
                                   const Vector& offset, const Vector& spacing, const Function& grid)
        : Explicit_cartesian_grid_domain_with_gradient(size_i, size_j, size_k, offset, spacing, grid, grad) {}

private:
    Gradient grad;
};

}  // namespace Isosurfacing
}  // namespace CGAL

#endif  // CGAL_EXPLICIT_CARTESIAN_GRID_DOMAIN_H
