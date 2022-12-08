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

#ifndef CGAL_BASE_DOMAIN_H
#define CGAL_BASE_DOMAIN_H

#include <CGAL/license/Isosurfacing_3.h>

#include <CGAL/Isosurfacing_3/internal/Cell_type.h>

#include <memory>

namespace CGAL {
namespace Isosurfacing {
namespace internal {

// A wrapper class to puzzle a domain together from different combinations of topology, geometry, function, and
// gradient.
template <class GeomTraits, typename Topology_, typename Geometry_, typename Function_, typename Gradient_>
class Base_domain {
public:
    typedef GeomTraits Geom_traits;
    typedef typename Geom_traits::FT FT;
    typedef typename Geom_traits::Point_3 Point;
    typedef typename Geom_traits::Vector_3 Vector;

    typedef std::shared_ptr<Topology_> Topology;
    typedef typename Topology_::Vertex_descriptor Vertex_descriptor;
    typedef typename Topology_::Edge_descriptor Edge_descriptor;
    typedef typename Topology_::Cell_descriptor Cell_descriptor;

    static constexpr Cell_type CELL_TYPE = Topology_::CELL_TYPE;
    static constexpr std::size_t VERTICES_PER_CELL = Topology_::VERTICES_PER_CELL;
    static constexpr std::size_t EDGES_PER_CELL = Topology_::EDGES_PER_CELL;

    typedef typename Topology_::Vertices_incident_to_edge Vertices_incident_to_edge;
    typedef typename Topology_::Cells_incident_to_edge Cells_incident_to_edge;
    typedef typename Topology_::Cell_vertices Cell_vertices;
    typedef typename Topology_::Cell_edges Cell_edges;

    typedef std::shared_ptr<Geometry_> Geometry;
    typedef std::shared_ptr<Function_> Function;
    typedef std::shared_ptr<Gradient_> Gradient;

public:
    // Create a base_domain from a topology, geometry, input function, and gradient
    Base_domain(const Topology& topo, const Geometry& geom, const Function& func, const Gradient& grad)
        : topo(topo), geom(geom), func(func), grad(grad) {}

    // Get the position of vertex v
    Point position(const Vertex_descriptor& v) const {
        return geom->operator()(v);
    }

    // Get the value of the function at vertex v
    FT value(const Vertex_descriptor& v) const {
        return func->operator()(v);
    }

    // Get the gradient at vertex v
    Vector gradient(const Point& p) const {
        return grad->operator()(p);
    }

    // Get a container with the two vertices incident to the edge e
    Vertices_incident_to_edge edge_vertices(const Edge_descriptor& e) const {
        return topo->edge_vertices(e);
    }

    // Get a container with all cells incident to the edge e
    Cells_incident_to_edge cells_incident_to_edge(const Edge_descriptor& e) const {
        return topo->cells_incident_to_edge(e);
    }

    // Get a container with all vertices of the cell c
    Cell_vertices cell_vertices(const Cell_descriptor& c) const {
        return topo->cell_vertices(c);
    }

    // Get a container with all edges of the cell c
    Cell_edges cell_edges(const Cell_descriptor& c) const {
        return topo->cell_edges(c);
    }

    // Iterate over all vertices v calling f(v) on every one
    template <typename Concurrency_tag, typename Functor>
    void iterate_vertices(Functor& f) const {
        topo->iterate_vertices(f, Concurrency_tag());
    }

    // Iterate over all edges e calling f(e) on every one
    template <typename Concurrency_tag, typename Functor>
    void iterate_edges(Functor& f) const {
        topo->iterate_edges(f, Concurrency_tag());
    }

    // Iterate over all cells c calling f(c) on every one
    template <typename Concurrency_tag, typename Functor>
    void iterate_cells(Functor& f) const {
        topo->iterate_cells(f, Concurrency_tag());
    }

private:
    const Topology topo;
    const Geometry geom;
    const Function func;
    const Gradient grad;
};

}  // namespace internal
}  // namespace Isosurfacing
}  // namespace CGAL

#endif  // CGAL_BASE_DOMAIN_H
