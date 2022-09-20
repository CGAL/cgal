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

#include <CGAL/Isosurfacing_3/internal/Cell_type.h>
#include <CGAL/license/Isosurfacing_3.h>

namespace CGAL {
namespace Isosurfacing {

template <class GeomTraits, typename Topology_, typename Geometry_, typename Function_, typename Gradient_>
class Base_domain {
public:
    typedef GeomTraits Geom_traits;
    typedef typename Geom_traits::FT FT;
    typedef typename Geom_traits::Point_3 Point;
    typedef typename Geom_traits::Vector_3 Vector;

    typedef Topology_ Topology;
    typedef typename Topology::Vertex_descriptor Vertex_descriptor;
    typedef typename Topology::Edge_descriptor Edge_descriptor;
    typedef typename Topology::Cell_descriptor Cell_descriptor;

    static constexpr Cell_type CELL_TYPE = Topology::CELL_TYPE;
    static constexpr std::size_t VERTICES_PER_CELL = Topology::VERTICES_PER_CELL;
    static constexpr std::size_t EDGES_PER_CELL = Topology::EDGES_PER_CELL;

    typedef typename Topology::Vertices_incident_to_edge Vertices_incident_to_edge;
    typedef typename Topology::Cells_incident_to_edge Cells_incident_to_edge;
    typedef typename Topology::Cell_vertices Cell_vertices;
    typedef typename Topology::Cell_edges Cell_edges;

    typedef Geometry_ Geometry;
    typedef Function_ Function;
    typedef Gradient_ Gradient;

public:
    Base_domain(const Topology& topo, const Geometry& geom, const Function& func, const Gradient& grad)
        : topo(&topo), geom(&geom), func(&func), grad(&grad) {}

    Point position(const Vertex_descriptor& v) const {
        return geom->operator()(v);
    }

    FT value(const Vertex_descriptor& v) const {
        return func->operator()(v);
    }

    Vector gradient(const Point& p) const {
        return grad->operator()(p);
    }

    Vertices_incident_to_edge edge_vertices(const Edge_descriptor& e) const {
        return topo->edge_vertices(e);
    }

    Cells_incident_to_edge cells_incident_to_edge(const Edge_descriptor& e) const {
        return topo->cells_incident_to_edge(e);
    }

    Cell_vertices cell_vertices(const Cell_descriptor& c) const {
        return topo->cell_vertices(c);
    }

    Cell_edges cell_edges(const Cell_descriptor& c) const {
        return topo->cell_edges(c);
    }

    template <typename Concurrency_tag, typename Functor>
    void iterate_vertices(Functor& f) const {
        topo->iterate_vertices(f, Concurrency_tag());
    }

    template <typename Concurrency_tag, typename Functor>
    void iterate_edges(Functor& f) const {
        topo->iterate_edges(f, Concurrency_tag());
    }

    template <typename Concurrency_tag, typename Functor>
    void iterate_cells(Functor& f) const {
        topo->iterate_cells(f, Concurrency_tag());
    }

private:
    const Topology* topo;
    const Geometry* geom;
    const Function* func;
    const Gradient* grad;
};

}  // namespace Isosurfacing
}  // namespace CGAL

#endif  // CGAL_BASE_DOMAIN_H
