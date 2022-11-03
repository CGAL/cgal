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

#ifndef CGAL_OCTREE_TOPOLOGY_H
#define CGAL_OCTREE_TOPOLOGY_H

#include <CGAL/Isosurfacing_3/internal/Cell_type.h>
#include <CGAL/Octree_wrapper.h>
#include <CGAL/license/Isosurfacing_3.h>
#include <CGAL/tags.h>

#ifdef CGAL_LINKED_WITH_TBB
#include <tbb/parallel_for.h>
#endif  // CGAL_LINKED_WITH_TBB

namespace CGAL {
namespace Isosurfacing {

template <class GeomTraits>  // TODO: should not be necessary
class Octree_topology {
public:
    typedef GeomTraits Geom_traits;
    typedef Octree_wrapper<Geom_traits> Octree_;
    typedef std::shared_ptr<Octree_wrapper<Geom_traits>> Octree;
    typedef typename Octree_::Vertex_handle Vertex_descriptor;
    typedef typename Octree_::Edge_handle Edge_descriptor;
    typedef typename Octree_::Voxel_handle Cell_descriptor;

    static constexpr Cell_type CELL_TYPE = CUBICAL_CELL;
    static constexpr std::size_t VERTICES_PER_CELL = 8;
    static constexpr std::size_t EDGES_PER_CELL = 12;

    typedef std::array<Vertex_descriptor, 2> Vertices_incident_to_edge;
    typedef std::array<Cell_descriptor, 4> Cells_incident_to_edge;  // TODO: not always 4
    typedef std::array<Vertex_descriptor, 8> Cell_vertices;
    typedef std::array<Edge_descriptor, 12> Cell_edges;

public:
    Octree_topology(const Octree& octree) : octree(octree) {}

    Vertices_incident_to_edge edge_vertices(const Edge_descriptor& e) const {
        return octree->edge_vertices(e);
    }

    Cells_incident_to_edge cells_incident_to_edge(const Edge_descriptor& e) const {
        return octree->edge_voxels(e);
    }

    Cell_vertices cell_vertices(const Cell_descriptor& c) const {
        return octree->voxel_vertices(c);
    }

    Cell_edges cell_edges(const Cell_descriptor& c) const {
        return octree->voxel_edges(c);
    }

    template <typename Functor>
    void iterate_vertices(Functor& f, Sequential_tag) const {
        for (const Vertex_descriptor& v : octree->leaf_vertices()) {
            f(v);
        }
    }

    template <typename Functor>
    void iterate_edges(Functor& f, Sequential_tag) const {
        for (const Edge_descriptor& e : octree->leaf_edges()) {
            f(e);
        }
    }

    template <typename Functor>
    void iterate_cells(Functor& f, Sequential_tag) const {
        for (const Cell_descriptor& v : octree->leaf_voxels()) {
            f(v);
        }
    }

#ifdef CGAL_LINKED_WITH_TBB
    template <typename Functor>
    void iterate_vertices(Functor& f, Parallel_tag) const {
        const auto& vertices = octree->leaf_vertices();

        auto iterator = [&](const tbb::blocked_range<std::size_t>& r) {
            for (std::size_t i = r.begin(); i != r.end(); i++) {
                f(vertices[i]);
            }
        };

        tbb::parallel_for(tbb::blocked_range<std::size_t>(0, vertices.size()), iterator);
    }

    template <typename Functor>
    void iterate_edges(Functor& f, Parallel_tag) const {
        const auto& edges = octree->leaf_edges();

        auto iterator = [&](const tbb::blocked_range<std::size_t>& r) {
            for (std::size_t i = r.begin(); i != r.end(); i++) {
                f(edges[i]);
            }
        };

        tbb::parallel_for(tbb::blocked_range<std::size_t>(0, edges.size()), iterator);
    }

    template <typename Functor>
    void iterate_cells(Functor& f, Parallel_tag) const {
        const auto& cells = octree->leaf_voxels();

        auto iterator = [&](const tbb::blocked_range<std::size_t>& r) {
            for (std::size_t i = r.begin(); i != r.end(); i++) {
                f(cells[i]);
            }
        };

        tbb::parallel_for(tbb::blocked_range<std::size_t>(0, cells.size()), iterator);
    }
#endif  // CGAL_LINKED_WITH_TBB

private:
    const Octree octree;
};

}  // namespace Isosurfacing
}  // namespace CGAL

#endif  // CGAL_OCTREE_TOPOLOGY_H
