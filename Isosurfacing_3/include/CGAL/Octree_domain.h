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

#ifndef CGAL_OCTREE_DOMAIN_H
#define CGAL_OCTREE_DOMAIN_H

#include <CGAL/Cell_type.h>
#include <CGAL/Default_gradients.h>
#include <CGAL/Octree_wrapper.h>

#ifdef CGAL_LINKED_WITH_TBB
#include <tbb/parallel_for.h>
#endif  // CGAL_LINKED_WITH_TBB

#include <array>

namespace CGAL {
namespace Isosurfacing {

template <typename GeomTraits, typename Gradient = Zero_gradient<GeomTraits>>
class Octree_domain {
public:
    typedef GeomTraits Geom_traits;
    typedef typename Geom_traits::FT FT;
    typedef typename Geom_traits::Point_3 Point;
    typedef typename Geom_traits::Vector_3 Vector;

    typedef Octree_wrapper<Geom_traits> Octree;
    typedef typename Octree::Vertex_handle Vertex_handle;
    typedef typename Octree::Edge_handle Edge_handle;
    typedef typename Octree::Voxel_handle Cell_handle;

    static constexpr Cell_type CELL_TYPE = CUBICAL_CELL;
    static constexpr std::size_t VERTICES_PER_CELL = 8;
    static constexpr std::size_t EDGES_PER_CELL = 12;

    typedef std::array<Vertex_handle, 2> Edge_vertices;
    typedef std::array<Cell_handle, 4> Cells_incident_to_edge;  // TODO: not alwayys 4
    typedef std::array<Vertex_handle, 8> Cell_vertices;
    typedef std::array<Edge_handle, 12> Cell_edges;

public:
    Octree_domain(const Octree& octree, const Gradient& grad = Gradient()) : octree_(&octree), grad(&grad) {}

    Point position(const Vertex_handle& v) const {
        return octree_->point(v);
    }

    Vector gradient(const Point& p) const {
        return grad->operator()(p);
    }

    FT value(const Vertex_handle& v) const {
        return octree_->vertex_value(v);
    }

    Edge_vertices edge_vertices(const Edge_handle& e_id) const {
        return octree_->edge_vertices(e_id);
    }

    Cells_incident_to_edge cells_incident_to_edge(const Edge_handle& e_id) const {
        return octree_->edge_voxels(e_id);
    }

    Cell_vertices cell_vertices(const Cell_handle& vox) const {
        return octree_->voxel_vertices(vox);
    }

    Cell_edges cell_edges(const Cell_handle& vox) const {
        return octree_->voxel_edges(vox);
    }

    template <typename Functor>
    void iterate_vertices(Functor& f, Sequential_tag = Sequential_tag()) const {
        for (const Vertex_handle& v : octree_->leaf_vertices()) {
            f(v);
        }
    }

    template <typename Functor>
    void iterate_edges(Functor& f, Sequential_tag = Sequential_tag()) const {
        for (const Edge_handle& e : octree_->leaf_edges()) {
            f(e);
        }
    }

    template <typename Functor>
    void iterate_cells(Functor& f, Sequential_tag = Sequential_tag()) const {
        for (const Cell_handle& v : octree_->leaf_voxels()) {
            f(v);
        }
    }

#ifdef CGAL_LINKED_WITH_TBB
    template <typename Functor>
    void iterate_vertices(Functor& f, Parallel_tag) const {
        const auto& vertices = octree_->leaf_vertices();

        auto iterator = [&](const tbb::blocked_range<std::size_t>& r) {
            for (std::size_t i = r.begin(); i != r.end(); i++) {
                f(vertices[i]);
            }
        };

        tbb::parallel_for(tbb::blocked_range<std::size_t>(0, vertices.size()), iterator);
    }

    template <typename Functor>
    void iterate_edges(Functor& f, Parallel_tag) const {
        const auto& edges = octree_->leaf_edges();

        auto iterator = [&](const tbb::blocked_range<std::size_t>& r) {
            for (std::size_t i = r.begin(); i != r.end(); i++) {
                f(edges[i]);
            }
        };

        tbb::parallel_for(tbb::blocked_range<std::size_t>(0, edges.size()), iterator);
    }

    template <typename Functor>
    void iterate_cells(Functor& f, Parallel_tag) const {
        const auto& cells = octree_->leaf_voxels();

        auto iterator = [&](const tbb::blocked_range<std::size_t>& r) {
            for (std::size_t i = r.begin(); i != r.end(); i++) {
                f(cells[i]);
            }
        };

        tbb::parallel_for(tbb::blocked_range<std::size_t>(0, cells.size()), iterator);
    }
#endif  // CGAL_LINKED_WITH_TBB

private:
    const Octree* octree_;

    const Gradient* grad;
};

}  // namespace Isosurfacing
}  // namespace CGAL

#endif  // CGAL_OCTREE_DOMAIN_H
