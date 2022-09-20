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

#ifndef CGAL_GRID_TOPOLOGY_H
#define CGAL_GRID_TOPOLOGY_H

#include <CGAL/Isosurfacing_3/internal/Cell_type.h>
#include <CGAL/Isosurfacing_3/internal/Tables.h>
#include <CGAL/license/Isosurfacing_3.h>
#include <CGAL/tags.h>

#include <array>

#ifdef CGAL_LINKED_WITH_TBB
#include <tbb/parallel_for.h>
#endif  // CGAL_LINKED_WITH_TBB

namespace CGAL {
namespace Isosurfacing {

class Grid_topology {
public:
    typedef std::array<std::size_t, 3> Vertex_descriptor;
    typedef std::array<std::size_t, 4> Edge_descriptor;
    typedef std::array<std::size_t, 3> Cell_descriptor;

    static constexpr Cell_type CELL_TYPE = CUBICAL_CELL;
    static constexpr std::size_t VERTICES_PER_CELL = 8;
    static constexpr std::size_t EDGES_PER_CELL = 12;

    typedef std::array<Vertex_descriptor, 2> Vertices_incident_to_edge;
    typedef std::array<Cell_descriptor, 4> Cells_incident_to_edge;
    typedef std::array<Vertex_descriptor, VERTICES_PER_CELL> Cell_vertices;
    typedef std::array<Edge_descriptor, EDGES_PER_CELL> Cell_edges;

public:
    Grid_topology(const std::size_t size_i, const std::size_t size_j, const std::size_t size_k)
        : size_i(size_i), size_j(size_j), size_k(size_k) {}

    Vertices_incident_to_edge edge_vertices(const Edge_descriptor& e) const {
        Vertices_incident_to_edge ev;
        ev[0] = {e[0], e[1], e[2]};
        ev[1] = {e[0], e[1], e[2]};
        ev[1][e[3]] += 1;
        return ev;
    }

    Cells_incident_to_edge cells_incident_to_edge(const Edge_descriptor& e) const {
        const int local = internal::Cube_table::edge_store_index[e[3]];
        auto neighbors = internal::Cube_table::edge_to_voxel_neighbor[local];

        Cells_incident_to_edge cite;
        for (std::size_t i = 0; i < cite.size(); i++) {
            for (std::size_t j = 0; j < cite[i].size(); j++) {
                cite[i][j] = e[j] + neighbors[i][j];
            }
        }
        return cite;
    }

    Cell_vertices cell_vertices(const Cell_descriptor& c) const {
        Cell_vertices cv;
        for (std::size_t i = 0; i < cv.size(); i++) {
            for (std::size_t j = 0; j < c.size(); j++) {
                cv[i][j] = c[j] + internal::Cube_table::local_vertex_position[i][j];
            }
        }
        return cv;
    }

    Cell_edges cell_edges(const Cell_descriptor& c) const {
        Cell_edges ce;
        for (std::size_t i = 0; i < ce.size(); i++) {
            for (std::size_t j = 0; j < c.size(); j++) {
                ce[i][j] = c[j] + internal::Cube_table::global_edge_id[i][j];
            }
            ce[i][3] = internal::Cube_table::global_edge_id[i][3];
        }
        return ce;
    }

    template <typename Functor>
    void iterate_vertices(Functor& f, Sequential_tag) const {
        for (std::size_t i = 0; i < size_i; i++) {
            for (std::size_t j = 0; j < size_j; j++) {
                for (std::size_t k = 0; k < size_k; k++) {
                    f({i, j, k});
                }
            }
        }
    }

    template <typename Functor>
    void iterate_edges(Functor& f, Sequential_tag) const {
        for (std::size_t i = 0; i < size_i - 1; i++) {
            for (std::size_t j = 0; j < size_j - 1; j++) {
                for (std::size_t k = 0; k < size_k - 1; k++) {
                    f({i, j, k, 0});
                    f({i, j, k, 1});
                    f({i, j, k, 2});
                }
            }
        }
    }

    template <typename Functor>
    void iterate_cells(Functor& f, Sequential_tag) const {
        for (std::size_t i = 0; i < size_i - 1; i++) {
            for (std::size_t j = 0; j < size_j - 1; j++) {
                for (std::size_t k = 0; k < size_k - 1; k++) {
                    f({i, j, k});
                }
            }
        }
    }

#ifdef CGAL_LINKED_WITH_TBB
    template <typename Functor>
    void iterate_vertices(Functor& f, Parallel_tag) const {
        auto iterator = [&f, size_i, size_j, size_k](const tbb::blocked_range<std::size_t>& r) {
            for (std::size_t i = r.begin(); i != r.end(); i++) {
                for (std::size_t j = 0; j < size_j; j++) {
                    for (std::size_t k = 0; k < size_k; k++) {
                        f({i, j, k});
                    }
                }
            }
        };

        tbb::parallel_for(tbb::blocked_range<std::size_t>(0, size_i), iterator);
    }

    template <typename Functor>
    void iterate_edges(Functor& f, Parallel_tag) const {
        auto iterator = [&f, size_i, size_j, size_k](const tbb::blocked_range<std::size_t>& r) {
            for (std::size_t i = r.begin(); i != r.end(); i++) {
                for (std::size_t j = 0; j < size_j - 1; j++) {
                    for (std::size_t k = 0; k < size_k - 1; k++) {
                        f({i, j, k, 0});
                        f({i, j, k, 1});
                        f({i, j, k, 2});
                    }
                }
            }
        };

        tbb::parallel_for(tbb::blocked_range<std::size_t>(0, size_i - 1), iterator);
    }

    template <typename Functor>
    void iterate_cells(Functor& f, Parallel_tag) const {
        auto iterator = [&f, size_i, size_j, size_k](const tbb::blocked_range<std::size_t>& r) {
            for (std::size_t i = r.begin(); i != r.end(); i++) {
                for (std::size_t j = 0; j < size_j - 1; j++) {
                    for (std::size_t k = 0; k < size_k - 1; k++) {
                        f({i, j, k});
                    }
                }
            }
        };

        tbb::parallel_for(tbb::blocked_range<std::size_t>(0, size_i - 1), iterator);
    }
#endif  // CGAL_LINKED_WITH_TBB

private:
    std::size_t size_i;
    std::size_t size_j;
    std::size_t size_k;
};

}  // namespace Isosurfacing
}  // namespace CGAL

#endif  // CGAL_GRID_TOPOLOGY_H
