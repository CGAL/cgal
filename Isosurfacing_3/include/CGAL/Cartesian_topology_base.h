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

#ifndef CGAL_CARTESIAN_TOPOLOGY_BASE_H
#define CGAL_CARTESIAN_TOPOLOGY_BASE_H

#include <CGAL/license/Isosurfacing_3.h>
#include <CGAL/Cell_type.h>
#include <CGAL/Isosurfacing_3/internal/Tables.h>
#include <CGAL/tags.h>

#include <array>

namespace CGAL {
namespace Isosurfacing {

class Cartesian_topology_base {
public:
    typedef std::array<std::size_t, 3> Vertex_handle;
    typedef std::array<std::size_t, 4> Edge_handle;
    typedef std::array<std::size_t, 3> Cell_handle;

    static constexpr Cell_type CELL_TYPE = CUBICAL_CELL;
    static constexpr std::size_t VERTICES_PER_CELL = 8;
    static constexpr std::size_t EDGES_PER_CELL = 12;

    typedef std::array<Vertex_handle, 2> Edge_vertices;
    typedef std::array<Cell_handle, 4> Cells_incident_to_edge;
    typedef std::array<Vertex_handle, VERTICES_PER_CELL> Cell_vertices;
    typedef std::array<Edge_handle, EDGES_PER_CELL> Cell_edges;

public:
    Edge_vertices edge_vertices(const Edge_handle& e) const {
        Edge_vertices ev;
        ev[0] = {e[0], e[1], e[2]};
        ev[1] = {e[0], e[1], e[2]};
        ev[1][e[3]] += 1;
        return ev;
    }

    Cells_incident_to_edge cells_incident_to_edge(const Edge_handle& e) const {
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

    Cell_vertices cell_vertices(const Cell_handle& v) const {
        Cell_vertices cv;
        for (std::size_t i = 0; i < cv.size(); i++) {
            for (std::size_t j = 0; j < v.size(); j++) {
                cv[i][j] = v[j] + internal::Cube_table::local_vertex_position[i][j];
            }
        }
        return cv;
    }

    Cell_edges cell_edges(const Cell_handle& v) const {
        Cell_edges ce;
        for (std::size_t i = 0; i < ce.size(); i++) {
            for (std::size_t j = 0; j < v.size(); j++) {
                ce[i][j] = v[j] + internal::Cube_table::global_edge_id[i][j];
            }
            ce[i][3] = internal::Cube_table::global_edge_id[i][3];
        }
        return ce;
    }

protected:
    template <typename Functor>
    void iterate_vertices_base(Functor& f, Sequential_tag, const std::size_t size_x, const std::size_t size_y,
                               const std::size_t size_z) const {

        for (std::size_t x = 0; x < size_x; x++) {
            for (std::size_t y = 0; y < size_y; y++) {
                for (std::size_t z = 0; z < size_z; z++) {
                    f({x, y, z});
                }
            }
        }
    }

    template <typename Functor>
    void iterate_edges_base(Functor& f, Sequential_tag, const std::size_t size_x, const std::size_t size_y,
                            const std::size_t size_z) const {

        for (std::size_t x = 0; x < size_x - 1; x++) {
            for (std::size_t y = 0; y < size_y - 1; y++) {
                for (std::size_t z = 0; z < size_z - 1; z++) {
                    f({x, y, z, 0});
                    f({x, y, z, 1});
                    f({x, y, z, 2});
                }
            }
        }
    }

    template <typename Functor>
    void iterate_cells_base(Functor& f, Sequential_tag, const std::size_t size_x, const std::size_t size_y,
                            const std::size_t size_z) const {

        for (std::size_t x = 0; x < size_x - 1; x++) {
            for (std::size_t y = 0; y < size_y - 1; y++) {
                for (std::size_t z = 0; z < size_z - 1; z++) {
                    f({x, y, z});
                }
            }
        }
    }
};

}  // namespace Isosurfacing
}  // namespace CGAL

#endif  // CGAL_CARTESIAN_TOPOLOGY_BASE_H
