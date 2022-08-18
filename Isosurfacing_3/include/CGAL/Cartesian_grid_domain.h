#ifndef CGAL_CARTESIAN_GRID_DOMAIN_H
#define CGAL_CARTESIAN_GRID_DOMAIN_H

#include <tbb/parallel_for.h>

#include "Cartesian_grid_3.h"
#include "Isosurfacing_3/internal/Tables.h"

namespace CGAL {
namespace Isosurfacing {

template <class GeomTraits>
class Cartesian_grid_domain {
public:
    typedef GeomTraits Geom_traits;
    typedef typename Geom_traits::FT FT;
    typedef typename Geom_traits::Point_3 Point;
    typedef typename Geom_traits::Vector_3 Vector;

    typedef std::array<std::size_t, 3> Vertex_handle;
    typedef std::array<std::size_t, 4> Edge_handle;
    typedef std::array<std::size_t, 3> Cell_handle;

    typedef std::array<Vertex_handle, 2> Edge_vertices;
    typedef std::array<Cell_handle, 4> Cells_incident_to_edge;
    typedef std::array<Vertex_handle, internal::Cube_table::N_VERTICES> Cell_vertices;
    typedef std::array<Edge_handle, internal::Cube_table::N_EDGES> Cell_edges;

public:
    Cartesian_grid_domain(const Cartesian_grid_3<Geom_traits>& grid) : grid(&grid) {}

    Point position(const Vertex_handle& v) const {
        const FT vx = grid->voxel_x();
        const FT vy = grid->voxel_y();
        const FT vz = grid->voxel_z();

        return Point(v[0] * vx + grid->offset_x(), v[1] * vy + grid->offset_y(), v[2] * vz + grid->offset_z());
    }

    Vector gradient(const Vertex_handle& v) const {
        const FT vx = grid->voxel_x();
        const FT vy = grid->voxel_y();
        const FT vz = grid->voxel_z();

        Vector g(v[0] * vx + grid->offset_x(), v[1] * vy + grid->offset_y(), v[2] * vz + grid->offset_z());
        return g / std::sqrt(g.squared_length());
    }

    FT value(const Vertex_handle& v) const {
        return grid->value(v[0], v[1], v[2]);
    }

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

    template <typename Functor>
    void iterate_vertices(Functor& f, Sequential_tag) const {
        const std::size_t size_x = grid->xdim();
        const std::size_t size_y = grid->ydim();
        const std::size_t size_z = grid->zdim();

        for (std::size_t x = 0; x < size_x - 1; x++) {
            for (std::size_t y = 0; y < size_y - 1; y++) {
                for (std::size_t z = 0; z < size_z - 1; z++) {
                    f({x, y, z});
                }
            }
        }
    }

#ifdef CGAL_LINKED_WITH_TBB
    template <typename Functor>
    void iterate_vertices(Functor& f, Parallel_tag) const {
        const std::size_t size_x = grid->xdim();
        const std::size_t size_y = grid->ydim();
        const std::size_t size_z = grid->zdim();

        auto iterator = [=, &f](const tbb::blocked_range<std::size_t>& r) {
            for (std::size_t x = r.begin(); x != r.end(); x++) {
                for (std::size_t y = 0; y < size_y - 1; y++) {
                    for (std::size_t z = 0; z < size_z - 1; z++) {
                        f({x, y, z});
                    }
                }
            }
        };

        tbb::parallel_for(tbb::blocked_range<std::size_t>(0, size_x - 1), iterator);
    }
#endif  // CGAL_LINKED_WITH_TBB

    template <typename Functor>
    void iterate_edges(Functor& f, Sequential_tag) const {
        const std::size_t size_x = grid->xdim();
        const std::size_t size_y = grid->ydim();
        const std::size_t size_z = grid->zdim();

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

#ifdef CGAL_LINKED_WITH_TBB
    template <typename Functor>
    void iterate_edges(Functor& f, Parallel_tag) const {
        const std::size_t size_x = grid->xdim();
        const std::size_t size_y = grid->ydim();
        const std::size_t size_z = grid->zdim();

        auto iterator = [=, &f](const tbb::blocked_range<std::size_t>& r) {
            for (std::size_t x = r.begin(); x != r.end(); x++) {
                for (std::size_t y = 0; y < size_y - 1; y++) {
                    for (std::size_t z = 0; z < size_z - 1; z++) {
                        f({x, y, z, 0});
                        f({x, y, z, 1});
                        f({x, y, z, 2});
                    }
                }
            }
        };

        tbb::parallel_for(tbb::blocked_range<std::size_t>(0, size_x - 1), iterator);
    }
#endif  // CGAL_LINKED_WITH_TBB

    template <typename Functor>
    void iterate_cells(Functor& f, Sequential_tag) const {
        const std::size_t size_x = grid->xdim();
        const std::size_t size_y = grid->ydim();
        const std::size_t size_z = grid->zdim();

        for (std::size_t x = 0; x < size_x - 1; x++) {
            for (std::size_t y = 0; y < size_y - 1; y++) {
                for (std::size_t z = 0; z < size_z - 1; z++) {
                    f({x, y, z});
                }
            }
        }
    }

#ifdef CGAL_LINKED_WITH_TBB
    template <typename Functor>
    void iterate_cells(Functor& f, Parallel_tag) const {
        const std::size_t size_x = grid->xdim();
        const std::size_t size_y = grid->ydim();
        const std::size_t size_z = grid->zdim();

        auto iterator = [=, &f](const tbb::blocked_range<std::size_t>& r) {
            for (std::size_t x = r.begin(); x != r.end(); x++) {
                for (std::size_t y = 0; y < size_y - 1; y++) {
                    for (std::size_t z = 0; z < size_z - 1; z++) {
                        f({x, y, z});
                    }
                }
            }
        };

        tbb::parallel_for(tbb::blocked_range<std::size_t>(0, size_x - 1), iterator);
    }
#endif  // CGAL_LINKED_WITH_TBB

private:
    const Cartesian_grid_3<Geom_traits>* grid;
};

}  // namespace Isosurfacing
}  // namespace CGAL

#endif  // CGAL_CARTESIAN_GRID_DOMAIN_H
