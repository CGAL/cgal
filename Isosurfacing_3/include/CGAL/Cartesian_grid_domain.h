#ifndef CGAL_CARTESIAN_GRID_DOMAIN_H
#define CGAL_CARTESIAN_GRID_DOMAIN_H

#include <CGAL/Cartesian_grid_3.h>
#include <CGAL/Cartesian_topology_base.h>
#include <CGAL/Isosurfacing_3/internal/Tables.h>

#ifdef CGAL_LINKED_WITH_TBB
#include <tbb/parallel_for.h>
#endif  // CGAL_LINKED_WITH_TBB

namespace CGAL {
namespace Isosurfacing {

template <class GeomTraits>
class Cartesian_grid_domain : public Cartesian_topology_base {
public:
    typedef GeomTraits Geom_traits;
    typedef typename Geom_traits::FT FT;
    typedef typename Geom_traits::Point_3 Point;
    typedef typename Geom_traits::Vector_3 Vector;

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

    template <typename Functor>
    void iterate_vertices(Functor& f, Sequential_tag tag = Sequential_tag()) const {
        iterate_vertices_base(f, tag, grid->xdim(), grid->ydim(), grid->zdim());
    }

    template <typename Functor>
    void iterate_edges(Functor& f, Sequential_tag tag = Sequential_tag()) const {
        iterate_edges_base(f, tag, grid->xdim(), grid->ydim(), grid->zdim());
    }

    template <typename Functor>
    void iterate_cells(Functor& f, Sequential_tag tag = Sequential_tag()) const {
        iterate_cells_base(f, tag, grid->xdim(), grid->ydim(), grid->zdim());
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
