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

#ifndef CGAL_CARTESIAN_GRID_DOMAIN_H
#define CGAL_CARTESIAN_GRID_DOMAIN_H

#include <CGAL/license/Isosurfacing_3.h>
#include <CGAL/Cartesian_grid_3.h>
#include <CGAL/Cartesian_topology_base.h>
#include <CGAL/Default_gradients.h>
#include <CGAL/Isosurfacing_3/internal/Tables.h>

#ifdef CGAL_LINKED_WITH_TBB
#include <tbb/parallel_for.h>
#endif  // CGAL_LINKED_WITH_TBB

namespace CGAL {
namespace Isosurfacing {

template <class GeomTraits, typename Gradient = Zero_gradient<GeomTraits>>
class Cartesian_grid_domain : public Cartesian_topology_base {
public:
    typedef GeomTraits Geom_traits;
    typedef typename Geom_traits::FT FT;
    typedef typename Geom_traits::Point_3 Point;
    typedef typename Geom_traits::Vector_3 Vector;
    typedef typename Geom_traits::Vector_3 Grid_spacing;
    typedef Cartesian_grid_3<Geom_traits> Grid;

public:
    Cartesian_grid_domain(const Grid& grid, const Gradient& grad = Gradient()) : grid(&grid), grad(&grad) {}

    Point position(const Vertex_handle& v) const {
        const Bbox_3& bbox = grid->get_bbox();
        const Vector& spacing = grid->get_spacing();

        return Point(v[0] * spacing.x() + bbox.xmin(), v[1] * spacing.y() + bbox.ymin(),
                     v[2] * spacing.z() + bbox.zmin());
    }

    Vector gradient(const Point& p) const {
        return grad->operator()(p);
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

        auto iterator = [&f, size_x, size_y, size_z](const tbb::blocked_range<std::size_t>& r) {
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

        auto iterator = [&f, size_x, size_y, size_z](const tbb::blocked_range<std::size_t>& r) {
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

        //#pragma omp parallel for
        // for (int x = 0; x < size_x - 1; x++) {
        //    for (std::size_t y = 0; y < size_y - 1; y++) {
        //        for (std::size_t z = 0; z < size_z - 1; z++) {
        //            f({(std::size_t)x, y, z});
        //        }
        //    }
        //}

        auto iterator = [&f, size_x, size_y, size_z](const tbb::blocked_range<std::size_t>& r) {
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
    const Grid* grid;

    const Gradient* grad;
};

}  // namespace Isosurfacing
}  // namespace CGAL

#endif  // CGAL_CARTESIAN_GRID_DOMAIN_H
