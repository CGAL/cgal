#ifndef CGAL_IMPLICIT_DOMAIN_H
#define CGAL_IMPLICIT_DOMAIN_H

#include <CGAL/Bbox_3.h>
#include <CGAL/Cartesian_topology_base.h>
#include <tbb/parallel_for.h>

namespace CGAL {
namespace Isosurfacing {

template <class GeomTraits, typename Function>
class Implicit_domain : public Cartesian_topology_base {
public:
    typedef GeomTraits Geom_traits;
    typedef typename Geom_traits::FT FT;
    typedef typename Geom_traits::Point_3 Point;
    typedef typename Geom_traits::Vector_3 Grid_spacing;

public:
    Implicit_domain(const Function& func, const Bbox_3& domain, const Grid_spacing& spacing)
        : func(&func), bbox(domain), spacing(spacing) {

        sizes[0] = domain.x_span() / spacing.x();
        sizes[1] = domain.y_span() / spacing.y();
        sizes[2] = domain.z_span() / spacing.z();
    }

    Point position(const Vertex_handle& v) const {
        return Point(v[0] * spacing.x() + bbox.xmin(), v[1] * spacing.y() + bbox.ymin(),
                     v[2] * spacing.z() + bbox.zmin());
    }

    FT value(const Vertex_handle& v) const {
        return func->operator()(position(v));
    }

    template <typename Functor>
    void iterate_vertices(Functor& f, Sequential_tag tag = Sequential_tag()) const {
        iterate_vertices_base(f, tag, sizes[0], sizes[1], sizes[2]);
    }

    template <typename Functor>
    void iterate_edges(Functor& f, Sequential_tag tag = Sequential_tag()) const {
        iterate_edges_base(f, tag, sizes[0], sizes[1], sizes[2]);
    }

    template <typename Functor>
    void iterate_cells(Functor& f, Sequential_tag tag = Sequential_tag()) const {
        iterate_cells_base(f, tag, sizes[0], sizes[1], sizes[2]);
    }

#ifdef CGAL_LINKED_WITH_TBB
    template <typename Functor>
    void iterate_vertices(Functor& f, Parallel_tag) const {
        const std::size_t size_x = sizes[0];
        const std::size_t size_y = sizes[1];
        const std::size_t size_z = sizes[2];

        auto iterator = [=, &f](const tbb::blocked_range<std::size_t>& r) {
            for (std::size_t x = r.begin(); x != r.end(); x++) {
                for (std::size_t y = 0; y < size_y; y++) {
                    for (std::size_t z = 0; z < size_z; z++) {
                        f({x, y, z});
                    }
                }
            }
        };

        tbb::parallel_for(tbb::blocked_range<std::size_t>(0, size_x), iterator);
    }

    template <typename Functor>
    void iterate_edges(Functor& f, Parallel_tag) const {
        const std::size_t size_x = sizes[0];
        const std::size_t size_y = sizes[1];
        const std::size_t size_z = sizes[2];

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
        const std::size_t size_x = sizes[0];
        const std::size_t size_y = sizes[1];
        const std::size_t size_z = sizes[2];

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
    const Function* func;

    Bbox_3 bbox;
    Grid_spacing spacing;

    std::array<std::size_t, 3> sizes;
};


template <typename Function, class GeomTraits = typename Function::Geom_traits>
Implicit_domain<GeomTraits, Function> create_implicit_domain(const Function& func, const Bbox_3& domain,
                                                             const typename GeomTraits::Vector_3& spacing) {
    return Implicit_domain<GeomTraits, Function>(func, domain, spacing);
}

}  // namespace Isosurfacing
}  // namespace CGAL

#endif  // CGAL_IMPLICIT_DOMAIN_H
