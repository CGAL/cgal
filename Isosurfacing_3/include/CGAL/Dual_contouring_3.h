#ifndef CGAL_DUAL_CONTOURING_3_H
#define CGAL_DUAL_CONTOURING_3_H

#include <CGAL/Cell_type.h>
#include <CGAL/Isosurfacing_3/internal/Dual_contouring_internal.h>
#include <CGAL/tags.h>

namespace CGAL {
namespace Isosurfacing {

template <typename Concurrency_tag = Sequential_tag, class Domain_, class PointRange, class PolygonRange,
          class Positioning = internal::Positioning::QEM_SVD<true>>
void make_quad_mesh_using_dual_contouring(const Domain_& domain, const typename Domain_::FT iso_value,
                                          PointRange& points, PolygonRange& polygons, const Positioning& positioning = Positioning()) {

    // static_assert(Domain_::CELL_TYPE & ANY_CELL);

    internal::Dual_contouring_position_functor<Domain_, Positioning> pos_func(domain, iso_value, positioning);
    domain.iterate_cells(pos_func, Concurrency_tag());

    internal::Dual_contouring_quads_functor<Domain_> quad_func(domain, iso_value);
    domain.iterate_edges(quad_func, Concurrency_tag());

    // write points and quads in ranges
    points.resize(pos_func.points_counter);
    for (const auto& vtop : pos_func.map_voxel_to_point) {
        points[pos_func.map_voxel_to_point_id[vtop.first]] = vtop.second;
    }

    polygons.reserve(quad_func.quads.size());
    for (const auto& q : quad_func.quads) {
        std::vector<std::size_t> vertex_ids;
        for (const auto& v_id : q.second) {
            vertex_ids.push_back(pos_func.map_voxel_to_point_id[v_id]);
        }
        polygons.push_back(vertex_ids);
    }
}

}  // namespace Isosurfacing
}  // namespace CGAL

#endif  // CGAL_DUAL_CONTOURING_3_H