#ifndef CGAL_TMC_3_H
#define CGAL_TMC_3_H

#include <CGAL/Cell_type.h>
#include <CGAL/Isosurfacing_3/internal/Tmc_internal.h>
#include <CGAL/tags.h>

namespace CGAL {
namespace Isosurfacing {

template <typename Concurrency_tag = Sequential_tag, class Domain_, class PointRange, class PolygonRange>
void make_triangle_mesh_using_tmc(const Domain_& domain, const typename Domain_::FT iso_value, PointRange& points,
                                  PolygonRange& polygons) {

    // static_assert(Domain_::CELL_TYPE & CUBICAL_CELL);

    internal::TMC_functor<Domain_, PointRange, PolygonRange> functor(domain, iso_value, points, polygons);
    domain.iterate_cells(functor, Concurrency_tag());
}

}  // namespace Isosurfacing
}  // namespace CGAL

#endif  // CGAL_TMC_3_H
