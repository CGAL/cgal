#ifndef CGAL_MARCHING_CUBES_H
#define CGAL_MARCHING_CUBES_H

#include "Isosurfacing_3/internal/Marching_cubes_3_internal.h"

namespace CGAL {
namespace Isosurfacing {

/**
 * \ingroup PkgIsosurfacing3Ref
 *
 * \brief Creates a polygon soup that represents an isosurface using the marching cubes algorithm.
 *
 * \details
 *
 * \tparam ConcurrencyTag determines if the algorithm is executed sequentially or in parallel.
 * \tparam Domain_ must be a model of `IsosurfacingDomain`.
 * \tparam PointRange is a model of the concept `RandomAccessContainer` and `BackInsertionSequence` whose value type can
 * be constructed from the point type of the polygon mesh. \tparam PolygonRange a model of the concept
 * `RandomAccessContainer` and `BackInsertionSequence` whose value type is itself a model of the concepts
 * `RandomAccessContainer` and `BackInsertionSequence` whose value type is `std::size_t`.
 *
 * \param domain the domain providing input data and its topology
 * \param iso_value value of the isosurface
 * \param points points making the polygons of the soup
 * \param polygons each element in the vector describes a polygon using the indices of the points in points
 */
template <typename ConcurrencyTag = Sequential_tag, class Domain_, class PointRange, class PolygonRange>
void make_triangle_mesh_using_marching_cubes(const Domain_& domain, const typename Domain_::FT iso_value,
                                             PointRange& points, PolygonRange& polygons) {

    internal::Marching_cubes_RG2<Domain_, PointRange, PolygonRange> functor(domain, iso_value, points, polygons);
    domain.iterate_cells(functor, ConcurrencyTag());
}

}  // namespace Isosurfacing
}  // namespace CGAL

#endif  // CGAL_MARCHING_CUBES_H
