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
//                 Daniel Zint

#ifndef CGAL_DUAL_CONTOURING_3_H
#define CGAL_DUAL_CONTOURING_3_H

#include <CGAL/Cell_type.h>
#include <CGAL/Isosurfacing_3/internal/Dual_contouring_internal.h>
#include <CGAL/license/Isosurfacing_3.h>
#include <CGAL/tags.h>

namespace CGAL {
namespace Isosurfacing {

/**
 * \ingroup PkgIsosurfacing3Ref
 *
 * \brief Creates an indexed face set that represents an isosurface using the Dual Contouring algorithm.
 *
 * \details
 *
 * \tparam ConcurrencyTag determines if the algorithm is executed sequentially or in parallel.
 *
 * \tparam Domain_ must be a model of `IsosurfacingDomain`.
 *
 * \tparam PointRange is a model of the concept `RandomAccessContainer` and `BackInsertionSequence` whose value type can
 * be constructed from the point type of the polygon mesh.
 * \tparam PolygonRange a model of the concept
 * `RandomAccessContainer` and `BackInsertionSequence` whose value type is itself a model of the concepts
 * `RandomAccessContainer` and `BackInsertionSequence` whose value type is `std::size_t`.
 *
 * \param domain the domain providing input data and its topology
 * \param iso_value value of the isosurface
 * \param points points making the polygons of the indexed face set
 * \param polygons each element in the vector describes a polygon using the indices of the points in points
 */
template <typename Concurrency_tag = Sequential_tag, class Domain_, class PointRange, class PolygonRange,
          class Positioning = internal::Positioning::QEM_SVD<true>>
void dual_contouring(const Domain_& domain, const typename Domain_::FT iso_value, PointRange& points,
                     PolygonRange& polygons, const Positioning& positioning = Positioning()) {

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