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

#include <CGAL/Isosurfacing_3/internal/Cell_type.h>
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
 * \tparam Positioning is a functor containing the operator() that takes `domain`, `iso_value`, `cell`, and `position`
 * as input and returns a boolean that is true if the isosurface intersects the cell.
 *
 * \param domain the domain providing input data and its topology
 * \param iso_value value of the isosurface
 * \param points points making the polygons of the indexed face set
 * \param polygons each element in the vector describes a polygon using the indices of the points in points
 * \param positioning the functor dealing with vertex positioning inside a voxel
 */
template <typename Concurrency_tag = Sequential_tag, class Domain_, class PointRange, class PolygonRange,
          class Positioning = internal::Positioning::QEM_SVD<true>>
void dual_contouring(const Domain_& domain, const typename Domain_::FT iso_value, PointRange& points,
                     PolygonRange& polygons, const Positioning& positioning = Positioning()) {

    // static_assert(Domain_::CELL_TYPE & ANY_CELL);

    internal::Dual_contouring_vertex_positioning<Domain_, Positioning> pos_func(domain, iso_value, positioning);
    domain.iterate_cells<Concurrency_tag>(pos_func);

    internal::Dual_contouring_face_generation<Domain_> face_generation(domain, iso_value);
    domain.iterate_edges<Concurrency_tag>(face_generation);

    // write points and faces in ranges
    points.resize(pos_func.points_counter);
    for (const auto& vtop : pos_func.map_voxel_to_point) {
        points[pos_func.map_voxel_to_point_id[vtop.first]] = vtop.second;
    }

    polygons.reserve(face_generation.faces.size());
    for (const auto& q : face_generation.faces) {
        std::vector<std::size_t> vertex_ids;
        for (const auto& v_id : q.second) {
            if (pos_func.map_voxel_to_point_id.count(v_id) > 0) {
                vertex_ids.push_back(pos_func.map_voxel_to_point_id[v_id]);
            }
        }
        if (vertex_ids.size() > 2) {
            polygons.push_back(vertex_ids);
        }
    }
}

}  // namespace Isosurfacing
}  // namespace CGAL

#endif  // CGAL_DUAL_CONTOURING_3_H