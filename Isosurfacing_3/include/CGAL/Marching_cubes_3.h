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

#ifndef CGAL_MARCHING_CUBES_3_H
#define CGAL_MARCHING_CUBES_3_H

#include <CGAL/Cell_type.h>
#include <CGAL/Isosurfacing_3/internal/Tmc_internal.h>
#include <CGAL/license/Isosurfacing_3.h>
#include <CGAL/tags.h>

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
 *
 * \tparam Domain_ must be a model of `IsosurfacingDomain`.
 *
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
template <typename Concurrency_tag = Sequential_tag, class Domain_, class PointRange, class TriangleRange>
void marching_cubes(const Domain_& domain, const typename Domain_::FT iso_value, PointRange& points,
                    TriangleRange& polygons, bool topologically_correct = true) {

    // static_assert(Domain_::CELL_TYPE & CUBICAL_CELL);

    if (topologically_correct) {
        internal::TMC_functor<Domain_, PointRange, PolygonRange> functor(domain, iso_value, points, polygons);
        domain.iterate_cells(functor, Concurrency_tag());
    } else {
        internal::Marching_cubes_functor<Domain_> functor(domain, iso_value);
        domain.iterate_cells(functor, Concurrency_tag());
        internal::to_indexed_face_set(functor.get_triangles(), points, polygons);
    }
}

}  // namespace Isosurfacing
}  // namespace CGAL

#endif  // CGAL_MARCHING_CUBES_3_H
