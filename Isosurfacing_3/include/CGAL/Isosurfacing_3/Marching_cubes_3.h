// Copyright (c) 2022-2023 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Julian Stahl

#ifndef CGAL_ISOSURFACING_3_MARCHING_CUBES_3_H
#define CGAL_ISOSURFACING_3_MARCHING_CUBES_3_H

#include <CGAL/license/Isosurfacing_3.h>

#include <CGAL/Isosurfacing_3/internal/Tmc_internal.h>
#include <CGAL/tags.h>

namespace CGAL {
namespace Isosurfacing {

/**
 * \ingroup PkgIsosurfacing3Ref
 *
 * \brief creates a triangular indexed face set that represents an isosurface using the Marching Cubes algorithm.
 *
 * @todo summary; citation; link to user manual.
 *
 * \tparam ConcurrencyTag enables sequential versus parallel algorithm.
 *                        Possible values are `Parallel_if_available_tag`, `Parallel_tag`, or `Sequential_tag`.
 * \tparam Domain_ must be a model of `IsosurfacingDomain`.
 * \tparam PointRange must be a model of the concept `RandomAccessContainer` and `BackInsertionSequence`
 *                    whose value type can be constructed from the point type of the domain.
 * \tparam PolygonRange must be a model of the concept `RandomAccessContainer` and `BackInsertionSequence`
 *                      whose value type is itself a model of the concepts `RandomAccessContainer`
 *                      and `BackInsertionSequence` whose value type is `std::size_t`.
 *
 * \param domain the domain providing input data and its topology
 * \param isovalue value of the isosurface
 * \param points points of the triangles in the created indexed face set
 * \param triangles each element in the vector describes a triangle using the indices of the points in `points`
 * \param topologically_correct whether the topologically correct variant of Marching Cubes should be used
 */
template <typename Concurrency_tag = Sequential_tag,
          typename Domain_,
          typename PointRange,
          typename TriangleRange>
void marching_cubes(const Domain_& domain,
                    const typename Domain_::FT isovalue,
                    PointRange& points,
                    TriangleRange& triangles,
                    bool topologically_correct = true)
{
  if(topologically_correct)
  {
    // run TMC and directly write the result to points and triangles
    internal::TMC_functor<Domain_, PointRange, TriangleRange> functor(domain, isovalue, points, triangles);
    domain.template iterate_cells<Concurrency_tag>(functor);
  }
  else
  {
    // run MC
    internal::Marching_cubes_3<Domain_> functor(domain, isovalue);
    domain.template iterate_cells<Concurrency_tag>(functor);

    // copy the result to points and triangles
    internal::to_indexed_face_set(functor.triangles(), points, triangles);
  }
}

} // namespace Isosurfacing
} // namespace CGAL

#endif // CGAL_ISOSURFACING_3_MARCHING_CUBES_3_H
