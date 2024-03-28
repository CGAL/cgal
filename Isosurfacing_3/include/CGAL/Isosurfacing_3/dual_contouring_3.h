// Copyright (c) 2022-2024 INRIA Sophia-Antipolis (France), GeometryFactory (France).
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
//                 Mael Rouxel-Labbé

#ifndef CGAL_ISOSURFACING_3_DUAL_CONTOURING_3_H
#define CGAL_ISOSURFACING_3_DUAL_CONTOURING_3_H

#include <CGAL/license/Isosurfacing_3.h>

#include <CGAL/Isosurfacing_3/internal/dual_contouring_functors.h>

#include <CGAL/tags.h>

namespace CGAL {
namespace Isosurfacing {

/**
 * \ingroup IS_Methods_grp
 *
 * \brief creates a polygon soup that discretizes an isosurface using the %Dual Contouring algorithm.
 *
 * \details The point placement strategy within each cell of the space partition is based on
 * Quadric Error Metrics ("QEM", or "QEF" in %Dual Contouring-related works).
 *
 * \tparam ConcurrencyTag enables sequential versus parallel algorithm.
 *                        Possible values are `Sequential_tag`, `Parallel_if_available_tag`, or `Parallel_tag`.
 * \tparam Domain must be a model of `IsosurfacingDomainWithGradient_3`.
 * \tparam PointRange must be a model of the concept `AssociativeContainer`
 *                    whose value type can be constructed from the point type of the domain.
 * \tparam PolygonRange must be a model of the concepts `RandomAccessContainer` and `BackInsertionSequence`
 *                      whose value type is itself a model of the concepts `RandomAccessContainer`
 *                      and `BackInsertionSequence` whose value type is `std::size_t`.
 * \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
 *
 * \param domain the domain providing the spacial partition and the values and gradient data
 * \param isovalue the value defining the isosurface
 * \param points the points of the polygons in the created polygon soup
 * \param polygons each element in the vector describes a polygon using the indices of the points in `points`
 * \param np optional \ref bgl_namedparameters "Named Parameters" described below
 *
 * \cgalNamedParamsBegin
 *   \cgalParamNBegin{constrain_to_cell}
 *     \cgalParamDescription{whether to constrain the vertex position to the geometrical space of its cell}
 *     \cgalParamType{Boolean}
 *     \cgalParamDefault{`false`}
 *     \cgalParamExtra{Constraining the vertex to its dual cell guarantees that the resulting
 *                     surface is a without self-intersections (non-manifoldness aside). Oppositely,
 *                     an unconstrained positioning strategy might produce better looking surfaces
 *                     near sharp features (ridges, corners), at the cost of possible self-intersections.}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{do_not_triangulate_faces}
 *     \cgalParamDescription{If `true`, the output will contain quadrilaterals.
 *                           If `false`, the output will contain triangles.}
 *     \cgalParamType{Boolean}
 *     \cgalParamDefault{`false` (faces are triangulated)}
 *     \cgalParamExtra{Triangulating faces is done by inserting the intersection between an edge and
 *                     the isosurface, and linking it to the dual points of the cells incident to the edge.
 *                     If `constrain_to_cell` is set to `false`, triangulation faces can result in additional
 *                     self-intersections. An alternative that has worse approximation but is less likely
 *                     to produce self-intersections is to use the function
 *                     `CGAL::Polygon_mesh_processing::triangulate_faces()`.}
 *   \cgalParamNEnd
 * \cgalNamedParamsEnd
 *
 * \sa `CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh()`
 */
template <typename ConcurrencyTag = CGAL::Sequential_tag,
          typename Domain,
          typename PointRange,
          typename PolygonRange,
          typename NamedParameters = parameters::Default_named_parameters>
void dual_contouring(const Domain& domain,
                     const typename Domain::Geom_traits::FT isovalue,
                     PointRange& points,
                     PolygonRange& polygons,
                     const NamedParameters& np = parameters::default_values())
{
  internal::Dual_contourer<ConcurrencyTag, Domain, internal::DC_Strategy::QEM> contourer;
  contourer(domain, isovalue, points, polygons, np);
}

} // namespace Isosurfacing
} // namespace CGAL

#endif // CGAL_ISOSURFACING_3_DUAL_CONTOURING_3_H
