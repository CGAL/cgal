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
//                 Daniel Zint

#ifndef CGAL_ISOSURFACING_3_DUAL_CONTOURING_3_H
#define CGAL_ISOSURFACING_3_DUAL_CONTOURING_3_H

#include <CGAL/license/Isosurfacing_3.h>

#include <CGAL/Isosurfacing_3/internal/dual_contouring_functors.h>

#include <CGAL/Container_helper.h>
#include <CGAL/tags.h>

namespace CGAL {
namespace Isosurfacing {

/**
 * \ingroup IS_Methods_grp
 *
 * \brief creates an indexed face set that represents an isosurface using the %Dual Contouring algorithm.
 *
 * \tparam ConcurrencyTag enables sequential versus parallel algorithm.
 *                        Possible values are `Sequential_tag`, `Parallel_if_available_tag`, or `Parallel_tag`.
 * \tparam Domain must be a model of `IsosurfacingDomainWithGradient_3`.
 * \tparam PointRange must be a model of the concept `AssociativeContainer`
 *                    whose value type can be constructed from the point type of the domain.
 * \tparam PolygonRange must be a model of the concepts `RandomAccessContainer` and `BackInsertionSequence`
 *                      whose value type is itself a model of the concepts `RandomAccessContainer`
 *                      and `BackInsertionSequence` whose value type is `std::size_t`.
 *
 * \param domain the domain providing input data and its topology
 * \param isovalue value of the isosurface
 * \param points points of the polygons in the created indexed face set
 * \param polygons each element in the vector describes a polygon using the indices of the points in `points`
 */
#ifdef DOXYGEN_RUNNING // Do not document Positioning
template <typename ConcurrencyTag = CGAL::Sequential_tag,
          typename Domain,
          typename PointRange,
          typename PolygonRange>
void dual_contouring(const Domain& domain,
                     const typename Domain::Geom_traits::FT isovalue,
                     PointRange& points,
                     PolygonRange& polygons)
#else
template <typename ConcurrencyTag = CGAL::Sequential_tag,
          typename Domain,
          typename PointRange,
          typename PolygonRange,
          typename Positioning = internal::Positioning::QEM_SVD<true> >
void dual_contouring(const Domain& domain,
                     const typename Domain::Geom_traits::FT isovalue,
                     PointRange& points,
                     PolygonRange& polygons,
                     const Positioning& positioning = Positioning())
#endif
{
  // create vertices in each relevant cell
  internal::Dual_contouring_vertex_positioning<Domain, Positioning> pos_func(domain, isovalue, positioning);
  domain.template iterate_cells<ConcurrencyTag>(pos_func);

  // connect vertices around an edge to form a face
  internal::Dual_contouring_face_generation<Domain> face_generation(domain, isovalue);
  domain.template iterate_edges<ConcurrencyTag>(face_generation);

  // copy vertices to point range
  CGAL::internal::resize(points, pos_func.points_counter);
  for(const auto& vtop : pos_func.map_voxel_to_point)
    points[pos_func.map_voxel_to_point_id[vtop.first]] = vtop.second;

  // copy faces to polygon range
  polygons.reserve(face_generation.faces.size());
  for(const auto& q : face_generation.faces)
  {
    std::vector<std::size_t> vertex_ids;
    for(const auto& v_id : q.second)
    {
      // ignore voxels that are outside the valid region and are not stored in the map
      if(pos_func.map_voxel_to_point_id.count(v_id) > 0)
        vertex_ids.push_back(pos_func.map_voxel_to_point_id[v_id]);
    }

    // ignore degenerated faces
    if(vertex_ids.size() > 2)
      polygons.push_back(vertex_ids);
  }
}

} // namespace Isosurfacing
} // namespace CGAL

#endif // CGAL_ISOSURFACING_3_DUAL_CONTOURING_3_H
