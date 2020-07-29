// Copyright (c) 2018, 2019 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Mael Rouxel-Labbé

#ifndef CGAL_POLYGON_MESH_PROCESSING_INTERNAL_SNAPPING_HELPER_H
#define CGAL_POLYGON_MESH_PROCESSING_INTERNAL_SNAPPING_HELPER_H

#include <CGAL/license/Polygon_mesh_processing/repair.h>

#include <CGAL/boost/graph/iterator.h>
#include <CGAL/number_utils.h>
#include <CGAL/Polygon_mesh_processing/internal/named_function_params.h>
#include <CGAL/Polygon_mesh_processing/internal/named_params_helper.h>

namespace CGAL {
namespace Polygon_mesh_processing {
namespace internal {

template <typename VertexRange,
          typename HalfedgeOutputIterator,
          typename PolygonMesh>
void vertices_as_halfedges(const VertexRange& vertex_range,
                           const PolygonMesh& pmesh,
                           HalfedgeOutputIterator out)
{
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor                vertex_descriptor;

  for(vertex_descriptor v : vertex_range)
    *out++ = halfedge(v, pmesh);
}

// Assigns at each vertex the 'tolerance' value as tolerance, but bounded by a percentage of the length of its shortest incident edge
template <typename HalfedgeRange,
          typename ToleranceMap,
          typename PolygonMesh,
          typename NamedParameters>
void assign_tolerance_with_local_edge_length_bound(const HalfedgeRange& halfedge_range,
                                                   ToleranceMap& tolerance_map,
                                                   const typename GetGeomTraits<PolygonMesh, NamedParameters>::type::FT tolerance,
                                                   PolygonMesh& mesh,
                                                   const NamedParameters& np)
{
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor                vertex_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor              halfedge_descriptor;

  typedef typename GetVertexPointMap<PolygonMesh, NamedParameters>::type              VPM;
  typedef typename GetGeomTraits<PolygonMesh, NamedParameters>::type                  GT;
  typedef typename GT::FT                                                             FT;

  using parameters::get_parameter;
  using parameters::choose_parameter;

  GT gt = choose_parameter<GT>(get_parameter(np, internal_np::geom_traits));
  VPM vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                             get_property_map(vertex_point, mesh));

  for(halfedge_descriptor hd : halfedge_range)
  {
    const vertex_descriptor vd = target(hd, mesh);
    CGAL::Halfedge_around_target_iterator<PolygonMesh> hit, hend;
    boost::tie(hit, hend) = CGAL::halfedges_around_target(vd, mesh);
    CGAL_assertion(hit != hend);

    FT sq_length = gt.compute_squared_distance_3_object()(get(vpm, source(*hit, mesh)),
                                                          get(vpm, target(*hit, mesh)));
    FT min_sq_dist = sq_length;
    ++hit;

    for(; hit!=hend; ++hit)
    {
      sq_length = gt.compute_squared_distance_3_object()(get(vpm, source(*hit, mesh)),
                                                         get(vpm, target(*hit, mesh)));

      if(sq_length < min_sq_dist)
        min_sq_dist = sq_length;
    }

#ifdef CGAL_PMP_SNAP_DEBUG_PP
    std::cout << "tolerance at vd: " << vd << " [" << get(vpm, vd) << "]: min of "
              << 0.9 * CGAL::approximate_sqrt(min_sq_dist) << " AND " << tolerance << std::endl;
#endif
    put(tolerance_map, vd, CGAL::min<FT>(0.9 * CGAL::approximate_sqrt(min_sq_dist), tolerance));
  }
}

template <typename HalfedgeRange,
          typename ToleranceMap,
          typename PolygonMesh>
void assign_tolerance_with_local_edge_length_bound(const HalfedgeRange& halfedge_range,
                                                   ToleranceMap& tolerance_map,
                                                   const typename GetGeomTraits<PolygonMesh>::type::FT tolerance,
                                                   PolygonMesh& mesh)
{
  return assign_tolerance_with_local_edge_length_bound(halfedge_range, tolerance_map, tolerance, mesh, CGAL::parameters::all_default());
}

template <typename GeomTraits>
bool is_collinear_with_tolerance(const typename GeomTraits::Point_3& p, // va == vb
                                 const typename GeomTraits::Point_3& pa,
                                 const typename GeomTraits::Point_3& pb,
                                 const GeomTraits& gt)
{
  typedef typename GeomTraits::FT                                                     FT;
  typedef typename GeomTraits::Vector_3                                               Vector_3;

  const Vector_3 va = gt.construct_vector_3_object()(p, pa);
  const Vector_3 vb = gt.construct_vector_3_object()(p, pb);
  const FT sp = gt.compute_scalar_product_3_object()(va, vb);

  // To avoid junctions like:
  //  ----va vb-----  from being locked since it could be needed in ------------- later
  //       | |                                                      ----va  vb---
  //                                                                    |   |
  if(sp < FT(0))
    return false;

  const FT sq_va_l = gt.compute_squared_distance_3_object()(p, pa);
  const FT sq_vb_l = gt.compute_squared_distance_3_object()(p, pb);
  const FT sq_cos = 0.99999228397046306; // CGAL::square(std::cos(1°));

  return (CGAL::square(sp) >= sq_va_l * sq_vb_l * sq_cos);
}

} // namespace internal
} // namespace Polygon_mesh_processing
} // namespace CGAL

#endif // CGAL_POLYGON_MESH_PROCESSING_INTERNAL_SNAPPING_HELPER_H
