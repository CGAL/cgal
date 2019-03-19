// Copyright (c) 2015, 2018 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
//
// Author(s)     :  Konstantinos Katrioplas,
//                  Mael Rouxel-Labbé

#ifndef CGAL_POLYGON_MESH_PROCESSING_SHAPE_PREDICATES_H
#define CGAL_POLYGON_MESH_PROCESSING_SHAPE_PREDICATES_H

#include <CGAL/license/Polygon_mesh_processing/repair.h>

#include <CGAL/Polygon_mesh_processing/internal/named_function_params.h>
#include <CGAL/Polygon_mesh_processing/internal/named_params_helper.h>

#include <CGAL/array.h>
#include <CGAL/boost/graph/iterator.h>
#include <CGAL/boost/graph/helpers.h>

#include <boost/graph/graph_traits.hpp>

#include <limits>
#include <map>
#include <utility>
#include <vector>

namespace CGAL {

namespace Polygon_mesh_processing {

/// \ingroup PMP_repairing_grp
/// checks whether an edge is degenerate.
/// An edge is considered degenerate if the geometric positions of its two extremities are identical.
///
/// @tparam PolygonMesh a model of `HalfedgeGraph`
/// @tparam NamedParameters a sequence of \ref pmp_namedparameters "Named Parameters"
///
/// @param e an edge of `pm`
/// @param pm polygon mesh containing `e`
/// @param np optional \ref pmp_namedparameters "Named Parameters" described below
///
/// \cgalNamedParamsBegin
///   \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `pm`.
///                                     The type of this map is model of `ReadWritePropertyMap`.
///                                     If this parameter is omitted, an internal property map for
///                                     `CGAL::vertex_point_t` should be available in `PolygonMesh`
/// \cgalParamEnd
///   \cgalParamBegin{geom_traits} a geometric traits class instance.
///                                The traits class must provide the nested type `Point_3`,
///                                and the nested functor `Equal_3` to check whether two points are identical.
///   \cgalParamEnd
/// \cgalNamedParamsEnd
///
/// \sa `degenerate_edges()`
///
/// \return `true` if the edge `e` is degenerate, `false` otherwise.
template <typename PolygonMesh, typename NamedParameters>
bool is_degenerate_edge(typename boost::graph_traits<PolygonMesh>::edge_descriptor e,
                        const PolygonMesh& pm,
                        const NamedParameters& np)
{
  using boost::get_param;
  using boost::choose_param;

  typedef typename GetVertexPointMap<PolygonMesh, NamedParameters>::const_type VertexPointMap;
  VertexPointMap vpmap = choose_param(get_param(np, internal_np::vertex_point),
                                      get_const_property_map(vertex_point, pm));

  typedef typename GetGeomTraits<PolygonMesh, NamedParameters>::type Traits;
  Traits traits = choose_param(get_param(np, internal_np::geom_traits), Traits());

  return traits.equal_3_object()(get(vpmap, source(e, pm)), get(vpmap, target(e, pm)));
}

template <typename PolygonMesh>
bool is_degenerate_edge(typename boost::graph_traits<PolygonMesh>::edge_descriptor e,
                        const PolygonMesh& pm)
{
  return is_degenerate_edge(e, pm, parameters::all_default());
}

/// \ingroup PMP_repairing_grp
/// checks whether a triangle face is degenerate.
/// A triangle face is considered degenerate if the geometric positions of its vertices are collinear.
///
/// @tparam TriangleMesh a model of `FaceGraph`
/// @tparam NamedParameters a sequence of \ref pmp_namedparameters "Named Parameters"
///
/// @param f a triangle face of `tm`
/// @param tm a triangle mesh containing `f`
/// @param np optional \ref pmp_namedparameters "Named Parameters" described below
///
/// \cgalNamedParamsBegin
///    \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `tm`.
///                                      The type of this map is model of `ReadWritePropertyMap`.
///                                      If this parameter is omitted, an internal property map for
///                                      `CGAL::vertex_point_t` should be available in `TriangleMesh`
/// \cgalParamEnd
///    \cgalParamBegin{geom_traits} a geometric traits class instance.
///                                 The traits class must provide the nested functor `Collinear_3`
///                                 to check whether three points are collinear.
///   \cgalParamEnd
/// \cgalNamedParamsEnd
///
/// \sa `degenerate_faces()`
///
/// \return `true` if the face `f` is degenerate, `false` otherwise.
template <typename TriangleMesh, typename NamedParameters>
bool is_degenerate_triangle_face(typename boost::graph_traits<TriangleMesh>::face_descriptor f,
                                 const TriangleMesh& tm,
                                 const NamedParameters& np)
{
  CGAL_precondition(CGAL::is_triangle(halfedge(f, tm), tm));

  using boost::get_param;
  using boost::choose_param;

  typedef typename GetVertexPointMap<TriangleMesh, NamedParameters>::const_type VertexPointMap;
  VertexPointMap vpmap = choose_param(get_param(np, internal_np::vertex_point),
                                      get_const_property_map(vertex_point, tm));

  typedef typename GetGeomTraits<TriangleMesh, NamedParameters>::type Traits;
  Traits traits = choose_param(get_param(np, internal_np::geom_traits), Traits());

  typename boost::graph_traits<TriangleMesh>::halfedge_descriptor h = halfedge(f, tm);

  return traits.collinear_3_object()(get(vpmap, source(h, tm)),
                                     get(vpmap, target(h, tm)),
                                     get(vpmap, target(next(h, tm), tm)));
}

template <typename TriangleMesh>
bool is_degenerate_triangle_face(typename boost::graph_traits<TriangleMesh>::face_descriptor f,
                                 const TriangleMesh& tm)
{
  return CGAL::Polygon_mesh_processing::is_degenerate_triangle_face(f, tm, parameters::all_default());
}

/// \ingroup PMP_repairing_grp
/// checks whether a triangle face is needle.
/// A triangle is said to be a <i>needle</i> if its longest edge is much longer than its shortest edge.
///
/// @tparam TriangleMesh a model of `FaceGraph`
/// @tparam NamedParameters a sequence of \ref pmp_namedparameters "Named Parameters"
///
/// @param f a triangle face of `tm`
/// @param tm triangle mesh containing `f`
/// @param threshold a bound on the ratio of the longest edge length and the shortest edge length
/// @param np optional \ref pmp_namedparameters "Named Parameters" described below
///
/// \cgalNamedParamsBegin
///    \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `tm`.
///                                      The type of this map is model of `ReadWritePropertyMap`.
///                                      If this parameter is omitted, an internal property map for
///                                      `CGAL::vertex_point_t` should be available in `TriangleMesh`
/// \cgalParamEnd
///    \cgalParamBegin{geom_traits} a geometric traits class instance.
///                                 The traits class must provide the nested type `FT` and
///                                 the nested functor `Compute_squared_distance_3`.
///   \cgalParamEnd
/// \cgalNamedParamsEnd
///
/// \return the shortest halfedge if the triangle face is a needle, and a null halfedge otherwise.
///         If the face contains degenerate edges, a halfedge corresponding to one of these edges is returned.
template <typename TriangleMesh, typename NamedParameters>
typename boost::graph_traits<TriangleMesh>::halfedge_descriptor
is_needle_triangle_face(typename boost::graph_traits<TriangleMesh>::face_descriptor f,
                        const TriangleMesh& tm,
                        const double threshold,
                        const NamedParameters& np)
{
  CGAL_precondition(threshold >= 1.);

  using boost::get_param;
  using boost::choose_param;

  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor   halfedge_descriptor;

  typedef typename GetVertexPointMap<TriangleMesh, NamedParameters>::const_type VertexPointMap;
  VertexPointMap vpmap = choose_param(get_param(np, internal_np::vertex_point),
                                      get_const_property_map(vertex_point, tm));

  typedef typename GetGeomTraits<TriangleMesh, NamedParameters>::type       Traits;
  Traits traits = choose_param(get_param(np, internal_np::geom_traits), Traits());

  typedef typename Traits::FT                                               FT;

  CGAL::Halfedge_around_face_iterator<TriangleMesh> hit, hend;
  boost::tie(hit, hend) = CGAL::halfedges_around_face(halfedge(f, tm), tm);
  CGAL_precondition(std::distance(hit, hend) == 3);

  const halfedge_descriptor h0 = *hit++;
  FT sq_length = traits.compute_squared_distance_3_object()(get(vpmap, source(h0, tm)),
                                                            get(vpmap, target(h0, tm)));

  FT min_sq_length = sq_length, max_sq_length = sq_length;
  halfedge_descriptor min_h = h0;

  for(; hit!=hend; ++hit)
  {
    const halfedge_descriptor h = *hit;
    sq_length = traits.compute_squared_distance_3_object()(get(vpmap, source(h, tm)),
                                                           get(vpmap, target(h, tm)));

    if(max_sq_length < sq_length)
      max_sq_length = sq_length;

    if(min_sq_length > sq_length)
    {
      min_h = h;
      min_sq_length = sq_length;
    }
  }

  if(min_sq_length == 0)
    return min_h;

  const FT sq_threshold = threshold * threshold;
  if(max_sq_length / min_sq_length >= sq_threshold)
  {
    CGAL_assertion(min_h != boost::graph_traits<TriangleMesh>::null_halfedge());
    return min_h;
  }
  else
    return boost::graph_traits<TriangleMesh>::null_halfedge();
}

template <typename TriangleMesh>
typename boost::graph_traits<TriangleMesh>::halfedge_descriptor
is_needle_triangle_face(typename boost::graph_traits<TriangleMesh>::face_descriptor f,
                        const TriangleMesh& tm,
                        const double threshold)
{
  return is_needle_triangle_face(f, tm, threshold, parameters::all_default());
}

/// \ingroup PMP_repairing_grp
/// checks whether a triangle face is a cap.
/// A triangle is said to be a <i>cap</i> if one of the its angles is close to `180` degrees.
///
/// @tparam TriangleMesh a model of `FaceGraph`
/// @tparam NamedParameters a sequence of \ref pmp_namedparameters "Named Parameters"
///
/// @param f a triangle face of `tm`
/// @param tm triangle mesh containing `f`
/// @param threshold the cosine of a minimum angle such that if `f` has an angle greater than this bound,
///                  it is a cap. The threshold is in range `[-1 0]` and corresponds to an angle
///                  between `90` and `180` degrees.
/// @param np optional \ref pmp_namedparameters "Named Parameters" described below
///
/// \cgalNamedParamsBegin
///    \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `tm`.
///                                      The type of this map is model of `ReadWritePropertyMap`.
///                                      If this parameter is omitted, an internal property map for
///                                      `CGAL::vertex_point_t` should be available in `TriangleMesh`
/// \cgalParamEnd
///    \cgalParamBegin{geom_traits} a geometric traits class instance.
///                                 The traits class must provide the nested type `Point_3` and
///                                 the nested functors `Compute_squared_distance_3`, `Construct_vector_3`,
///                                 and `Compute_scalar_product_3`.
///   \cgalParamEnd
/// \cgalNamedParamsEnd
///
/// \return the halfedge opposite of the largest angle if the face is a cap, and a null halfedge otherwise.
template <typename TriangleMesh, typename NamedParameters>
typename boost::graph_traits<TriangleMesh>::halfedge_descriptor
is_cap_triangle_face(typename boost::graph_traits<TriangleMesh>::face_descriptor f,
                     const TriangleMesh& tm,
                     const double threshold,
                     const NamedParameters& np)
{
  CGAL_precondition(CGAL::is_triangle_mesh(tm));
  CGAL_precondition(threshold >= -1.);
  CGAL_precondition(threshold <= 0.);

  using boost::get_param;
  using boost::choose_param;

  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor     vertex_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor   halfedge_descriptor;

  typedef typename GetVertexPointMap<TriangleMesh, NamedParameters>::const_type VertexPointMap;
  VertexPointMap vpmap = choose_param(get_param(np, internal_np::vertex_point),
                                      get_const_property_map(vertex_point, tm));

  typedef typename GetGeomTraits<TriangleMesh, NamedParameters>::type       Traits;
  Traits traits = choose_param(get_param(np, internal_np::geom_traits), Traits());

  typedef typename Traits::FT                                               FT;
  typedef typename Traits::Vector_3                                         Vector_3;

  const FT sq_threshold = threshold * threshold;
  const halfedge_descriptor h0 = halfedge(f, tm);

  cpp11::array<FT, 3> sq_lengths;
  int pos = 0;
  for(halfedge_descriptor h : halfedges_around_face(h0, tm))
  {
    const FT sq_d = traits.compute_squared_distance_3_object()(get(vpmap, source(h, tm)),
                                                               get(vpmap, target(h, tm)));

    // If even one edge is degenerate, it cannot be a cap
    if(sq_d == 0)
      return boost::graph_traits<TriangleMesh>::null_halfedge();

    sq_lengths[pos++] = sq_d;
  }

  pos = 0;
  for(halfedge_descriptor h : halfedges_around_face(h0, tm))
  {
    const vertex_descriptor v0 = source(h, tm);
    const vertex_descriptor v1 = target(h, tm);
    const vertex_descriptor v2 = target(next(h, tm), tm);
    const Vector_3 a = traits.construct_vector_3_object()(get(vpmap, v1), get(vpmap, v2));
    const Vector_3 b = traits.construct_vector_3_object()(get(vpmap, v1), get(vpmap, v0));
    const FT dot_ab = traits.compute_scalar_product_3_object()(a, b);
    const bool neg_sp = (dot_ab <= 0);
    const FT sq_a = sq_lengths[(pos+1)%3];
    const FT sq_b = sq_lengths[pos];
    const FT sq_cos = dot_ab * dot_ab / (sq_a * sq_b);

    if(neg_sp && sq_cos >= sq_threshold)
      return prev(h, tm);

    ++pos;
  }
  return boost::graph_traits<TriangleMesh>::null_halfedge();
}

template <typename TriangleMesh>
typename boost::graph_traits<TriangleMesh>::halfedge_descriptor
is_cap_triangle_face(typename boost::graph_traits<TriangleMesh>::face_descriptor f,
                     const TriangleMesh& tm,
                     const double threshold)
{
  return is_cap_triangle_face(f, tm, threshold, parameters::all_default());
}

} } // end namespaces CGAL and PMP

#endif // CGAL_POLYGON_MESH_PROCESSING_SHAPE_PREDICATES_H
