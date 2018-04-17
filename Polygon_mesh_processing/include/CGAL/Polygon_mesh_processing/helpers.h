// Copyright (c) 2015 GeometryFactory (France).
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
// Author(s)     :  Konstantinos Katrioplas

#ifndef CGAL_POLYGON_MESH_PROCESSING_HELPERS_H
#define CGAL_POLYGON_MESH_PROCESSING_HELPERS_H

#include <CGAL/Polygon_mesh_processing/internal/named_function_params.h>
#include <CGAL/Polygon_mesh_processing/internal/named_params_helper.h>


namespace CGAL {

namespace Polygon_mesh_processing {

namespace internal {

template <typename PolygonMesh>
void merge_identical_points(PolygonMesh& mesh,
                            typename boost::graph_traits<PolygonMesh>::vertex_descriptor v_keep,
                            typename boost::graph_traits<PolygonMesh>::vertex_descriptor v_rm)
{
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;
  halfedge_descriptor h = halfedge(v_rm, mesh);
  halfedge_descriptor start = h;

  do{
    set_target(h, v_keep, mesh);
    h = opposite(next(h, mesh), mesh);
  } while( h != start );

  remove_vertex(v_rm, mesh);
}
} // end internal

/// \ingroup PMP_repairing_grp
/// checks whether a vertex is non-manifold.
///
/// @tparam PolygonMesh a model of `FaceListGraph` and `MutableFaceGraph`
///
/// @param v the vertex
/// @param tm triangle mesh containing v
///
/// \return true if the vertrex is non-manifold
template <typename PolygonMesh>
bool is_non_manifold_vertex(typename boost::graph_traits<PolygonMesh>::vertex_descriptor v,
                            const PolygonMesh& tm)
{
  CGAL_assertion(CGAL::is_triangle_mesh(tm));

  typedef boost::graph_traits<PolygonMesh> GT;
  typedef typename GT::halfedge_descriptor halfedge_descriptor;

  boost::unordered_set<halfedge_descriptor> halfedges_handled;

  BOOST_FOREACH(halfedge_descriptor h, halfedges_around_target(v, tm))
    halfedges_handled.insert(h);

  BOOST_FOREACH(halfedge_descriptor h, halfedges(tm))
  {
    if(v == target(h, tm))
    {
      if(halfedges_handled.count(h) == 0)
        return true;
    }
  }
  return false;
}

/// \ingroup PMP_repairing_grp
/// checks whether an edge is degenerate.
/// An edge is considered degenerate if the points of its vertices are identical.
///
/// @tparam PolygonMesh a model of `HalfedgeGraph`
/// @tparam NamedParameters a sequence of \ref pmp_namedparameters "Named Parameters"
///
/// @param e the edge
/// @param pm polygon mesh containing e
/// @param np optional \ref pmp_namedparameters "Named Parameters" described below
///
/// \cgalNamedParamsBegin
///    \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `pmesh`. The type of this map is model of `ReadWritePropertyMap`.
/// If this parameter is omitted, an internal property map for
/// `CGAL::vertex_point_t` should be available in `PolygonMesh`
/// \cgalParamEnd
///    \cgalParamBegin{geom_traits} a geometric traits class instance.
///       The traits class must provide the nested type `Point_3`,
///       and the nested functor :
///         - `Equal_3` to check whether 2 points are identical
///   \cgalParamEnd
/// \cgalNamedParamsEnd
///
/// \return true if the edge is degenerate
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

  if ( traits.equal_3_object()(get(vpmap, target(e, pm)), get(vpmap, source(e, pm))) )
    return true;
  return false;
}

template <typename PolygonMesh>
bool is_degenerate_edge(typename boost::graph_traits<PolygonMesh>::edge_descriptor e,
                        const PolygonMesh& pm)
{
  return is_degenerate_edge(e, pm, parameters::all_default());
}

/// \ingroup PMP_repairing_grp
/// checks whether a triangle face is degenerate.
/// A triangle face is degenerate if its points are collinear.
///
/// @tparam TriangleMesh a model of `FaceGraph`
/// @tparam NamedParameters a sequence of \ref pmp_namedparameters "Named Parameters"
///
/// @param f the triangle face
/// @param tm triangle mesh containing f
/// @param np optional \ref pmp_namedparameters "Named Parameters" described below
///
/// \cgalNamedParamsBegin
///    \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `pmesh`. The type of this map is model of `ReadWritePropertyMap`.
/// If this parameter is omitted, an internal property map for
/// `CGAL::vertex_point_t` should be available in `TriangleMesh`
/// \cgalParamEnd
///    \cgalParamBegin{geom_traits} a geometric traits class instance.
///       The traits class must provide the nested type `Point_3`,
///       and the nested functor :
///         - `Collinear_3` to check whether 3 points are collinear
///   \cgalParamEnd
/// \cgalNamedParamsEnd
///
/// \return true if the triangle face is degenerate
template <typename TriangleMesh, typename NamedParameters>
bool is_degenerate_triangle_face(typename boost::graph_traits<TriangleMesh>::face_descriptor f,
                                 const TriangleMesh& tm,
                                 const NamedParameters& np)
{
  CGAL_assertion(CGAL::is_triangle_mesh(tm));

  using boost::get_param;
  using boost::choose_param;

  typedef typename GetVertexPointMap<TriangleMesh, NamedParameters>::const_type VertexPointMap;
  VertexPointMap vpmap = choose_param(get_param(np, internal_np::vertex_point),
                                      get_const_property_map(vertex_point, tm));
  typedef typename GetGeomTraits<TriangleMesh, NamedParameters>::type Traits;
  Traits traits = choose_param(get_param(np, internal_np::geom_traits), Traits());

  typename boost::graph_traits<TriangleMesh>::halfedge_descriptor hd = halfedge(f,tm);
  const typename Traits::Point_3& p1 = get(vpmap, target( hd, tm) );
  const typename Traits::Point_3& p2 = get(vpmap, target(next(hd, tm), tm) );
  const typename Traits::Point_3& p3 = get(vpmap, source( hd, tm) );
  return traits.collinear_3_object()(p1, p2, p3);

}

template <typename TriangleMesh>
bool is_degenerate_triangle_face(typename boost::graph_traits<TriangleMesh>::face_descriptor f,
                                 const TriangleMesh& tm)
{
  return CGAL::Polygon_mesh_processing::is_degenerate_triangle_face(f, tm, parameters::all_default());
}

/// \ingroup PMP_repairing_grp
/// checks whether a triangle face is needle.
/// A triangle is needle if its longest edge is much longer than the shortest one.
///
/// @tparam TriangleMesh a model of `FaceGraph`
/// @tparam NamedParameters a sequence of \ref pmp_namedparameters "Named Parameters"
///
/// @param f the triangle face
/// @param tm triangle mesh containing f
/// @param threshold the cosine of an angle of f.
///        The threshold is in range [0 1] and corresponds to
///        angles between 0 and 90 degrees.
/// @param np optional \ref pmp_namedparameters "Named Parameters" described below
///
/// \cgalNamedParamsBegin
///    \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `pmesh`. The type of this map is model of `ReadWritePropertyMap`.
/// If this parameter is omitted, an internal property map for
/// `CGAL::vertex_point_t` should be available in `TriangleMesh`
/// \cgalParamEnd
///    \cgalParamBegin{geom_traits} a geometric traits class instance.
///       The traits class must provide the nested type `Point_3`.
///   \cgalParamEnd
/// \cgalNamedParamsEnd
///
/// \return true if the triangle face is a needle
template <typename TriangleMesh, typename NamedParameters>
bool is_needle_triangle_face(typename boost::graph_traits<TriangleMesh>::face_descriptor f,
                             const TriangleMesh& tm,
                             const double threshold,
                             const NamedParameters& np)
{
  CGAL_assertion(CGAL::is_triangle_mesh(tm));
  CGAL_assertion(threshold >= 0);
  CGAL_assertion(threshold <= 1);

  using boost::get_param;
  using boost::choose_param;

  typedef typename GetVertexPointMap<TriangleMesh, NamedParameters>::const_type VertexPointMap;
  VertexPointMap vpmap = choose_param(get_param(np, internal_np::vertex_point),
                                      get_const_property_map(vertex_point, tm));
  typedef typename GetGeomTraits<TriangleMesh, NamedParameters>::type Traits;
  typedef typename Traits::FT FT;
  typedef boost::graph_traits<TriangleMesh> GT;
  typedef typename GT::halfedge_descriptor halfedge_descriptor;
  typedef typename GT::vertex_descriptor vertex_descriptor;
  typedef typename boost::property_traits<VertexPointMap>::value_type Point_type;
  typedef typename Kernel_traits<Point_type>::Kernel::Vector_3 Vector;

  BOOST_FOREACH(halfedge_descriptor h, halfedges_around_face(halfedge(f, tm), tm))
  {
    vertex_descriptor v0 = source(h, tm);
    vertex_descriptor v1 = target(h, tm);
    vertex_descriptor v2 = target(next(h, tm), tm);
    Vector a = get(vpmap, v0) - get (vpmap, v1);
    Vector b = get(vpmap, v2) - get(vpmap, v1);
    FT aa = a.squared_length();
    FT bb = b.squared_length();
    FT squared_dot_ab = ((a*b)*(a*b)) / (aa * bb);

    if(squared_dot_ab > threshold * threshold)
      return true;
  }
  return false;

}

template <typename TriangleMesh>
bool is_needle_triangle_face(typename boost::graph_traits<TriangleMesh>::face_descriptor f,
                             const TriangleMesh& tm,
                             const double threshold)
{
  return is_needle_triangle_face(f, tm, threshold, parameters::all_default());
}

/// \ingroup PMP_repairing_grp
/// checks whether a triangle face is a cap.
/// A triangle is a cap if it has an angle very close to 180 degrees.
///
/// @tparam TriangleMesh a model of `FaceGraph`
/// @tparam NamedParameters a sequence of \ref pmp_namedparameters "Named Parameters"
///
/// @param f the triangle face
/// @param tm triangle mesh containing f
/// @param threshold the cosine of an angle of f.
///        The threshold is in range [-1 0] and corresponds to
///        angles between 90 and 180 degrees.
/// @param np optional \ref pmp_namedparameters "Named Parameters" described below
///
/// \cgalNamedParamsBegin
///    \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `pmesh`. The type of this map is model of `ReadWritePropertyMap`.
/// If this parameter is omitted, an internal property map for
/// `CGAL::vertex_point_t` should be available in `TriangleMesh`
/// \cgalParamEnd
///    \cgalParamBegin{geom_traits} a geometric traits class instance.
///       The traits class must provide the nested type `Point_3`
///   \cgalParamEnd
/// \cgalNamedParamsEnd
///
/// \return true if the triangle face is a cap
template <typename TriangleMesh, typename NamedParameters>
bool is_cap_triangle_face(typename boost::graph_traits<TriangleMesh>::face_descriptor f,
                             const TriangleMesh& tm,
                             const double threshold,
                             const NamedParameters& np)
{
  CGAL_assertion(CGAL::is_triangle_mesh(tm));
  CGAL_assertion(threshold >= -1);
  CGAL_assertion(threshold <= 0);

  using boost::get_param;
  using boost::choose_param;

  typedef typename GetVertexPointMap<TriangleMesh, NamedParameters>::const_type VertexPointMap;
  VertexPointMap vpmap = choose_param(get_param(np, internal_np::vertex_point),
                                      get_const_property_map(vertex_point, tm));
  typedef typename GetGeomTraits<TriangleMesh, NamedParameters>::type Traits;
  typedef typename Traits::FT FT;
  typedef boost::graph_traits<TriangleMesh> GT;
  typedef typename GT::halfedge_descriptor halfedge_descriptor;
  typedef typename GT::vertex_descriptor vertex_descriptor;
  typedef typename boost::property_traits<VertexPointMap>::value_type Point_type;
  typedef typename Kernel_traits<Point_type>::Kernel::Vector_3 Vector;

  BOOST_FOREACH(halfedge_descriptor h, halfedges_around_face(halfedge(f, tm), tm))
  {
    vertex_descriptor v0 = source(h, tm);
    vertex_descriptor v1 = target(h, tm);
    vertex_descriptor v2 = target(next(h, tm), tm);
    Vector a = get(vpmap, v0) - get (vpmap, v1);
    Vector b = get(vpmap, v2) - get(vpmap, v1);
    FT aa = a.squared_length();
    FT bb = b.squared_length();
    FT squared_dot_ab = ((a*b)*(a*b)) / (aa * bb);

    if(squared_dot_ab > threshold * threshold)
      return true;
  }
  return false;
}

template <typename TriangleMesh>
bool is_cap_triangle_face(typename boost::graph_traits<TriangleMesh>::face_descriptor f,
                             const TriangleMesh& tm,
                             const double threshold)
{
  return is_cap_triangle_face(f, tm, threshold, parameters::all_default());
}




} } // end namespaces CGAL and PMP



#endif // CGAL_POLYGON_MESH_PROCESSING_HELPERS_H

