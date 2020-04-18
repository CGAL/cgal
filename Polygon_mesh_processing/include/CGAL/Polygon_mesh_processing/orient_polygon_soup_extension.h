// Copyright (c) 2019 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Sebastien Loriot and Maxime Gimeno

#ifndef CGAL_ORIENT_POLYGON_SOUP_EXTENSION_H
#define CGAL_ORIENT_POLYGON_SOUP_EXTENSION_H

#include <CGAL/license/Polygon_mesh_processing/repair.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_traits.h>

#include <CGAL/Polygon_mesh_processing/shape_predicates.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>

#include <boost/iterator/counting_iterator.hpp>
#include <boost/iterator/filter_iterator.hpp>

#ifdef CGAL_LINKED_WITH_TBB
#include <tbb/parallel_for.h>
#endif // CGAL_LINKED_WITH_TBB


namespace CGAL {

namespace Polygon_mesh_processing {

/*!
 * \ingroup PMP_orientation_grp
 * duplicates each point \a p at which the intersection
 * of an infinitesimally small ball centered at \a p
 * with the polygons incident to it is not a topological disk.
 *
 * @tparam PointRange a model of the concepts `RandomAccessContainer`
 * and `BackInsertionSequence` whose `value_type` is the point type.
 * @tparam PolygonRange a model of the concept `RandomAccessContainer`
 * whose `value_type` is a model of the concept `RandomAccessContainer`
 * whose `value_type` is `std::size_t`, and is also a model of `BackInsertionSequence`.
 *
 * @param points points of the soup of polygons. Some additional points might be pushed back to resolve
 *               non-manifoldness or non-orientability issues.
 * @param polygons each element in the vector describes a polygon using the indices of the points in `points`.
 *                 If needed the order of the indices of a polygon might be reversed.
 * @return `false` if some points were duplicated, thus producing a self-intersecting surface mesh.
 * @return `true` otherwise.
 * @sa `orient_polygon_soup()`
 */
template <class PointRange, class PolygonRange>
bool
duplicate_non_manifold_edges_in_polygon_soup(PointRange& points,
                                            PolygonRange& polygons)
{
  std::size_t inital_nb_pts = points.size();
  typedef CGAL::Polygon_mesh_processing::internal::
    Polygon_soup_orienter<PointRange, PolygonRange> Orienter;

  Orienter orienter(points, polygons);
  orienter.fill_edge_map();
  // make edges to duplicate
  for(std::size_t i1=0;i1<points.size();++i1)
    for(const typename Orienter::Internal_map_type::value_type& i2_and_pids : orienter.edges[i1])
      if (i2_and_pids.second.size() > 1)
        orienter.set_edge_marked(i1,i2_and_pids.first,orienter.marked_edges);
  orienter.duplicate_singular_vertices();

  return inital_nb_pts==points.size();
}

/*!
 * \ingroup PMP_orientation_grp
 *
 * orients each triangle of a triangle soup using the orientation of its
 * closest non degenerate triangle in `tm_ref`.
 *
 * \tparam Concurrency_tag enables sequential versus parallel orientation.
                        Possible values are `Sequential_tag` (the default),
                        `Parallel_if_available_tag`, and `Parallel_tag`.
 * \tparam PointRange a model of the concepts `RandomAccessContainer`
 * and `BackInsertionSequence` whose value type is the point type.
 * @tparam TriangleRange a model of the concept `RandomAccessContainer`
 * whose `value_type` is a model of the concept `RandomAccessContainer`
 * whose `value_type` is `std::size_t`and of size 3.
 * @tparam TriangleMesh a model of `FaceListGraph` and `MutableFaceGraph` .
 *
 * \param tm_ref the reference triangle_mesh.
 * \param points the points of the soup.
 * \param triangles the triangles of the soup.
 * @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
 *
 * \cgalNamedParamsBegin
 *   \cgalParamNBegin{vertex_point_map}
 *     \cgalParamDescription{a property map associating points to the vertices of `tm_ref`}
 *     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<TriangleMesh>::%vertex_descriptor`
 *                    as key type and `%Point_3` as value type}
 *     \cgalParamDefault{`boost::get(CGAL::vertex_point, tm_ref)`}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{geom_traits}
 *     \cgalParamDescription{an instance of a geometric traits class}
 *     \cgalParamType{The traits class must provide the nested functor `Collinear_3` to check whether three points are collinear. }
 *     \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
 *     \cgalParamExtra{The geometric traits class must be compatible with the vertex point type.}
 *   \cgalParamNEnd
 * \cgalNamedParamsEnd
 *
 * \attention The types of points in `PointRange`, `geom_traits` and `vertex_point_map` must be the same.
 */

template <class Concurrency_tag = Sequential_tag, class PointRange, class TriangleRange,
          class TriangleMesh, class NamedParameters>
void
orient_triangle_soup_with_reference_triangle_mesh(
    const TriangleMesh& tm_ref,
    PointRange& points,
    TriangleRange& triangles,
    const NamedParameters& np)
{
  namespace PMP = CGAL::Polygon_mesh_processing;

  typedef boost::graph_traits<TriangleMesh> GrT;
  typedef typename GrT::face_descriptor face_descriptor;
  typedef typename PointRange::value_type Point_3;
  typedef typename GetGeomTraits<TriangleMesh, NamedParameters>::type K;
  typedef typename
  GetVertexPointMap<TriangleMesh, NamedParameters>::const_type Vpm;

  Vpm vpm = parameters::choose_parameter(parameters::get_parameter(np, internal_np::vertex_point),
                                         get_const_property_map(CGAL::vertex_point, tm_ref));


  typedef std::function<bool(face_descriptor)> Face_predicate;
  Face_predicate is_not_deg =
      [&tm_ref, np](face_descriptor f)
  {
    return !PMP::is_degenerate_triangle_face(f, tm_ref, np);
  };

  // build a tree filtering degenerate faces
  typedef CGAL::AABB_face_graph_triangle_primitive<TriangleMesh, Vpm> Primitive;
  typedef CGAL::AABB_traits<K, Primitive> Tree_traits;

  boost::filter_iterator<Face_predicate, typename GrT::face_iterator>
      begin(is_not_deg, faces(tm_ref).begin(), faces(tm_ref).end()),
      end(is_not_deg, faces(tm_ref).end(), faces(tm_ref).end());

  CGAL::AABB_tree<Tree_traits> tree(begin, end, tm_ref, vpm);

  // now orient the faces
  tree.build();
  tree.accelerate_distance_queries();
  auto process_facet =
      [&points, &tree, &tm_ref, &triangles](std::size_t fid) {
    const Point_3& p0 = points[triangles[fid][0]];
    const Point_3& p1 = points[triangles[fid][1]];
    const Point_3& p2 = points[triangles[fid][2]];
    const Point_3 mid = CGAL::centroid(p0, p1, p2);
    std::pair<Point_3, face_descriptor> pt_and_f =
        tree.closest_point_and_primitive(mid);
    auto face_ref_normal = PMP::compute_face_normal(pt_and_f.second, tm_ref);
    if(face_ref_normal * cross_product(p1-p0, p2-p0) < 0) {
      std::swap(triangles[fid][1], triangles[fid][2]);
    }
  };

#if !defined(CGAL_LINKED_WITH_TBB)
  CGAL_static_assertion_msg (!(boost::is_convertible<Concurrency_tag, CGAL::Parallel_tag>::value),
                             "Parallel_tag is enabled but TBB is unavailable.");
#else
  if (boost::is_convertible<Concurrency_tag,CGAL::Parallel_tag>::value)
    tbb::parallel_for(std::size_t(0), triangles.size(), std::size_t(1), process_facet);
  else
#endif
    std::for_each(
          boost::counting_iterator<std::size_t> (0),
          boost::counting_iterator<std::size_t> (triangles.size()),
          process_facet);
}


template <class Concurrency_tag = Sequential_tag, class PointRange, class TriangleRange,
          class TriangleMesh>
void
orient_triangle_soup_with_reference_triangle_mesh(
  const TriangleMesh& tm_ref,
  PointRange& points,
  TriangleRange& triangles)
{
  orient_triangle_soup_with_reference_triangle_mesh<Concurrency_tag>(tm_ref, points, triangles, CGAL::parameters::all_default());
}

}}//end namespace CGAL::Polygon_mesh_processing
#endif // CGAL_ORIENT_POLYGON_SOUP_EXTENSION_H
