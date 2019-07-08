// Copyright (c) 2019 GeometryFactory (France).
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
// Author(s)     : Sebastien Loriot and Maxime Gimeno

#ifndef CGAL_ORIENT_POLYGON_SOUP_EXTENSION_H
#define CGAL_ORIENT_POLYGON_SOUP_EXTENSION_H

#include <CGAL/license/Polygon_mesh_processing/orientation.h>

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
 * Duplicate each point \a p at which the intersection
 * of an infinitesimally small ball centered at \a p
 * with the polygons incident to it is not a topological disk.
 *
 * @tparam PointRange a model of the concepts `RandomAccessContainer`
 * and `BackInsertionSequence` whose value type is the point type.
 * @tparam PolygonRange a model of the concept `RandomAccessContainer`
 * whose value_type is a model of the concept `RandomAccessContainer`
 * whose value_type is `std::size_t`.
 *
 * @param points points of the soup of polygons. Some additional points might be pushed back to resolve
 *               non-manifoldness or non-orientability issues.
 * @param polygons each element in the vector describes a polygon using the index of the points in `points`.
 *                 If needed the order of the indices of a polygon might be reversed.
 * @return `true`  if the orientation operation succeded.
 * @return `false` if some points were duplicated, thus producing a self-intersecting polyhedron.
 * @sa `orient_polygon_soup()`
 */
template <class PointRange, class PolygonRange>
bool
duplicate_incompatible_edges_in_polygon_soup(PointRange& points,
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
 * Orient each triangle of a triangle soup using the orientation of its
 * closest non degenerate triangle in `tm_ref`.
 * \tparam Concurrency_tag enables sequential versus parallel orientation.
                        Possible values are `Sequential_tag` (the default) and
                        `Parallel_tag`.
 * \tparam PointRange a model of the concepts `RandomAccessContainer`
 * and `BackInsertionSequence` whose value type is the point type.
 * @tparam TriangleRange a model of the concept `RandomAccessContainer`
 * whose value_type is a model of the concept `RandomAccessContainer`
 * whose value_type is `std::size_t`and of size 3.
 * @tparam TriangleMesh a model of `FaceListGraph` and `MutableFaceGraph` .
 *
 * \param tm_ref the reference TriangleMesh.
 * \param points the points of the soup.
 * \param triangles the triangles of the soup.
 */
template <class Concurrency_tag, class PointRange, class TriangleRange,
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
  typedef typename Kernel_traits<Point_3>::Kernel K;
  typedef typename Polygon_mesh_processing::
      GetVertexPointMap<TriangleMesh, NamedParameters>::const_type Vpm;

  Vpm vpm = choose_param(get_param(np, internal_np::vertex_point),
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
      const auto& p0 = points[triangles[fid][0]];
      const auto& p1 = points[triangles[fid][1]];
      const auto& p2 = points[triangles[fid][2]];
      const typename K::Point_3 mid = CGAL::centroid(p0, p1, p2);
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


template <class Concurrency_tag, class PointRange, class TriangleRange,
          class TriangleMesh>
void
orient_triangle_soup_with_reference_triangle_mesh(
  const TriangleMesh& tm_ref,
  PointRange& points,
  TriangleRange& triangles)
{
  orient_triangle_soup_with_reference_triangle_mesh<Concurrency_tag>(tm_ref, points, triangles, CGAL::Polygon_mesh_processing::parameters::all_default());
}

} }//end namespace CGAL::Polygon_mesh_processing
#endif // CGAL_ORIENT_POLYGON_SOUP_EXTENSION_H
