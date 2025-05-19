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
// Author(s)     : Sebastien Loriot
//                 Maxime Gimeno
//                 Mael Rouxel-Labbé

#ifndef CGAL_ORIENT_POLYGON_SOUP_EXTENSION_H
#define CGAL_ORIENT_POLYGON_SOUP_EXTENSION_H

#include <CGAL/license/Polygon_mesh_processing/combinatorial_repair.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_triangle_primitive_3.h>
#include <CGAL/AABB_traits_3.h>

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
 *
 * duplicates each point \a p at which the intersection
 * of an infinitesimally small ball centered at \a p
 * with the polygons incident to it is not a topological disk.
 *
 * @tparam PointRange a model of the concepts `RandomAccessContainer`
 * and `BackInsertionSequence` whose `value_type` is the point type
 * @tparam PolygonRange a model of the concept `RandomAccessContainer`
 * whose `value_type` is a model of the concept `RandomAccessContainer`
 * whose `value_type` is `std::size_t`, and is also a model of `BackInsertionSequence`
 *
 * @param points points of the soup of polygons. Some additional points might be pushed back to resolve
 *               non-manifoldness or non-orientability issues.
 * @param polygons each element in the vector describes a polygon using the indices of the points in `points`.
 *                 If needed the order of the indices of a polygon might be reversed.
 *
 * @return `false` if some points were duplicated, thus producing a self-intersecting surface mesh.
 * @return `true` otherwise.
 *
 * @sa `orient_polygon_soup()`
 * @sa `duplicate_non_manifold_vertices()`
 */
template <class PointRange, class PolygonRange>
bool
duplicate_non_manifold_edges_in_polygon_soup(PointRange& points,
                                             PolygonRange& polygons)
{
  std::size_t inital_nb_pts = points.size();
  typedef CGAL::Polygon_mesh_processing::internal::
    Polygon_soup_orienter<PointRange, PolygonRange> Orienter;

  Default_orientation_visitor visitor;
  Orienter orienter(points, polygons, visitor);
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
 * closest non degenerate triangle in a triangle soup.
 *
 * \tparam Concurrency_tag enables sequential versus parallel orientation.
 *                         Possible values are `Sequential_tag` (the default),
 *                         `Parallel_if_available_tag`, and `Parallel_tag`.
 * \tparam ReferencePointRange a model of the concept `RandomAccessContainer`
 *                             whose `value_type` is the point type
 * \tparam ReferenceTriangleRange a model of the concept `RandomAccessContainer`
 *                                whose `value_type` is a model of the concept `RandomAccessContainer`
 *                                whose `value_type` is `std::size_t` and is of size 3
 * \tparam PointRange a model of the concept `RandomAccessContainer` whose `value_type` is the point type
 * \tparam TriangleRange a model of the concept `RandomAccessContainer`
 *                       whose `value_type` is a model of the concept `RandomAccessContainer`
 *                       whose `value_type` is `std::size_t`and is of size 3.
 * \tparam NamedParameters1 a sequence of \ref bgl_namedparameters "Named Parameters"
 * \tparam NamedParameters2 a sequence of \ref bgl_namedparameters "Named Parameters"
 *
 * \param ref_points the points of the reference soup
 * \param ref_faces triples of indices of points in `ref_points` defining the triangles of the reference soup
 * \param points the points of the soup to be oriented
 * \param faces triples of indices of points in `points` defining the triangles of the soup
 * \param np1 an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
 *
 * \cgalNamedParamsBegin
 *   \cgalParamNBegin{point_map}
 *     \cgalParamDescription{a property map associating points to the elements of the point set `ref_points`}
 *     \cgalParamType{a model of `ReadablePropertyMap` whose key type is the value type
 *                    of the iterator of `ReferencePointRange` and whose value type is `geom_traits::Point_3`}
 *     \cgalParamDefault{`CGAL::Identity_property_map<geom_traits::Point_3>`}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{geom_traits}
 *     \cgalParamDescription{an instance of a geometric traits class}
 *     \cgalParamType{a class model of `Kernel`}
 *     \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
 *   \cgalParamNEnd
 * \cgalNamedParamsEnd
 *
 * \param np2 an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
 *
 * \cgalNamedParamsBegin
 *   \cgalParamNBegin{point_map}
 *     \cgalParamDescription{a property map associating points to the elements of the point set `points`}
 *     \cgalParamType{a model of `ReadablePropertyMap` whose key type is the value type
 *                    of the iterator of `PointRange` and whose value type is `geom_traits::Point_3`}
 *     \cgalParamDefault{`CGAL::Identity_property_map<geom_traits::Point_3>`}
 *   \cgalParamNEnd
 * \cgalNamedParamsEnd
 *
 * \attention The types of points in `ReferencePointRange`, `PointRange`, and `geom_traits` must be the same.
 */
template <class Concurrency_tag = CGAL::Sequential_tag,
          class ReferencePointRange, class ReferenceFaceRange, class PointRange, class FaceRange,
          class NamedParameters1 = parameters::Default_named_parameters,
          class NamedParameters2 = parameters::Default_named_parameters>
void orient_triangle_soup_with_reference_triangle_soup(const ReferencePointRange& ref_points,
                                                       const ReferenceFaceRange& ref_faces,
                                                       const PointRange& points,
                                                       FaceRange& faces,
                                                       const NamedParameters1& np1 = parameters::default_values(),
                                                       const NamedParameters2& np2 = parameters::default_values())
{
  using parameters::choose_parameter;
  using parameters::get_parameter;

  typedef Point_set_processing_3_np_helper<ReferencePointRange, NamedParameters1> NP_helper1;
  typedef typename NP_helper1::Const_point_map PointMap1;

  typedef Point_set_processing_3_np_helper<ReferencePointRange, NamedParameters2> NP_helper2;
  typedef typename NP_helper2::Const_point_map PointMap2;
  typedef typename boost::property_traits<PointMap2>::reference PM2_Point_ref;

  typedef typename boost::property_traits<PointMap1>::value_type Point_3;
  static_assert(std::is_same<Point_3, typename boost::property_traits<PointMap2>::value_type>::value);

  typedef typename CGAL::Kernel_traits<Point_3>::Kernel K;
  typedef typename K::Triangle_3 Triangle;
  typedef typename K::Vector_3 Vector;

  PointMap1 point_map1 = NP_helper1::get_const_point_map(ref_points, np1);
  PointMap2 point_map2 = NP_helper2::get_const_point_map(points, np2);

  K k = choose_parameter<K>(get_parameter(np1, internal_np::geom_traits));

  typename K::Construct_centroid_3 centroid = k.construct_centroid_3_object();
  typename K::Construct_vector_3 vector = k.construct_vector_3_object();
  typename K::Is_degenerate_3 is_degenerate = k.is_degenerate_3_object();
  typename K::Compute_scalar_product_3 scalar_product = k.compute_scalar_product_3_object();
  typename K::Construct_cross_product_vector_3 cross_product = k.construct_cross_product_vector_3_object();

  // build a tree filtering degenerate faces
  std::vector<Triangle> ref_triangles;
  ref_triangles.reserve(ref_faces.size());
  for(const auto& f : ref_faces)
  {
    Triangle tr(get(point_map1, ref_points[f[0]]),
                get(point_map1, ref_points[f[1]]),
                get(point_map1, ref_points[f[2]]));

    if(!is_degenerate(tr))
      ref_triangles.emplace_back(tr);
  }

  typedef typename std::vector<Triangle>::const_iterator Iterator;
  typedef CGAL::AABB_triangle_primitive_3<K, Iterator> Primitive;
  typedef CGAL::AABB_traits_3<K, Primitive> Tree_traits;

  CGAL::AABB_tree<Tree_traits> tree(ref_triangles.begin(), ref_triangles.end());

  // now orient the faces
  tree.build();
  tree.accelerate_distance_queries();

  auto process_facet = [&](const std::size_t fid)
  {
    PM2_Point_ref p0 = get(point_map2, points[faces[fid][0]]);
    PM2_Point_ref p1 = get(point_map2, points[faces[fid][1]]);
    PM2_Point_ref p2 = get(point_map2, points[faces[fid][2]]);
    const Point_3 mid = centroid(p0, p1, p2);

    auto pt_and_ref_tr = tree.closest_point_and_primitive(mid);
    const Triangle& ref_tr = *(pt_and_ref_tr.second);
    Vector ref_n = cross_product(vector(ref_tr[0], ref_tr[1]),
                                 vector(ref_tr[0], ref_tr[2]));
    if(is_negative(scalar_product(ref_n, cross_product(vector(p0, p1), vector(p0, p2)))))
      std::swap(faces[fid][1], faces[fid][2]);
  };

#if !defined(CGAL_LINKED_WITH_TBB)
  static_assert (!std::is_convertible<Concurrency_tag, CGAL::Parallel_tag>::value,
                 "Parallel_tag is enabled but TBB is unavailable.");
#else
  if(std::is_convertible<Concurrency_tag,CGAL::Parallel_tag>::value)
    tbb::parallel_for(std::size_t(0), faces.size(), std::size_t(1), process_facet);
  else
#endif
    std::for_each(boost::counting_iterator<std::size_t>(0),
                  boost::counting_iterator<std::size_t>(faces.size()),
                  process_facet);
}

/*!
 * \ingroup PMP_orientation_grp
 *
 * orients each triangle of a triangle soup using the orientation of its
 * closest non degenerate triangle in `tm_ref`.
 *
 * \tparam Concurrency_tag enables sequential versus parallel orientation.
 *                         Possible values are `Sequential_tag` (the default),
 *                         `Parallel_if_available_tag`, and `Parallel_tag`.
 * \tparam PointRange a model of the concept `RandomAccessContainer` whose value type is the point type
 * \tparam TriangleRange a model of the concept `RandomAccessContainer`
 *                       whose `value_type` is a model of the concept `RandomAccessContainer`
 *                       whose `value_type` is `std::size_t`and of size 3
 * \tparam TriangleMesh a model of `FaceListGraph`
 *
 * \param tm_ref the reference triangle_mesh
 * \param points the points of the soup
 * \param triangles the triangles of the soup
 * \param np1 an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
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
 *     \cgalParamType{a class model of `Kernel`}
 *     \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
 *   \cgalParamNEnd
 * \cgalNamedParamsEnd
 *
 * \param np2 an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
 *
 * \cgalNamedParamsBegin
 *   \cgalParamNBegin{point_map}
 *     \cgalParamDescription{a property map associating points to the elements of the point set `points`}
 *     \cgalParamType{a model of `ReadablePropertyMap` whose key type is the value type
 *                    of the iterator of `PointRange` and whose value type is `geom_traits::Point_3`}
 *     \cgalParamDefault{`CGAL::Identity_property_map<geom_traits::Point_3>`}
 *   \cgalParamNEnd
 * \cgalNamedParamsEnd
 *
 * \attention The types of points in `PointRange`, `geom_traits`, and `vertex_point_map` must be the same.
 *
 * \sa `orient_polygon_soup()`
 */
template <class Concurrency_tag = Sequential_tag,
          class PointRange, class TriangleRange, class TriangleMesh,
          class NamedParameters1 = parameters::Default_named_parameters,
          class NamedParameters2 = parameters::Default_named_parameters>
void orient_triangle_soup_with_reference_triangle_mesh(const TriangleMesh& tm_ref,
                                                       const PointRange& points,
                                                       TriangleRange& triangles,
                                                       const NamedParameters1& np1 = parameters::default_values(),
                                                       const NamedParameters2& np2 = parameters::default_values())
{
  namespace PMP = CGAL::Polygon_mesh_processing;

  using parameters::choose_parameter;
  using parameters::get_parameter;

  typedef boost::graph_traits<TriangleMesh> GrT;
  typedef typename GrT::face_descriptor face_descriptor;
  typedef typename GetGeomTraits<TriangleMesh, NamedParameters1>::type K;
  typedef typename K::Vector_3 Vector;

  typedef typename GetVertexPointMap<TriangleMesh, NamedParameters1>::const_type VPM;
  typedef typename boost::property_traits<VPM>::value_type Point_3;

  VPM vpm = choose_parameter(get_parameter(np1, internal_np::vertex_point),
                             get_const_property_map(CGAL::vertex_point, tm_ref));

  typedef Point_set_processing_3_np_helper<PointRange, NamedParameters2> NP_helper;
  typedef typename NP_helper::Const_point_map PointMap;
  typedef typename boost::property_traits<PointMap>::reference PM2_Point_ref;

  PointMap point_map = NP_helper::get_const_point_map(points, np2);

  static_assert(std::is_same<Point_3, typename boost::property_traits<PointMap>::value_type>::value);

  K k = choose_parameter<K>(get_parameter(np1, internal_np::geom_traits));

  typename K::Construct_centroid_3 centroid = k.construct_centroid_3_object();
  typename K::Construct_vector_3 vector = k.construct_vector_3_object();
  typename K::Compute_scalar_product_3 scalar_product = k.compute_scalar_product_3_object();
  typename K::Construct_cross_product_vector_3 cross_product = k.construct_cross_product_vector_3_object();

  // build a tree filtering degenerate faces
  typedef std::function<bool(face_descriptor)> Face_predicate;
  Face_predicate is_not_deg = [&](face_descriptor f)
  {
    return !PMP::is_degenerate_triangle_face(f, tm_ref, np1);
  };

  typedef CGAL::AABB_face_graph_triangle_primitive<TriangleMesh, VPM> Primitive;
  typedef CGAL::AABB_traits_3<K, Primitive> Tree_traits;

  boost::filter_iterator<Face_predicate, typename GrT::face_iterator>
      begin(is_not_deg, faces(tm_ref).begin(), faces(tm_ref).end()),
      end(is_not_deg, faces(tm_ref).end(), faces(tm_ref).end());

  CGAL::AABB_tree<Tree_traits> tree(begin, end, tm_ref, vpm);

  // now orient the faces
  tree.build();
  tree.accelerate_distance_queries();
  auto process_facet = [&](const std::size_t fid)
  {
    PM2_Point_ref p0 = get(point_map, points[triangles[fid][0]]);
    PM2_Point_ref p1 = get(point_map, points[triangles[fid][1]]);
    PM2_Point_ref p2 = get(point_map, points[triangles[fid][2]]);
    const Point_3 mid = centroid(p0, p1, p2);

    std::pair<Point_3, face_descriptor> pt_and_f = tree.closest_point_and_primitive(mid);
    Vector face_ref_normal = PMP::compute_face_normal(pt_and_f.second, tm_ref,
                                                      CGAL::parameters::vertex_point_map(vpm));
    if(is_negative(scalar_product(face_ref_normal, cross_product(vector(p0,p1), vector(p0, p2)))))
      std::swap(triangles[fid][1], triangles[fid][2]);
  };

#if !defined(CGAL_LINKED_WITH_TBB)
  static_assert (!std::is_convertible<Concurrency_tag, CGAL::Parallel_tag>::value,
                 "Parallel_tag is enabled but TBB is unavailable.");
#else
  if (std::is_convertible<Concurrency_tag,CGAL::Parallel_tag>::value)
    tbb::parallel_for(std::size_t(0), triangles.size(), std::size_t(1), process_facet);
  else
#endif
    std::for_each(boost::counting_iterator<std::size_t> (0),
                  boost::counting_iterator<std::size_t> (triangles.size()),
                  process_facet);
}

} // namespace Polygon_mesh_processing
} // namespace CGAL

#endif // CGAL_ORIENT_POLYGON_SOUP_EXTENSION_H
