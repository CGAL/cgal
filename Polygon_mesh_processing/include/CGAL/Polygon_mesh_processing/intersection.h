// Copyright (c) 2016 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Maxime Gimeno and Sebastien Loriot

#ifndef CGAL_POLYGON_MESH_PROCESSING_INTERSECTION_H
#define CGAL_POLYGON_MESH_PROCESSING_INTERSECTION_H

#include <CGAL/license/Polygon_mesh_processing/corefinement.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/Polygon_mesh_processing/internal/soup_and_face_graph_tree_helper.h>
#include <CGAL/AABB_segment_primitive_3.h>

#include <CGAL/boost/iterator/counting_iterator.hpp>
#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/Polygon_mesh_processing/bbox.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/Side_of_triangle_mesh.h>

#include <boost/range/has_range_iterator.hpp>
#include <boost/iterator/function_output_iterator.hpp>

#include <exception>
#include <iterator>
#include <utility>
#include <vector>
#include <type_traits>

namespace CGAL {
namespace Polygon_mesh_processing {
namespace internal {

template<class TM,
         class Kernel,
         class Box,
         class OutputIterator,
         class VertexPointMap1,
         class VertexPointMap2>
struct Intersect_faces
{
  // typedefs
  typedef typename Kernel::Segment_3    Segment;
  typedef typename Kernel::Triangle_3   Triangle;

  typedef typename boost::graph_traits<TM>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::property_map<TM, boost::vertex_point_t>::const_type Ppmap;

  // members
  const TM& m_tm1;
  const VertexPointMap1 m_vpmap1;
  const TM& m_tm2;
  const VertexPointMap2 m_vpmap2;
  mutable OutputIterator  m_iterator;

  typename Kernel::Construct_triangle_3 triangle_functor;
  typename Kernel::Do_intersect_3       do_intersect_3_functor;

  Intersect_faces(const TM& tm1, const TM& tm2,
                   OutputIterator it,
                   VertexPointMap1 vpmap1, VertexPointMap2 vpmap2,
                   const Kernel& kernel)
    :
      m_tm1(tm1),
      m_vpmap1(vpmap1),
      m_tm2(tm2),
      m_vpmap2(vpmap2),
      m_iterator(it),
      triangle_functor(kernel.construct_triangle_3_object()),
      do_intersect_3_functor(kernel.do_intersect_3_object())
  { }

  void operator()(const Box* b, const Box* c) const
  {
    halfedge_descriptor h = halfedge(b->info(), m_tm1);
    halfedge_descriptor g = halfedge(c->info(), m_tm2);


    // check for geometric intersection
    Triangle t1 = triangle_functor( get(m_vpmap1, target(h,m_tm1)),
                                    get(m_vpmap1, target(next(h,m_tm1),m_tm1)),
                                    get(m_vpmap1, source(h,m_tm1)));

    Triangle t2 = triangle_functor( get(m_vpmap2, target(g,m_tm2)),
                                    get(m_vpmap2, target(next(g,m_tm2),m_tm2)),
                                    get(m_vpmap2, source(g,m_tm2)));
    if(do_intersect_3_functor(t1, t2)){
      *m_iterator++ = std::make_pair(b->info(), c->info());
    }
  } // end operator ()
}; // end struct Intersect_faces

template<class TM,
         class Kernel,
         class Box,
         class OutputIterator,
         class Polyline,
         class VertexPointMap>
struct Intersect_face_polyline
{
  // wrapper to check whether anything is inserted to output iterator

  // typedefs
  typedef typename Kernel::Segment_3    Segment;
  typedef typename Kernel::Triangle_3   Triangle;

  typedef typename boost::graph_traits<TM>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<TM>::face_descriptor face_descriptor;
  typedef typename boost::property_map<TM, boost::vertex_point_t>::const_type Ppmap;
  typedef typename boost::property_traits<Ppmap>::value_type Point;


  // members
  const TM& m_tm;
  const std::vector<face_descriptor>& faces;
  const VertexPointMap m_vpmap;
  const Polyline& polyline;
  mutable OutputIterator  m_iterator;

  typename Kernel::Construct_segment_3  segment_functor;
  typename Kernel::Construct_triangle_3 triangle_functor;
  typename Kernel::Do_intersect_3       do_intersect_3_functor;

  Intersect_face_polyline(const TM& tm,
                           const std::vector<face_descriptor>& faces,
                           const Polyline& polyline,
                           OutputIterator it,
                           VertexPointMap vpmap,
                           const Kernel& kernel)
    :
      m_tm(tm),
      faces(faces),
      m_vpmap(vpmap),
      polyline(polyline),
      m_iterator(it),
      segment_functor(kernel.construct_segment_3_object()),
      triangle_functor(kernel.construct_triangle_3_object()),
      do_intersect_3_functor(kernel.do_intersect_3_object())
  { }

  void operator()(const Box* b, const Box* c) const
  {
    halfedge_descriptor h = halfedge(faces[b->info()], m_tm);


    // check for geometric intersection
    Triangle t = triangle_functor( get(m_vpmap, target(h,m_tm)),
                                   get(m_vpmap, target(next(h,m_tm),m_tm)),
                                   get(m_vpmap, source(h,m_tm)));

    Segment s = segment_functor(polyline[c->info()], polyline[c->info() + 1]);
    if(do_intersect_3_functor(t, s)){
      *m_iterator++ = std::make_pair(b->info(), c->info());
    }
  } // end operator ()
}; // end struct Intersect_face_polyline

template<class TM,
         class Kernel,
         class Box,
         class PolylineRange,
         class OutputIterator,
         class VertexPointMap>
struct Intersect_face_polylines
{
  // wrapper to check whether anything is inserted to output iterator

  // typedefs
  typedef typename Kernel::Segment_3    Segment;
  typedef typename Kernel::Triangle_3   Triangle;

  typedef typename boost::graph_traits<TM>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<TM>::face_descriptor face_descriptor;
  typedef typename boost::property_map<TM, boost::vertex_point_t>::const_type Ppmap;
  typedef typename boost::property_traits<Ppmap>::value_type Point;

  // members
  const TM& m_tm;
  const std::vector<face_descriptor>& faces;
  const VertexPointMap m_vpmap;
  const PolylineRange& polylines;
  mutable OutputIterator  m_iterator;

  typename Kernel::Construct_segment_3  segment_functor;
  typename Kernel::Construct_triangle_3 triangle_functor;
  typename Kernel::Do_intersect_3       do_intersect_3_functor;

  Intersect_face_polylines(const TM& tm,
                           const std::vector<face_descriptor>& faces,
                           const PolylineRange& polylines,
                           OutputIterator it,
                           VertexPointMap vpmap,
                           const Kernel& kernel)
    :
      m_tm(tm),
      faces(faces),
      m_vpmap(vpmap),
      polylines(polylines),
      m_iterator(it),
      segment_functor(kernel.construct_segment_3_object()),
      triangle_functor(kernel.construct_triangle_3_object()),
      do_intersect_3_functor(kernel.do_intersect_3_object())
  { }

  void operator()(const Box* b, const Box* c) const
  {
    halfedge_descriptor h = halfedge(faces[b->info().second], m_tm);


    // check for geometric intersection
    Triangle t = triangle_functor( get(m_vpmap, target(h,m_tm)),
                                   get(m_vpmap, target(next(h,m_tm),m_tm)),
                                   get(m_vpmap, source(h,m_tm)));

    Segment s = segment_functor(polylines[c->info().first][c->info().second], polylines[c->info().first][c->info().second + 1]);
    if(do_intersect_3_functor(t, s)){
      *m_iterator++ = std::make_pair(b->info().second, c->info());
    }
  } // end operator ()
}; // end struct Intersect_face_polylines


template<class Polyline,
         class Kernel,
         class Box,
         class OutputIterator>
struct Intersect_polylines
{
  // typedefs
  typedef typename Kernel::Segment_3    Segment;
  typedef typename Kernel::Point_3 Point;


  // members
  const Polyline& polyline1;
  const Polyline& polyline2;
  mutable OutputIterator  m_iterator;

  typename Kernel::Construct_segment_3  segment_functor;
  typename Kernel::Do_intersect_3       do_intersect_3_functor;

  Intersect_polylines(const Polyline& polyline1,
                      const Polyline& polyline2,
                      OutputIterator it,
                      const Kernel& kernel)
    :
      polyline1(polyline1),
      polyline2(polyline2),
      m_iterator(it),
      segment_functor(kernel.construct_segment_3_object()),
      do_intersect_3_functor(kernel.do_intersect_3_object())
  { }

  void operator()(const Box* b, const Box* c) const
  {


    // check for geometric intersection

    Segment s1 = segment_functor(polyline1[b->info()], polyline1[b->info() + 1]);
    Segment s2 = segment_functor(polyline2[c->info()], polyline2[c->info() + 1]);
    if(do_intersect_3_functor(s1, s2)){
      *m_iterator++ = std::make_pair(b->info(), c->info());
    }
  } // end operator ()
}; // end struct Intersect_polylines

template<class PolylineRange,
         class Kernel,
         class Box,
         class OutputIterator>
struct Intersect_polyline_ranges
{
  // typedefs
  typedef typename Kernel::Segment_3    Segment;
  typedef typename Kernel::Point_3 Point;


  // members
  const PolylineRange& polyline1;
  const PolylineRange& polyline2;
  mutable OutputIterator  m_iterator;

  typename Kernel::Construct_segment_3  segment_functor;
  typename Kernel::Do_intersect_3       do_intersect_3_functor;

  Intersect_polyline_ranges(const PolylineRange& polyline1,
                            const PolylineRange& polyline2,
                            OutputIterator it,
                            const Kernel& kernel)
    :
      polyline1(polyline1),
      polyline2(polyline2),
      m_iterator(it),
      segment_functor(kernel.construct_segment_3_object()),
      do_intersect_3_functor(kernel.do_intersect_3_object())
  { }

  void operator()(const Box* b, const Box* c) const
  {


    // check for geometric intersection

    Segment s1 = segment_functor(polyline1[b->info().first][b->info().second], polyline1[b->info().first][b->info().second + 1]);
    Segment s2 = segment_functor(polyline2[c->info().first][c->info().second], polyline2[c->info().first][c->info().second + 1]);

    if(do_intersect_3_functor(s1, s2)){
      *m_iterator++ = std::make_pair(b->info(), c->info());
    }
  } // end operator ()
}; // end struct Intersect_polyline_ranges

struct Throw_at_first_output {
  class Throw_at_first_output_exception: public std::exception
  { };

  template<class T>
  void operator()(const T& /* t */) const {
    throw Throw_at_first_output_exception();
  }
};

// Note this is not officially documented
/* \ingroup PMP_intersection_grp
 *
 * reports all the pairs of faces intersecting between two triangulated surface meshes.
 *
 * This function depends on the package \ref PkgBoxIntersectionD.
 *
 * @pre \link CGAL::is_triangle_mesh `CGAL::is_triangle_mesh(tm1)` \endlink
 * @pre \link CGAL::is_triangle_mesh `CGAL::is_triangle_mesh(tm2)` \endlink
 *
 * \tparam TriangleMesh a model of `FaceListGraph`
 * \tparam FaceRange range of `boost::graph_traits<TriangleMesh>::%face_descriptor`,
 *  model of `RandomAccessRange`.
 * \tparam OutputIterator a model of `OutputIterator` holding objects of type
 *   `std::pair<boost::graph_traits<TriangleMesh>::%face_descriptor, boost::graph_traits<TriangleMesh>::%face_descriptor>`
 * \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
 *
 * \param face_range1 the range of faces of `tm1` to check for intersections.
 * \param face_range2 the range of faces of `tm2` to check for intersections.
 * \param tm1 the first triangulated surface mesh.
 * \param tm2 the second triangulated surface mesh.
 * \param out output iterator to be filled with all pairs of faces that intersect.
 *  First and second element in the pairs correspond to faces of `tm1` and `tm2` respectively
 * \param np1 an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
 * \param np2 an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
 *
 * \cgalNamedParamsBegin
 *   \cgalParamNBegin{vertex_point_map}
 *     \cgalParamDescription{a property map associating points to the vertices of `tm1` (`tm2`)}
 *     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<TriangleMesh>::%vertex_descriptor`
 *                    as key type and `%Point_3` as value type}
 *     \cgalParamDefault{`boost::get(CGAL::vertex_point, tm1 (tm2))`}
 *     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
 *                     should be available for the vertices of `tm1` (`tm2`)}
 *     \cgalParamExtra{Both vertex point maps must have the same value type}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{geom_traits}
 *     \cgalParamDescription{an instance of a geometric traits class}
 *     \cgalParamType{a class model of `PMPSelfIntersectionTraits`}
 *     \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
 *     \cgalParamExtra{The geometric traits class must be compatible with the vertex point type.}
 *     \cgalParamExtra{np1 only}
 *   \cgalParamNEnd
 * \cgalNamedParamsEnd
 *
 * \return `out`
 */
template <class TriangleMesh,
          class FaceRange,
          class OutputIterator,
          class NamedParameters1,
          class NamedParameters2>
OutputIterator
compute_face_face_intersection(const FaceRange& face_range1,
                               const FaceRange& face_range2,
                               const TriangleMesh& tm1,
                               const TriangleMesh& tm2,
                               OutputIterator out,
                               const NamedParameters1& np1,
                               const NamedParameters2& np2)
{
  using parameters::get_parameter;
  using parameters::choose_parameter;

  CGAL_precondition(CGAL::is_triangle_mesh(tm1));
  CGAL_precondition(CGAL::is_triangle_mesh(tm2));

  using Concurrency_tag = typename internal_np::Lookup_named_param_def <
                                          internal_np::concurrency_tag_t,
                                          NamedParameters1,
                                          Sequential_tag
                                        > ::type;

  using GT = typename GetGeomTraits<TriangleMesh, NamedParameters1>::type;

  using VPM1 = typename GetVertexPointMap<TriangleMesh, NamedParameters1>::const_type;
  using VPM2 = typename GetVertexPointMap<TriangleMesh, NamedParameters2>::const_type;

  using face_descriptor = typename boost::graph_traits<TriangleMesh>::face_descriptor;

  using AABB_tree_helper_1 = internal::AABB_tree_graph_helper<TriangleMesh, GT, VPM1>;
  using AABB_tree_helper_2 = internal::AABB_tree_graph_helper<TriangleMesh, GT, VPM2>;
  using Tree_1 = typename AABB_tree_helper_1::Tree;
  using Tree_2 = typename AABB_tree_helper_2::Tree;
  AABB_tree_helper_1 helper_1;
  AABB_tree_helper_2 helper_2;

  // compute intersections filtered out by boxes
  CGAL::Bbox_3 b1 = CGAL::Polygon_mesh_processing::bbox(tm1, np1),
               b2 = CGAL::Polygon_mesh_processing::bbox(tm2, np2);

  if(!CGAL::do_overlap(b1, b2))
    return out;
  Bbox_3 bb((std::max)(b1.xmin(),b2.xmin()), (std::max)(b1.ymin(),b2.ymin()), (std::max)(b1.zmin(),b2.zmin()),
            (std::min)(b1.xmax(),b2.xmax()), (std::min)(b1.ymax(),b2.ymax()), (std::min)(b1.zmax(),b2.zmax()));

  VPM1 vpm1 = choose_parameter(get_parameter(np1, internal_np::vertex_point),
                               get_const_property_map(boost::vertex_point, tm1));
  VPM2 vpm2 = choose_parameter(get_parameter(np2, internal_np::vertex_point),
                               get_const_property_map(boost::vertex_point, tm2));
  static_assert(
      (std::is_same<
       typename boost::property_traits<VPM1>::value_type,
       typename boost::property_traits<VPM2>::value_type
       >::value) );

  std::vector<face_descriptor> face_to_test1;
  std::vector<face_descriptor> face_to_test2;

  for(face_descriptor f : face_range1)
    if(do_overlap(Polygon_mesh_processing::face_bbox(f, tm1), bb))
      face_to_test1.emplace_back(f);
  for(face_descriptor f : face_range2)
    if(do_overlap(Polygon_mesh_processing::face_bbox(f, tm2), bb))
      face_to_test2.emplace_back(f);

  Tree_1 tree1;
  Tree_2 tree2;
  helper_1.template build<Concurrency_tag>(tree1, tm1, vpm1);
  helper_2.template build<Concurrency_tag>(tree2, tm2, vpm2);

  CGAL::AABB_trees::all_pairs_of_intersecting_primitives(tree1, tree2, out, parameters::concurrency_tag(Concurrency_tag()));
  return out;
}

// Note this is not officially documented
/* \ingroup PMP_intersection_grp
 *
 * reports all the pairs of segments and faces intersecting between
 * a triangulated surface mesh and a polyline.
 *
 * This function depends on the package \ref PkgBoxIntersectionD.
 *
 * \attention If a polyline vertex intersects a face, the intersection will
 * be reported twice (or more if it is on a vertex, edge, or point).
 *
 * @pre \link CGAL::is_triangle_mesh `CGAL::is_triangle_mesh(tm)` \endlink
 *
 * \tparam TriangleMesh a model of `FaceListGraph`
 * \tparam FaceRange range of `boost::graph_traits<TriangleMesh>::%face_descriptor`,
 *  model of `RandomAccessRange`.
 * \tparam Polyline a `RandomAccessRange` of points. The point type of the range must be
 * the same as the value type of the vertex point map.
 * \tparam OutputIterator a model of `OutputIterator` holding objects of type
 *   `std::pair<std::size_t, std::size_t>`.This `OutputIterator` will hold the position of the
 *  elements in their respective range. In the case of the polyline, this position is the index of
 * the segment intersecting the face (which is the index  of the first point of the
 * segment following the range order.)
 * \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
 *
 * \param face_range the range of faces of `tm` to check for intersections.
 * \param polyline the polyline to check for intersections.
 * \param tm the triangulated surface mesh to check for intersections.
 * \param out output iterator to be filled with all pairs of face-segment that intersect
 * \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
 *
 * \cgalNamedParamsBegin
 *   \cgalParamNBegin{vertex_point_map}
 *     \cgalParamDescription{a property map associating points to the vertices of `tm`}
 *     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<TriangleMesh>::%vertex_descriptor` as key type and `%Point_3` as value type}
 *     \cgalParamDefault{`boost::get(CGAL::vertex_point, tm)`}
 *     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
 *                     should be available for the vertices of `tm`.}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{geom_traits}
 *     \cgalParamDescription{an instance of a geometric traits class}
 *     \cgalParamType{a class model of `PMPSelfIntersectionTraits`}
 *     \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
 *     \cgalParamExtra{The geometric traits class must be compatible with the vertex point type.}
 *   \cgalParamNEnd
 * \cgalNamedParamsEnd
 *
 * \return `out`
 */
template <class TriangleMesh,
          class FaceRange,
          class Polyline,
          class OutputIterator,
          class NamedParameters>
OutputIterator
compute_face_polyline_intersection(const FaceRange& face_range,
                                   const Polyline& polyline,
                                   const TriangleMesh& tm,
                                   OutputIterator out,
                                   const NamedParameters& np)
{
  using parameters::get_parameter;
  using parameters::choose_parameter;

  CGAL_precondition(CGAL::is_triangle_mesh(tm));

  using Concurrency_tag = typename internal_np::Lookup_named_param_def <
                                          internal_np::concurrency_tag_t,
                                          NamedParameters,
                                          Sequential_tag
                                        >::type;

  using GT = typename GetGeomTraits<TriangleMesh, NamedParameters>::type;

  using VPM = typename GetVertexPointMap<TriangleMesh, NamedParameters>::const_type;
  using face_descriptor = typename boost::graph_traits<TriangleMesh>::face_descriptor;
  using Segment_3 = typename GT::Segment_3;

  using AABB_tree_helper = internal::AABB_tree_graph_helper<TriangleMesh, GT, VPM>;
  using Face_tree = typename AABB_tree_helper::Tree;

  using Segment_primitive = CGAL::AABB_segment_primitive_3<GT, typename std::vector<Segment_3>::const_iterator>;
  using Segment_traits = CGAL::AABB_traits_3<GT, Segment_primitive>;
  using Segment_tree = CGAL::AABB_tree<Segment_traits>;

  VPM vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                             get_const_property_map(boost::vertex_point, tm));

  Bbox_3 b1 = Polygon_mesh_processing::bbox(tm, np),
         b2 = bbox_3(polyline.begin(), polyline.end());

  if(!CGAL::do_overlap(b1, b2))
    return out;
  Bbox_3 bb((std::max)(b1.xmin(),b2.xmin()), (std::max)(b1.ymin(),b2.ymin()), (std::max)(b1.zmin(),b2.zmin()),
            (std::min)(b1.xmax(),b2.xmax()), (std::min)(b1.ymax(),b2.ymax()), (std::min)(b1.zmax(),b2.zmax()));

  std::vector<Segment_3> segments;
  if(polyline.size() < 2)
    return out;

  auto pit = std::begin(polyline);
  auto prev = pit++;
  for(; pit != std::end(polyline); ++pit, ++prev)
    segments.emplace_back(*prev, *pit);

  std::vector<face_descriptor> face_to_test;
  for(face_descriptor f : face_range)
    if(do_overlap(Polygon_mesh_processing::face_bbox(f, tm), bb))
      face_to_test.emplace_back(f);

  AABB_tree_helper helper;
  Face_tree face_tree;

  helper.template build<Concurrency_tag>(face_tree, tm, vpm);
  Segment_tree segment_tree(segments.begin(), segments.end());

  auto out_converter = boost::make_function_output_iterator(
    [&](const auto& pair)    {
      *out++ = std::make_pair(pair.first,
                    static_cast<std::size_t>(pair.second - segments.begin()));
  });
  CGAL::AABB_trees::all_pairs_of_intersecting_primitives(face_tree, segment_tree, out_converter, parameters::concurrency_tag(Concurrency_tag()));
  return out;
}

// Note this is not officially documented
/* \ingroup PMP_intersection_grp
 *
 * reports all the pairs of segments and faces intersecting between
 * a triangulated surface mesh and a range of polylines.
 *
 * This function depends on the package \ref PkgBoxIntersectionD.
 *
 * \attention If a polyline vertex intersects a face, the intersection will
 * be reported twice (even more if it is on a vertex, edge, or point).
 *
 * @pre \link CGAL::is_triangle_mesh `CGAL::is_triangle_mesh(tm)` \endlink
 *
 * \tparam TriangleMesh a model of `FaceListGraph`
 * \tparam FaceRange range of `boost::graph_traits<TriangleMesh>::%face_descriptor`,
 *  model of `RandomAccessRange`.
 * \tparam PolylineRange a `RandomAccessRange` of `RandomAccessRange` of points. The point type of the range must be
 * the same as the value type of the vertex point map.
 * \tparam OutputIterator a model of `OutputIterator` holding objects of type
 *   `std::pair<std::size_t, std::pair<std::size_t, std::size_t> >`.
 * Each pair holds the index of the face and a pair containing the index of the polyline in the range and the index of
 * the first point of the segment in the polyline.
 * \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
 *
 * \param face_range the range of `tm` faces to check for intersections.
 * \param polyline_range the range of polylines to check for intersections.
 * \param tm the triangulated surface mesh to check for intersections.
 * \param out output iterator to be filled with all pairs of face-segment that intersect
 * \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
 *
 * \cgalNamedParamsBegin
 *   \cgalParamNBegin{vertex_point_map}
 *     \cgalParamDescription{a property map associating points to the vertices of `tm`}
 *     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<TriangleMesh>::%vertex_descriptor`
 *                    as key type and `%Point_3` as value type}
 *     \cgalParamDefault{`boost::get(CGAL::vertex_point, tm)`}
 *     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
 *                     should be available for the vertices of `tm`.}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{geom_traits}
 *     \cgalParamDescription{an instance of a geometric traits class}
 *     \cgalParamType{a class model of `PMPSelfIntersectionTraits`}
 *     \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
 *     \cgalParamExtra{The geometric traits class must be compatible with the vertex point type.}
 *   \cgalParamNEnd
 * \cgalNamedParamsEnd
 * \return `out`
 */
template <class TriangleMesh,
          class FaceRange,
          class PolylineRange,
          class OutputIterator,
          class NamedParameters>
OutputIterator
compute_face_polylines_intersection(const FaceRange& face_range,
                                    const PolylineRange& polyline_range,
                                    const TriangleMesh& tm,
                                    OutputIterator out,
                                    const NamedParameters& np)
{
  using parameters::get_parameter;
  using parameters::choose_parameter;

  CGAL_precondition(CGAL::is_triangle_mesh(tm));

  using Concurrency_tag = typename internal_np::Lookup_named_param_def <
                                          internal_np::concurrency_tag_t,
                                          NamedParameters,
                                          Sequential_tag
                                        >::type;

  using GT = typename GetGeomTraits<TriangleMesh, NamedParameters>::type;

  using VPM = typename GetVertexPointMap<TriangleMesh, NamedParameters>::const_type;
  using face_descriptor = typename boost::graph_traits<TriangleMesh>::face_descriptor;
  using Segment_3 = typename GT::Segment_3;

  using AABB_tree_helper = internal::AABB_tree_graph_helper<TriangleMesh, GT, VPM>;
  using Face_tree = typename AABB_tree_helper::Tree;

  using Segment_primitive = CGAL::AABB_segment_primitive_3<GT, typename std::vector<Segment_3>::const_iterator>;
  using Segment_traits = CGAL::AABB_traits_3<GT, Segment_primitive>;
  using Segment_tree = CGAL::AABB_tree<Segment_traits>;

  VPM vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                             get_const_property_map(boost::vertex_point, tm));

  Bbox_3 b1 = Polygon_mesh_processing::bbox(tm, np);
  auto polyline_it = polyline_range.begin();
  Bbox_3 b2 = bbox_3(polyline_it->begin(), polyline_it->end());
  for(++polyline_it; polyline_it != polyline_range.end(); ++polyline_it)
    b2 += bbox_3(polyline_it->begin(), polyline_it->end());

  if(!CGAL::do_overlap(b1, b2))
    return out;
  Bbox_3 bb((std::max)(b1.xmin(),b2.xmin()), (std::max)(b1.ymin(),b2.ymin()), (std::max)(b1.zmin(),b2.zmin()),
            (std::min)(b1.xmax(),b2.xmax()), (std::min)(b1.ymax(),b2.ymax()), (std::min)(b1.zmax(),b2.zmax()));

  // We store all polylines in one range and we store the offset of the beginning of each polyline
  std::vector<Segment_3> segments;
  std::vector<std::size_t> polyline_begin;
  for(std::size_t pi = 0; pi < polyline_range.size(); ++pi)
  {
    const auto& polyline = polyline_range[pi];

    if(polyline.size() < 2)
      continue;

    polyline_begin.push_back(segments.size());

    auto pit = std::begin(polyline);
    auto prev = pit++;
    for(; pit != std::end(polyline); ++pit, ++prev)
      segments.emplace_back(*prev, *pit);
  }

  std::vector<face_descriptor> face_to_test;
  for(face_descriptor f : face_range)
    if(do_overlap(Polygon_mesh_processing::face_bbox(f, tm), bb))
      face_to_test.emplace_back(f);

  AABB_tree_helper helper;
  Face_tree face_tree;

  helper.template build<Concurrency_tag>(face_tree, tm, vpm);
  Segment_tree segment_tree(segments.begin(), segments.end());

  auto out_converter =
  boost::make_function_output_iterator(
    [&](const auto& pair)
    {
      std::size_t index = pair.second - segments.begin();

      // Find the polyline containing this segment.
      auto it = std::upper_bound(polyline_begin.begin(),
                                 polyline_begin.end(),
                                 index);

      std::size_t polyline_index = static_cast<std::size_t>(it - polyline_begin.begin() - 1);
      std::size_t segment_index = index - polyline_begin[polyline_index];

      *out++ = std::make_pair(pair.first,
                    std::make_pair(polyline_index, segment_index));
  });
  CGAL::AABB_trees::all_pairs_of_intersecting_primitives(face_tree, segment_tree, out_converter, parameters::concurrency_tag(Concurrency_tag()));
  return out;
}

// Note this is not officially documented
/* \ingroup PMP_intersection_grp
 *
 * detects and records intersections between two ranges of polylines.
 *
 * This function depends on the package \ref PkgBoxIntersectionD.
 *
 * \attention If a polyline vertex intersects another polyline, the intersection will
 * be reported twice (even more if it is on a vertex).
 *
 * \tparam PolylineRange a `RandomAccessRange` of `RandomAccessRange` of points.
 * \tparam OutputIterator a model of `OutputIterator` holding objects of type
 *   `std::pair<std::pair<std::size_t, std::size_t>, std::pair<std::size_t, std::size_t> >`.
 * Each pair holds the index of the face and a pair containing the index of the polyline in the range and the index of
 * the first point of the segment in the polyline.
 * \tparam Kernel a model of `Kernel`
 *
 * \param polylines1 the first range of polylines to check for intersections.
 * \param polylines2 the second range of polylines to check for intersections.
 * \param out output iterator to be filled with all pairs of segments that intersect
 * \param K an instance of `Kernel`
 *
 * \return `out`
 */
template < class PolylineRange,
           class OutputIterator,
           class Kernel>
OutputIterator
compute_polylines_polylines_intersection(const PolylineRange& polylines1,
                                         const PolylineRange& polylines2,
                                         OutputIterator out,
                                         const Kernel& /*K*/)
{
  using parameters::get_parameter;

  using GT = Kernel;
  using Segment_3 = typename GT::Segment_3;

  using Segment_primitive = CGAL::AABB_segment_primitive_3<GT, typename std::vector<Segment_3>::const_iterator>;
  using Segment_traits = CGAL::AABB_traits_3<GT, Segment_primitive>;
  using Segment_tree = CGAL::AABB_tree<Segment_traits>;

  std::vector<Segment_3> segments1, segments2;
  std::vector<std::size_t> polyline_begin1, polyline_begin2;

  for(std::size_t pi = 0; pi < polylines1.size(); ++pi)
  {
    const auto& polyline = polylines1[pi];
    polyline_begin1.push_back(segments1.size());

    auto pit = std::begin(polyline);
    auto prev = pit++;
    for(; pit != std::end(polyline); ++pit, ++prev)
      segments1.emplace_back(*prev, *pit);
  }

  for(std::size_t pi = 0; pi < polylines2.size(); ++pi)
  {
    const auto& polyline = polylines2[pi];
    polyline_begin2.push_back(segments2.size());

    auto pit = std::begin(polyline);
    auto prev = pit++;
    for(; pit != std::end(polyline); ++pit, ++prev)
      segments2.emplace_back(*prev, *pit);
  }

  if(segments1.empty() || segments2.empty())
    return out;

  Segment_tree tree1(segments1.begin(), segments1.end());
  Segment_tree tree2(segments2.begin(), segments2.end());

  auto convert_out =
    boost::make_function_output_iterator(
      [&](const auto& pair)
      {
        const std::size_t offset1 = static_cast<std::size_t>(pair.first - segments1.cbegin());
        const std::size_t offset2 = static_cast<std::size_t>(pair.second - segments2.cbegin());

        const auto it1 = std::upper_bound(polyline_begin1.begin(),
                           polyline_begin1.end(),
                           offset1);

        const auto it2 = std::upper_bound(polyline_begin2.begin(),
                           polyline_begin2.end(),
                           offset2);

        const std::size_t pi1 = static_cast<std::size_t>(it1 - polyline_begin1.begin() - 1);
        const std::size_t pi2 = static_cast<std::size_t>(it2 - polyline_begin2.begin() - 1);

        *out++ = std::make_pair(
          std::make_pair(pi1, offset1 - polyline_begin1[pi1]),
          std::make_pair(pi2, offset2 - polyline_begin2[pi2]));
      });

  CGAL::AABB_trees::all_pairs_of_intersecting_primitives(
    tree1,
    tree2,
    convert_out);

  return out;
}

// Note this is not officially documented
/* \ingroup PMP_intersection_grp
 *
 * detects and records intersections between two polylines.
 *
 * This function depends on the package \ref PkgBoxIntersectionD.
 *
 * \attention If a polyline vertex intersects another polyline, the intersection will
 * be reported twice (even more if it is on a vertex).
 *
 * \tparam Polyline a `RandomAccessRange` of points.
 * \tparam OutputIterator a model of `OutputIterator` holding objects of type
 *   `std::pair<std::size_t, std::size_t>`. This OutputIterator will hold the position of the
 *  elements in their respective range. This position is the index of the segment that holds the
 * intersection, so it is the index of the first point of the segment following the range order.
 * \tparam Kernel a model of `Kernel`
 *
 * \param polyline1 the first polyline to check for intersections.
 * \param polyline2 the second polyline to check for intersections.
 * \param out output iterator to be filled with all pairs of segments that intersect
 * \param K an instance of `Kernel`
 *
 * \return `out`
 */
template < class Polyline,
           class OutputIterator,
           class Kernel>
OutputIterator
compute_polyline_polyline_intersection(const Polyline& polyline1,
                                       const Polyline& polyline2,
                                       OutputIterator out,
                                       const Kernel& K)
{
  std::array<Polyline, 1> polylines1{{polyline1}};
  std::array<Polyline, 1> polylines2{{polyline2}};

  auto polyline_out =
    boost::make_function_output_iterator(
      [&](const auto& pair)
      {
        *out++ = std::make_pair(pair.first.second,
                                pair.second.second);
      });

  compute_polylines_polylines_intersection(polylines1, polylines2, polyline_out, K);
  return out;
}

// Note this is not officially documented
/* \ingroup PMP_intersection_grp
 *
 * reports all the pairs of faces intersecting between two triangulated surface meshes.
 *
 * This function depends on the package \ref PkgBoxIntersectionD.
 *
 * @pre \link CGAL::is_triangle_mesh `CGAL::is_triangle_mesh(tm1)` \endlink
 * @pre \link CGAL::is_triangle_mesh `CGAL::is_triangle_mesh(tm2)` \endlink
 *
 * \tparam TriangleMesh a model of `FaceListGraph`
 * \tparam OutputIterator a model of `OutputIterator` holding objects of type
 *   `std::pair<boost::graph_traits<TriangleMesh>::%face_descriptor, boost::graph_traits<TriangleMesh>::%face_descriptor>`
 * \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
 *
 * \param tm1 the first triangulated surface mesh to check for intersections
 * \param tm2 the second triangulated surface mesh to check for intersections
 * \param out output iterator to be filled with all pairs of faces that intersect
 * \param np1 an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
 * \param np2 an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
 *
 * \cgalNamedParamsBegin
 *   \cgalParamNBegin{vertex_point_map}
 *     \cgalParamDescription{a property map associating points to the vertices of `tm1` (`tm2`)}
 *     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<TriangleMesh>::%vertex_descriptor`
 *                    as key type and `%Point_3` as value type}
 *     \cgalParamDefault{`boost::get(CGAL::vertex_point, tm1 (tm2))`}
 *     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
 *                     should be available for the vertices of `tm1` (`tm2`)}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{geom_traits}
 *     \cgalParamDescription{an instance of a geometric traits class}
 *     \cgalParamType{a class model of `PMPSelfIntersectionTraits`}
 *     \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
 *     \cgalParamExtra{The geometric traits class must be compatible with the vertex point type.}
 *     \cgalParamExtra{np1 only}
 *   \cgalParamNEnd
 * \cgalNamedParamsEnd
 *
 * \return `out`
 */
template <class TriangleMesh,
          class OutputIterator,
          class NamedParameters1,
          class NamedParameters2>
OutputIterator
compute_face_face_intersection(const TriangleMesh& tm1,
                               const TriangleMesh& tm2,
                               OutputIterator out,
                               const NamedParameters1& np1,
                               const NamedParameters2& np2)
{
  return compute_face_face_intersection(faces(tm1), faces(tm2),
                                        tm1, tm2, out, np1, np2);
}

// Note this is not officially documented
/* \ingroup PMP_intersection_grp
 *
 * detects and records intersections between a triangulated surface mesh and a polyline.
 *
 * This function depends on the package \ref PkgBoxIntersectionD.
 *
 * \attention If a polyline vertex intersects a face or another polyline, the intersection will
 * be reported twice (even more if it is on a vertex, edge, or point).
 *
 * @pre \link CGAL::is_triangle_mesh `CGAL::is_triangle_mesh(tm)` \endlink
 *
 * \tparam TriangleMesh a model of `FaceListGraph`
 * \tparam Polyline a `RandomAccessRange` of points. The point type of the range must be the
 * same as the value type of the vertex point map.
 * \cgalDescribePolylineType
 * \tparam OutputIterator a model of `OutputIterator` holding objects of type
 *   `std::pair<std::size_t, std::size_t>`. This OutputIterator will hold the position of the
 * elements in their respective range. In the case of the polyline, this position is the index
 * of the segment that holds the intersection, so it is the index of the first point of the
 * segment following the range order.
 * \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
 *
 * \param tm the triangulated surface mesh to check for intersections.
 * \param polyline the polyline to check for intersections.
 * \param out output iterator to be filled with all pairs of face-segment that intersect
 * \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
 *
 * \cgalNamedParamsBegin
 *   \cgalParamNBegin{vertex_point_map}
 *     \cgalParamDescription{a property map associating points to the vertices of `tm`}
 *     \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<TriangleMesh>::%vertex_descriptor` as key type and `%Point_3` as value type}
 *     \cgalParamDefault{`boost::get(CGAL::vertex_point, tm)`}
 *     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
 *                     should be available for the vertices of `tm`.}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{geom_traits}
 *     \cgalParamDescription{an instance of a geometric traits class}
 *     \cgalParamType{a class model of `PMPSelfIntersectionTraits`}
 *     \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
 *     \cgalParamExtra{The geometric traits class must be compatible with the vertex point type.}
 *   \cgalParamNEnd
 * \cgalNamedParamsEnd
 *
 * \return `out`
 */
template <class TriangleMesh,
          class Polyline,
          class OutputIterator,
          class NamedParameters>
OutputIterator
compute_face_polyline_intersection(const TriangleMesh& tm,
                                   const Polyline& polyline,
                                   OutputIterator out,
                                   const NamedParameters& np)
{
  return compute_face_polyline_intersection(faces(tm), polyline, tm, out, np);
}

// functions to check for overlap of meshes
template <class GT, class TriangleMesh, class VPM>
void get_one_point_per_cc(TriangleMesh& tm,
                          const VPM& vpm,
                          std::vector<typename GT::Point_3>& points_of_interest)
{
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor face_descriptor;
  std::unordered_map<face_descriptor, int> fid_map;
  int id = 0;
  for(face_descriptor fd : faces(tm))
  {
    fid_map.insert(std::make_pair(fd,id++));
  }
  boost::associative_property_map< std::unordered_map<face_descriptor, int> >
      fid_pmap(fid_map);
  std::unordered_map<face_descriptor, int> fcc_map;

  int nb_cc = Polygon_mesh_processing::connected_components(tm,
                                                            boost::make_assoc_property_map(fcc_map),
                                                            parameters::face_index_map(fid_pmap));
  std::vector<bool> is_cc_treated(nb_cc, false);
  points_of_interest.resize(nb_cc);
  int cc_treated = 0;
  for(face_descriptor fd : faces(tm))
  {
    int cc=fcc_map[fd];
    if(!is_cc_treated[cc])
    {
      points_of_interest[cc]=get(vpm, target(halfedge(fd, tm),tm));
      is_cc_treated[cc] = true;
      if(++cc_treated == nb_cc)
        break;
    }
  }
}

//this assumes the meshes does not intersect
template <class TriangleMesh, class VPM, class GT, class AABB_tree>
bool is_mesh2_in_mesh1_impl(const AABB_tree& tree1,
                            const std::vector<typename GT::Point_3>& points_of_interest2,
                            const GT& gt)
{
  //for each CC, take a point on it and test bounded side
  Side_of_triangle_mesh<TriangleMesh, GT, VPM> sotm(tree1, gt);
  for(const typename GT::Point_3& p : points_of_interest2)
  {
    if(sotm(p) == CGAL::ON_BOUNDED_SIDE) // sufficient as we know meshes do not intersect
    {
      return true;
    }
  }
  return false;
}

template <class TriangleMesh, class VPM1, class VPM2, class GT>
bool is_mesh2_in_mesh1(const TriangleMesh& tm1,
                       const TriangleMesh& tm2,
                       const VPM1& vpm1,
                       const VPM2& vpm2,
                       const GT& gt)
{
  typedef CGAL::AABB_face_graph_triangle_primitive<TriangleMesh, VPM1> Primitive;
  typedef CGAL::AABB_traits_3<GT, Primitive> Traits;
  typedef CGAL::AABB_tree<Traits> AABBTree;

  AABBTree tree1(faces(tm1).begin(), faces(tm1).end(), tm1, vpm1);
  std::vector<typename GT::Point_3> points_of_interest2;
  get_one_point_per_cc<GT>(tm2, vpm2, points_of_interest2);

  return is_mesh2_in_mesh1_impl<TriangleMesh, VPM1>(tree1, points_of_interest2, gt);
}


}// namespace internal

/**
 * \ingroup PMP_intersection_grp
 *
 * returns `true` if there exists a segment of a polyline of `polylines1`
 * and a segment of a polyline of `polylines2` which intersect, and `false` otherwise.
 *
 * This function depends on the package \ref PkgBoxIntersectionD.
 *
 * \tparam PolylineRange a `RandomAccessRange` of `RandomAccessRange` of points.
 *         The point type must be from a 3D point from a \cgal Kernel.
 * \cgalDescribePolylineType
 *
 * @param polylines1 the first range of polylines to check for intersections.
 * @param polylines2 the second range of polylines to check for intersections.
 *
 */
template <class PolylineRange>
bool do_intersect(const PolylineRange& polylines1,
                  const PolylineRange& polylines2
#ifndef DOXYGEN_RUNNING
                  , const std::enable_if_t<
                      boost::has_range_iterator<
                      typename boost::mpl::eval_if<
                        boost::has_range_iterator<PolylineRange>,
                        boost::range_value<PolylineRange>,
                        std::false_type >::type
                    >::value
                   >* = 0//end enable_if
#endif
    )
{
  typedef typename boost::range_value<PolylineRange>::type Polyline;
  typedef typename boost::range_value<Polyline>::type Point;
  typedef typename CGAL::Kernel_traits<Point>::Kernel K;
  try
  {
    typedef boost::function_output_iterator<internal::Throw_at_first_output> OutputIterator;
    internal::compute_polylines_polylines_intersection(polylines1, polylines2, OutputIterator(), K());
  }
  catch( internal::Throw_at_first_output::Throw_at_first_output_exception& )
  { return true; }

  return false;
}

/**
 * \ingroup PMP_intersection_grp
 *
 * returns `true` if there exists a segment of `polyline1` and a segment of `polyline2` which intersect,
 * and `false` otherwise.
 *
 * This function depends on the package \ref PkgBoxIntersectionD.
 *
 * \tparam Polyline a `RandomAccessRange` of points.
 *         The point type must be from a 3D point type from \cgal Kernel.
 * \cgalDescribePolylineType
 *
 * @param polyline1 the first polyline to check for intersections.
 * @param polyline2 the second polyline to check for intersections.
 *
 */
template <class Polyline>
bool do_intersect(const Polyline& polyline1,
                  const Polyline& polyline2
#ifndef DOXYGEN_RUNNING
                , const std::enable_if_t<
                    boost::has_range_const_iterator<Polyline>::value
                  >* = 0,
                  const std::enable_if_t<
                    !boost::has_range_iterator<
                      typename boost::mpl::eval_if<
                        boost::has_range_iterator<Polyline>,
                        boost::range_value<Polyline>,
                        std::false_type
                      >::type
                    >::value
                  >* = 0//end enable_if
#endif
                 )
{
  typedef typename boost::range_value<Polyline>::type Point;
  typedef typename CGAL::Kernel_traits<Point>::Kernel K;
  try
  {
    typedef boost::function_output_iterator<internal::Throw_at_first_output> OutputIterator;
    internal::compute_polyline_polyline_intersection(polyline1, polyline2, OutputIterator(), K());
  }
  catch( internal::Throw_at_first_output::Throw_at_first_output_exception& )
  { return true; }

  return false;
}

/**
 * \ingroup PMP_intersection_grp
 *
 * \brief returns `true` if there exists a face of `tm1` and a face of `tm2` which intersect, and `false` otherwise.
 *
 * If `do_overlap_test_of_bounded_sides` is set to `true`, the overlap of bounded sides are tested as well.
 * In that case, the meshes must be closed.
 *
 * This function depends on the package \ref PkgBoxIntersectionD.
 *
 * @pre \link CGAL::is_triangle_mesh `CGAL::is_triangle_mesh(tm1)` \endlink
 * @pre \link CGAL::is_triangle_mesh `CGAL::is_triangle_mesh(tm2)` \endlink
 * @pre `!do_overlap_test_of_bounded_sides` || \link CGAL::is_closed `CGAL::is_closed(tm1)` \endlink
 * @pre `!do_overlap_test_of_bounded_sides` || \link CGAL::is_closed `CGAL::is_closed(tm2)` \endlink
 *
 * @tparam TriangleMesh a model of `FaceListGraph`
 * @tparam NamedParameters1 a sequence of \ref bgl_namedparameters "Named Parameters" for `tm1`
 * @tparam NamedParameters2 a sequence of \ref bgl_namedparameters "Named Parameters" for `tm2`
 *
 * @param tm1 the first triangulated surface mesh to check for intersections
 * @param tm2 the second triangulated surface mesh to check for intersections
 * @param np1 an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
 * @param np2 an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
 *
 * \cgalNamedParamsBegin
 *   \cgalParamNBegin{vertex_point_map}
 *     \cgalParamDescription{a property map associating points to the vertices of `tm1` (`tm2`)}
 *     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<TriangleMesh>::%vertex_descriptor`
 *                    as key type and `%Point_3` as value type}
 *     \cgalParamDefault{`boost::get(CGAL::vertex_point, tm1 (tm2))`}
 *     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
 *                     should be available for the vertices of `tm1` (`tm2`)}
 *     \cgalParamExtra{Both vertex point maps must have the same value type}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{geom_traits}
 *     \cgalParamDescription{an instance of a geometric traits class}
 *     \cgalParamType{a class model of `PMPSelfIntersectionTraits`}
 *     \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
 *     \cgalParamExtra{The geometric traits class must be compatible with the vertex point type.}
 *     \cgalParamExtra{np1 only}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{do_overlap_test_of_bounded_sides}
 *     \cgalParamDescription{If `true`, also tests the overlap of the bounded sides of `tm1` and `tm2`.
 *                           If `false`, only the intersection of surface triangles is tested.}
 *     \cgalParamType{Boolean}
 *     \cgalParamDefault{`false`}
 *   \cgalParamNEnd
 * \cgalNamedParamsEnd
 *
 * \see `intersecting_meshes()`
 */
template <class TriangleMesh,
          class CGAL_NP_TEMPLATE_PARAMETERS_1,
          class CGAL_NP_TEMPLATE_PARAMETERS_2>
bool do_intersect(const TriangleMesh& tm1,
                  const TriangleMesh& tm2,
                  const CGAL_NP_CLASS_1& np1 = parameters::default_values(),
                  const CGAL_NP_CLASS_2& np2 = parameters::default_values()
#ifndef DOXYGEN_RUNNING
                  , const std::enable_if_t<
                            !boost::has_range_const_iterator<TriangleMesh>::value
                  >* = 0
#endif
                  )
{
  using parameters::choose_parameter;
  using parameters::get_parameter;

  bool test_overlap =  choose_parameter(get_parameter(np1, internal_np::overlap_test),false) ||
                       choose_parameter(get_parameter(np2, internal_np::overlap_test),false);

  CGAL_precondition(CGAL::is_triangle_mesh(tm1));
  CGAL_precondition(CGAL::is_triangle_mesh(tm2));
  CGAL_precondition(!test_overlap || CGAL::is_closed(tm1));
  CGAL_precondition(!test_overlap || CGAL::is_closed(tm2));

  try
  {
    typedef boost::function_output_iterator<internal::Throw_at_first_output> OutputIterator;
    internal::compute_face_face_intersection(tm1,tm2, OutputIterator(), np1, np2);
  }
  catch( internal::Throw_at_first_output::Throw_at_first_output_exception& )
  { return true; }

  if (test_overlap)
  {
    typedef typename GetVertexPointMap<TriangleMesh, CGAL_NP_CLASS_1>::const_type VertexPointMap1;
    typedef typename GetVertexPointMap<TriangleMesh, CGAL_NP_CLASS_2>::const_type VertexPointMap2;
    VertexPointMap1 vpm1 = choose_parameter(get_parameter(np1, internal_np::vertex_point),
                                        get_const_property_map(boost::vertex_point, tm1));
    VertexPointMap2 vpm2 = choose_parameter(get_parameter(np2, internal_np::vertex_point),
                                        get_const_property_map(boost::vertex_point, tm2));
    typedef typename GetGeomTraits<TriangleMesh, CGAL_NP_CLASS_1>::type GeomTraits;
    GeomTraits gt = choose_parameter<GeomTraits>(get_parameter(np1, internal_np::geom_traits));

    return internal::is_mesh2_in_mesh1(tm1, tm2, vpm1, vpm2, gt) ||
           internal::is_mesh2_in_mesh1(tm2, tm1, vpm2, vpm1, gt);
  }
  return false;
}

/**
 * \ingroup PMP_intersection_grp
 *
 * returns `true` if there exists a face of `tm` and a segment of a polyline of `polylines` which intersect,
 * and `false` otherwise.
 *
 * This function depends on the package \ref PkgBoxIntersectionD.
 *
 * @pre \link CGAL::is_triangle_mesh `CGAL::is_triangle_mesh(tm)` \endlink
 *
 * \tparam TriangleMesh a model of `FaceListGraph`
 * \tparam PolylineRange a `RandomAccessRange` of `RandomAccessRange` of points. The point type of the range must be the
 *  same as the value type of the vertex point map.
 * \cgalDescribePolylineType
 * @tparam NamedParameters a sequence of \ref bgl_namedparameters
 *
 * @param tm the triangulated surface mesh to check for intersections
 * @param polylines the range of polylines to check for intersections.
 * @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
 *
 * \cgalNamedParamsBegin
 *   \cgalParamNBegin{vertex_point_map}
 *     \cgalParamDescription{a property map associating points to the vertices of `tm`}
 *     \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<TriangleMesh>::%vertex_descriptor` as key type and `%Point_3` as value type}
 *     \cgalParamDefault{`boost::get(CGAL::vertex_point, tm)`}
 *     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
 *                     should be available for the vertices of `tm`.}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{geom_traits}
 *     \cgalParamDescription{an instance of a geometric traits class}
 *     \cgalParamType{a class model of `PMPSelfIntersectionTraits`}
 *     \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
 *     \cgalParamExtra{The geometric traits class must be compatible with the vertex point type.}
 *   \cgalParamNEnd
 * \cgalNamedParamsEnd
 *
 */
template <class TriangleMesh,
          class PolylineRange,
          class NamedParameters = parameters::Default_named_parameters>
bool do_intersect(const TriangleMesh& tm,
                  const PolylineRange& polylines,
                  const NamedParameters& np = parameters::default_values()
#ifndef DOXYGEN_RUNNING
                , const std::enable_if_t<
                    boost::has_range_iterator<
                      typename boost::mpl::eval_if<
                        boost::has_range_iterator<PolylineRange>,
                        boost::range_value<PolylineRange>,
                        std::false_type
                      >::type
                    >::value
                  >* = 0//end enable_if
#endif
                 )
{
  CGAL_precondition(CGAL::is_triangle_mesh(tm));
  try
  {
    typedef boost::function_output_iterator<internal::Throw_at_first_output> OutputIterator;
    internal::compute_face_polylines_intersection(faces(tm), polylines, tm, OutputIterator(), np);
  }
  catch( internal::Throw_at_first_output::Throw_at_first_output_exception& )
  { return true; }

  return false;
}

/**
 * \ingroup PMP_intersection_grp
 *
 * returns `true` if there exists a face of `tm` and a segment of `polyline` which intersect, and `false` otherwise.
 *
 * This function depends on the package \ref PkgBoxIntersectionD.
 *
 * @pre \link CGAL::is_triangle_mesh `CGAL::is_triangle_mesh(tm)` \endlink
 *
 * \tparam TriangleMesh a model of `FaceListGraph`
 * \tparam Polyline a `RandomAccessRange` of points. The point type of the range must be the
 *  same as the value type of the vertex point map.
 * \cgalDescribePolylineType
 * @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
 *
 * @param tm the triangulated surface mesh to check for intersections
 * @param polyline the polyline to check for intersections.
 * @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
 *
 * \cgalNamedParamsBegin
 *   \cgalParamNBegin{vertex_point_map}
 *     \cgalParamDescription{a property map associating points to the vertices of `tm`}
 *     \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<TriangleMesh>::%vertex_descriptor` as key type and `%Point_3` as value type}
 *     \cgalParamDefault{`boost::get(CGAL::vertex_point, tm)`}
 *     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
 *                     should be available for the vertices of `tm`.}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{geom_traits}
 *     \cgalParamDescription{an instance of a geometric traits class}
 *     \cgalParamType{a class model of `PMPSelfIntersectionTraits`}
 *     \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
 *     \cgalParamExtra{The geometric traits class must be compatible with the vertex point type.}
 *   \cgalParamNEnd
 * \cgalNamedParamsEnd
 *
 */
template <class TriangleMesh,
          class Polyline,
          class CGAL_NP_TEMPLATE_PARAMETERS>
bool do_intersect(const TriangleMesh& tm,
                  const Polyline& polyline,
                  const CGAL_NP_CLASS& np = parameters::default_values()
#ifndef DOXYGEN_RUNNING
                , const std::enable_if_t<!(
                      std::is_same_v<TriangleMesh, Polyline> || // Added to please MSVC 2015
                      !boost::has_range_iterator<Polyline>::value || // not a range
                      boost::has_range_iterator<
                        typename boost::mpl::eval_if<
                          boost::has_range_iterator<Polyline>,
                          boost::range_value<Polyline>,
                          std::false_type>::type
                        >::value
                    )
                  >* = 0
#endif
                 )
{
  CGAL_precondition(CGAL::is_triangle_mesh(tm));
  try
  {
    typedef boost::function_output_iterator<internal::Throw_at_first_output> OutputIterator;
    internal::compute_face_polyline_intersection(tm,polyline, OutputIterator(), np);
  }
  catch( internal::Throw_at_first_output::Throw_at_first_output_exception& )
  { return true; }

  return false;
}

namespace internal {

template<class TriangleMeshRange,
         class GT,
         typename OutputIterator,
         class NamedParametersRange>
struct Mesh_callback
{
  typedef typename boost::range_value<TriangleMeshRange>::type TriangleMesh;
  typedef typename boost::range_value<NamedParametersRange>::type NamedParameters;
  typedef typename GetVertexPointMap<TriangleMesh, NamedParameters>::const_type VPM;
  typedef CGAL::AABB_face_graph_triangle_primitive<TriangleMesh, VPM> Primitive;
  typedef CGAL::AABB_traits_3<GT, Primitive> Traits;
  typedef CGAL::AABB_tree<Traits> AABBTree;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor face_descriptor;

  // fill them in the operator and test for inclusion with
  // Side_of_triangle_mesh.
  const TriangleMeshRange& meshes;
  OutputIterator m_iterator;
  const bool report_overlap;
  const NamedParametersRange& nps;
  std::vector<AABBTree*> trees;
  GT gt;

  std::vector<std::vector<CGAL::Point_3<GT> > > points_of_interest;

  Mesh_callback(const TriangleMeshRange& meshes,
                OutputIterator iterator,
                const bool report_overlap,
                const GT& gt,
                const NamedParametersRange& nps)
    : meshes(meshes), m_iterator(iterator),
      report_overlap(report_overlap), nps(nps), gt(gt)
  {
    std::size_t size = std::distance(meshes.begin(), meshes.end());
    trees = std::vector<AABBTree*>(size, nullptr);
    points_of_interest.resize(size);
  }

  ~Mesh_callback()
  {
    for(AABBTree* tree : trees)
    {
      delete tree;
    }
  }

  template<class TriangleMesh,
           class VPM>
  bool is_mesh2_in_mesh1(const TriangleMesh& tm1,
                         const TriangleMesh& tm2,
                         const std::size_t mesh_id_1,
                         const std::size_t mesh_id_2,
                         const VPM& vpm1,
                         const VPM& vpm2)
  {
    //test if tm2 is included in tm1

    //get AABB_tree for tm1
    if(!trees[mesh_id_1])
    {
      trees[mesh_id_1] = new AABBTree(faces(tm1).begin(),
                                      faces(tm1).end(),
                                      tm1, vpm1);
    }
    //get a face-index map for tm2
    if(points_of_interest[mesh_id_2].size() == 0)
      get_one_point_per_cc<GT>(tm2, vpm2, points_of_interest[mesh_id_2]);

    //test if tm2 is included in tm1:
    return is_mesh2_in_mesh1_impl<TriangleMesh, VPM, GT>(
      *trees[mesh_id_1],  points_of_interest[mesh_id_2], gt);
  }

  template<class Mesh_box>
  void operator()(const Mesh_box* b1, const Mesh_box* b2)
  {
    using parameters::choose_parameter;
    using parameters::get_parameter;

    std::size_t mesh_id_1 = std::distance(meshes.begin(), b1->info());
    std::size_t mesh_id_2 = std::distance(meshes.begin(), b2->info());


    VPM vpm1 = choose_parameter(get_parameter(*(nps.begin() + mesh_id_1), internal_np::vertex_point),
                                get_const_property_map(CGAL::vertex_point, *b1->info()));

    VPM vpm2 = choose_parameter(get_parameter(*(nps.begin() + mesh_id_2), internal_np::vertex_point),
                                get_const_property_map(CGAL::vertex_point, *b2->info()));

    //surfacic test
    if(Polygon_mesh_processing::do_intersect(*b1->info(),
                                             *b2->info(),
                                             parameters::vertex_point_map(vpm1)
                                             .geom_traits(gt),
                                             parameters::vertex_point_map(vpm2)
                                             .geom_traits(gt)))
    {
      *m_iterator++ = std::make_pair(mesh_id_1, mesh_id_2);
    }
    //volumic test
    else if(report_overlap)
    {
      if(!CGAL::do_overlap(b1->bbox(), b2->bbox()))
        return;
      if(is_mesh2_in_mesh1(*b1->info(), *b2->info(), mesh_id_1, mesh_id_2, vpm1, vpm2))
        *m_iterator++ = std::make_pair(mesh_id_1, mesh_id_2);
      else if(is_mesh2_in_mesh1(*b2->info(), *b1->info(), mesh_id_2, mesh_id_1, vpm2, vpm1))
        *m_iterator++ = std::make_pair(mesh_id_2, mesh_id_1);
    }
  }
};
}//end internal

namespace internal {
  template <class Iterator>
  struct AABB_indexed_bbox_primitive
  {
    using Id = std::size_t;
    using Datum = CGAL::Bbox_3;

    AABB_indexed_bbox_primitive() = default;

    AABB_indexed_bbox_primitive(Iterator it)
      : m_it(it)
    {}

    Id id() const
    {
      return m_it->index;
    }

    const Datum& datum() const
    {
      return m_it->bbox;
    }

  private:
    Iterator m_it;
  };
}

/*!
 * \ingroup PMP_intersection_grp
 *
 * detects and reports all the pairs of meshes intersecting in a range of triangulated surface meshes.
 * A pair of meshes intersecting is put in the output iterator `out` as a `std::pair<std::size_t, std::size_t>`,
 * each index referring to the index of the triangle mesh in the input range.
 * If `do_overlap_test_of_bounded_sides` is `true`, the overlap of bounded sides are tested as well. In that case, the meshes must be closed.
 *
 * This function depends on the package \ref PkgBoxIntersectionD.
 *
 * \tparam TriangleMeshRange a model of `RandomAccessRange` of triangulated surface meshes model of `FaceListGraph`.
 * \tparam OutputIterator an output iterator in which `std::pair<std::size_t, std::size_t>` can be put.
 * \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters" for the algorithm
 * \tparam NamedParametersRange a range of \ref bgl_namedparameters "Named Parameters" for the meshes.
 *
 * \param range the range of triangulated surface meshes to be checked for intersections.
 * \param out output iterator used to collect pairs of intersecting meshes.
 * \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
 *
 * \cgalNamedParamsBegin
 *   \cgalParamNBegin{geom_traits}
 *     \cgalParamDescription{an instance of a geometric traits class}
 *     \cgalParamType{a class model of `PMPSelfIntersectionTraits`}
 *     \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`,
 *                       where `Point` is the value type of the vertex point map of the meshes}
 *   \cgalParamNEnd
 *     \cgalParamNBegin{concurrency_tag}
 *       \cgalParamDescription{a tag indicating if the task should be done using one or several threads.}
 *       \cgalParamType{Either `CGAL::Sequential_tag`, or `CGAL::Parallel_tag`, or `CGAL::Parallel_if_available_tag`}
 *       \cgalParamDefault{`CGAL::Sequential_tag`}
 *       \cgalParamExtra{`np1` only}
 *    \cgalParamNEnd
 *
 *   \cgalParamNBegin{do_overlap_test_of_bounded_sides}
 *     \cgalParamDescription{If `true`, reports also overlap of bounded sides of meshes.
 *                           If `false`, only the intersection of surface triangles are tested.}
 *     \cgalParamType{Boolean}
 *     \cgalParamDefault{`false`}
 *   \cgalParamNEnd
 * \cgalNamedParamsEnd
 *
 * \param nps an optional range of sequences of \ref bgl_namedparameters "Named Parameters" among the ones listed below
 *
 * \cgalNamedParamsBegin
 *   \cgalParamNBegin{vertex_point_map}
 *     \cgalParamDescription{a property map associating points to the vertices of a mesh `tm`}
 *     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<TriangleMesh>::%vertex_descriptor`
 *                    as key type and `%Point_3` as value type}
 *     \cgalParamDefault{`boost::get(CGAL::vertex_point, tm)`}
 *     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
 *                     should be available for the vertices of `tm`.}
 *     \cgalParamExtra{All vertex point maps must have the same value type}
 *   \cgalParamNEnd
 * \cgalNamedParamsEnd
 *
 * \see `do_intersect()`
 */
template <class TriangleMeshRange,
          class OutputIterator,
          class NamedParameters,
          class NamedParametersRange>
OutputIterator
intersecting_meshes(const TriangleMeshRange& range,
                    OutputIterator out,
                    const NamedParameters& /*np*/,
                    const NamedParametersRange& nps)
{
  using parameters::choose_parameter;
  using parameters::get_parameter;

  using Concurrency_tag = typename internal_np::Lookup_named_param_def <
                                          internal_np::concurrency_tag_t,
                                          NamedParameters,
                                          Sequential_tag
                                        > ::type;

  using TriangleMeshIterator = typename TriangleMeshRange::const_iterator;
  using TriangleMesh = typename std::iterator_traits<TriangleMeshIterator>::value_type;
  using GT = typename GetGeomTraits<TriangleMesh, NamedParameters>::type;

  struct Indexed_bbox
  {
    Bbox_3 bbox;
    std::size_t index;
  };

  std::vector<Indexed_bbox> indexed_bboxes;
  indexed_bboxes.reserve(range.size());
  for(std::size_t i = 0; i < range.size(); ++i)
    indexed_bboxes.push_back({ Polygon_mesh_processing::bbox(range[i], nps[i]), i });

  if(indexed_bboxes.empty())
    return out;

  // // AABB_traits expects a bbox map associating the primitive ID
  // // with its bounding box.
  using Bbox_map = boost::vector_property_map<Bbox_3>;

  // The primitive ID is the index of the mesh.
  struct AABB_indexed_bbox_primitive
  {
    using Id = std::size_t;
    using Datum = Bbox_3;
    using Point = typename GT::Point_3;

    AABB_indexed_bbox_primitive() = default;
    AABB_indexed_bbox_primitive(typename std::vector<Indexed_bbox>::const_iterator it, const Bbox_map&)
      : m_it(it)
    {}

    Id id() const { return m_it->index; }
    const Datum& datum() const { return m_it->bbox; }
    Point reference_point() const {
      const Bbox_3& b = m_it->bbox;
      return Point(b.xmin(), b.ymin(), b.zmin());
    }

  private:
    typename std::vector<Indexed_bbox>::const_iterator m_it;
  };

  using Primitive = AABB_indexed_bbox_primitive;
  using Traits = CGAL::AABB_traits_3<GT, Primitive, Bbox_map>;
  using Tree = CGAL::AABB_tree<Traits>;

  boost::vector_property_map<Bbox_3> bbox_map(indexed_bboxes.size());
  for(const auto& ib : indexed_bboxes)
    put(bbox_map, ib.index, ib.bbox);
  Tree tree(indexed_bboxes.begin(), indexed_bboxes.end(), bbox_map);

  std::vector<std::pair<std::size_t, std::size_t>> candidates;
  CGAL::AABB_trees::all_pairs_of_intersecting_primitives<Concurrency_tag>(tree, std::back_inserter(candidates));

  for(const auto& p : candidates)
  {
    const TriangleMesh &mesh1 = range[p.first];
    const TriangleMesh &mesh2 = range[p.second];

    if(CGAL::Polygon_mesh_processing::do_intersect(mesh1, mesh2, nps[p.first], nps[p.second]))
      *out++ = p;
  }
  return out;
}

template <class TriangleMeshRange, class NamedParameters, class OutputIterator>
OutputIterator intersecting_meshes(const TriangleMeshRange& range,
                                         OutputIterator out,
                                   const NamedParameters& np)
{
  std::vector<parameters::Default_named_parameters> nps(
    std::distance(range.begin(), range.end()), parameters::default_values());
  return intersecting_meshes(range, out, np, nps);
}

template <class TriangleMeshRange, class OutputIterator>
OutputIterator intersecting_meshes(const TriangleMeshRange& range,
                                         OutputIterator out)
{
  return intersecting_meshes(range, out, parameters::default_values());
}

} // namespace Polygon_mesh_processing
} // namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_POLYGON_MESH_PROCESSING_INTERSECTION_H
