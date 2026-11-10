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
 * \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters" for the algorithm
 *
 * @param polylines1 the first range of polylines to check for intersections.
 * @param polylines2 the second range of polylines to check for intersections.
 * @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
 *
 * \cgalNamedParamsBegin
 *     \cgalParamNBegin{concurrency_tag}
 *       \cgalParamDescription{a tag indicating if the task should be done using one or several threads.}
 *       \cgalParamType{Either `CGAL::Sequential_tag`, or `CGAL::Parallel_tag`, or `CGAL::Parallel_if_available_tag`}
 *       \cgalParamDefault{`CGAL::Sequential_tag`}
 *       \cgalParamExtra{`np1` only}
 *    \cgalParamNEnd
 * \cgalNamedParamsEnd
 *
 */
template <class PolylineRange, class NamedParameters = parameters::Default_named_parameters>
bool do_intersect(const PolylineRange& polylines1,
                  const PolylineRange& polylines2,
                  const NamedParameters& np  = parameters::default_values()
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
  CGAL_USE(np);
  using Concurrency_tag = typename internal_np::Lookup_named_param_def <
                                          internal_np::concurrency_tag_t,
                                          NamedParameters,
                                          Sequential_tag
                                        > ::type;

  using Polyline = typename boost::range_value<PolylineRange>::type;
  using Point = typename boost::range_value<Polyline>::type;
  using GT = typename CGAL::Kernel_traits<Point>::Kernel;
  using Segment_3 = typename GT::Segment_3;

  using Segment_primitive = CGAL::AABB_segment_primitive_3<GT, typename std::vector<Segment_3>::const_iterator>;
  using Segment_traits = CGAL::AABB_traits_3<GT, Segment_primitive>;
  using Segment_tree = CGAL::AABB_tree<Segment_traits>;

  Bbox_3 b1;
  for(const auto& polyline : polylines1)
    b1 += bbox_3(polyline.begin(), polyline.end());
  Bbox_3 b2;
  for(const auto& polyline : polylines2)
    b2 += bbox_3(polyline.begin(), polyline.end());

  if(!CGAL::do_overlap(b1, b2))
    return false;

  Bbox_3 bb((std::max)(b1.xmin(),b2.xmin()), (std::max)(b1.ymin(),b2.ymin()), (std::max)(b1.zmin(),b2.zmin()),
            (std::min)(b1.xmax(),b2.xmax()), (std::min)(b1.ymax(),b2.ymax()), (std::min)(b1.zmax(),b2.zmax()));

  // We store all polylines in one filtered range
  std::vector<Segment_3> segments1;
  for (const auto& polyline : polylines1)
  {
    if(polyline.size() < 2) continue;

    auto pit = std::begin(polyline);
    auto prev = pit++;
    for(; pit != std::end(polyline); ++pit, ++prev)
      if (do_overlap(bb, prev->bbox()+pit->bbox()))
        segments1.emplace_back(*prev, *pit);
  }
  std::vector<Segment_3> segments2;
  for (const auto& polyline : polylines2)
  {
    if(polyline.size() < 2) continue;

    auto pit = std::begin(polyline);
    auto prev = pit++;
    for(; pit != std::end(polyline); ++pit, ++prev)
      if (do_overlap(bb, prev->bbox()+pit->bbox()))
        segments2.emplace_back(*prev, *pit);
  }

  Segment_tree tree1(segments1.begin(), segments1.end());
  Segment_tree tree2(segments2.begin(), segments2.end());

  return AABB_trees::do_intersect(tree1, tree2, parameters::concurrency_tag(Concurrency_tag()));
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
 * \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters" for the algorithm
 *
 * @param polyline1 the first polyline to check for intersections.
 * @param polyline2 the second polyline to check for intersections.
 * @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
 *
 * \cgalNamedParamsBegin
 *     \cgalParamNBegin{concurrency_tag}
 *       \cgalParamDescription{a tag indicating if the task should be done using one or several threads.}
 *       \cgalParamType{Either `CGAL::Sequential_tag`, or `CGAL::Parallel_tag`, or `CGAL::Parallel_if_available_tag`}
 *       \cgalParamDefault{`CGAL::Sequential_tag`}
 *       \cgalParamExtra{`np1` only}
 *    \cgalParamNEnd
 * \cgalNamedParamsEnd
 *
 */
template <class Polyline, class NamedParameters = parameters::Default_named_parameters>
bool do_intersect(const Polyline& polyline1,
                  const Polyline& polyline2,
                  const NamedParameters& np = parameters::default_values()
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
  auto polylines1 = make_range(&polyline1,&polyline1+1);
  auto polylines2 = make_range(&polyline2,&polyline2+1);

  return do_intersect(polylines1, polylines2, np);
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
  *  \cgalParamNBegin{concurrency_tag}
 *     \cgalParamDescription{a tag indicating if the task should be done using one or several threads.}
 *     \cgalParamType{Either `CGAL::Sequential_tag`, or `CGAL::Parallel_tag`, or `CGAL::Parallel_if_available_tag`}
 *     \cgalParamDefault{`CGAL::Sequential_tag`}
 *     \cgalParamExtra{`np1` only}
 *   \cgalParamNEnd
 *   \cgalParamNBegin{vertex_point_map}
 *     \cgalParamDescription{a property map associating points to the vertices of `tm1` (`tm2`)}
 *     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<TriangleMesh>::%vertex_descriptor`
 *                    as key type and `%Point_3` as value type}
 *     \cgalParamDefault{`boost::get(CGAL::vertex_point, tm1 (tm2))`}
 *     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
 *                     should be available for the vertices of `tm1` (`tm2`)}
 *     \cgalParamExtra{Both vertex point maps must have the same value type}
 *   \cgalParamNEnd
 *   \cgalParamNBegin{geom_traits}
 *     \cgalParamDescription{an instance of a geometric traits class}
 *     \cgalParamType{a class model of `PMPSelfIntersectionTraits`}
 *     \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
 *     \cgalParamExtra{The geometric traits class must be compatible with the vertex point type.}
 *     \cgalParamExtra{np1 only}
 *   \cgalParamNEnd
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
  CGAL_precondition(CGAL::is_triangle_mesh(tm1));
  CGAL_precondition(CGAL::is_triangle_mesh(tm2));

  using parameters::choose_parameter;
  using parameters::get_parameter;

  bool test_overlap =  choose_parameter(get_parameter(np1, internal_np::overlap_test),false) ||
                       choose_parameter(get_parameter(np2, internal_np::overlap_test),false);
  CGAL_precondition(!test_overlap || CGAL::is_closed(tm1));
  CGAL_precondition(!test_overlap || CGAL::is_closed(tm2));

  using Concurrency_tag = typename internal_np::Lookup_named_param_def <
                                          internal_np::concurrency_tag_t,
                                          CGAL_NP_CLASS_1,
                                          Sequential_tag
                                        > ::type;

  using GT = typename GetGeomTraits<TriangleMesh, CGAL_NP_CLASS_1>::type;

  using VPM1 = typename GetVertexPointMap<TriangleMesh, CGAL_NP_CLASS_1>::const_type;
  using VPM2 = typename GetVertexPointMap<TriangleMesh, CGAL_NP_CLASS_2>::const_type;

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
    return false;

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

// TODO: can be done in parallel
  std::vector<face_descriptor> faces_to_test1;
  std::vector<face_descriptor> faces_to_test2;

  for(face_descriptor f : faces(tm1))
    if(do_overlap(Polygon_mesh_processing::face_bbox(f, tm1, parameters::vertex_point_map(vpm1)), bb))
      faces_to_test1.emplace_back(f);
  for(face_descriptor f : faces(tm2)  )
    if(do_overlap(Polygon_mesh_processing::face_bbox(f, tm2, parameters::vertex_point_map(vpm2)), bb))
      faces_to_test2.emplace_back(f);
//---

  Tree_1 tree1;
  Tree_2 tree2;
  helper_1.template build<Concurrency_tag>(faces_to_test1, tree1, tm1, vpm1);
  helper_2.template build<Concurrency_tag>(faces_to_test2, tree2, tm2, vpm2);

  bool inter = AABB_trees::do_intersect(tree1, tree2, parameters::concurrency_tag(Concurrency_tag()));

  if (inter) return true;

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
 *   \cgalParamNBegin{concurrency_tag}
 *     \cgalParamDescription{a tag indicating if the task should be done using one or several threads.}
 *     \cgalParamType{Either `CGAL::Sequential_tag`, or `CGAL::Parallel_tag`, or `CGAL::Parallel_if_available_tag`}
 *     \cgalParamDefault{`CGAL::Sequential_tag`}
 *     \cgalParamExtra{`np1` only}
 *   \cgalParamNEnd
 *   \cgalParamNBegin{vertex_point_map}
 *     \cgalParamDescription{a property map associating points to the vertices of `tm`}
 *     \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<TriangleMesh>::%vertex_descriptor` as key type and `%Point_3` as value type}
 *     \cgalParamDefault{`boost::get(CGAL::vertex_point, tm)`}
 *     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
 *                     should be available for the vertices of `tm`.}
 *   \cgalParamNEnd
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
  if (faces(tm).size()==0 || polylines.empty())
    return false;

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

  Bbox_3 b1;
  for (face_descriptor f : faces(tm))
   b1 += Polygon_mesh_processing::face_bbox(f, tm, np);
  Bbox_3 b2;
  for(const auto& polyline : polylines)
    b2 += bbox_3(polyline.begin(), polyline.end());

  if(!CGAL::do_overlap(b1, b2))
    return false;

  Bbox_3 bb((std::max)(b1.xmin(),b2.xmin()), (std::max)(b1.ymin(),b2.ymin()), (std::max)(b1.zmin(),b2.zmin()),
            (std::min)(b1.xmax(),b2.xmax()), (std::min)(b1.ymax(),b2.ymax()), (std::min)(b1.zmax(),b2.zmax()));

  // We store all polylines in one filtered range
  std::vector<Segment_3> segments;
  for (const auto& polyline : polylines)
  {
    if(polyline.size() < 2) continue;

    auto pit = std::begin(polyline);
    auto prev = pit++;
    for(; pit != std::end(polyline); ++pit, ++prev)
      if (do_overlap(bb, prev->bbox()+pit->bbox()))
        segments.emplace_back(*prev, *pit);
  }

  std::vector<face_descriptor> faces_to_test;
  for(face_descriptor f : faces(tm))
    if(do_overlap(Polygon_mesh_processing::face_bbox(f, tm, parameters::vertex_point_map(vpm)), bb))
      faces_to_test.emplace_back(f);

  AABB_tree_helper helper;
  Face_tree face_tree;

  helper.template build<Concurrency_tag>(faces_to_test, face_tree, tm, vpm);
  Segment_tree segment_tree(segments.begin(), segments.end());

  return CGAL::AABB_trees::do_intersect(face_tree, segment_tree, parameters::concurrency_tag(Concurrency_tag()));
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
 *    \cgalParamNBegin{concurrency_tag}
 *      \cgalParamDescription{a tag indicating if the task should be done using one or several threads.}
 *      \cgalParamType{Either `CGAL::Sequential_tag`, or `CGAL::Parallel_tag`, or `CGAL::Parallel_if_available_tag`}
 *      \cgalParamDefault{`CGAL::Sequential_tag`}
 *      \cgalParamExtra{`np1` only}
 *   \cgalParamNEnd
 *   \cgalParamNBegin{vertex_point_map}
 *     \cgalParamDescription{a property map associating points to the vertices of `tm`}
 *     \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<TriangleMesh>::%vertex_descriptor` as key type and `%Point_3` as value type}
 *     \cgalParamDefault{`boost::get(CGAL::vertex_point, tm)`}
 *     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
 *                     should be available for the vertices of `tm`.}
 *   \cgalParamNEnd
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
  auto polylines = make_range(&polyline,&polyline+1);
  return do_intersect(tm, polylines, np);
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
 *   \cgalParamNBegin{concurrency_tag}
 *     \cgalParamDescription{a tag indicating if the task should be done using one or several threads.}
 *     \cgalParamType{Either `CGAL::Sequential_tag`, or `CGAL::Parallel_tag`, or `CGAL::Parallel_if_available_tag`}
 *     \cgalParamDefault{`CGAL::Sequential_tag`}
 *     \cgalParamExtra{`np1` only}
 *   \cgalParamNEnd
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
                    const NamedParameters& np,
                    const NamedParametersRange& nps)
{
  using parameters::choose_parameter;
  using parameters::get_parameter;

  using Concurrency_tag = typename internal_np::Lookup_named_param_def <
                                          internal_np::concurrency_tag_t,
                                          NamedParameters,
                                          Sequential_tag
                                        > ::type;

  bool report_overlap = choose_parameter(get_parameter(np, internal_np::overlap_test),false);

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

  // The primitive ID is the index of the mesh.
  struct AABB_indexed_bbox_primitive
  {
    using Id = std::size_t;
    struct Datum : public Bbox_3
    {
      Datum() = default;
      Datum(const Bbox_3& bbox) : Bbox_3(bbox) {}

      const Datum& bbox() const { return *this; }
    };
    using Point = typename GT::Point_3;

    AABB_indexed_bbox_primitive() = default;
    AABB_indexed_bbox_primitive(typename std::vector<Indexed_bbox>::const_iterator it)
      : m_it(it)
    {}

    Id id() const { return m_it->index; }
    Datum datum() const { return Datum(m_it->bbox); }
    Point reference_point() const {
      const Bbox_3& b = m_it->bbox;
      return Point(b.xmin(), b.ymin(), b.zmin());
    }

  private:
    typename std::vector<Indexed_bbox>::const_iterator m_it;
  };

  using Primitive = AABB_indexed_bbox_primitive;
  using Traits = CGAL::AABB_traits_3<GT, Primitive>;
  using Tree = CGAL::AABB_tree<Traits>;

  Tree tree(indexed_bboxes.begin(), indexed_bboxes.end());

  std::vector<std::pair<std::size_t, std::size_t>> candidates;
  CGAL::AABB_trees::all_pairs_of_intersecting_primitives<Concurrency_tag>(tree, std::back_inserter(candidates));

  for(const auto& p : candidates)
  {
    const TriangleMesh &mesh1 = range[p.first];
    const TriangleMesh &mesh2 = range[p.second];

    if(CGAL::Polygon_mesh_processing::do_intersect(mesh1, mesh2, nps[p.first].do_overlap_test_of_bounded_sides(report_overlap),
                                                                 nps[p.second].do_overlap_test_of_bounded_sides(report_overlap)))
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
