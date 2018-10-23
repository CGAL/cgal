// Copyright (c) 2018 GeometryFactory (France).
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
// Author(s)     : Sebastien Loriot
//                 Mael Rouxel-Labbé

#ifndef CGAL_POLYGON_MESH_PROCESSING_INTERNAL_SNAP_H
#define CGAL_POLYGON_MESH_PROCESSING_INTERNAL_SNAP_H

#include <CGAL/license/Polygon_mesh_processing/repair.h>

#include <CGAL/Polygon_mesh_processing/internal/named_function_params.h>
#include <CGAL/Polygon_mesh_processing/internal/named_params_helper.h>
#include <CGAL/Polygon_mesh_processing/border.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_halfedge_graph_segment_primitive.h>
#include <CGAL/assertions.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/boost/graph/iterator.h>
#include <CGAL/box_intersection_d.h>
#include <CGAL/circulator.h>
#include <CGAL/Dynamic_property_map.h>
#include <CGAL/number_utils.h>
#include <CGAL/unordered.h>

#include <boost/bimap.hpp>
#include <boost/bimap/multiset_of.hpp>
#include <boost/bimap/set_of.hpp>
#include <boost/foreach.hpp>
#include <boost/function.hpp>
#include <boost/functional.hpp>
#include <boost/type_traits/is_same.hpp>

#include <iostream>
#include <iterator>
#include <limits>
#include <map>
#include <utility>
#include <vector>

namespace CGAL {

namespace Polygon_mesh_processing {

namespace internal {

// Assigns at each vertex the length of its shortest incident edge as 'epsilon' value
template <typename HalfedgeRange,
          typename ToleranceMap,
          typename PolygonMesh,
          typename SourceNamedParameters>
void compute_tolerance_at_vertices(const HalfedgeRange& hrange,
                                   ToleranceMap& tol_pmap,
                                   PolygonMesh& mesh,
                                   const SourceNamedParameters& snp)
{
  using boost::get_param;
  using boost::choose_param;

  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor                vertex_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor              halfedge_descriptor;

  typedef typename GetVertexPointMap<PolygonMesh, SourceNamedParameters>::type        SVPM;
  typedef typename GetGeomTraits<PolygonMesh, SourceNamedParameters>::type            GT;
  typedef typename GT::FT                                                             FT;

  GT gt = choose_param(get_param(snp, internal_np::geom_traits), GT());

  SVPM svpm = choose_param(get_param(snp, internal_np::vertex_point),
                           get_property_map(vertex_point, mesh));

  BOOST_FOREACH(halfedge_descriptor hd, hrange)
  {
    const vertex_descriptor vd = target(hd, mesh);
    CGAL::Halfedge_around_target_iterator<PolygonMesh> hit, hend;
    boost::tie(hit, hend) = CGAL::halfedges_around_target(vd, mesh);
    CGAL_assertion(hit != hend);

    FT sq_length = gt.compute_squared_distance_3_object()(get(svpm, source(*hit, mesh)),
                                                          get(svpm, target(*hit, mesh)));
    FT min_sq_dist = sq_length;
    ++hit;

    for(; hit!=hend; ++hit)
    {
      sq_length = gt.compute_squared_distance_3_object()(get(svpm, source(*hit, mesh)),
                                                         get(svpm, target(*hit, mesh)));

      if(sq_length < min_sq_dist)
        min_sq_dist = sq_length;
    }

    put(tol_pmap, vd, CGAL::approximate_sqrt(min_sq_dist / 4.));
  }
}

template <typename PolygonMesh, typename GeomTraits,
          typename VertexCorrespondenceMap,
          typename SVPM, typename TVPM,
          typename ToleranceMap,
          typename Box>
struct Vertex_proximity_report
{
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor       vertex_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::edge_descriptor         edge_descriptor;

  typedef typename GeomTraits::FT                                            FT;
  typedef typename boost::property_traits<SVPM>::value_type                  Point;

  typedef VertexCorrespondenceMap                                            Vertex_correspondence_map;
  typedef typename Vertex_correspondence_map::left_iterator                  VCM_left_iterator;
  typedef typename Vertex_correspondence_map::value_type                     VCM_value_type;

  Vertex_proximity_report(Vertex_correspondence_map& vertex_map,
                          const SVPM& svpm, const TVPM& tvpm,
                          const ToleranceMap& tol_pmap,
                          const GeomTraits& gt)
    :
      m_vertex_map(vertex_map),
      svpm(svpm), tvpm(tvpm),
      tol_pmap(tol_pmap),
      gt(gt)
  { }

  void operator()(const Box& a, const Box& b)
  {
    vertex_descriptor va = a.info();
    vertex_descriptor vb = b.info();

    if(va == vb)
      return;

    const Point& sp = get(svpm, va);
    const Point& tp = get(tvpm, vb);
    FT tol = get(tol_pmap, va);

    const FT sq_dist = gt.compute_squared_distance_3_object()(sp, tp);
    CGAL::Comparison_result res = CGAL::compare(sq_dist, tol * tol);

    if(res == CGAL::LARGER)
      return;

    std::pair<VCM_left_iterator, bool> it = m_vertex_map.left.insert(std::make_pair(va, vb));
    const vertex_descriptor vb2 = it.first->second;

    // If there are multiple candidates for a source vertex, keep the closest target vertex
    if(vb2 != vb)
    {
      typename CGAL::cpp11::unordered_map<vertex_descriptor, FT>::iterator dist_it =
          m_sq_distance_to_snapped_point.find(va);
      CGAL_assertion(dist_it != m_sq_distance_to_snapped_point.end());

      const FT sq_dist_to_prev_best = dist_it->second;
      if(CGAL::compare(sq_dist, sq_dist_to_prev_best) == CGAL::SMALLER)
      {
        VCM_left_iterator hint = it.first;
        ++hint;
        m_vertex_map.left.erase(it.first);
        m_vertex_map.left.insert(hint, std::make_pair(va, vb));
        dist_it->second = sq_dist;
      }
    }
    else
    {
      m_sq_distance_to_snapped_point[va] = sq_dist;
    }
  }

private:
  Vertex_correspondence_map& m_vertex_map;
  CGAL::cpp11::unordered_map<vertex_descriptor/*source*/, FT> m_sq_distance_to_snapped_point;

  const SVPM& svpm;
  const TVPM& tvpm;
  const ToleranceMap& tol_pmap;
  const GeomTraits& gt;
};

// This is the function if you know what you're doing with the ranges
//
// \ingroup PMP_repairing_grp
//
// Attempts to snap the vertices in `source_hrange` onto the vertices in `target_hrange`.
// A vertex of the source range is only snapped to a vertex of the target mesh
// if its distance to the target mesh vertex is smaller than a user-chosen bound.
// If any source vertex can be snapped onto multiple vertices of the target
// range, the closest one is chosen.
// If multiple vertices within the source range are to be snapped to the same target vertex,
// then the snapping is not performed for these vertices.
//
// @tparam PolygonMesh a model of `FaceListGraph` and `MutableFaceGraph`
// @tparam SourceHalfedgeRange a model of `Range` with value type `boost::graph_traits<PolygonMesh>::%halfedge_descriptor`
// @tparam TargetHalfedgeRange a model of `Range` with value type `boost::graph_traits<PolygonMesh>::%halfedge_descriptor`
// @tparam ToleranceMap a model of `ReadablePropertyMap` with key type `boost::graph_traits<PolygonMesh>::%vertex_descriptor`
//                      and value type `GetGeomTraits<PolygonMesh, SourceNamedParameters>::type::FT`
// @tparam SourceNamedParameters a sequence of \ref pmp_namedparameters "Named Parameters"
// @tparam TargetNamedParameters a sequence of \ref pmp_namedparameters "Named Parameters"
//
// @param source_hrange a range of vertices of the source mesh whose positions can be changed.
//                      the vertices must be border vertices of `smesh`.
// @param smesh the source mesh whose border vertices might be moved
// @param target_hrange a range of vertices of the target mesh which are potential new positions
//                      for the vertices in the source range
// @param tmesh the target mesh to which the vertices in `target_hrange` belong
// @param tol_pmap a tolerance map associating to each vertex of the source range a tolerance value:
//               potential projection targets are sought in a sphere centered at the vertex and
//               whose radius is the tolerance value.
// @param snp optional \ref pmp_namedparameters "Named Parameters" related to the source mesh,
//            amongst those described below:
//
// \cgalNamedParamsBegin
//    \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of the source mesh.
//                                      The type of this map is model of `ReadWritePropertyMap`.
//                                      If this parameter is omitted, an internal property map for
//                                      `CGAL::vertex_point_t` must be available in `PolygonMesh`
//    \cgalParamEnd
//    \cgalParamBegin{geom_traits} a geometric traits class instance.
//       The traits class must provide the nested types `Point_3` and `Vector_3`,
//       and the nested functors :
//         - `Construct_bbox_3` to construct a bounding box of a point,
//         - `Compute_squared_distance_3` to compute the distance between two points,
//
//       and, for each functor `Foo`, a function `Foo foo_object()`
//   \cgalParamEnd
// \cgalNamedParamsEnd
//
// @param tnp optional \ref pmp_namedparameters "Named Parameters" related to the target mesh,
//            amongst those described below:
//
// \cgalNamedParamsBegin
//    \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of the target mesh.
//                                      The type of this map is model of `ReadablePropertyMap`.
//                                      If this parameter is omitted, an internal property map for
//                                      `CGAL::vertex_point_t` must be available in `PolygonMesh`
//    \cgalParamEnd
// \cgalNamedParamsEnd
//
// @return the number of snapped vertices
//
template <typename PolygonMesh,
          typename SourceHalfedgeRange, typename TargetHalfedgeRange,
          typename ToleranceMap,
          typename SourceNamedParameters, typename TargetNamedParameters>
std::size_t snap_vertex_range_onto_vertex_range(const SourceHalfedgeRange& source_hrange,
                                                PolygonMesh& smesh,
                                                const TargetHalfedgeRange& target_hrange,
                                                const PolygonMesh& tmesh,
                                                const ToleranceMap& tol_pmap,
                                                const SourceNamedParameters& snp,
                                                const TargetNamedParameters& tnp)
{
  using boost::get_param;
  using boost::choose_param;

  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor      vertex_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor    halfedge_descriptor;

  typedef CGAL::Box_intersection_d::Box_with_info_d<double, 3, vertex_descriptor>     Box;

  typedef typename GetVertexPointMap<PolygonMesh, SourceNamedParameters>::type        SVPM;
  typedef typename GetVertexPointMap<PolygonMesh, TargetNamedParameters>::const_type  TVPM;
  typedef typename boost::property_traits<TVPM>::value_type                           Point;

  typedef typename GetGeomTraits<PolygonMesh, SourceNamedParameters>::type            GT;

  if(is_empty_range(source_hrange.begin(), source_hrange.end()) ||
     is_empty_range(target_hrange.begin(), target_hrange.end()))
    return 0;

  CGAL_static_assertion((boost::is_same<Point, typename GT::Point_3>::value));

  GT gt = choose_param(get_param(snp, internal_np::geom_traits), GT());

  SVPM svpm = choose_param(get_param(snp, internal_np::vertex_point),
                           get_property_map(vertex_point, smesh));

  TVPM tvpm = choose_param(get_param(tnp, internal_np::vertex_point),
                           get_const_property_map(vertex_point, tmesh));

  // Try to snap vertices
  std::vector<Box> boxes;
  BOOST_FOREACH(halfedge_descriptor hd, source_hrange)
  {
    const vertex_descriptor vd = target(hd, smesh);
    const double eps = CGAL::to_double(get(tol_pmap, vd));
    const Bbox_3 pb = gt.construct_bbox_3_object()(get(svpm, vd));
    const Bbox_3 b(pb.xmin() - eps, pb.ymin() - eps, pb.zmin() - eps, pb.xmax() + eps, pb.ymax() + eps, pb.zmax() + eps);
    boxes.push_back(Box(b, vd));
  }

  std::vector<Box> target_boxes;
  BOOST_FOREACH(halfedge_descriptor hd, target_hrange)
  {
    const vertex_descriptor vd = target(hd, tmesh);
    const Point& p = get(tvpm, vd);
    target_boxes.push_back(Box(gt.construct_bbox_3_object()(p), vd));
  }

  // the correspondence map, multiset of targets because the mapping is not necessarily surjective
  typedef boost::bimap<boost::bimaps::set_of<vertex_descriptor /*source*/>,
                       boost::bimaps::multiset_of<vertex_descriptor /*target*/> >  Vertex_correspondence_map;
  typedef typename Vertex_correspondence_map::right_iterator                       VCM_right_it;
  typedef Vertex_proximity_report<PolygonMesh, GT, Vertex_correspondence_map,
                                  SVPM, TVPM, ToleranceMap, Box>                   Reporter;

  Vertex_correspondence_map vertex_map;
  Reporter vpr(vertex_map, svpm, tvpm, tol_pmap, gt);

  // Shenanigans to pass a reference as callback(which is copied by value by 'box_intersection_d')
  boost::function<void(const Box&, const Box&)> callback(boost::ref(vpr));

  CGAL::box_intersection_d(boxes.begin(), boxes.end(),
                           target_boxes.begin(), target_boxes.end(),
                           callback);

  if(vertex_map.empty())
    return 0;

  const bool is_same_mesh = (&smesh == &tmesh);

  std::size_t counter = 0;

  // Now, move the source vertices when the mapping is surjective
  VCM_right_it vmc_it = vertex_map.right.begin();
  VCM_right_it last = --(vertex_map.right.end());
  VCM_right_it end = vertex_map.right.end();
  for(; vmc_it!=end;)
  {
    const vertex_descriptor vs = vmc_it->second;
    const vertex_descriptor vt = vmc_it->first;
    CGAL_assertion(vs != vt);

    // Check that the next iterator is not also the same target vertex, otherwise that means that
    // we have multiple source vertices snapping to the same target vertex
    // In this case, ignore all those mappings.

    if(vmc_it == last)
    {
      ++counter;
      put(svpm, vs, get(tvpm, vt));
      return counter;
    }

    bool skipped = false;
    VCM_right_it next_it = vmc_it;
    ++next_it;

    // As long as we are on the same target, ignore sources
    while(vt == next_it->first)
    {
      skipped = true;

      if(next_it == last)
        return counter;

      ++next_it;
    }

    if(skipped)
    {
      vmc_it = next_it;
    }
    else // a single target, thus move the source to the target
    {
      ++counter;
      ++vmc_it;
      put(svpm, vs, get(tvpm, vt));

      // if both ranges are from the same mesh, we don't want to move 'vs' on 'vt' and then 'vt'
      // somewhere else. Thus, if 'vs' gets moved to 'vt', lock 'vt' in place by removing from
      // the list of positions to change.
      //
      // @todo On paper, we could have that the lock of 'vt' means that another vertex 'vt2' can
      // now move because 'vt' and 'vt2' were moving to the same vertex...
      if(is_same_mesh)
        vertex_map.right.erase(vt);
    }
  }

  return counter;
}

template <typename PolygonMesh, typename SourceHalfedgeRange, typename TargetHalfedgeRange,
          typename SourceNamedParameters, typename TargetNamedParameters>
std::size_t snap_vertex_range_onto_vertex_range(const SourceHalfedgeRange& source_hrange,
                                                PolygonMesh& smesh,
                                                const TargetHalfedgeRange& target_hrange,
                                                const PolygonMesh& tmesh,
                                                const SourceNamedParameters& snp,
                                                const TargetNamedParameters& tnp)
{
  typedef typename GetGeomTraits<PolygonMesh, SourceNamedParameters>::type            GT;
  typedef typename GT::FT                                                             FT;
  typedef CGAL::dynamic_vertex_property_t<FT>                                         Vertex_property_tag;
  typedef typename boost::property_map<PolygonMesh, Vertex_property_tag>::type        Tolerance_map;

  Tolerance_map tol_pmap = get(Vertex_property_tag(), smesh);
  compute_tolerance_at_vertices(source_hrange, tol_pmap, smesh, snp);

  return snap_vertex_range_onto_vertex_range(source_hrange, smesh, target_hrange, tmesh, tol_pmap, snp, tnp);
}

template <typename PolygonMesh, typename SourceHalfedgeRange, typename TargetHalfedgeRange, typename ToleranceMap>
std::size_t snap_vertex_range_onto_vertex_range(const SourceHalfedgeRange& source_hrange,
                                                PolygonMesh& smesh,
                                                const TargetHalfedgeRange& target_hrange,
                                                const PolygonMesh& tmesh,
                                                const ToleranceMap& tol_pmap)
{
  return snap_vertex_range_onto_vertex_range(source_hrange, smesh, target_hrange, tmesh, tol_pmap,
                                             CGAL::parameters::all_default(),
                                             CGAL::parameters::all_default());
}

template <typename PolygonMesh, typename SourceHalfedgeRange, typename TargetHalfedgeRange>
std::size_t snap_vertex_range_onto_vertex_range(const SourceHalfedgeRange& source_hrange,
                                                PolygonMesh& smesh,
                                                const TargetHalfedgeRange& target_hrange,
                                                const PolygonMesh& tmesh)
{
  return snap_vertex_range_onto_vertex_range(source_hrange, smesh, target_hrange, tmesh,
                                             CGAL::parameters::all_default(),
                                             CGAL::parameters::all_default());
}

template <typename PolygonMesh>
std::size_t snap_vertex_range_onto_vertex_range(PolygonMesh& smesh, const PolygonMesh& tmesh)
{
  return snap_vertex_range_onto_vertex_range(halfedges(smesh), smesh, halfedges(tmesh), tmesh);
}

template <typename PolygonMesh, typename HalfedgeRange, typename ToleranceMap>
std::size_t snap_border_vertices_onto_vertex_range(PolygonMesh& smesh,
                                                   const HalfedgeRange& target_hrange,
                                                   const PolygonMesh& tmesh,
                                                   const ToleranceMap& tol_pmap)
{
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor    halfedge_descriptor;

  std::vector<halfedge_descriptor> border_vertices;
  border_halfedges(smesh, std::back_inserter(border_vertices));

  return snap_vertex_range_onto_vertex_range(border_vertices, smesh, target_hrange, tmesh, tol_pmap);
}

template <typename PolygonMesh, typename HalfedgeRange>
std::size_t snap_border_vertices_onto_vertex_range(PolygonMesh& smesh,
                                                   const HalfedgeRange& target_hrange,
                                                   const PolygonMesh& tmesh)
{
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor    halfedge_descriptor;

  std::vector<halfedge_descriptor> border_vertices;
  border_halfedges(smesh, std::back_inserter(border_vertices));

  return snap_vertex_range_onto_vertex_range(border_vertices, smesh, target_hrange, tmesh);
}

// \ingroup PMP_repairing_grp
//
// Attempts to snap the border vertices of the source mesh onto the vertices of the target mesh.
//
// A vertex of the source range is only snapped to a vertex of the target mesh
// if its distance to the target mesh vertex is smaller than a user-chosen bound.
// If any source vertex can be snapped onto multiple vertices of the target
// range, the closest one is chosen.
// If multiple vertices within the source range are to be snapped to the same target vertex,
// then the snapping is not performed for these vertices.
//
// @tparam PolygonMesh a model of `FaceListGraph` and `MutableFaceGraph`
// @tparam ToleranceMap a model of `ReadablePropertyMap` with key type `boost::graph_traits<PolygonMesh>::%vertex_descriptor`
//                      and value type the number type associated with the traits of the mesh.
//
// @param smesh the source mesh whose border vertices might be moved
// @param tmesh the target mesh whose vertices are potential projection targets
// @param tol_pmap a tolerance map associating to each vertex of the source range a tolerance value:
//               potential projection targets are sought in a sphere centered at the vertex and
//               whose radius is the tolerance value.
//
// @pre `smesh` and `tmesh` are different meshes
//
// \return the number of snapped vertices
//
template <typename PolygonMesh, typename ToleranceMap>
std::size_t snap_border_vertices_onto_vertex_range(PolygonMesh& smesh,
                                                   const PolygonMesh& tmesh,
                                                   const ToleranceMap& tol_pmap)
{
  return snap_border_vertices_onto_vertex_range(smesh, halfedges(tmesh), tmesh, tol_pmap);
}

template <typename PolygonMesh>
std::size_t snap_border_vertices_onto_vertex_range(PolygonMesh& smesh, const PolygonMesh& tmesh)
{
  return snap_border_vertices_onto_vertex_range(smesh, halfedges(tmesh), tmesh);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

// Adapted from <CGAL/internal/AABB_tree/AABB_traversal_traits.h>
template <typename AABBTraits>
class Projection_traits
{
  typedef typename AABBTraits::FT                               FT;
  typedef typename AABBTraits::Point_3                          Point;
  typedef typename AABBTraits::Primitive                        Primitive;
  typedef typename AABBTraits::Bounding_box                     Bounding_box;
  typedef typename AABBTraits::Primitive::Id                    Primitive_id;
  typedef typename AABBTraits::Point_and_primitive_id           Point_and_primitive_id;
  typedef typename AABBTraits::Object_and_primitive_id          Object_and_primitive_id;

  typedef CGAL::AABB_node<AABBTraits>                           Node;

public:
  explicit Projection_traits(const AABBTraits& traits)
    : m_traits(traits),
      m_closest_point_initialized(false)
  {}

  bool go_further() const { return true; }

  void intersection(const Point& query, const Primitive& primitive)
  {
    // skip the primitive if one of its endpoints is the query
    typename Primitive::Datum s = primitive.datum(m_traits.shared_data());
    if(s[0] == query || s[1] == query)
      return;

    if(!m_closest_point_initialized)
    {
      m_closest_point_initialized = true;
      m_closest_primitive = primitive.id();
      m_closest_point = primitive.reference_point(m_traits.shared_data());
    }

    Point new_closest_point = m_traits.closest_point_object()(query, primitive, m_closest_point);
    if(new_closest_point != m_closest_point)
    {
      m_closest_primitive = primitive.id();
      m_closest_point = new_closest_point; // this effectively shrinks the sphere
    }
  }

  bool do_intersect(const Point& query, const Node& node) const
  {
    return !m_closest_point_initialized ||
           m_traits.compare_distance_object()(query, node.bbox(), m_closest_point) == CGAL::SMALLER;
  }

  Point closest_point() const { return m_closest_point; }
  typename Primitive::Id closest_primitive_id() const { return m_closest_primitive; }
  bool closest_point_initialized() const { return m_closest_point_initialized; }

private:
  Point m_closest_point;
  typename Primitive::Id m_closest_primitive;
  const AABBTraits& m_traits;
  bool m_closest_point_initialized;
};

template <typename PolygonMesh, typename GeomTraits, typename VPM>
struct Compare_points_along_edge
{
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor    halfedge_descriptor;
  typedef typename GeomTraits::Point_3                                      Point;

  Compare_points_along_edge(halfedge_descriptor support,
                            const GeomTraits& gt, const VPM& vpm, const PolygonMesh& pm)
    : gt(gt), m_src(get(vpm, target(opposite(support, pm), pm)))
  { }

  bool operator()(const Point& p1, const Point& p2) const {
    return gt.less_distance_to_point_3_object()(m_src, p1, p2);
  }

  const GeomTraits& gt;
  Point m_src;
};

// Try to merge vertices in 'hrange' (vd = target(hd, pm)) onto the edges in the 'erange', non-conformingly
template <typename HalfedgeRange, typename TriangleMesh,
          typename ToleranceMap, typename NamedParameters>
std::size_t snap_vertex_range_onto_vertex_range_non_conforming(const HalfedgeRange source_hrange,
                                                               const HalfedgeRange target_hrange,
                                                               TriangleMesh& pm,
                                                               const ToleranceMap& tol_pmap,
                                                               const NamedParameters& np)
{
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor     vertex_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor   halfedge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::edge_descriptor       edge_descriptor;

  typedef typename GetVertexPointMap<TriangleMesh, NamedParameters>::type   VPM;

  typedef typename GetGeomTraits<TriangleMesh, NamedParameters>::type       GT;
  typedef typename GT::FT                                                   FT;
  typedef typename GT::Point_3                                              Point;

  typedef Compare_points_along_edge<TriangleMesh, GT, VPM>                  Point_along_edge_comparer;

  typedef CGAL::AABB_halfedge_graph_segment_primitive<TriangleMesh, VPM>    Primitive;
  typedef CGAL::AABB_traits<GT, Primitive>                                  AABB_Traits;
  typedef CGAL::AABB_tree<AABB_Traits>                                      AABB_tree;

  using boost::get_param;
  using boost::choose_param;

  // Need to have a triangle mesh to ensure that the refinement point of a border edge 'e' has visibility
  // of the third point of the face 'f' incident to 'e'
  CGAL_precondition(is_triangle_mesh(pm));

  std::size_t snapped_n = 0;

  VPM vpm = choose_param(get_param(np, internal_np::vertex_point), get_property_map(vertex_point, pm));
  const GT gt = choose_param(get_param(np, internal_np::geom_traits), GT());

  // start by snapping vertices together to simplify things
  snapped_n = snap_vertex_range_onto_vertex_range(source_hrange, pm, target_hrange, pm, tol_pmap, np, np);

  std::map<Point, unsigned int> multiplicity;
  BOOST_FOREACH(halfedge_descriptor hd, source_hrange)
  {
    vertex_descriptor vd = target(hd, pm);
    ++multiplicity.insert(std::make_pair(get(vpm, vd), 0)).first->second;
  }

  // Since we're inserting primitives one by one, we can't pass this shared data in the constructor of the tree
  AABB_Traits aabb_traits;
  aabb_traits.set_shared_data(pm, vpm);

  // Fill the AABB-tree
  AABB_tree aabb_tree(aabb_traits);
  BOOST_FOREACH(halfedge_descriptor hd, target_hrange)
  {
    CGAL_precondition(is_border(edge(hd, pm), pm));
    aabb_tree.insert(Primitive(edge(hd, pm), pm, vpm));
  }

  // Collect border points that can be projected onto a border edge
  std::vector<std::pair<Point, halfedge_descriptor> > edges_to_split;
  BOOST_FOREACH(halfedge_descriptor hd, source_hrange)
  {
    const vertex_descriptor vd = target(hd, pm);
    const Point& query = get(vpm, vd);
    if(multiplicity[query] > 1)
      continue; // skip points that are already attached to another border

    // use the current halfedge as hint
    Projection_traits<AABB_Traits> traversal_traits(aabb_tree.traits());
    aabb_tree.traversal(query, traversal_traits);

    if(!traversal_traits.closest_point_initialized())
      continue;

    const FT sq_eps = CGAL::square(get(tol_pmap, vd));

    const Point& closest_p = traversal_traits.closest_point();
    if(gt.compute_squared_distance_3_object()(query, closest_p) < sq_eps)
    {
      edge_descriptor closest = traversal_traits.closest_primitive_id();
      CGAL_assertion(get(vpm, source(closest, pm)) != query &&
                     get(vpm, target(closest, pm)) != query);

      halfedge_descriptor clos_hd = halfedge(closest, pm);
      CGAL_assertion(is_border(edge(clos_hd, pm), pm));

      if(!is_border(clos_hd, pm))
        clos_hd = opposite(clos_hd, pm);

      const Point new_pos = gt.construct_midpoint_3_object()(query, closest_p);
      put(vpm, vd, new_pos);

      edges_to_split.push_back(std::make_pair(new_pos, clos_hd));
    }
  }

  // Sort points falling on the same edge and split the edge
  std::map<halfedge_descriptor, std::vector<Point> > points_per_edge;
  typedef std::pair<Point, halfedge_descriptor> Pair_type;
  BOOST_FOREACH(const Pair_type& p, edges_to_split)
    points_per_edge[p.second].push_back(p.first);

  typedef std::pair<const halfedge_descriptor, std::vector<Point> > Map_type;
  BOOST_FOREACH(Map_type& mt, points_per_edge)
  {
    const halfedge_descriptor mt_hd = mt.first;
    const halfedge_descriptor mt_hd_opp = opposite(mt_hd, pm);
    CGAL_assertion(!is_border(mt_hd_opp, pm));

    std::cout <<" --" << std::endl;
    std::cout << "hd: " << get(vpm, source(mt_hd, pm)) << " -- " << get(vpm, target(mt_hd, pm)) << std::endl;

    if(mt.second.size() > 1)
    {
      /// \todo Sorting the projected positions is too simple (for example, this doesn't work
      ///       for a zigzaging polyline). Rather, points should be sorted following
      ///       the order along the matching border. Note that this requires identifying
      ///       matching polylines because a sequence of zigzaging points might not all project
      ///       onto the same halfedge.
      Point_along_edge_comparer cmp(mt_hd, gt, vpm, pm);
      std::sort(mt.second.begin(), mt.second.end(), cmp);
    }

    BOOST_FOREACH(const Point& p, mt.second)
    {
      halfedge_descriptor res = CGAL::Euler::split_edge(mt_hd, pm);
      put(vpm, target(res, pm), p);
      ++snapped_n;

      CGAL::Euler::split_face(mt_hd_opp, prev(prev(mt_hd_opp, pm), pm), pm);
    }
  }

  return snapped_n;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

template <typename TriangleMesh>
std::size_t snap_border_vertices_non_conforming(TriangleMesh& pm)
{
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor        halfedge_descriptor;

  typedef typename GetGeomTraits<TriangleMesh>::type                             GT;
  typedef typename GT::FT                                                        FT;

  typedef CGAL::dynamic_vertex_property_t<FT>                                    Vertex_property_tag;
  typedef typename boost::property_map<TriangleMesh, Vertex_property_tag>::type  Tolerance_map;

  std::vector<halfedge_descriptor> border_vertices;
  border_halfedges(pm, std::back_inserter(border_vertices));

  Tolerance_map tol_pmap = get(Vertex_property_tag(), pm);
  compute_tolerance_at_vertices(border_vertices, tol_pmap, pm, CGAL::parameters::all_default());

  return snap_vertex_range_onto_vertex_range_non_conforming(border_vertices, border_vertices, pm, tol_pmap, CGAL::parameters::all_default());
}

template <typename TriangleMesh, typename ToleranceMap, typename NamedParameters>
std::size_t snap_border_vertices_non_conforming(TriangleMesh& pm,
                                                const ToleranceMap& tol_pmap,
                                                const NamedParameters& np)
{
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor        halfedge_descriptor;

  std::vector<halfedge_descriptor> border_vertices;
  border_halfedges(pm, std::back_inserter(border_vertices));

  return snap_vertex_range_onto_vertex_range_non_conforming(border_vertices, border_vertices, pm, tol_pmap, np);
}

// Attempt to snap the vertices of the border of 'pm' onto border edges of 'pm'.
template <typename TriangleMesh, typename ToleranceMap>
std::size_t snap_border_vertices_non_conforming(TriangleMesh& pm,
                                                const ToleranceMap& tol_pmap)
{
  return snap_border_vertices_non_conforming(pm, tol_pmap, CGAL::parameters::all_default());
}

} // end namespace internal

} // end namespace Polygon_mesh_processing

} // end namespace CGAL

#endif // CGAL_POLYGON_MESH_PROCESSING_INTERNAL_SNAP_H
