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

#include <CGAL/assertions.h>
#include <CGAL/box_intersection_d.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/boost/graph/iterator.h>
#include <CGAL/Bbox_3.h>
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
#include <fstream>
#include <limits>
#include <utility>
#include <vector>

namespace CGAL {

namespace Polygon_mesh_processing {

namespace internal {

// Assigns at each vertex the length of its shortest incident edge as 'epsilon' value
template <typename ToleranceMap,
          typename PolygonMesh,
          typename SourceNamedParameters>
void compute_tolerance_at_vertices(ToleranceMap& tol_vm,
                                   PolygonMesh& smesh,
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
                           get_property_map(vertex_point, smesh));

  std::vector<halfedge_descriptor> border_vertices;
  border_halfedges(smesh, std::back_inserter(border_vertices));

  BOOST_FOREACH(halfedge_descriptor hd, border_vertices)
  {
    CGAL_assertion(CGAL::is_border(hd, smesh));
    vertex_descriptor vd = target(hd, smesh);

    CGAL::Halfedge_around_target_iterator<PolygonMesh> hit, hend;
    boost::tie(hit, hend) = CGAL::halfedges_around_target(vd, smesh);
    CGAL_assertion(hit != hend);

    FT sq_length = gt.compute_squared_distance_3_object()(get(svpm, source(*hit, smesh)),
                                                          get(svpm, target(*hit, smesh)));
    FT min_sq_dist = sq_length;
    ++hit;

    for(; hit!=hend; ++hit)
    {
      sq_length = gt.compute_squared_distance_3_object()(get(svpm, source(*hit, smesh)),
                                                         get(svpm, target(*hit, smesh)));

      if(sq_length < min_sq_dist)
        min_sq_dist = sq_length;
    }

    put(tol_vm, vd, CGAL::approximate_sqrt(min_sq_dist));
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

  typedef typename GeomTraits::FT                                            FT;
  typedef typename boost::property_traits<SVPM>::value_type                  Point;

  typedef VertexCorrespondenceMap                                            Vertex_correspondence_map;
  typedef typename Vertex_correspondence_map::left_iterator                  VCM_left_iterator;
  typedef typename Vertex_correspondence_map::value_type                     VCM_value_type;

  Vertex_proximity_report(Vertex_correspondence_map& vertex_map,
                          const SVPM& svpm, const TVPM& tvpm,
                          const ToleranceMap& tol_vm,
                          const GeomTraits& gt)
    :
      m_vertex_map(vertex_map),
      svpm(svpm), tvpm(tvpm),
      tol_vm(tol_vm),
      gt(gt)
  { }

  void operator()(const Box& a, const Box& b)
  {
    vertex_descriptor va = a.info();
    vertex_descriptor vb = b.info();

    const Point& sp = get(svpm, va);
    const Point& tp = get(tvpm, vb);
    FT tol = get(tol_vm, va);

    const FT sq_dist = gt.compute_squared_distance_3_object()(sp, tp);
    CGAL::Comparison_result res = CGAL::compare(sq_dist, tol * tol);

    if(res != CGAL::LARGER)
    {
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
  }

private:
  Vertex_correspondence_map& m_vertex_map;
  CGAL::cpp11::unordered_map<vertex_descriptor/*source*/, FT> m_sq_distance_to_snapped_point;

  const SVPM& svpm;
  const TVPM& tvpm;
  const ToleranceMap& tol_vm;
  const GeomTraits& gt;
};

// \ingroup PMP_repairing_grp
//
// Attempts to snap the border vertices of the source mesh onto the target mesh.
// A border vertex of the source mesh is snapped to a vertex of the target mesh
// if its distance to the target mesh vertex is smaller than a user-chosen bound.
// If any source target vertex can be projected onto multiple vertices of the target
// mesh, the snapping process is aborted and no vertex is snapped.
//
// @tparam PolygonMesh a model of `FaceListGraph` and `MutableFaceGraph`
// @tparam SourceNamedParameters a sequence of \ref pmp_namedparameters "Named Parameters"
// @tparam TargetNamedParameters a sequence of \ref pmp_namedparameters "Named Parameters"
//
// @param smesh the source mesh whose border vertices might be moved.
// @param tmesh the target mesh whose vertices are taken as potential new position for the border
//              vertices of the source mesh.
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
//         - `Construct_vector_3` to construct a vector from three coordinates,
//         - `Construct_translated_point_3` to translate a point with a vector
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
// \return the number of snapped vertices
//
template <typename PolygonMesh, typename ToleranceMap,
          typename SourceNamedParameters, typename TargetNamedParameters>
std::size_t snap_border(PolygonMesh& smesh,
                        const PolygonMesh& tmesh,
                        const ToleranceMap& tol_vm,
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
  typedef typename GT::FT                                                             FT;

  CGAL_static_assertion((boost::is_same<Point, typename GT::Point_3>::value));

  GT gt = choose_param(get_param(snp, internal_np::geom_traits), GT());

  SVPM svpm = choose_param(get_param(snp, internal_np::vertex_point),
                           get_property_map(vertex_point, smesh));

  TVPM tvpm = choose_param(get_param(tnp, internal_np::vertex_point),
                           get_const_property_map(vertex_point, tmesh));

  // Extract border vertices
  std::vector<halfedge_descriptor> border_vertices;
  border_halfedges(smesh, std::back_inserter(border_vertices));

  // Try to snap vertices
  std::vector<Box> boxes;
  BOOST_FOREACH(halfedge_descriptor hd, border_vertices)
  {
    vertex_descriptor vd = target(hd, smesh);
    const FT eps = get(tol_vm, vd);

    const typename GT::Vector_3 eps_v = gt.construct_vector_3_object()(eps, eps, eps);
    const typename GT::Vector_3 opp_eps_v = gt.construct_vector_3_object()(-eps, -eps, -eps);
    const Point& pm = gt.construct_translated_point_3_object()(get(svpm, vd), opp_eps_v);
    const Point& pM = gt.construct_translated_point_3_object()(get(svpm, vd), eps_v);

    Bbox_3 bbox = gt.construct_bbox_3_object()(pm) + gt.construct_bbox_3_object()(pM);
    boxes.push_back(Box(bbox, vd));
  }

  std::vector<Box> target_boxes;
  BOOST_FOREACH(vertex_descriptor vd, vertices(tmesh))
  {
    const Point& p = get(tvpm, vd);
    target_boxes.push_back(Box(gt.construct_bbox_3_object()(p), vd));
  }

  // the correspondence map, multiset of targets because the mapping is not necessarily surjective
  typedef boost::bimap<boost::bimaps::set_of<vertex_descriptor /*source*/>,
                       boost::bimaps::multiset_of<vertex_descriptor /*target*/> >  Vertex_correspondence_map;

  typedef typename Vertex_correspondence_map::right_iterator                       VCM_right_it;

  Vertex_correspondence_map vertex_map;

  // Shenanigans to pass a reference as callback (which is copied by value by 'box_intersection_d')
  typedef Vertex_proximity_report<PolygonMesh, GT, Vertex_correspondence_map, SVPM, TVPM, ToleranceMap, Box>     Reporter;
  Reporter vpr(vertex_map, svpm, tvpm, tol_vm, gt);
  boost::function<void(const Box&, const Box&)> callback(boost::ref(vpr));

  CGAL::box_intersection_d(boxes.begin(), boxes.end(),
                           target_boxes.begin(), target_boxes.end(),
                           callback);

  if(vertex_map.empty())
    return 0;

  std::size_t counter = 0;

  // Now, move the source vertices when the mapping is surjective
  VCM_right_it vmc_it = vertex_map.right.begin();
  VCM_right_it last = --(vertex_map.right.end());
  VCM_right_it end = vertex_map.right.end();
  for(; vmc_it!=end;)
  {
    const vertex_descriptor vs = vmc_it->second;
    const vertex_descriptor vt = vmc_it->first;

    // Check that the next iterator is not also the same target vertex, otherwise that means that
    // we have multiple source vertices projecting to the same target vertex
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
    }
  }

  return counter;
}

template <typename PolygonMesh,
          typename SourceNamedParameters, typename TargetNamedParameters>
std::size_t snap_border(PolygonMesh& smesh,
                        const PolygonMesh& tmesh,
                        const double epsilon,
                        const SourceNamedParameters& snp,
                        const TargetNamedParameters& tnp)
{
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor    vertex_descriptor;
  typedef typename GetGeomTraits<PolygonMesh>::type::FT                   FT;

  Constant_property_map<vertex_descriptor, FT> constant_tolerance_map(epsilon);

  return snap_border(smesh, tmesh, constant_tolerance_map, snp, tnp);
}

template <typename PolygonMesh>
std::size_t snap_border(PolygonMesh& smesh,
                        const PolygonMesh& tmesh,
                        const double epsilon)
{
  return snap_border(smesh, tmesh, epsilon,
                     CGAL::parameters::all_default(), CGAL::parameters::all_default());
}

template <typename PolygonMesh,
          typename SourceNamedParameters, typename TargetNamedParameters>
std::size_t snap_border(PolygonMesh& smesh,
                        const PolygonMesh& tmesh,
                        const SourceNamedParameters& snp,
                        const TargetNamedParameters& tnp)
{
  typedef typename GetGeomTraits<PolygonMesh, SourceNamedParameters>::type            GT;
  typedef typename GT::FT                                                             FT;
  typedef CGAL::dynamic_vertex_property_t<FT>                                         Vertex_property_tag;
  typedef typename boost::property_map<PolygonMesh, Vertex_property_tag>::type        Tolerance_map;

  Tolerance_map tol_vm = get(Vertex_property_tag(), smesh);
  compute_tolerance_at_vertices(tol_vm, smesh, snp);

  return snap_border(smesh, tmesh, tol_vm, snp, tnp);
}

template <typename PolygonMesh>
std::size_t snap_border(PolygonMesh& smesh,
                        const PolygonMesh& tmesh)
{
  return snap_border(smesh, tmesh, CGAL::parameters::all_default(), CGAL::parameters::all_default());
}

} // end namespace internal

} // end namespace Polygon_mesh_processing

} // end namespace CGAL

#endif // CGAL_POLYGON_MESH_PROCESSING_INTERNAL_SNAP_H
