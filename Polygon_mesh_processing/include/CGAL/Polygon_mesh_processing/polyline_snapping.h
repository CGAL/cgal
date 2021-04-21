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
// Author(s)     : Simon Giraudot

#ifndef CGAL_POLYGON_MESH_PROCESSING_POLYLINE_SNAPPING_H
#define CGAL_POLYGON_MESH_PROCESSING_POLYLINE_SNAPPING_H

#include <CGAL/license/Polygon_mesh_processing/corefinement.h>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/box_intersection_d.h>
#include <CGAL/Iterator_range.h>
#include <CGAL/tuple.h>

#include <CGAL/boost/iterator/counting_iterator.hpp>
#include <boost/range/const_iterator.hpp>
#include <boost/range/value_type.hpp>
#include <boost/foreach.hpp>
#include <boost/iterator/function_output_iterator.hpp>
#include <boost/graph/properties.hpp>

#include <cmath>

/// \cond SKIP_IN_MANUAL

namespace CGAL {
namespace Polygon_mesh_processing{
namespace internal {

template <typename Exact_kernel>
bool exact_snapping (const typename Exact_kernel::Segment_3& s0,
                     const typename Exact_kernel::Segment_3& s1,
                     typename Exact_kernel::Point_3& result)
{
  typedef typename Exact_kernel::Point_3 Point_3;
  typedef typename Exact_kernel::Segment_3 Segment_3;
  typedef typename Exact_kernel::Vector_3 Vector_3;
  typedef typename Exact_kernel::Plane_3 Plane_3;

  CGAL_assertion (s0.source() != s0.target());
  CGAL_assertion (s1.source() != s1.target());

  Vector_3 v0 = s0.to_vector();
  Vector_3 v1 = s1.to_vector();
  Vector_3 normal = CGAL::cross_product (v0, v1);

  // Collinear segments
  if (normal == CGAL::NULL_VECTOR)
  {
    return false;
  }

  Plane_3 plane0 (s0.source(), normal);
  Plane_3 plane1 (s1.source(), normal);

  Segment_3 s0proj (plane1.projection (s0.source()),
                    plane1.projection (s0.target()));
  Segment_3 s1proj (plane0.projection (s1.source()),
                    plane0.projection (s1.target()));

  CGAL_assertion (s0proj.source() != s0proj.target());
  CGAL_assertion (s1proj.source() != s1proj.target());

  if (!CGAL::do_intersect (s0, s1proj) ||
      !CGAL::do_intersect (s1, s0proj))
    return false;

  Point_3 p0, p1;

  typename std::result_of<typename Exact_kernel::Intersect_3(Segment_3, Segment_3)>::type
    intersection0 = intersection(s0, s1proj);
  if (!intersection0)
    return false;

  if (const Point_3* p = boost::get<Point_3>(&*intersection0))
    p0 = *p;
  else // Coplanar segments
  {
    return false;
  }

  typename std::result_of<typename Exact_kernel::Intersect_3(Segment_3, Segment_3)>::type
    intersection1 = intersection(s1, s0proj);
  if (!intersection1)
    return false;

  if (const Point_3* p = boost::get<Point_3>(&*intersection1))
    p1 = *p;
  else // Coplanar segments
  {
    return false;
  }

  result = CGAL::midpoint (p0, p1);

  return true;
}

template <class Kernel,
          bool Has_exact_constructions=
          !std::is_floating_point<typename Kernel::FT>::value >
class Exact_snapping;

template <class Kernel>
class Exact_snapping<Kernel,false>
{
//typedefs
  typedef Kernel Input_kernel;
  typedef CGAL::Exact_predicates_exact_constructions_kernel Exact_kernel;
  typedef CGAL::Cartesian_converter<Exact_kernel, Input_kernel> Exact_to_inexact;
  typedef CGAL::Cartesian_converter<Input_kernel, Exact_kernel> Inexact_to_exact;

  typedef typename Input_kernel::Point_3 Inexact_point_3;
  typedef typename Input_kernel::Line_3 Inexact_line_3;
  typedef typename Input_kernel::Segment_3 Inexact_segment_3;
  typedef typename Exact_kernel::Point_3 Exact_point_3;
  typedef typename Exact_kernel::Line_3 Exact_line_3;
  typedef typename Exact_kernel::Segment_3 Exact_segment_3;

  Exact_kernel ek;
  Exact_to_inexact to_inexact;
  Inexact_to_exact to_exact;

public:

  std::pair<Inexact_point_3, bool> intersection (const Inexact_segment_3& s0, const Inexact_segment_3& s1) const
  {
    Exact_segment_3 exact_s0 = to_exact(s0);
    Exact_segment_3 exact_s1 = to_exact(s1);
    Exact_point_3 result;

    if (exact_snapping<Exact_kernel> (exact_s0, exact_s1, result))
      return std::make_pair (to_inexact (result), true);

    return std::make_pair (Inexact_point_3(), false);
  }

};

template <class Kernel>
class Exact_snapping<Kernel,true>
{
//typedefs
  typedef typename Kernel::Point_3 Point_3;
  typedef typename Kernel::Line_3 Line_3;
  typedef typename Kernel::Segment_3 Segment_3;
public:

  Exact_snapping()
  {

  }

  std::pair<Point_3, bool> intersection (const Segment_3& s0, const Segment_3& s1) const
  {
    Point_3 result;
    bool out = exact_snapping<Kernel> (s0, s1, result);
    return std::make_pair (result, out);
  }

};

}// namespace internal


// Note this is not officially documented
/*!
  Computes snapping points among a set of polylines.

  The output is given as a set of `std::tuple` with the following
  elements in the following order:

  - index `i` of the first polyline
  - index `ii` of the segment in polyline `i`
  - index `j` of the second polyline
  - index `jj` of the segment in polyline `j`
  - barycentric coordinate `c_i` of snapping point `p_i` with respect
    to points `p[i][ii]` and `p[i][ii+1]`
  - barycentric coordinate `c_j` of snapping point `p_j` with respect
    to points `p[j][jj]` and `p[j][jj+1]`

  \tparam PolylineRange a `RandomAccessRange` of `RandomAccessRange`
  of `Kernel::Point_3`.
  \tparam OutputIterator a model of `OutputIterator` that accepts
  objects of type `std::tuple<std::size_t, std::size_t, std::size_t, std::size_t, FT, FT>`
  \tparam Kernel a model of `Kernel`

  \param polylines the input polyline range.
  \param tolerance the minimum distance between two polylines such
  that the snapping points are computed
  \param output output iterator

 */
template <class PolylineRange, class OutputIterator, class Kernel>
void polyline_snapping (const PolylineRange& polylines,
                        double tolerance,
                        OutputIterator output,
                        const Kernel&)
{
  typedef typename CGAL::Box_intersection_d::Box_with_info_d<double, 3, std::pair<std::size_t, std::size_t> > Box;
  typedef typename Kernel::FT FT;
  typedef typename Kernel::Point_3 Point_3;
  typedef typename Kernel::Segment_3 Segment_3;

  std::size_t size = 0;
  for (std::size_t i = 0; i < polylines.size(); ++ i)
    size += polylines[i].size() - 1;

  std::vector<Box> boxes;
  boxes.reserve (size);
  for (std::size_t i = 0; i < polylines.size(); ++ i)
    for (std::size_t j = 0; j < polylines[i].size()-1; ++ j)
    {
      const Point_3& p1 = polylines[i][j];
      const Point_3& p2 = polylines[i][j+1];
      CGAL::Bbox_3 bbox = p1.bbox() + p2.bbox();
      bbox = CGAL::Bbox_3 (bbox.xmin() - tolerance / 2.,
                           bbox.ymin() - tolerance / 2.,
                           bbox.zmin() - tolerance / 2.,
                           bbox.xmax() + tolerance / 2.,
                           bbox.ymax() + tolerance / 2.,
                           bbox.zmax() + tolerance / 2.);

      boxes.push_back(Box(bbox, std::make_pair (i, j)));
    }

  internal::Exact_snapping<Kernel> exact_snapping;
  std::ptrdiff_t cutoff = 2000;

  std::size_t idx_start = 0;
  for (std::size_t i = 0; i < polylines.size() - 1; ++ i)
  {
    std::vector<const Box*> box1_ptr
      (boost::make_counting_iterator<const Box*>(&boxes[idx_start]),
       boost::make_counting_iterator<const Box*>(&boxes[idx_start + polylines[i].size() - 1]));
    std::vector<const Box*> box2_ptr
      (boost::make_counting_iterator<const Box*>(&boxes[idx_start + polylines[i].size() - 1]),
       boost::make_counting_iterator<const Box*>(&boxes.back()));
    idx_start += polylines[i].size() - 1;

    CGAL::box_intersection_d
      (box1_ptr.begin(), box1_ptr.end(),
       box2_ptr.begin(), box2_ptr.end(),
       [&](const Box* b, const Box* c)
       {
         Segment_3 s1 (polylines[b->info().first][b->info().second],
                       polylines[b->info().first][b->info().second + 1]);
         Segment_3 s2 (polylines[c->info().first][c->info().second],
                       polylines[c->info().first][c->info().second + 1]);
         if (CGAL::squared_distance (s1, s2) > tolerance * tolerance)
           return;

         Point_3 new_point;
         bool okay;
         std::tie (new_point, okay) = exact_snapping.intersection (s1, s2);

         if (!okay)
           return;

         FT coord1 = CGAL::approximate_sqrt (CGAL::squared_distance (s1.source(), s1.supporting_line().projection(new_point))
                                             / s1.squared_length());
         FT coord2 = CGAL::approximate_sqrt (CGAL::squared_distance (s2.source(), s2.supporting_line().projection(new_point))
                                             / s2.squared_length());

         *(output ++) = std::make_tuple (b->info().first, b->info().second,
                                         c->info().first, c->info().second,
                                         coord1, coord2);
       },
       cutoff);
  }
}

template <class PolylineRange, class Kernel>
void polyline_snapping (PolylineRange& polylines,
                        double tolerance,
                        const Kernel&)
{
  using FT = typename Kernel::FT;
  using Point_3 = typename Kernel::Point_3;
  using Polyline = typename PolylineRange::value_type;

  using Points = std::vector<std::pair<FT, Point_3> >;
  using Segment_points = std::vector<Points>;
  using Polyline_points = std::vector<Segment_points>;
  using Snapping_point = std::tuple<std::size_t, std::size_t, std::size_t, std::size_t, FT, FT>;

  Polyline_points polyline_points (polylines.size());
  for (std::size_t i = 0; i < polylines.size(); ++ i)
    polyline_points[i].resize (polylines[i].size() - 1);

  polyline_snapping
    (polylines, tolerance,
     boost::make_function_output_iterator
     ([&](const Snapping_point& s)
      {
        Point_3 a = CGAL::barycenter (polylines[get<0>(s)][get<1>(s)],     (1. - get<4>(s)),
                                      polylines[get<0>(s)][get<1>(s) + 1],       get<4>(s));
        Point_3 b = CGAL::barycenter (polylines[get<2>(s)][get<3>(s)],     (1. - get<5>(s)),
                                      polylines[get<2>(s)][get<3>(s) + 1],       get<5>(s));

        Point_3 p = CGAL::midpoint(a,b);

        polyline_points[get<0>(s)][get<1>(s)].emplace_back(get<4>(s), p);
        polyline_points[get<2>(s)][get<3>(s)].emplace_back(get<5>(s), p);
      }),
     Kernel());

  for (std::size_t npoly = 0; npoly < polylines.size(); ++ npoly)
  {
    Polyline& input = polylines[npoly];
    Polyline polyline;

    Segment_points& seg_points = polyline_points[npoly];

    for (std::size_t i = 0; i < input.size() - 1; ++ i)
    {
      polyline.push_back(input[i]);

      Points& points = seg_points[i];
      if (points.size() > 1)
        std::sort (points.begin(), points.end(),
                   [](const std::pair<FT, Point_3>& a, const std::pair<FT, Point_3>& b) -> bool
                   { return a.first < b.first; });
      for (const auto& p : points)
        polyline.push_back(p.second);
    }
    polyline.push_back(input.back());

    input.swap(polyline);
  }
}

template <class Graph, class PointMap>
void polyline_snapping (Graph& graph, PointMap point_map, double tolerance)
{
  using Point_3 = typename boost::property_traits<PointMap>::value_type;
  using Kernel = typename Kernel_traits<Point_3>::Kernel;
  using Segment_3 = typename Kernel::Segment_3;
  using FT = typename Kernel::FT;

  using vertex_descriptor = typename boost::graph_traits<Graph>::vertex_descriptor;
  using edge_descriptor = typename boost::graph_traits<Graph>::edge_descriptor;

  using Box = Box_intersection_d::Box_with_info_d<double, 3, edge_descriptor>;

  std::vector<Box> boxes;
  boxes.reserve (num_edges(graph));

  for (edge_descriptor ed : make_range(edges(graph).first, edges(graph).second))
  {
    vertex_descriptor vs = source (ed, graph);
    vertex_descriptor vt = target (ed, graph);

    const Point_3& ps = get (point_map, vs);
    const Point_3& pt = get (point_map, vt);

    Bbox_3 bbox = ps.bbox() + pt.bbox();
    bbox = Bbox_3 (bbox.xmin() - tolerance / 2.,
                   bbox.ymin() - tolerance / 2.,
                   bbox.zmin() - tolerance / 2.,
                   bbox.xmax() + tolerance / 2.,
                   bbox.ymax() + tolerance / 2.,
                   bbox.zmax() + tolerance / 2.);

    boxes.emplace_back(bbox, ed);
  }

  internal::Exact_snapping<Kernel> exact_snapping;
  std::ptrdiff_t cutoff = 2000;

  using Vertex_position = std::pair<FT, vertex_descriptor>;
  std::map<edge_descriptor, std::vector<Vertex_position>> new_vertices;

  // Compute intersections
  box_self_intersection_d
    (boxes.begin(), boxes.end(),
     [&](const Box& a, const Box& b)
     {
       edge_descriptor ea = a.info();
       edge_descriptor eb = b.info();

       vertex_descriptor vas = source (ea, graph);
       vertex_descriptor vat = target (ea, graph);
       vertex_descriptor vbs = source (eb, graph);
       vertex_descriptor vbt = target (eb, graph);

       // Skip adjacent edges
       if (vas == vbs || vat == vbs || vas == vbt || vat == vbt)
         return;

       // Skip segments too far apart
       Segment_3 sa (get (point_map, vas), get (point_map, vat));
       Segment_3 sb (get (point_map, vbs), get (point_map, vbt));
       if (squared_distance (sa, sb) > tolerance * tolerance)
         return;

       // Compute exact snapping point
       Point_3 new_point;
       bool okay;
       std::tie (new_point, okay) = exact_snapping.intersection (sa, sb);
       if (!okay)
         return;

       // Create new vertex
       vertex_descriptor vd = add_vertex(graph);
       put (point_map, vd, new_point);

       // Compute barycentric coord along each segment
       FT ba = approximate_sqrt (squared_distance (sa.source(), sa.supporting_line().projection(new_point))
                                 / sa.squared_length());
       FT bb = approximate_sqrt (squared_distance (sb.source(), sb.supporting_line().projection(new_point))
                                 / sb.squared_length());

       // Store relative position of new vertex on edges
       auto ma = new_vertices.insert (std::make_pair(ea, std::vector<Vertex_position>()));
       auto mb = new_vertices.insert (std::make_pair(eb, std::vector<Vertex_position>()));
       ma.first->second.emplace_back (ba, vd);
       mb.first->second.emplace_back (bb, vd);
     },
     cutoff);

  // Refine edges
  for (auto& m : new_vertices)
  {
    edge_descriptor ed = m.first;
    std::vector<Vertex_position>& vertices = m.second;
    vertices.emplace_back (FT(0), source(ed, graph));
    vertices.emplace_back (FT(1), target(ed, graph));

    // Sort vertices along the edge
    std::sort (vertices.begin(), vertices.end(),
               [](const Vertex_position& a, const Vertex_position& b) -> bool
               { return a.first < b.first; });

    // Remove current edge
    remove_edge(ed, graph);

    // Add edges after subdivision
    for (std::size_t i = 0; i < vertices.size() - 1; ++ i)
      add_edge (vertices[i].second, vertices[i+1].second, graph);
  }

#if 0
  // Detect edges smaller than threshold
  using vedge = std::pair<vertex_descriptor, vertex_descriptor>;
  std::vector<vedge> vedges;
  for (edge_descriptor ed : make_range(edges(graph).first, edges(graph).second))
  {
    vertex_descriptor vs = source (ed, graph);
    vertex_descriptor vt = target (ed, graph);
    vedges.emplace_back (vs, vt);
  }

  // Map deleted vertex -> valid vertex
  std::map<vertex_descriptor, vertex_descriptor> map_v2v;
  // Map vertex -> barycenter weight
  std::map<vertex_descriptor, std::size_t> map_v2b;
  for (vertex_descriptor vd : make_range(vertices(graph).first, vertices(graph).second))
  {
    map_v2v.insert (std::make_pair(vd, vd));
    map_v2b.insert (std::make_pair(vd, 1));
  }

  // Collapse edges
  for (vedge ve : vedges)
  {
    vertex_descriptor vs = map_v2v[ve.first];
    vertex_descriptor vt = map_v2v[ve.second];

    if (vs == vt)
      continue;

    const Point_3& ps = get (point_map, vs);
    const Point_3& pt = get (point_map, vt);

    if (squared_distance(ps, ps) > tolerance * tolerance)
      continue;

    std::size_t& bs = map_v2b[vs];
    std::size_t& bt = map_v2b[vt];
    Point_3 bary = barycenter (ps, bs, pt, bt);
    bs += bt;

    put (point_map, vs, bary);

    map_v2v[vt] = vs;
  }
#endif

}


} } //end of namespace CGAL::Polygon_mesh_processing

/// \endcond

#endif // CGAL_POLYGON_MESH_PROCESSING_POLYLINE_SNAPPING_H
