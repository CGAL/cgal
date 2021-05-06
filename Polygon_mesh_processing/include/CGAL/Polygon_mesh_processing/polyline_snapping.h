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

template <typename Kernel>
struct Snapping_result
{
  using Point_3 = typename Kernel::Point_3;
  using Segment_3 = typename Kernel::Segment_3;
  enum Type { EMPTY, POINT, SEGMENT };

  Type type;
  boost::variant<Point_3, Segment_3> first;
  boost::variant<Point_3, Segment_3> second;

  Snapping_result() : type (EMPTY) { }
  Snapping_result (const Point_3& first, const Point_3& second)
    : type (POINT), first (first), second (second) { }
  Snapping_result (const Segment_3& first, const Segment_3& second)
    : type (SEGMENT), first (first), second (second) { }
};

template <typename Exact_kernel>
Snapping_result<Exact_kernel>
exact_snapping (const typename Exact_kernel::Segment_3& s0,
                const typename Exact_kernel::Segment_3& s1,
                const typename Exact_kernel::FT& tolerance)
{
  using Point_3 = typename Exact_kernel::Point_3;
  using Point_2 = typename Exact_kernel::Point_2;
  using Segment_3 = typename Exact_kernel::Segment_3;
  using Segment_2 = typename Exact_kernel::Segment_2;
  using Line_3 = typename Exact_kernel::Line_3;
  using Line_2 = typename Exact_kernel::Line_2;
  using Vector_3 = typename Exact_kernel::Vector_3;
  using Vector_2 = typename Exact_kernel::Vector_2;
  using Plane_3 =  typename Exact_kernel::Plane_3;
  using FT =  typename Exact_kernel::FT;
  using Result = Snapping_result<Exact_kernel>;

  CGAL_assertion (s0.source() != s0.target());
  CGAL_assertion (s1.source() != s1.target());

  Vector_3 v0 = s0.to_vector();
  Vector_3 v1 = s1.to_vector();
  Vector_3 normal = CGAL::cross_product (v0, v1);

  // Collinear segments
  if (normal == CGAL::NULL_VECTOR)
  {
    Line_3 l0 = s0.supporting_line();
    Line_3 l1 = s1.supporting_line();
    Vector_3 v0 = s0.to_vector();
    v0 = v0 / approximate_sqrt(v0 * v0);
    Vector_3 v1 = s1.to_vector();
    v1 = v1 / approximate_sqrt(v1 * v1);

    std::vector<std::pair<FT, FT>> points_along_lines;

    auto create_projection = [&](const Point_3& p)
                             {
                               Point_3 p0 = l0.projection(p);
                               Vector_3 vec0 (s0.source(), p0);
                               FT pos0 = vec0 * v0;

                               Point_3 p1 = l1.projection(p);
                               Vector_3 vec1 (s1.source(), p1);
                               FT pos1 = vec1 * v1;

                               points_along_lines.emplace_back(pos0, pos1);
                             };
    create_projection (s0.source());
    create_projection (s0.target());
    create_projection (s1.source());
    create_projection (s1.target());

    std::sort (points_along_lines.begin(), points_along_lines.end(),
               [](const auto& a, const auto& b) -> bool
               {
                 return a.first < b.first;
               });

    Segment_3 seg0 (s0.source() + points_along_lines[1].first * v0,
                    s0.source() + points_along_lines[2].first * v0);
    Segment_3 seg1 (s1.source() + points_along_lines[1].second * v1,
                    s1.source() + points_along_lines[2].second * v1);


    return Result(seg0, seg1);
  }

  // Coplanar segments
  if (coplanar (s0.source(), s0.target(), s1.source(), s1.target()))
  {
    Plane_3 plane
      (s0.source(), s0.target(),
       collinear(s0.source(), s0.target(), s1.target()) ? s1.source() : s1.target());
    Vector_3 base1 = plane.base1();
    base1 = base1 / approximate_sqrt(base1 * base1);
    Vector_3 base2 = plane.base2();
    base2 = base2 / approximate_sqrt(base2 * base2);

    auto to_2d = [&](const Point_3& p) -> Point_2
                 {
                   Vector_3 v (s0.source(), p);
                   return Point_2 (v * base1, v * base2);
                 };
    auto to_3d = [&](const Point_2& p) -> Point_3
                 {
                   return s0.source() + p.x() * base1 + p.y() * base2;
                 };

    auto coplanar_snapping
      = [&](const Segment_3& seg, const Segment_3& ref) -> std::pair<Point_3, typename Result::Type>
      {
        Segment_2 seg2 (to_2d (seg.source()), to_2d (seg.target()));
        Segment_2 ref2 (to_2d (ref.source()), to_2d (ref.target()));

        auto seg2_i_ref2 = intersection(seg2, ref2);
        if (seg2_i_ref2)
        {
          const Point_2* point = boost::get<Point_2>(&*seg2_i_ref2);
          CGAL_assertion (point != nullptr);
          return std::make_pair(to_3d(*point), Result::POINT);
        }

        Vector_2 vec2 = ref2.to_vector();
        Vector_2 ortho = vec2.perpendicular(CLOCKWISE);
        ortho = ortho / approximate_sqrt(ortho * ortho);

        Segment_2 sc (ref2.source() + ortho * tolerance,
                      ref2.target() + ortho * tolerance);
        Segment_2 scc (ref2.source() - ortho * tolerance,
                       ref2.target() - ortho * tolerance);

        Orientation source_sc = orientation (sc.source(), sc.target(), seg2.source());
        Orientation source_scc = orientation (scc.source(), scc.target(), seg2.source());
        Orientation target_sc = orientation (sc.source(), sc.target(), seg2.target());
        Orientation target_scc = orientation (scc.source(), scc.target(), seg2.target());

        bool source_in = (source_sc != source_scc);
        bool target_in = (target_sc != target_scc);

        // Whole segment is inside tolerance
        if (source_in && target_in)
          return std::make_pair(Point_3(), Result::SEGMENT);


        if (!source_in && !target_in)
        {
          Line_2 lc = sc.supporting_line();
          auto inter = intersection (seg2, lc);
          // Segment intersect both borders
          if (inter)
            return std::make_pair(Point_3(), Result::SEGMENT);

          // Segment does not intersect at all
          return std::make_pair(Point_3(), Result::EMPTY);
        }

        // One vertex only is in tolerance
        if (source_in)
          return std::make_pair(seg.source(), Result::POINT);
        // else
        return std::make_pair(seg.target(), Result::POINT);
      };

    auto crop0 = coplanar_snapping (s0, s1);
    auto crop1 = coplanar_snapping (s1, s0);

    if (crop0.second == Result::POINT && crop1.second == Result::POINT)
      return Result(crop0.first, crop1.first);

    if (crop0.second == Result::POINT && crop1.second == Result::SEGMENT)
      return Result (crop0.first, s1.supporting_line().projection(crop0.first));

    if (crop0.second == Result::SEGMENT && crop1.second == Result::POINT)
      return Result (crop1.first, s0.supporting_line().projection(crop1.first));

    if (crop0.second == Result::SEGMENT && crop1.second == Result::SEGMENT)
    {
      // Note: there might be very specific cases where these segments should be cropped
      return Result (s0, s1);
    }

    return Result();
  }

  // General case (segments are not coplanar/collinear)

  Plane_3 plane0 (s0.source(), normal);
  Plane_3 plane1 (s1.source(), normal);

  Segment_3 s0proj (plane1.projection (s0.source()),
                    plane1.projection (s0.target()));
  Segment_3 s1proj (plane0.projection (s1.source()),
                    plane0.projection (s1.target()));

  CGAL_assertion (s0proj.source() != s0proj.target());
  CGAL_assertion (s1proj.source() != s1proj.target());

  typename std::result_of<typename Exact_kernel::Intersect_3(Segment_3, Segment_3)>::type
    intersection0 = intersection(s0, s1proj);
  if (!intersection0)
    return Result();

  typename std::result_of<typename Exact_kernel::Intersect_3(Segment_3, Segment_3)>::type
    intersection1 = intersection(s1, s0proj);
  CGAL_assertion (bool(intersection1));

  if (const Point_3* p0 = boost::get<Point_3>(&*intersection0))
  {
    const Point_3* p1 = boost::get<Point_3>(&*intersection1);
    CGAL_assertion(p1 != nullptr);
    return Result (*p0, *p1);
  }

  return Result();
}

template <class Kernel,
          bool Has_exact_constructions=
          !std::is_floating_point<typename Kernel::FT>::value >
class Exact_snapping;

template <class Kernel>
class Exact_snapping<Kernel,false>
{
//typedefs
  using Exact_kernel = Exact_predicates_exact_constructions_kernel;
  using Exact_FT = typename Exact_kernel::FT;
  using Exact_segment_3 = typename Exact_kernel::Segment_3;
  using Exact_point_3 = typename Exact_kernel::Point_3;
  using Exact_result = Snapping_result<Exact_kernel>;

  using Inexact_kernel = Kernel;
  using Inexact_FT = typename Inexact_kernel::FT;
  using Inexact_segment_3 = typename Inexact_kernel::Segment_3;
  using Inexact_point_3 = typename Inexact_kernel::Point_3;
  using Inexact_result = Snapping_result<Inexact_kernel>;

  using Exact_to_inexact = Cartesian_converter<Exact_kernel, Inexact_kernel>;
  using Inexact_to_exact = Cartesian_converter<Inexact_kernel, Exact_kernel>;

  Exact_to_inexact to_inexact;
  Inexact_to_exact to_exact;

public:

  Inexact_result intersection (const Inexact_segment_3& s0, const Inexact_segment_3& s1,
                               const Inexact_FT& tolerance) const
  {
    Exact_result result = exact_snapping<Exact_kernel> (to_exact(s0), to_exact(s1), to_exact(tolerance));
    if (result.type == Exact_result::EMPTY)
      return Inexact_result();
    if (result.type == Exact_result::POINT)
      return Inexact_result (to_inexact (get<const Exact_point_3&>(result.first)),
                             to_inexact (get<const Exact_point_3&>(result.second)));
    // else Segment
    return Inexact_result (to_inexact (get<const Exact_segment_3&>(result.first)),
                           to_inexact (get<const Exact_segment_3&>(result.second)));
  }

};

template <class Kernel>
class Exact_snapping<Kernel,true>
{
  using FT = typename Kernel::FT;
  using Segment_3 = typename Kernel::Segment_3;
  using Point_3 = typename Kernel::Point_3;
  using Result = Snapping_result<Kernel>;

public:

  Result intersection (const Segment_3& s0, const Segment_3& s1, const FT& tolerance) const
  {
    return exact_snapping<Kernel>(s0,s1,tolerance);
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
         std::tie (new_point, okay) = exact_snapping.intersection (s1, s2, tolerance);

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
  using Result = internal::Snapping_result<Kernel>;

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

  std::map<Point_3, vertex_descriptor> map_p2v;

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
       Result snapping = exact_snapping.intersection (sa, sb, tolerance);
       if (snapping.type == Result::EMPTY)
         return;

       if (snapping.type == Result::POINT)
       {
         const Point_3& a = get<const Point_3&>(snapping.first);
         const Point_3& b = get<const Point_3&>(snapping.second);
         Point_3 new_point = midpoint(a,b);

         // Create new vertex
         vertex_descriptor vd = add_vertex(graph);
         put (point_map, vd, new_point);

         // Compute barycentric coord along each segment
         FT ba = approximate_sqrt (squared_distance (sa.source(), a) / sa.squared_length());
         FT bb = approximate_sqrt (squared_distance (sb.source(), b) / sb.squared_length());

         // Store relative position of new vertex on edges
         auto ma = new_vertices.insert (std::make_pair(ea, std::vector<Vertex_position>()));
         auto mb = new_vertices.insert (std::make_pair(eb, std::vector<Vertex_position>()));
         ma.first->second.emplace_back (ba, vd);
         mb.first->second.emplace_back (bb, vd);
       }
       else // snapping.type == Result::SEGMENT)
       {
         std::cerr << "Warning: coplanar segment snapping not handled" << std::endl;
       }
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


}


} } //end of namespace CGAL::Polygon_mesh_processing

/// \endcond

#endif // CGAL_POLYGON_MESH_PROCESSING_POLYLINE_SNAPPING_H
