// Copyright (c) 2023 University of Genova (Italy).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Claudio Mancinelli

#ifndef CGAL_POLYGON_MESH_PROCESSING_BSURF_LOCALLY_SHORTEST_PATH_H
#define CGAL_POLYGON_MESH_PROCESSING_BSURF_LOCALLY_SHORTEST_PATH_H

// #include <CGAL/license/Polygon_mesh_processing/bsurf.h>

#include <CGAL/Named_function_parameters.h>
#include <CGAL/Polygon_mesh_processing/locate.h>
#include <CGAL/boost/graph/Dual.h>
#include <CGAL/boost/graph/named_params_helper.h>

#include <CGAL/boost/graph/dijkstra_shortest_paths.h>
#include <boost/graph/graph_traits.hpp>

namespace CGAL {
namespace Polygon_mesh_processing {

template <class FT, class TriangleMesh, class EdgeLocationRange>
void locally_shortest_path(const Face_location<TriangleMesh, FT> &src,
                           const Face_location<TriangleMesh, FT> &tgt,
                           const TriangleMesh &tmesh,
                           EdgeLocationRange &edge_locations);

template <class TriangleMesh, class FT>
using Edge_location =
    std::pair<typename boost::graph_traits<TriangleMesh>::edge_descriptor,
              std::array<FT, 2>>;

template <class TriangleMesh, class FT>
using Bezier_segment = std::array<Face_location<TriangleMesh, FT>, 4>;


template <typename FT, typename TriangleMesh,
          typename NamedParameters = parameters::Default_named_parameters>
#ifdef DOXYGEN_RUNNING
Point construct_point(
    const Edge_location<TriangleMesh, FT> &loc,
#else
typename internal::Location_traits<TriangleMesh, NamedParameters>::Point
construct_point(const Edge_location<TriangleMesh, FT> &loc,
#endif
    const TriangleMesh &tm,
    const NamedParameters &np = parameters::default_values()) {
  typedef typename boost::graph_traits<TriangleMesh>::edge_descriptor
      edge_descriptor;
  typedef
      typename GetGeomTraits<TriangleMesh, NamedParameters>::type Geom_traits;

  typedef typename GetVertexPointMap<TriangleMesh, NamedParameters>::const_type
      VertexPointMap;
  typedef typename boost::property_traits<VertexPointMap>::value_type Point;
  typedef typename boost::property_traits<VertexPointMap>::reference
      Point_reference;

  using parameters::choose_parameter;
  using parameters::get_parameter;

  CGAL_precondition(CGAL::is_triangle_mesh(tm));

  VertexPointMap vpm = parameters::choose_parameter(
      parameters::get_parameter(np, internal_np::vertex_point),
      get_const_property_map(boost::vertex_point, tm));
  Geom_traits gt = choose_parameter<Geom_traits>(
      get_parameter(np, internal_np::geom_traits));

  edge_descriptor ed = loc.first;
  const Point_reference p0 = get(vpm, source(ed, tm));
  const Point_reference p1 = get(vpm, target(ed, tm));

  internal::Barycentric_point_constructor<Geom_traits, Point> bp_constructor;
  return bp_constructor(p0, loc.second[0], p1, loc.second[1], gt);
}

namespace internal {

template <class K, class TriangleMesh, class VertexPointMap>
struct Locally_shortest_path_imp
{
  using face_descriptor =
      typename boost::graph_traits<TriangleMesh>::face_descriptor;
  using vertex_descriptor =
      typename boost::graph_traits<TriangleMesh>::vertex_descriptor;
  using halfedge_descriptor =
      typename boost::graph_traits<TriangleMesh>::halfedge_descriptor;

  using Point_2 = typename K::Point_2;
  using Point_3 = typename K::Point_3;
  using Vector_2 = typename K::Vector_2;
  using Vector_3 = typename K::Vector_3;
  using FT = typename K::FT;

// #ifdef CGAL_DEBUG_BSURF
  static
  void dump_path(const std::vector<halfedge_descriptor>& path,
                 const std::vector<FT>& lerps,
                 const Face_location<TriangleMesh,FT>& src,
                 const Face_location<TriangleMesh,FT>& tgt,
                 const TriangleMesh& mesh)
  {
    static int i = -1;
    std::cout << "dump current path in path_"+std::to_string(i)+".polylines.txt\n";
    std::ofstream out("path_"+std::to_string(++i)+".polylines.txt");
    out << path.size()+2 << " " << construct_point(src, mesh);
    for(std::size_t i=0; i<path.size(); ++i)
      out << " " << construct_point(Edge_location<TriangleMesh, FT>(path[i], make_array(lerps[i], 1.-lerps[i])), mesh);
    out << " " << construct_point(tgt, mesh) << "\n";
  }
// #endif

  // TODO: recode using CGAL code?
  static
  Vector_2
  intersect_circles(const Vector_2 &c2, FT R2,
                    const Vector_2 &c1, FT R1)
  {
    auto R = (c2 - c1).squared_length();
    assert(R > 0);
    auto invR = FT(1) / R;
    Vector_2 result = c1+c2;

    result = result + (c2 - c1) * ((R1 - R2) * invR);
    auto A = 2 * (R1 + R2) * invR;
    auto B = (R1 - R2) * invR;
    auto s = A - B * B - 1;
    assert(s >= 0);
    result = result + Vector_2(c2.y() - c1.y(), c1.x() - c2.x()) * sqrt(s);
    return result / 2.;
  }

  static std::array<Vector_2, 3>
  init_flat_triangles(face_descriptor f,
                      const VertexPointMap &vpm, const TriangleMesh &mesh,
                      const Face_location<TriangleMesh, FT> &src)
  {
    halfedge_descriptor h = halfedge(f, mesh);
    std::array<vertex_descriptor, 3> triangle_vertices = make_array(
        source(h, mesh), target(h, mesh), target(next(h, mesh), mesh));

    std::array<Vector_2, 3> tr2d;
    tr2d[0] = Vector_2(0, 0);
    tr2d[1] =
        Vector_2(0, sqrt(squared_distance(get(vpm, triangle_vertices[0]),
                                          get(vpm, triangle_vertices[1]))));
    auto rx = squared_distance(get(vpm, triangle_vertices[0]),
                               get(vpm, triangle_vertices[2]));
    auto ry = squared_distance(get(vpm, triangle_vertices[1]),
                               get(vpm, triangle_vertices[2]));
    tr2d[2] = intersect_circles(tr2d[0], rx, tr2d[1], ry);

    return tr2d;
  }

  static std::array<Vector_2, 2>
  init_source_triangle(halfedge_descriptor hopp,
                       const VertexPointMap &vpm,
                       const TriangleMesh &mesh,
                       Face_location<TriangleMesh, FT> src)
  {
    halfedge_descriptor h = opposite(hopp, mesh);
    std::array<vertex_descriptor, 3> triangle_vertices = make_array(
        source(h, mesh), target(h, mesh), target(next(h, mesh), mesh));

    std::array<Vector_2, 3> tr2d;
    tr2d[0] = Vector_2(0, 0);
    tr2d[1] =
        Vector_2(0, sqrt(squared_distance(get(vpm, triangle_vertices[0]),
                                          get(vpm, triangle_vertices[1]))));
    FT rx = squared_distance(get(vpm, triangle_vertices[0]),
                               get(vpm, triangle_vertices[2]));
    FT ry = squared_distance(get(vpm, triangle_vertices[1]),
                               get(vpm, triangle_vertices[2]));
    tr2d[2] = intersect_circles(tr2d[0], rx, tr2d[1], ry);


    halfedge_descriptor href = halfedge(src.first, mesh);
    if (href!=h)
    {
      if (href==next(h, mesh))
      {
        std::array<FT, 3> tmp = CGAL::make_array(src.second[2], src.second[0], src.second[1]);
        src.second = tmp;
      }
      else
      {
        CGAL_assertion(next(href, mesh)==h);
        std::array<FT, 3> tmp = CGAL::make_array(src.second[1], src.second[2], src.second[0]);
        src.second = tmp;
      }
    }

    auto point_coords = tr2d[0] * src.second[0] + tr2d[1] * src.second[1] +
                        tr2d[2] * src.second[2];
#ifdef CGAL_DEBUG_BSURF
    std::cout << "4 " << tr2d[0] - point_coords << " 0 "
              << tr2d[1] - point_coords << " 0 "
              << tr2d[2] - point_coords << " 0 "
              << tr2d[0] - point_coords << " 0\n";
#endif
    return make_array(tr2d[0] - point_coords,
                      tr2d[1] - point_coords);
  }

  static std::array<Vector_2, 2>
  init_target_triangle(halfedge_descriptor h,
                       const std::array<Vector_2, 2>& flat_tid,
                       const VertexPointMap &vpm, const TriangleMesh &mesh,
                       Face_location<TriangleMesh, FT> tgt)
  {
    std::array<vertex_descriptor, 3> triangle_vertices = make_array(
        source(h, mesh), target(h, mesh), target(next(h, mesh), mesh));

    std::array<Vector_2, 3> tr2d;
    tr2d[0] = flat_tid[1];
    tr2d[1] = flat_tid[0];
    FT rx = squared_distance(get(vpm, triangle_vertices[0]),
                             get(vpm, triangle_vertices[2]));
    FT ry = squared_distance(get(vpm, triangle_vertices[1]),
                             get(vpm, triangle_vertices[2]));
    tr2d[2] = intersect_circles(tr2d[0], rx, tr2d[1], ry);


    halfedge_descriptor href = halfedge(tgt.first, mesh);
    if (href!=h)
    {
      if (href==next(h, mesh))
      {
        std::array<FT, 3> tmp = CGAL::make_array(tgt.second[2], tgt.second[0], tgt.second[1]);
        tgt.second = tmp;
      }
      else
      {
        CGAL_assertion(next(href, mesh)==h);
        std::array<FT, 3> tmp = CGAL::make_array(tgt.second[1], tgt.second[2], tgt.second[0]);
        tgt.second = tmp;
      }
    }

    auto point_coords = tr2d[0] * tgt.second[0] + tr2d[1] * tgt.second[1] +
                        tr2d[2] * tgt.second[2];

#ifdef CGAL_DEBUG_BSURF
    std::cout << "4 " << tr2d[0]  << " 0 "
              << tr2d[1]  << " 0 "
              << tr2d[2] << " 0 "
              << tr2d[0]  << " 0\n";
#endif

    return make_array(point_coords,
                      point_coords);
  }

  // static
  // std::array<Vector_2, 3>
  // unfold_face(halfedge_descriptor h,
  //             const VertexPointMap &vpm, const TriangleMesh &mesh,
  //             const std::array<Vector_2, 3>& flat_tid)
  // {
  //   halfedge_descriptor h_opp = opposite(h_curr, mesh);



  //   vertex_descriptor v = target(next(h_curr,mesh),mesh);
  //   vertex_descriptor a = target(h_curr,mesh);
  //   vertex_descriptor b = source(h_curr, mesh);
  //   FT r0 = squared_distance(get(vpm,v), get(vpm,a));
  //   FT r1 = squared_distance(get(vpm,v), get(vpm,b));

  //   Vector_2 v2 = intersect_circles(flat_tid[0], r0, flat_tid[1], r1);


  //   std::array<Vector_2, 2> res;
  //   if(next(h_curr, mesh) == h_next_opp)
  //   {
  //      res[0]=flat_tid[0];
  //      //res[2]=flat_tid[1];
  //      res[1]=v2;
  //   }
  //   else
  //   {
  //     CGAL_assertion(prev(h_curr, mesh) == h_next_opp);
  //     res[0]=v2;
  //     res[1]=flat_tid[1];
  //     //res[2]=flat_tid[0];
  //   }

  //   return res;
  // }
static
  std::array<Vector_2, 2>
  unfold_face(halfedge_descriptor h_curr, halfedge_descriptor h_next,
              const VertexPointMap &vpm, const TriangleMesh &mesh,
              const std::array<Vector_2, 2>& flat_tid)
  {
    halfedge_descriptor h_next_opp = opposite(h_next, mesh);
    CGAL_assertion(face( h_curr, mesh) == face(h_next_opp, mesh));


    vertex_descriptor v = target(next(h_curr,mesh),mesh);
    vertex_descriptor a = target(h_curr,mesh);
    vertex_descriptor b = source(h_curr, mesh);
    FT r0 = squared_distance(get(vpm,v), get(vpm,a));
    FT r1 = squared_distance(get(vpm,v), get(vpm,b));

    Vector_2 v2 = intersect_circles(flat_tid[1], r1, flat_tid[0], r0);


    std::array<Vector_2, 2> res;
    if(next(h_curr, mesh) == h_next_opp)
    {
       res[0]=flat_tid[0];
       //res[2]=flat_tid[1];
       res[1]=v2;
#ifdef CGAL_DEBUG_BSURF
       std::cout << "4 " << res[0]  << " 0 "
                 << res[1]  << " 0 "
                 << flat_tid[1] << " 0 "
                 << res[0]  << " 0\n";
#endif
    }
    else
    {
      CGAL_assertion(prev(h_curr, mesh) == h_next_opp);
      res[0]=v2;
      res[1]=flat_tid[1];
      //res[2]=flat_tid[0];
#ifdef CGAL_DEBUG_BSURF
      std::cout << "4 " << res[0]  << " 0 "
                << res[1]  << " 0 "
                << flat_tid[0] << " 0 "
                << res[0]  << " 0\n";
#endif
    }

    return res;
  }
  static
  std::vector< std::array<Vector_2, 2>>
  unfold_strip(const std::vector<halfedge_descriptor>& initial_path,
              const Face_location<TriangleMesh, FT>& src,
              const Face_location<TriangleMesh, FT>& tgt,
              const VertexPointMap &vpm, const TriangleMesh &mesh)
  {
    std::size_t s=initial_path.size();
    std::vector<std::array<Vector_2, 2>> result(s+1);
    result[0]=init_source_triangle(initial_path[0], vpm, mesh, src);
    for(std::size_t i=1;i<s;++i)
      result[i]=unfold_face(initial_path[i-1],initial_path[i], vpm, mesh, result[i-1]);

    result[s]=init_target_triangle(initial_path.back(), result[s-1], vpm, mesh, tgt);

    return result;
  }
  struct funnel_point {
    int   face = 0;
    Vector_2 pos;
  };
  //TODO: reimplement using CGAL
  static
  FT intersect_segments(const Vector_2 &start1, const Vector_2 &end1,
                        const Vector_2 &start2, const Vector_2 &end2)
  {
    if (end1 == start2) return 0;
    if (end2 == start1) return 1;
    if (start2 == start1) return 0;
    if (end2 == end1) return 1;
    auto a   = end1 - start1;    // direction of line a
    auto b   = start2 - end2;    // direction of line b, reversed
    auto d   = start2 - start1;  // right-hand side
    auto det = a.x() * b.y() - a.y() * b.x();
    assert(det);
    return (a.x() * d.y() - a.y() * d.x()) / det;
  }
  static
  int max_curvature_point(const std::vector<funnel_point> &path) {
    // Among vertices around which the path curves, find the vertex
    // with maximum angle. We are going to fix that vertex. Actually, max_index is
    // the index of the first face containing that vertex.
    std::size_t max_index = -1;
    FT max_angle = 0.;
    for (std::size_t i = 1; i < path.size()-1; ++i) {
      Vector_2 pos   = path[i].pos;
      Vector_2 prev  = path[i - 1].pos;
      Vector_2 next  = path[i + 1].pos;
      Vector_2 v0    = pos - prev;
      v0 = v0 / sqrt(v0.squared_length());
      Vector_2 v1    = next - pos;
      v1 = v1 / sqrt(v1.squared_length());
      FT angle = 1 - scalar_product(v0, v1);
      if (angle > max_angle) {
        max_index = path[i].face;
        max_angle = angle;
      }
    }

std::cout << "funnels ("<< max_index << ")";
for (auto f : path)
  std::cout << " " << f.pos << " |";
std::cout << "\n";

    return max_index;
  }

  static
  std::vector<FT>
  funnel(const std::vector< std::array<Vector_2, 2>>& portals, std::size_t& max_index)
  {
    // Find straight path.
    Vector_2 start(NULL_VECTOR);
    int apex_index  = 0;
    int left_index  = 0;
    int right_index = 0;
    Vector_2 apex        = start;
    Vector_2 left_bound  = portals[0][0];
    Vector_2 right_bound = portals[0][1];

    // Add start point.
    std::vector<funnel_point> points = std::vector<funnel_point>{{apex_index, apex}};
    points.reserve(portals.size());
    // @Speed: is this slower than an inlined function?
    auto area = [](const Vector_2 a, const Vector_2 b, const Vector_2 c) {
      return determinant(b - a, c - a);
    };

    for (std::size_t i = 0; i < portals.size(); ++i) 
    {
      auto left = portals[i][0], right = portals[i][1];
      // Update right vertex.
      if (area(apex, right_bound, right) <= 0) {
        if (apex == right_bound || area(apex, left_bound, right) > 0) {
          // Tighten the funnel.
          right_bound = right;
          right_index = i;
        } else {
          // Right over left, insert left to path and restart scan from
          // portal left point.
          if (left_bound != apex) {
            points.push_back({left_index, left_bound});
            // Make current left the new apex.
            apex       = left_bound;
            apex_index = left_index;
            // Reset portal
            left_bound  = apex;
            right_bound = apex;
            left_index  = apex_index;
            right_index = apex_index;
            // Restart scan
            i = apex_index;
            continue;
          }
        }
      }

      // Update left vertex.
      if (area(apex, left_bound, left) >= 0) {
        if (apex == left_bound || area(apex, right_bound, left) < 0) {
          // Tighten the funnel.
          left_bound = left;
          left_index = i;
        } else {
          if (right_bound != apex) {
            points.push_back({right_index, right_bound});
            // Make current right the new apex.
            apex       = right_bound;
            apex_index = right_index;
            // Reset portal
            left_bound  = apex;
            right_bound = apex;
            left_index  = apex_index;
            right_index = apex_index;
            // Restart scan
            i = apex_index;
            continue;
          }
        }
      }
    }

    // This happens when we got an apex on the last edge of the strip
    if (points.back().pos != portals.back()[0]) {
      points.push_back({(int)portals.size() - 1, portals.back()[0]});
    }
    assert(points.back().pos == portals.back()[0]);
    assert(points.back().pos == portals.back()[1]);

    std::vector<double> lerps;
    lerps.reserve(portals.size());
    for (std::size_t i = 0; i < points.size() - 1; i++) {
      auto a = points[i].pos;
      auto b = points[i + 1].pos;
      for (auto k = points[i].face; k < points[i + 1].face; ++k) {
        auto portal = portals[k];
        //      assert(cross(b - a, portal.second - portal.first) > 0);

std::cout << "i=" << i << "\n";
std::cout << "a=" << a << " b=" << b << " portal[0]=" << portal[0] << " portal[1]=" << portal[1] << "\n";

        FT s = intersect_segments(a, b, portal[0], portal[1]);

std::cout << "s=" << s << "\n";

        lerps.push_back(std::clamp(s, 0.0, 1.0));
      }
    }

    auto index = 1;
    std::cout << "setting funnel_point indices\n";
    for (std::size_t i = 0; i < portals.size(); ++i) 
    {
      std::cout << "  i=" << i << " index = " << index << "\n";
      std::cout << "  portals[i][0]=" << portals[i][0] << " portals[i][1]=" << portals[i][1] << "\n";
      std::cout << "  points[index].pos = " <<points[index].pos << "\n";
      if ((portals[i][0] == points[index].pos) ||
          (portals[i][1] == points[index].pos)) 
      {
        std::cout << "  setting point["<<index<<"].face="<< i << "\n"; 
        points[index].face = i;
        index += 1;
      }
    }
    max_index = max_curvature_point(points);
    // assert(lerps.size() == portals.size() - 1);
    return lerps;
  }
  static
  void straighten_path(std::vector< std::array<Vector_2, 2>>& portals,
                       std::vector<FT>& lerps,
                       std::vector<halfedge_descriptor>& path,
                       const Face_location<TriangleMesh, FT>& src,
                       const Face_location<TriangleMesh, FT>& tgt,
                       const VertexPointMap &vpm, const TriangleMesh &mesh, std::size_t index)
  {
// #ifdef CGAL_DEBUG_BSURF
    dump_path(path, lerps, src, tgt, mesh);
// #endif

    vertex_descriptor vertex=boost::graph_traits<TriangleMesh>::null_vertex();

    // TODO: use a while loop breaking when no apex vertices not already visited are available
    for (std::size_t i = 0; i < portals.size() * 2 && index != std::size_t(-1); i++)
    {
// #ifdef CGAL_DEBUG_BSURF
      std::cout << "Improving path " << path.size() << " hedges\n";
      std::cout << "src = " << construct_point(src, mesh) << "\n";
      std::cout << "tgt = " << construct_point(tgt, mesh) << "\n";
// #endif

      vertex_descriptor new_vertex=boost::graph_traits<TriangleMesh>::null_vertex();
      halfedge_descriptor h_curr       = path[index];
      halfedge_descriptor h_next       = path[index + 1];
      bool is_target = false;
      if (lerps[index] == 0) {
        new_vertex = target(h_curr,mesh);
        is_target = true;
      } else if (lerps[index] == 1) {
        new_vertex = source(h_curr,mesh);
      }
      if (new_vertex == vertex) break;
      vertex = new_vertex;

// #ifdef CGAL_DEBUG_BSURF
      std::cout << "  Current strip with Apex: " << get(vpm, new_vertex) <<  "\n";
      for (auto h : path)
      {
        std::cout << "  4 " << get(vpm, source(h, mesh))
                  << "  " << get(vpm, target(h, mesh))
                  << "  " << get(vpm, target(next(h, mesh), mesh))
                  << "  " << get(vpm, source(h, mesh)) << std::endl;

      }
// #endif

      // if I hit the source vertex v of h_curr, then h_next has v as source, thus we turn ccw around v in path
      // Similarly, if I hit the target vertex v of h_curr, then  h_next has v as target, thus we turn cw around v in path

      if ( !(!is_target || opposite(next(h_curr, mesh), mesh)==h_next) )
        std::cout <<  edge(h_curr, mesh) << " |  " << edge(opposite(next(h_curr, mesh), mesh), mesh) << " vs " << edge(h_next, mesh) << "\n";

      CGAL_assertion(!is_target || opposite(next(h_curr, mesh), mesh)==h_next);
      CGAL_assertion(is_target || opposite(prev(h_curr, mesh), mesh)==h_next);

      std::size_t curr_index = index+1;
      std::vector<halfedge_descriptor> new_hedges;
      if (is_target)
      {
        face_descriptor target_face;

        while (target(path[curr_index], mesh) == new_vertex)
        {
          if(curr_index==path.size()-1)
          {
            target_face=tgt.first;
            curr_index=path.size();
            break;
          }
          ++curr_index;
        }
        if (curr_index != path.size())
          target_face = face(opposite(path[curr_index], mesh), mesh);

        halfedge_descriptor h_loop=opposite(prev(opposite(h_curr, mesh), mesh), mesh);
        do {
          new_hedges.push_back(h_loop);
          h_loop=opposite(prev(h_loop,mesh), mesh);
        }
        while(target_face!=face(h_loop, mesh));
        new_hedges.push_back(h_loop);
      }
      else
      {
        face_descriptor target_face;

std::cout << "index = " << index << "\n";
std::cout << "path.size() = " << path.size() << "\n";

        while (source(path[curr_index], mesh) == new_vertex)
        {
          if(curr_index==path.size()-1)
          {
            target_face=tgt.first;
            curr_index=path.size();
            break;
          }
          ++curr_index;
        }
        if (curr_index != path.size())
          target_face=face(opposite(path[curr_index], mesh), mesh);

        halfedge_descriptor h_loop=opposite(next(opposite(h_curr, mesh), mesh), mesh); // skip the face before h_curr (as we won't remove it from path)
        do {
          new_hedges.push_back(h_loop);
          h_loop=opposite(next(h_loop,mesh), mesh);
        }
        while(target_face!=face(h_loop, mesh));
        new_hedges.push_back(h_loop);
      }

      // replace the halfedges incident to the apex vertex with the opposite part of the ring
      std::vector<halfedge_descriptor> new_path(path.begin(),path.begin()+index);
      new_path.insert(new_path.end(), new_hedges.begin(), new_hedges.end());
      new_path.insert(new_path.end(), path.begin()+curr_index, path.end());
      path.swap(new_path);

      portals=unfold_strip(path,src,tgt,vpm,mesh);
      lerps=funnel(portals,index);
// #ifdef CGAL_DEBUG_BSURF
      dump_path(path, lerps, src, tgt, mesh);
// #endif
    }

#ifdef CGAL_DEBUG_BSURF
    std::cout << "  Final strip\n";
    for (auto h : path)
    {
      std::cout << "  4 " << get(vpm, source(h, mesh))
                << "  " << get(vpm, target(h, mesh))
                << "  " << get(vpm, target(next(h, mesh), mesh))
                << "  " << get(vpm, source(h, mesh)) << "\n";

    }
#endif

  }

// TODO:  starting from here, we can move that to a de Casteljau Impl struct

  template <class EdgeLocationRange>
  static
  std::vector<Point_3>
  get_positions(const EdgeLocationRange& edge_locations,
                const TriangleMesh& mesh,
                const Face_location<TriangleMesh, FT>& src,
                const Face_location<TriangleMesh, FT>& tgt)
  {
    std::vector<Point_3> result;
    result.reserve(edge_locations.size()+2);
    result.push_back(construct_point(src,mesh));
    for(auto& e: edge_locations)
        result.push_back(construct_point(e,mesh));

    result.push_back(construct_point(tgt,mesh));
//TODO: we must guarantee that result is sorted and unique (rounding issue?)
    return result;
  }

  template <class EdgeLocationRange>
  static
  std::vector<FT>
  path_parameters(const EdgeLocationRange& edge_locations,
                  const TriangleMesh& mesh,
                  const Face_location<TriangleMesh, FT>& src,
                  const Face_location<TriangleMesh, FT>& tgt)
  {
    std::vector<Point_3> pos=get_positions(edge_locations,mesh,src,tgt);
    FT L=0.;
    std::vector<FT> result(pos.size());
    for(std::size_t i=0;i<pos.size();++i)
    {
      if(i) L+=sqrt(squared_distance(pos[i],pos[i-1]));
      result[i]=L;
    }

    for(auto& t:result) t/=L;

    return result;
  }

  template <class EdgeLocationRange>
  static
  Face_location<TriangleMesh, FT>
  eval_point_on_geodesic(const EdgeLocationRange& edge_locations,
                         const TriangleMesh& mesh,
                         const Face_location<TriangleMesh, FT>& src,
                         const Face_location<TriangleMesh, FT>& tgt,
                         const std::vector<FT>& parameters,/// edge length parameterization of the path from src to tgt through edge_locations
                         const FT& t)
  {
    if (t==0) return src;
    if (t==1) return tgt;

    if(src.first==tgt.first)
    {
      std::array<FT,3> bary;
      bary[0]=(1-t)*src.second[0]+t*tgt.second[0];
      bary[1]=(1-t)*src.second[1]+t*tgt.second[1];
      bary[2]=(1-t)*src.second[2]+t*tgt.second[2];
      return {src.first,bary};
    }

    std::size_t i = 0;
    for (; i < parameters.size() - 1; i++)
    {
      if (parameters[i + 1] >= t) break;
    }
    FT t_low = parameters[i];
    FT t_high = parameters[i + 1];
    CGAL_assertion(t_high!=t_low);
    FT alpha = (t - t_low) / (t_high - t_low);
    std::array<FT,3> bary_low;
    std::array<FT,3> bary_high;

    // warning there is an offset of the index: parameters contains one extra element (src) at 0
    // while edge_locations does not
    face_descriptor curr_tid = i==0?src.first:face(halfedge(edge_locations[i-1].first,mesh),mesh);
    halfedge_descriptor h_face = halfedge(curr_tid, mesh);
    auto edge_barycentric_coordinate =
      [&mesh, h_face](halfedge_descriptor h_edge,
                     const std::array<FT,2>& bary_edge)
    {
      std::array<FT,3> bary_edge_in_face;
      if (h_face!=h_edge)
      {
        if (h_face==next(h_edge, mesh))
        {
          bary_edge_in_face[0]=bary_edge[1];
          bary_edge_in_face[1]=0;
          bary_edge_in_face[2]=bary_edge[0];
        }
        else
        {
          bary_edge_in_face[0]=0;
          bary_edge_in_face[1]=bary_edge[0];
          bary_edge_in_face[2]=bary_edge[1];
        }
      }
      else
      {
        bary_edge_in_face[0]=bary_edge[0];
        bary_edge_in_face[1]=bary_edge[1];
        bary_edge_in_face[2]=0;
      }

      return bary_edge_in_face;
    };

    if(i==0)
      bary_low=src.second;
    else
    {
      halfedge_descriptor h_low = halfedge(edge_locations[i-1].first, mesh);
      bary_low = edge_barycentric_coordinate(h_low, edge_locations[i-1].second);
    }

    if(i==parameters.size()-2)
      bary_high=tgt.second;
    else
    {
      halfedge_descriptor h_high = opposite(halfedge(edge_locations[i].first, mesh), mesh);
      CGAL_assertion(face(h_high,mesh)==curr_tid);
      std::array<FT,2> edge_bary_high=edge_locations[i].second;
      std::swap(edge_bary_high[0],edge_bary_high[1]);
      bary_high = edge_barycentric_coordinate(h_high, edge_bary_high);
    }

    std::array<FT,3> bary;
    bary[0]=(1-alpha)*bary_low[0]+alpha*bary_high[0];
    bary[1]=(1-alpha)*bary_low[1]+alpha*bary_high[1];
    bary[2]=(1-alpha)*bary_low[2]+alpha*bary_high[2];

    return {curr_tid,bary};
  }

  static
  Face_location<TriangleMesh, FT>
  geodesic_lerp(const TriangleMesh &mesh,
                const Face_location<TriangleMesh, FT>& src,
                const Face_location<TriangleMesh, FT>& tgt,const FT& t)
  {
    std::vector<Edge_location<TriangleMesh, FT>> edge_locations;
    locally_shortest_path<FT>(src,tgt,mesh, edge_locations);
    std::vector<FT> parameters=path_parameters(edge_locations,mesh,src,tgt);
    Face_location<TriangleMesh, FT> point =
      eval_point_on_geodesic(edge_locations,mesh,src,tgt,parameters,t);
    return point;
  }


  static
  std::pair<Bezier_segment<TriangleMesh, FT>,Bezier_segment<TriangleMesh,FT>>
  subdivide_bezier_polygon(const TriangleMesh& mesh,
                           const Bezier_segment<TriangleMesh,FT>& polygon,
                           const FT& t)
  {
     Face_location<TriangleMesh, FT> Q0 = geodesic_lerp(mesh, polygon[0], polygon[1], t);
     Face_location<TriangleMesh, FT> Q1 = geodesic_lerp(mesh, polygon[1], polygon[2], t);
     Face_location<TriangleMesh, FT> Q2 = geodesic_lerp(mesh, polygon[2], polygon[3], t);
     Face_location<TriangleMesh, FT> R0 = geodesic_lerp(mesh, Q0, Q1, t);
     Face_location<TriangleMesh, FT> R1 = geodesic_lerp(mesh, Q1, Q2, t);
     Face_location<TriangleMesh, FT> S  = geodesic_lerp(mesh, R0, R1, t);

    return {{polygon[0], Q0, R0, S}, {S, R1, Q2, polygon[3]}};
  }
};

} // namespace internal

template <class FT, class TriangleMesh, class EdgeLocationRange>
void locally_shortest_path(const Face_location<TriangleMesh, FT> &src,
                           const Face_location<TriangleMesh, FT> &tgt,
                           const TriangleMesh &tmesh,
                           EdgeLocationRange &edge_locations)
{
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor
      face_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor
      halfedge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::edge_descriptor
      edge_descriptor;

  //TODO replace with named parameter
  using VPM = typename boost::property_map<TriangleMesh, CGAL::vertex_point_t>::const_type;
  using K =  typename Kernel_traits<typename boost::property_traits<VPM>::value_type>::type;
  using Impl = internal::Locally_shortest_path_imp<K, TriangleMesh, VPM>;
  VPM vpm = get(CGAL::vertex_point, tmesh);

  typename boost::property_map<
      TriangleMesh, CGAL::dynamic_face_property_t<face_descriptor>>::const_type
      predecessor_map =
          get(CGAL::dynamic_face_property_t<face_descriptor>(), tmesh);
  typename boost::property_map<TriangleMesh,
                               CGAL::dynamic_face_property_t<FT>>::const_type
      distance_map = get(CGAL::dynamic_face_property_t<FT>(), tmesh);
  typename boost::property_map<
      TriangleMesh, CGAL::dynamic_edge_property_t<FT>>::const_type weight_map =
      get(CGAL::dynamic_edge_property_t<FT>(), tmesh);


//TODO: handle boundary edges
  Dual dual(tmesh);

  // TODO: fill the weight map using something better than euclidean distance
  for (edge_descriptor ed : edges(tmesh))
  {
    halfedge_descriptor h=halfedge(ed, tmesh), hopp=opposite(h, tmesh);
    put(weight_map, ed,
        sqrt(squared_distance(
              centroid(get(vpm, source(h, tmesh)), get(vpm, target(h, tmesh)), get(vpm, target(next(h, tmesh), tmesh))),
              centroid(get(vpm, source(hopp, tmesh)), get(vpm, target(hopp, tmesh)), get(vpm, target(next(hopp, tmesh), tmesh)))
              )));
  }

  // TODO try stopping dijkstra as soon tgt is out of the queue.
  boost::dijkstra_shortest_paths(dual, src.first,
                                 boost::distance_map(distance_map)
                                     .predecessor_map(predecessor_map)
                                     .weight_map(weight_map));

  std::vector<halfedge_descriptor> initial_path;

  auto common_halfedge = [](face_descriptor f1, face_descriptor f2,
                            const TriangleMesh &tmesh) {
    halfedge_descriptor h = halfedge(f1, tmesh);
    for (int i = 0; i < 3; ++i) {
      if (face(opposite(h, tmesh), tmesh) == f2)
        return h;
      h = next(h, tmesh);
    }
    CGAL_assertion(!"faces do no share a common edge");
    return halfedge_descriptor();
  };

  face_descriptor current_face = tgt.first;
  while (true) {
    face_descriptor prev = get(predecessor_map, current_face);
    halfedge_descriptor h = common_halfedge(current_face, prev, tmesh);
    initial_path.push_back(h);
    if (prev == src.first)
      break;
    current_face = prev;
  }
  std::reverse(initial_path.begin(), initial_path.end());

  std::vector< std::array<typename K::Vector_2, 2>> portals=Impl::unfold_strip(initial_path,src,tgt,vpm,tmesh);
  std::cout << "portals.size() " << portals.size() << "\n";
  std::cout << "initial_path.size() " << initial_path.size() << "\n";
  std::size_t max_index=0;
  std::vector<FT> lerps=Impl::funnel(portals,max_index);
  Impl::straighten_path(portals,lerps,initial_path,src,tgt,vpm,tmesh,max_index);
  CGAL_assertion(lerps.size()==initial_path.size());

  //TODO: tmp for testing
  edge_locations.reserve(initial_path.size());
  for(std::size_t i=0; i<initial_path.size(); ++i)
  {
    edge_locations.emplace_back(initial_path[i], make_array(lerps[i], 1.-lerps[i]));
  }
}

template <class TriangleMesh, class FT>
std::vector<Face_location<TriangleMesh, FT>>
recursive_de_Casteljau(const TriangleMesh &mesh,
                       const Bezier_segment<TriangleMesh, FT>& control_points,
                       const int num_subdiv)
{
  //TODO replace with named parameter
  using VPM = typename boost::property_map<TriangleMesh, CGAL::vertex_point_t>::const_type;
  using K =  typename Kernel_traits<typename boost::property_traits<VPM>::value_type>::type;
  using Impl = internal::Locally_shortest_path_imp<K, TriangleMesh, VPM>;

  std::vector<Bezier_segment<TriangleMesh, FT>> segments(1,control_points);
  std::vector<Bezier_segment<TriangleMesh, FT>> result;
  for (auto subdivision = 0; subdivision < num_subdiv; subdivision++)
  {
    result.clear();
    result.reserve(segments.size() * 2);
    for (std::size_t i = 0; i < segments.size(); ++i)
    {
      auto [split0, split1] = Impl::subdivide_bezier_polygon(mesh, segments[i], 0.5);
      result.push_back(split0);
      result.push_back(split1);
    }
    std::swap(segments, result);
  }

  // nasty trick to build the vector from a pair of iterators
  // using the fact that data in array and vector are contiguous
  return {(Face_location<TriangleMesh, FT>*)segments.data(),
          (Face_location<TriangleMesh, FT>*)segments.data() + segments.size() * 4};
}

} // namespace Polygon_mesh_processing
} // namespace CGAL

#endif
