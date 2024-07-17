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

#include <boost/graph/dijkstra_shortest_paths.hpp>
//TODO: split to avoid redundant linking
#include <Eigen/Dense>
#include <boost/graph/graph_traits.hpp>

#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Projection_traits_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>

#ifdef CGAL_DEBUG_BSURF
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#endif

namespace CGAL {
namespace Polygon_mesh_processing {

#ifndef CGAL_BSURF_USE_DIJKSTRA_SP
template <class FT>
struct Dual_geodesic_solver;
#endif

template <class FT, class TriangleMesh, class EdgeLocationRange>
void locally_shortest_path(Face_location<TriangleMesh, FT> src,
                           Face_location<TriangleMesh, FT> tgt,
                           const TriangleMesh &tmesh,
                           EdgeLocationRange &edge_locations
#ifndef CGAL_BSURF_USE_DIJKSTRA_SP
                           , const Dual_geodesic_solver<FT>& solver = Dual_geodesic_solver<FT>()
#endif
);


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

#ifdef CGAL_DEBUG_BSURF
  static
  void dump_path(const std::vector<halfedge_descriptor>& path,
                 const std::vector<FT>& lerps,
                 const Face_location<TriangleMesh,FT>& src,
                 const Face_location<TriangleMesh,FT>& tgt,
                 const TriangleMesh& mesh)
  {
    static int i = -1;
    std::cout << "dump current path in path_"+std::to_string(i+1)+".polylines.txt\n";
    std::ofstream out("path_"+std::to_string(++i)+".polylines.txt");
    out << path.size()+2 << " " << construct_point(src, mesh);
    for(std::size_t i=0; i<path.size(); ++i)
      out << " " << construct_point(Edge_location<TriangleMesh, FT>(path[i], make_array(lerps[i], 1.-lerps[i])), mesh);
    out << " " << construct_point(tgt, mesh) << "\n";
  }
#endif

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
  init_flat_triangle( const halfedge_descriptor& h,
                      const VertexPointMap &vpm, const TriangleMesh &mesh)
  {
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

  //case k=0 assume that flat_tid has been flattened putting x-axis aligned with h_edge
  static
  std::array<Vector_2, 3>
  unfold_face(const halfedge_descriptor& h_edge,
              const VertexPointMap &vpm, const TriangleMesh &mesh,
              const std::array<Vector_2, 3>& flat_tid,const int k=0)
  {
    halfedge_descriptor h_opp = opposite(h_edge, mesh);


    vertex_descriptor v = target(next(h_opp,mesh),mesh);
    vertex_descriptor a = source(h_edge, mesh);
    vertex_descriptor b = target(h_edge,mesh);
    FT r0 = squared_distance(get(vpm,v), get(vpm,a));
    FT r1 = squared_distance(get(vpm,v), get(vpm,b));

    Vector_2 v2 = intersect_circles(flat_tid[(k+1)%3], r1, flat_tid[k], r0);

    halfedge_descriptor h_ref_opp = halfedge(face(h_opp,mesh),mesh);

    std::array<Vector_2, 3> res;

    if(h_ref_opp==h_opp)
    {
      res[0]=flat_tid[(k+1)%3];
      res[1]=flat_tid[k];
      res[2]=v2;
    } else if(next(h_ref_opp,mesh)==h_opp)
    {
       res[0]=v2;
       res[1]=flat_tid[(k+1)%3];
       res[2]=flat_tid[k];
    }else
    {
       assert(prev(h_ref_opp,mesh)==h_opp);
       res[0]=flat_tid[k];
       res[1]=v2;
       res[2]=flat_tid[(k+1)%3];
    }



    return res;
  }

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
    std::vector<std::array<Vector_2, 2>> result(s+1); // Vector_2 should be Point_2 as they are funnel endpoints in 2D (after flattening)
    result[0]=init_source_triangle(initial_path[0], vpm, mesh, src);
#ifdef CGAL_DEBUG_BSURF
    std::cout << "unfolding faces\n";
#endif
    for(std::size_t i=1;i<s;++i)
      result[i]=unfold_face(initial_path[i-1],initial_path[i], vpm, mesh, result[i-1]);

    result[s]=init_target_triangle(initial_path.back(), result[s-1], vpm, mesh, tgt);
#ifdef CGAL_DEBUG_BSURF
    std::cout << "done\n";
#endif

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
    if (start2 == start1)
    {
      if (end2 == end1) return 1; // the portail and path segment coincide
      return 0;
    }
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

#ifdef CGAL_DEBUG_BSURF
    std::cout << "funnels ("<< max_index << ")";
    for (auto f : path)
      std::cout << " " << f.pos << " |";
    std::cout << "\n";
#endif

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
    auto area = [](const Vector_2 a, const Vector_2 b, const Vector_2 c) { // TODO replace with orientation predicate
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

#ifdef CGAL_DEBUG_BSURF
        std::cout << "i=" << i << "\n";
        std::cout << "a=" << a << " b=" << b << " portal[0]=" << portal[0] << " portal[1]=" << portal[1] << "\n";
#endif
        FT s = intersect_segments(a, b, portal[0], portal[1]);

#ifdef CGAL_DEBUG_BSURF
        std::cout << "s=" << s << "\n";
#endif
        lerps.push_back(std::clamp(s, 0.0, 1.0));
      }
    }

    auto index = 1;
#ifdef CGAL_DEBUG_BSURF
    std::cout << "setting funnel_point indices\n";
#endif
    for (std::size_t i = 0; i < portals.size(); ++i)
    {
#ifdef CGAL_DEBUG_BSURF
      std::cout << "  i=" << i << " index = " << index << "\n";
      std::cout << "  portals[i][0]=" << portals[i][0] << " portals[i][1]=" << portals[i][1] << "\n";
      std::cout << "  points[index].pos = " <<points[index].pos << "\n";
#endif
      if ((portals[i][0] == points[index].pos) ||
          (portals[i][1] == points[index].pos))
      {
#ifdef CGAL_DEBUG_BSURF
        std::cout << "  setting point["<<index<<"].face="<< i << "\n";
#endif
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
                       Face_location<TriangleMesh, FT>& src,
                       Face_location<TriangleMesh, FT>& tgt,
                       const VertexPointMap &vpm, const TriangleMesh &mesh, std::size_t index)
  {
#ifdef CGAL_DEBUG_BSURF
    dump_path(path, lerps, src, tgt, mesh);
#endif

    typedef boost::graph_traits<TriangleMesh> BGT;
    vertex_descriptor vertex=BGT::null_vertex();

    // TODO: use a while loop breaking when no apex vertices not already visited are available
    for (std::size_t i = 0; i < portals.size() * 2 && index != std::size_t(-1); ++i)
    {
#ifdef CGAL_DEBUG_BSURF
      std::cout << "Improving path " << path.size() << " hedges\n";
      std::cout << "src = " << construct_point(src, mesh) << "\n";
      std::cout << "tgt = " << construct_point(tgt, mesh) << "\n";
#endif

      vertex_descriptor new_vertex=BGT::null_vertex();
      halfedge_descriptor h_curr       = path[index];

      bool is_target = false;
      if (lerps[index] == 0) {
        new_vertex = target(h_curr,mesh);
        is_target = true;
      } else if (lerps[index] == 1) {
        new_vertex = source(h_curr,mesh);
      }
      if (new_vertex == vertex) break;
      vertex = new_vertex;

#ifdef CGAL_DEBUG_BSURF
      std::cout << "  Current strip with Apex: " << get(vpm, new_vertex) <<  "\n";
      std::cout << "  is_target? " << is_target << "\n";
      std::cout << "  -- path --\n";
      for (auto h : path)
      {
        std::cout << "  4 " << get(vpm, source(h, mesh))
                  << "  " << get(vpm, target(h, mesh))
                  << "  " << get(vpm, target(next(h, mesh), mesh))
                  << "  " << get(vpm, source(h, mesh)) << std::endl;
        std::cout << edge(h, mesh) << std::endl;
      }
      std::cout << "  ----------\n";
#endif

      // if I hit the source vertex v of h_curr, then h_next has v as source, thus we turn ccw around v in path
      // Similarly, if I hit the target vertex v of h_curr, then  h_next has v as target, thus we turn cw around v in path
      std::size_t curr_index = index+1;

      std::vector<halfedge_descriptor> new_hedges;

      if (is_target)
      {
        face_descriptor target_face;

        if (curr_index == path.size())
        {
          target_face=tgt.first;
        }
        else
        {
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
        }

        halfedge_descriptor h_loop=opposite(prev(opposite(h_curr, mesh), mesh), mesh);

        while(target_face!=face(h_loop, mesh))
        {
          new_hedges.push_back(h_loop);
          h_loop=opposite(prev(h_loop,mesh), mesh);
        }
        new_hedges.push_back(h_loop);
      }
      else
      {
        face_descriptor target_face;

        if (curr_index == path.size())
        {
          target_face=tgt.first;
        }
        else
        {
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
        }
        halfedge_descriptor h_loop=opposite(next(opposite(h_curr, mesh), mesh), mesh); // skip the face before h_curr (as we won't remove it from path)

        while(target_face!=face(h_loop, mesh))
        {
          new_hedges.push_back(h_loop);
          h_loop=opposite(next(h_loop,mesh), mesh);
        }
        new_hedges.push_back(h_loop);
      }

      // replace the halfedges incident to the apex vertex with the opposite part of the ring
      std::vector<halfedge_descriptor> new_path(path.begin(),path.begin()+index);
      new_path.insert(new_path.end(), new_hedges.begin(), new_hedges.end());
      new_path.insert(new_path.end(), path.begin()+curr_index, path.end());
      path.swap(new_path);

      strip_path(mesh, src, tgt, path);

#ifdef CGAL_DEBUG_BSURF
      std::cout << "  -- new path --\n";
      for (auto h : path)
      {
        std::cout << "  4 " << get(vpm, source(h, mesh))
                  << "  " << get(vpm, target(h, mesh))
                  << "  " << get(vpm, target(next(h, mesh), mesh))
                  << "  " << get(vpm, source(h, mesh)) << std::endl;
        std::cout << face(h, mesh) << std::endl;
        std::cout << edge(h, mesh) << std::endl;
      }
      std::cout << "  ----------\n";
#endif

      portals=unfold_strip(path,src,tgt,vpm,mesh);
      lerps=funnel(portals,index);
      CGAL_assertion(lerps.size()==path.size());

#ifdef CGAL_DEBUG_BSURF
      dump_path(path, lerps, src, tgt, mesh);
#endif
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



  #if 1
  //:::::::::::::::::::::Parallel Transport::::::::::::::::::::::::::::
  static
  Vector_3 face_normal(const VertexPointMap &vpm,
                       const TriangleMesh &mesh,
                       const face_descriptor& f)
  {  halfedge_descriptor h=halfedge(f,mesh);
     Point_3 p0=get(vpm,source(h,mesh));
     Point_3 p1=get(vpm,target(h,mesh));
     Point_3 p2=get(vpm,target(next(h,mesh),mesh));
     Vector_3 n = triangle_normal(p0,p1,p2);
     n=n/sqrt(n.squared_length());

     return n;
  }
  static
  std::vector<Vector_3> polar_basis(const VertexPointMap &vpm,
                                    const TriangleMesh &mesh,
                                    const face_descriptor& f,
                                    const Vector_3& normal)
  {
     halfedge_descriptor h=halfedge(f,mesh);
     Point_3 p0=get(vpm,source(h,mesh));
     Point_3 p1=get(vpm,target(h,mesh));
     Vector_3 e=p1-p0;
     e=e/sqrt(e.squared_length());
     Vector_3 e_perp=cross_product(normal,e);

     return std::vector<Vector_3>{e,e_perp,normal};

  }


  // //:::::::::::::::::::::Straightest Geodesic::::::::::::::::::::::::::::
  static
  std::tuple<bool, std::array<double,3>> point_in_triangle(const VertexPointMap &vpm,
                                             const TriangleMesh &mesh,
                                             const face_descriptor& face,
                                             const Point_3 point,
                                             float tol=1e-5) {
  // http://www.r-5.org/files/books/computers/algo-list/realtime-3d/Christer_Ericson-Real-Time_Collision_Detection-EN.pdf
  // pag.48
  std::array<double,3> b = make_array(0.,0.,0.);
  halfedge_descriptor h=halfedge(face,mesh);
  Point_3 v0 = get(vpm,source(h,mesh));
  Point_3 v1 = get(vpm,target(h,mesh));
  Point_3 v2 = get(vpm,target(next(h,mesh),mesh));

  Vector_3 u = v1 - v0, v = v2 - v0, w = point - v0;
  double d00 = u.squared_length(), d01 = u*v, d11 = v.squared_length(), d20 = w*u,
       d21 = w*v, d = d00 * d11 - d01 * d01;

  if (d == 0)
    return {false, make_array(0.,0.,0.)};

  b[2] = (d00 * d21 - d01 * d20) / d;
  assert(!isnan(b[2]));
  b[1] = (d11 * d20 - d01 * d21) / d;
  assert(!isnan(b[1]));
  b[0] = 1 - b[1] - b[2];
  assert(!isnan(b[0]));

  for (auto i = 0; i < 3; ++i) {
    if (b[i] < -tol || b[i] > 1.0 + tol)
      return {false, make_array(0.,0.,0.)};
  }

  return {true, b};
}
  static
  Eigen::Matrix3d rot_matrix (const FT& angle, const Vector_3& axis)
  {
    Eigen::Matrix3d result;
    double c=std::cos(angle);
    double s=std::sin(angle);
    result <<c + (1 - c) * std::pow(axis.x(),2),(1 - c) * axis.x() * axis.y() + s * axis.z(),
            (1 - c) * axis.x() * axis.z() - s * axis.y(),
            (1 - c) * axis.x() * axis.y() - s * axis.z(),c + (1 - c) * axis.y() * axis.y(),
            (1 - c) * axis.y() * axis.z() + s * axis.x(),(1 - c) * axis.x() * axis.z() + s * axis.y(),
            (1 - c) * axis.y() * axis.z() - s * axis.x(),c + (1 - c) * axis.z() * axis.z();

    return result;
  }
  static
  Vector_3 rotate_vector(const Vector_3& v,const Vector_3& axis, const FT& angle)
  {
    Eigen::Matrix3d rot=rot_matrix(angle,axis);
    Eigen::Vector3d wrap_v;
    wrap_v<<v.x(),v.y(),v.z();
    Eigen::Vector3d rotated=rot*wrap_v;

    return Vector_3{rotated(0),rotated(1),rotated(2)};
  }
  static
  void parallel_transport_through_flattening(Vector_3& v,
                                             const VertexPointMap &vpm,
                                             const TriangleMesh &mesh,
                                             const face_descriptor& from,
                                             const face_descriptor& to)
  {
      halfedge_descriptor h_ref = halfedge(from, mesh);
      halfedge_descriptor h=h_ref;
      for (int i = 0; i < 3; ++i) {
        if (face(opposite(h_ref, mesh), mesh) == to)
         {
           h=h_ref;
           break;
         }
        h_ref = next(h_ref, mesh);
      }

    std::array<Vector_2,3> flat_from = init_flat_triangle(h,vpm,mesh);
     std::array<Vector_2,3> flat_to = unfold_face(h,vpm,mesh,flat_from);
    Vector_2 bary = Vector_2{0.333, 0.333};
    Vector_2 c0 = 0.33*(flat_from[0]+flat_from[1]+flat_from[2]);
    Vector_2 c1 = 0.33*(flat_to[0]+flat_to[1]+flat_to[2]);
    Vector_2 e0 = flat_from[0] - c0;
    Vector_2 e1 = flat_to[0] - c1;

    Vector_2 w = c1 - c0;
    FT phi_ij = angle(e0, w);
    if (e0.x()*w.y()-e0.y()*w.x() < 0)
      phi_ij = 2 * M_PI - phi_ij;
    w *= -1;
    FT phi_ji = angle(e1, w);

    if (e1.x()*w.y()-e1.y()*w.x() < 0)
      phi_ji = 2 * M_PI - phi_ji;

    Vector_3 n_from= face_normal(vpm, mesh, from);
    std::vector<Vector_3> e_from = polar_basis(vpm, mesh, from);
    double teta = angle(e_from, v);
    if (cross_product(e_from, v)*n_from < 0)
      teta = 2 * M_PI - teta;

    std::vector<Vector_3> e_to = polar_basis(vpm, mesh, from ,to);
    Vector_3 n_to = face_normal(vpm, mesh, to);
    double rot = teta + phi_ji + M_PI - phi_ij;
    e_to =e_to*sqrt(v.squared_length());
    v = rotate_vector(e_to, n_to, rot);

  }
  static
  Eigen::Matrix3d transformation_matrix(const Vector_2 &x1, const Vector_2 &y1,
                                       const Vector_2 &O1)
  {
    Eigen::Matrix3d T = Eigen::Matrix3d::Zero();
    T << x1.x(), y1.x(), O1.x(), x1.y(), y1.y(), O1.y(), 0, 0, 1;

    return T;
  }

  //h_ref is the reference halfedge of the face we are in, h_edge is the halfedge along which we want to unfold
  static
  Vector_2 compute_new_dir(const halfedge_descriptor h_ref,
                           const halfedge_descriptor h_edge,
                           const std::array<FT,3>& prev_coords,
                           const std::array<FT,3>& curr_coords,
                           const VertexPointMap& vpm,
                           const TriangleMesh& mesh)
  {
    int k=0;

    if(h_edge==prev(h_ref,mesh))
      k=2;
    else if(h_edge==next(h_ref,mesh))
      k=1;
    else
      assert(h_edge==h_ref);

    std::array<Vector_2,3> flat_curr = init_flat_triangle(h_ref,vpm,mesh);
    std::array<Vector_2,3> flat_prev = unfold_face(h_edge,vpm,mesh,flat_curr,k);

    Vector_2 prev_flat_p=prev_coords[0]*flat_prev[0]+prev_coords[1]*flat_prev[1]+prev_coords[2]*flat_prev[2];
    Vector_2 curr_flat_p=curr_coords[0]*flat_curr[0]+curr_coords[1]*flat_curr[1]+curr_coords[2]*flat_curr[2];

#ifdef CGAL_DEBUG_BSURF
    std::cout<<"   k is "<<k<<"\n";
    std::cout<<"   Flat_prev: ";
    std::cout << "  " << flat_prev[0] << " " << flat_prev[1] << " " << flat_prev[2] << "\n";
    std::cout<<"   Flat_curr: ";
    std::cout << "  " << flat_curr[0] << " " << flat_curr[1] << " " << flat_curr[2] << "\n";
    //Unfold face take care of the orientation so flat_t1[1] - flat_t1[0] is our old
    //reference frame
    std::cout<<"   h_ref_curr"<<std::endl;
    std::cout << "  " << source(h_ref,mesh) << " " << target(h_ref,mesh) << " " << target(next(h_ref,mesh),mesh) << "\n";
    halfedge_descriptor h_ref_prev=halfedge(face(opposite(h_edge,mesh),mesh),mesh);
    std::cout<<"   h_ref_prev"<<std::endl;
    std::cout << "  " << source(h_ref_prev,mesh) << " " << target(h_ref_prev,mesh) << " " << target(next(h_ref_prev,mesh),mesh) << "\n";
    std::cout<<"     prev_Flat_p "<<prev_flat_p.x()<<" "<< prev_flat_p.y()<<" )\n";
    std::cout<<"     curr_Flat_p "<<curr_flat_p.x()<<" "<< curr_flat_p.y()<<" )\n";
#endif

    return curr_flat_p-prev_flat_p;

  }

  //h_edge is the halfedge along which we want to unfold
  // TODO: use interval to control the error?
  static
  Vector_2 parallel_transport_f2f(const halfedge_descriptor h_edge, // face(opposite(h_edge)) = previous face face(h_edge)==current_face
                                  const Vector_2& prev_dir,
                                  const VertexPointMap& vpm,
                                  const TriangleMesh& mesh)
  {
    halfedge_descriptor h_ref = halfedge(face(h_edge,mesh), mesh);
    int k=0;

    if(h_edge==prev(h_ref,mesh))
      k=2;
    else if(h_edge==next(h_ref,mesh))
      k=1;
    else
      assert(h_edge==h_ref);

    std::array<Vector_2,3> flat_curr = init_flat_triangle(h_ref,vpm,mesh);
    std::array<Vector_2,3> flat_prev = unfold_face(h_edge,vpm,mesh,flat_curr,k);

    Vector_2 prev_origin = flat_prev[0];
    Vector_2 prev_y = flat_prev[1]-flat_prev[0];
    Vector_2 prev_x = Vector_2(prev_y.y(), - prev_y.x());

    Eigen::Matrix3d T= transformation_matrix(prev_x, prev_y, prev_origin);
    return switch_reference_frame_vector(T, prev_dir);
  }

  static Vector_2 switch_reference_frame_vector(const Eigen::Matrix3d &T, const Vector_2 &p) {
    Eigen::Vector3d V;
    V << p.x(), p.y(), 0;
    Eigen::Vector3d t = T * V;
    return Vector_2(FT(t(0)), FT(t(1)));
  }

  static
  void parallel_transport_through_flattening(Vector_2& v,
                                             const VertexPointMap &vpm,
                                             const TriangleMesh &mesh,
                                             const face_descriptor& from,
                                             const face_descriptor& to)
  {


  }

  Vector_3 parallel_transport_along_path(const std::vector<Edge_location<TriangleMesh, FT>>& edge_locations,
                                         const VertexPointMap &vpm,
                                         const TriangleMesh& mesh,
                                         const Face_location<TriangleMesh, FT>& src,
                                         const Face_location<TriangleMesh, FT>& tgt,
                                         const Vector_3& v)
  {
    if(edge_locations.size()==0)
    return v;
    Vector_3 w=v;
    face_descriptor f_start=src.first;
    for(std::size_t i=0;i<edge_locations.size()-1;++i)
      parallel_transport_through_flattening(w,vpm,mesh,f_start,face(halfedge(edge_locations[i],mesh),mesh));

    return w;
  }
  // static
  // std::tuple<Vector_3,face_descriptor>
  // polthier_condition_at_vert(const TriangleMesh& mesh,
  //                            const VertexPointMap &vpm
  //                            const vertex_descriptor& vid,
  //                            const face_descriptor& tid,
  //                            const int kv,
  //                            const Vector_3 &dir)
  // {
  //   FT total_angle=get_total_angle(vid,mesh,vpm);
  //   FT theta = 0.5 * total_angle;
  //   halfedge_descriptor h=halfedge(tid,mesh);
  //   while(target(h,mesh)!=vid)
  //     h=next(h,mesh);

  //   Point_3 vert = get(vpm,vid);
  //   Point_3 vert_adj=get(vpm,source(h,mesh));
  //   Vector_3 v = vert - vert_adj;

  //   FT acc = angle(v, dir);
  //   FT prev_angle = acc;
  //   face_descriptor curr_tid = tid;

  //   while (acc < theta) {
  //     h=prev(opposite(h.mesh),mesh);
  //     prev_angle = acc;
  //     Point_3 next_vert_adj=get(vpm,source(h,mesh));
  //     acc += angle(vert_adj - vert,next_vert_adj - vert);
  //     vert_adj=next_vert_adj;
  //   }
  //   auto offset = theta - prev;
  //   Point_3 prev_vert_adj=get(vpm,target(next(h,mesh,mesh)));

  //   Vector_3 result=Vector_3{prev_vert_adj.x()-vert.x(),prev_vert_adj.y()-vert.y(),prev_vert_adj.z()-vert.z()};
  //   face_descriptor f=face(h,mesh);
  //   Vector_3 n=face_normal(vpm,mesh,f);
  //   result=rotate_vector(result,n,offset);
  //   return {result,f};
  // }


  static
  std::vector<Vector_3>
  get_3D_basis_at_point(const Face_location<TriangleMesh, FT>& p,
                        const TriangleMesh& mesh,
                        const VertexPointMap &vpm)
  {
    halfedge_descriptor h=halfedge(p.first,mesh);

    Point_3 p0=get(vpm, source(h, mesh));
    Point_3 p1=get(vpm, target(h, mesh));
    Point_3 p2= get(vpm,target(next(h, mesh), mesh));

    Vector_3 e=p1-p0;
    e/=sqrt(e.squared_length());
    Vector_3 n= triangle_normal(p0,p1,p2);
    Vector_3 e1= cross_product(n,e);

    return {e,e1,n};
  }


  static
  std::tuple<bool, int>
  point_is_edge(const Face_location<TriangleMesh, FT>& p,
                const FT& tol=1e-5)
  {
    auto bary=p.second;
    if (bary[0] > tol && bary[1] > tol && bary[2] <= tol)
      return {true, 0};
    if (bary[1] > tol && bary[2] > tol && bary[0] <= tol)
      return {true, 1};
    if (bary[2] > tol && bary[0] > tol && bary[1] <= tol)
      return {true, 2};

    return {false, -1};
  }

//TODO get rid of tol or at least get it dependent on the input?
//
  static
  std::tuple<bool,int>
  point_is_vert(const Face_location<TriangleMesh, FT>& p,const FT& tol=1e-5)
  {
     auto bary=p.second;
     if (bary[0] > tol && bary[1] <= tol && bary[2] <= tol)
         return {true, 0};
     if (bary[1] > tol && bary[0] <= tol && bary[2] <= tol)
        return {true, 1};
     if (bary[2] > tol && bary[0] <= tol && bary[1] <= tol)
        return {true, 2};

   return {false, -1};
  }

  static
  //https://rootllama.wordpress.com/2014/06/20/ray-line-segment-intersection-test-in-2d/
  std::pair<FT, FT>
  intersect(const Vector_2 &direction, const Vector_2 &left,
            const Vector_2 &right,const Vector_2& origin)
  {
    Vector_2 v1 = origin-left;
    Vector_2 v2 = right - left;
    Vector_2 v3(-direction.y(), direction.x());
    double t0 = (v2.x()*v1.y()-v2.y()*v1.x()) / (v2*v3);
    double t1 = v1*v3/ ( v2*v3 );
    return std::make_pair(t0, t1);
  };

  static
  std::tuple<int,double>
  segment_in_tri(const Vector_2& p,
                 const std::array<Vector_2, 3>& tri,
                 const Vector_2& dir,const int offset)
  {
    //rotated the triangle in order to test intersection at meaningful edges before
    std::array<Vector_2, 3> rotated_tri=tri;
    // TODO rotated_tri.begin() + offset ?
    for(int k=0;k<offset;++k)
      std::rotate(rotated_tri.begin(), rotated_tri.begin() + 1, rotated_tri.end());
#ifdef CGAL_DEBUG_BSURF
    std::cout << "  offset " << offset << "\n";
    std::cout << "  tri " << tri[0] << " " << tri[1] << " " << tri[2] << "\n";
    std::cout << "  p " << p << "\n";
    std::cout << "  p+dir " << p+dir << "\n";
    std::cout << "  " << rotated_tri[0] << " " << rotated_tri[1] << " " << rotated_tri[2] << "\n";
#endif
    for (auto i = 0; i < 3; ++i)
    {
      auto [t0, t1] = intersect(dir, rotated_tri[(i+1)%3], rotated_tri[(i+2) % 3],p);
#ifdef CGAL_DEBUG_BSURF
      std::cout  << "  t0/t1 " << t0 << " / " << t1 << "\n";
#endif
      //TODO: replace intersection with CGAL code
      if (t0 > 0 && t1 >= -1e-4 && t1 <= 1 + 1e-4)
      {
#ifdef CGAL_DEBUG_BSURF
          int h=((i+1)%3 + offset)%3;
          Vector_2 intersection_point=(1-t1)*tri[h]+t1*tri[(h+1)%3];
          std::cout<<"  2D intersection point: "<< intersection_point << std::endl;
#endif
          return {((i+1)%3 + offset)%3, std::clamp(t1, 0., 1.)}; //return the offset w.r.t h_ref
      }
    }
    CGAL_assertion(false);
    return {-1,-1};
  }

  // return an angle in degrees
  static
  FT get_total_angle(const vertex_descriptor& vid,
                     const TriangleMesh& mesh,
                     const VertexPointMap &vpm)
  {
    halfedge_descriptor h_ref= halfedge(vid,mesh);
    halfedge_descriptor h_start=h_ref;
    halfedge_descriptor h_next = next(h_start,mesh);

    FT theta=approximate_angle(get(vpm,source(h_start,mesh))-get(vpm,target(h_start,mesh)),
                   get(vpm,target(h_next,mesh))-get(vpm,target(h_start,mesh)));
    h_start=opposite(h_next,mesh);
    h_next=next(h_start,mesh);
    while(h_start!=h_ref)
    {
      theta+=approximate_angle(get(vpm,source(h_start,mesh))-get(vpm,target(h_start,mesh)),
                   get(vpm,target(h_next,mesh))-get(vpm,target(h_start,mesh)));
      h_start=opposite(h_next,mesh);
      h_next=next(h_start,mesh);
    }

    return theta;
  }

  static
  std::tuple<Vector_2,face_descriptor, halfedge_descriptor>
  polthier_condition_at_vert(const TriangleMesh& mesh,
                             const VertexPointMap &vpm,
                             const vertex_descriptor& vid,
                             const face_descriptor& tid,
                             const FT& init_angle)
  {
     auto get_vid_offset=[&mesh](const halfedge_descriptor& h_ref,const vertex_descriptor& vid)
    {
      if(source(h_ref,mesh)==vid)
        return 0;

      if(target(h_ref,mesh)==vid)
        return 1;

      if(target(next(h_ref,mesh),mesh)==vid)
        return 2;

      std::cerr<<"Error! Halfedges are in different faces"<<std::endl;

      CGAL_assertion(false);

      return -1;
    };
    //TODO use interval for robustness + snap if ambiguous
    FT total_angle=get_total_angle(vid,mesh,vpm);
    FT theta = 0.5 * total_angle; //in degrees
    halfedge_descriptor h=halfedge(tid,mesh);

    while(target(h,mesh)!=vid)
      h=next(h,mesh);


    Point_3 vert = get(vpm,vid);
    Point_3 vert_adj=get(vpm,source(h,mesh));
    FT acc = init_angle;

    FT prev_angle = acc;
#ifdef CGAL_DEBUG_BSURF
    std::cout<<"initial h "<< edge(h,mesh)<<std::endl;
    std::cout<<"acc "<< acc<<std::endl;
    std::cout<<"theta "<< theta<<std::endl;
#endif
    while (acc < theta) {
      h=prev(opposite(h,mesh),mesh);
      prev_angle = acc;
      Point_3 next_vert_adj=get(vpm,source(h,mesh));
      acc += approximate_angle(vert_adj - vert,next_vert_adj - vert);
      vert_adj=next_vert_adj;
    }
    auto offset = theta - prev_angle;
    //moving to radiants
    offset*=CGAL_PI/180.;
    Point_3 prev_vert_adj=get(vpm,target(next(h,mesh),mesh));

    FT l = sqrt(squared_distance(prev_vert_adj,vert));
    FT phi = approximate_angle(vert - prev_vert_adj, vert_adj - prev_vert_adj);
    phi*=CGAL_PI/180.;
    FT x = l * std::sin(offset) / std::sin(M_PI - phi - offset);
    FT alpha = x / sqrt(squared_distance(vert_adj, prev_vert_adj));
    halfedge_descriptor prev_h=prev(h,mesh);

    std::array<Vector_2,3> flat_tid = init_flat_triangle(halfedge(face(h,mesh),mesh),vpm,mesh);
    int kv=get_vid_offset(halfedge(face(h,mesh),mesh),vid);


    Vector_2 q = (1 - alpha) * flat_tid[(kv+1)%3] + alpha * flat_tid[(kv+2)%3];
    Vector_2 new_dir = q - flat_tid[kv];
#ifdef CGAL_DEBUG_BSURF
    Vector_2 flat_p=flat_tid[kv];
    std::cout<<"final h "<< edge(h,mesh)<<std::endl;
    std::cout<<" prev_vert_adj "<< prev_vert_adj<<std::endl;
    std::cout<<" vert_adj "<< vert_adj<<std::endl;
    std::cout<< "------ outgoing tid"<<std::endl;
    std::cout << "  " << flat_tid[0] << " " << flat_tid[1] << " " << flat_tid[2] << "\n";

    std::cout<< "------ flat_p"<<std::endl;
    std::cout << "  " << flat_p<< "\n";
    std::cout<< "------ outgoing dir"<<std::endl;
    std::cout << "  " << new_dir[0]+flat_p[0] << " " << new_dir[1]+flat_p[1] << "\n";
    std::cout<<"outgoing face "<<face(prev_h,mesh)<<std::endl;
    std::cout<< "h_curr "<< edge(h,mesh)<<std::endl;
#endif
    return {new_dir,face(prev_h,mesh), h};
  }

/////////////////////////////////////////
/////////////////////////////////////////
/////////////////////////////////////////
/////////////////////////////////////////

  // dir is defined wrt the halfedge of start_loc.first that is used as y axis
  static
  std::vector<Face_location<TriangleMesh, FT>>
  straightest_geodesic(const Face_location<TriangleMesh, FT>& start_loc,
                       const TriangleMesh& mesh,
                       const VertexPointMap &vpm,
                       Vector_2 dir,const FT& len)
  {
#ifdef CGAL_DEBUG_BSURF
    {
      std::ofstream log("/tmp/log.txt");
      log << std::setprecision(17);
      log << start_loc.first << " " << start_loc.second[0]  << " " << start_loc.second[1]  << " " << start_loc.second[2] << "\n";
      log << "dir = " << dir << "\n";
      log << "len = " << len << "\n";
    }
#endif

    auto get_halfedge_offset=[&mesh](const halfedge_descriptor& h_ref,const halfedge_descriptor& h_curr)
    {
      if( h_ref == h_curr) return 0;
      if( next(h_ref, mesh) == h_curr ) return 1;
      if( prev(h_ref, mesh) == h_curr ) return 2;

      std::cout<<"Error! Halfedges are in different faces"<<std::endl;

      CGAL_assertion(false);
      return -1;
    };

#ifdef CGAL_DEBUG_BSURF
    auto fn = compute_face_normal(start_loc.first, mesh);
    double theta= std::acos(dir.y()/sqrt(dir*dir));
    if(dir.x()<0)
      theta*=-1;
    auto dir3 = rotate_vector(mesh.point(target(halfedge(start_loc.first, mesh), mesh)) -
                              mesh.point(source(halfedge(start_loc.first, mesh), mesh)),
                              fn, theta);

    std::cout << "direction " << construct_point(start_loc, mesh) << " " << construct_point(start_loc, mesh)+dir3 << "\n";
#endif

    auto get_vid_offset=[&mesh](const halfedge_descriptor& h_ref,const vertex_descriptor& vid)
    {
      if(source(h_ref,mesh)==vid)
        return 0;

      if(target(h_ref,mesh)==vid)
        return 1;

      if(target(next(h_ref,mesh),mesh)==vid)
        return 2;

      std::cerr<<"Error! Halfedges are in different faces"<<std::endl;

      CGAL_assertion(false);

      return -1;
    };

    auto get_halfedge=[&mesh](const int k,const halfedge_descriptor& h_ref)
    {
      switch(k)
      {
        case 0:
         return h_ref;

        case 1:
         return next(h_ref,mesh);

        default:
         return prev(h_ref,mesh);
      }
    };

    auto edge_barycentric_coordinate =
      [&mesh](const halfedge_descriptor h_edge,
              const halfedge_descriptor h_face,
              const std::array<FT,2>& bary_edge)
    {
      std::array<FT,3> bary_edge_in_face=make_array(0.,0.,0.);
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
          assert(h_face==prev(h_edge,mesh));
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

    std::vector<Face_location<TriangleMesh, FT>> result;
    FT accumulated=0.;
    face_descriptor curr_tid=start_loc.first;
    std::array<Vector_2, 3> curr_flat_tid=init_flat_triangle(halfedge(curr_tid,mesh),vpm,mesh);
    Vector_2 flat_p= start_loc.second[0]*curr_flat_tid[0]+start_loc.second[1]*curr_flat_tid[1]+start_loc.second[2]*curr_flat_tid[2];
    Face_location<TriangleMesh, FT> curr_p=start_loc;

    descriptor_variant<TriangleMesh> var_des = get_descriptor_from_location(start_loc,mesh);
    switch(var_des.index())
    {
      case 0:
      {
        vertex_descriptor src_vertex = std::get<vertex_descriptor>(var_des);
        halfedge_descriptor href = halfedge(curr_tid, mesh);
        int src_index = source(href, mesh)==src_vertex?0:(target(href,mesh)==src_vertex?1:2);
        // first get the direction in the basis of the triangle to get its angle with an edge incident to the vertex

        Vector_2 v02 = curr_flat_tid[(src_index+2)%3] - curr_flat_tid[src_index];
        Vector_2 v01 = curr_flat_tid[(src_index+1)%3] - curr_flat_tid[src_index];
        double unnorm_cos_theta_1 = v02 * dir;
        double unnorm_cos_theta_2 = v02 * v01;

        // 2 orientation tests to check if the angle is larger or smaller than Pi
        // 2 orientation tests are needed as curr_flat_tid orientation could be inverted during flattening
        Orientation dir_ori=orientation(v02,dir), v01_ori=orientation(v02,v01);

        std::cout << "dir_ori!=v01_ori ? " << (dir_ori!=v01_ori) << "\n";

        // check if dir is  the cone centered at curr_flat_tid[src_vertex] (TODO: should be a predicate!)
        if ( dir_ori!=v01_ori ||
            unnorm_cos_theta_2 / std::sqrt(v01*v01) > unnorm_cos_theta_1 / std::sqrt(dir*dir)) // > because cos is decreasing from 0 -> Pi, should be < for the angle`
        {
          //TODO: should be made robust with snapping and predicates or something à la robust construction traits

          // compute all the angles around the vertex from where the path starts
          std::vector<double> angles;
          double acc_angle=0;
          halfedge_descriptor hloop=href;
          while (target(hloop, mesh)!=src_vertex) hloop=next(hloop, mesh);
          halfedge_descriptor hstart=hloop;
          Vector_3 n(get(vpm,target(hloop, mesh)), get(vpm,source(hloop, mesh)));
          n /= std::sqrt(n*n);
          do{
            Vector_3 nv(get(vpm,target(hloop, mesh)), get(vpm,target(next(hloop, mesh), mesh)));
            nv /= std::sqrt(nv*nv);
            angles.push_back( std::acos(n * nv) );
            acc_angle+=angles.back();
            n = nv;
            hloop=opposite(next(hloop, mesh), mesh);
          }while(hloop!=hstart);

          std::cout << std::acos(unnorm_cos_theta_2 / std::sqrt((v01*v01) * (v02*v02))) << "  vs  " << angles[0] << "\n";

          CGAL_assertion( std::acos(unnorm_cos_theta_2 / std::sqrt((v01*v01) * (v02*v02))) == angles[0]); // will not be true because of rounding

          // normalise the angle to bring them in [0,1]
          for (double& angle : angles)
            angle /= acc_angle;

          // compute and normalise the target angle
          double target_theta = std::acos(unnorm_cos_theta_1 / std::sqrt((dir*dir) * (v02*v02)));
          if ( dir_ori!=v01_ori ) target_theta = 2 * CGAL_PI - target_theta;
          target_theta /= acc_angle;

          double curr_angle = 0;
          int ia = 0;
          do{
            curr_angle+=angles[ia];
            if (curr_angle >= target_theta) break;
            hloop=opposite(next(hloop, mesh), mesh);
            ++ia;
          }while(hloop!=hstart);

          CGAL_assertion(hloop!=hstart);

          // angle in target face wrt hloop
          double delta = (target_theta - curr_angle + angles[ia]) * acc_angle;

          // using the law of the sinus we have:
          CGAL_assertion(target(hloop, mesh) == src_vertex);
          Vector_3 vloop(get(vpm, source(hloop, mesh)),get(vpm, target(hloop, mesh)));
          Vector_3 vprev(get(vpm, source(hloop, mesh)),get(vpm, target(next(hloop, mesh), mesh)));
          double dvloop = std::sqrt(vloop*vloop);
          double dvprev = std::sqrt(vprev*vprev);
          double theta = std::acos( vloop*vprev / dvloop / dvprev );
          double d_opp_delta =  dvloop * sin(delta) / sin(CGAL_PI-delta-theta);

          double alpha = d_opp_delta/dvprev;

          halfedge_descriptor href = halfedge(face(hloop,mesh),mesh);
          curr_flat_tid=init_flat_triangle(href,vpm,mesh);
          int src_id = source(href, mesh)==src_vertex?0:(target(href,mesh)==src_vertex?1:2);

#ifdef CGAL_DEBUG_BSURF
          std::cout << "new_face= " << face(href, mesh) << "\n";
          std::array<Vector_3, 3> debug_pt;
          debug_pt[0]=mesh.point(source(href,mesh)) - Point_3(0.,0.,0.);
          debug_pt[1]=mesh.point(target(href,mesh)) - Point_3(0.,0.,0.);
          debug_pt[2]=mesh.point(target(next(href,mesh),mesh)) - Point_3(0.,0.,0.);
          std::cout << "direction inter pt: " << (alpha) * debug_pt[(src_id+1)%3] + (1-alpha) * debug_pt[(src_id+2)%3] << "\n";
#endif
          // point on the opposite edge of src_vertex in face(href, mesh)
          Vector_2 ip = (alpha) * curr_flat_tid[(src_id+1)%3] + (1-alpha) * curr_flat_tid[(src_id+2)%3];
          dir = ip - curr_flat_tid[src_id];
          curr_tid=face(href, mesh);
          flat_p = curr_flat_tid[src_id];
          curr_p.second=CGAL::make_array(0.,0.,0.);
          curr_p.second[src_id]=1.;
          curr_p.first=curr_tid;
        }
      }
      break;
      case 1:
      {
        halfedge_descriptor h = std::get<halfedge_descriptor>(var_des);
        halfedge_descriptor href = halfedge(face(h,mesh), mesh);
        int opp_id = h==href?2:(next(href,mesh)==h?0:1);
        Vector_2 opp = curr_flat_tid[opp_id];

        // check if the line starts in the current face (TODO should be a predicate)
        if ( (opp-flat_p) * dir < 0 )
        {
          // take the same point but in the opposite face
          halfedge_descriptor h_start=h;
          h=opposite(h, mesh);
          curr_tid=face(h, mesh);
          curr_p.first=curr_tid;
          href=halfedge(curr_tid, mesh);
          std::array<FT, 2> bary2 = CGAL::make_array(start_loc.second[(opp_id+2)%3], start_loc.second[(opp_id+1)%3]); // swap done here
          int opp_id = h==href?2:(next(href,mesh)==h?0:1);
          curr_p.second[0]=FT(0);
          curr_p.second[(opp_id+1)%3]=bary2[0];
          curr_p.second[(opp_id+2)%3]=bary2[1];


          // dir must also be updated:
          // unfold the new triangle into the basis of the input one
          // compute the angle between the new triangle ref halfedge and dir
          // compute the new dir thanks to that angle.
          halfedge_descriptor start_href=halfedge(start_loc.first, mesh);
          int k=h_start==prev(start_href,mesh)?2 :( h_start==next(start_href,mesh)?1:0);
          std::array<Vector_2,3> new_flat_tid_in_curr_basis=unfold_face(h_start,vpm, mesh,curr_flat_tid, k);

          Vector_2 ybase = new_flat_tid_in_curr_basis[1]-new_flat_tid_in_curr_basis[0];
          double cos_theta = ybase * dir / std::sqrt(ybase*ybase) / std::sqrt(dir*dir);
          double theta = std::acos(cos_theta);
          dir = Vector_2(std::cos(CGAL_PI/2.-theta), std::sin(CGAL_PI/2.-theta));

          curr_flat_tid=init_flat_triangle(halfedge(curr_tid,mesh),vpm,mesh);
          flat_p= curr_p.second[0]*curr_flat_tid[0]+curr_p.second[1]*curr_flat_tid[1]+curr_p.second[2]*curr_flat_tid[2];
        }
      }
      break;
      default:
      break;
    }

    Face_location<TriangleMesh, FT> prev_p;
    Vector_2 curr_dir=dir;
    halfedge_descriptor h_ref=halfedge(curr_tid,mesh);
    halfedge_descriptor h_curr=h_ref;



    result.push_back(curr_p);
#ifdef CGAL_DEBUG_BSURF
  std::cout << "p= " << construct_point(curr_p,mesh) << ")\n";
#endif

    //TODO not sure why we need that
    auto [is_vert, kv] = point_is_vert(curr_p);
    auto [is_edge, ke] = point_is_edge(curr_p);
    if (is_vert)
      h_curr=get_halfedge(kv,h_ref);
    else if (is_edge)
      h_curr=get_halfedge(ke,h_ref);

#ifdef CGAL_DEBUG_BSURF
    std::cout << "Accumulated loop starts\n";
#endif
    while (accumulated < len)
    {
#ifdef CGAL_DEBUG_BSURF
      std::cout << "--->" << accumulated << " vs " << len << "\n";
#endif
      int curr_offset=get_halfedge_offset(h_ref,h_curr);

      auto [k, t1] = segment_in_tri(flat_p, curr_flat_tid, curr_dir,curr_offset);
#ifdef CGAL_DEBUG_BSURF
      std::cout << "  t1 = " << t1 << std::endl;
#endif
      CGAL_assertion(k!=-1);
      std::array<FT,3> new_bary=make_array(0.,0.,0.);
      Face_location<TriangleMesh, FT> point_on_edge;

      new_bary[k] = 1 - t1;
      new_bary[(k + 1) % 3] = t1;
      point_on_edge.first=curr_tid;
      point_on_edge.second=new_bary;
      std::tie(is_vert, kv) = point_is_vert(point_on_edge);
#ifdef CGAL_DEBUG_BSURF
      std::cout << "  is_vert? " << is_vert << "\n";
#endif
      if (is_vert)
      {
        vertex_descriptor vid = target(get_halfedge(kv, prev(h_ref,mesh)), mesh);
#ifdef CGAL_DEBUG_BSURF
        std::cout<< "hit vertex "<< vid <<std::endl;
        std::cout<< "href "<< edge(h_ref,mesh) <<std::endl;
        std::cout<< "kv "<< kv <<std::endl;
#endif
        accumulated +=
            sqrt(squared_distance(construct_point(curr_p,mesh), get(vpm,vid)));

  //   Point_3 vert = get(vpm,vid);
  //   Point_3 vert_adj=get(vpm,source(h,mesh));
  //   Vector_3 v = vert - vert_adj;

        //TODO add a 2D version of approximate_angle in CGAL
        // FT init_angle = approximate_angle(curr_flat_tid[(kv+2)%3]-curr_flat_tid[kv], curr_dir);
        auto tmp =curr_flat_tid[kv]-curr_flat_tid[(kv+2)%3];
        FT init_angle = approximate_angle(Vector_3(tmp.x(), tmp.y(), 0), Vector_3(curr_dir.x(), curr_dir.y(), 0));
        std::tie(curr_dir, curr_tid, h_curr) =
            polthier_condition_at_vert(mesh,vpm,vid,curr_tid, init_angle);

        h_ref=halfedge(curr_tid,mesh);
        kv=get_vid_offset(h_ref,vid);
        CGAL_assertion(kv!=-1);
        std::array<FT,3> tmp_bary=make_array(0.,0.,0.);
        tmp_bary[kv]=1;

        curr_flat_tid=init_flat_triangle(h_ref,vpm,mesh);
        prev_p = curr_p;

        curr_p=Face_location<TriangleMesh,FT>(curr_tid,tmp_bary);
        int k=get_vid_offset(h_ref,target(h_curr,mesh));
        flat_p=curr_flat_tid[k];
      }
      else
      {
        h_curr=opposite(get_halfedge(k, h_ref),mesh);

        face_descriptor adj = face(h_curr,mesh);
        std::array<FT,2> curr_alpha=make_array(t1,1-t1); //reversed because will switch face
        new_bary=edge_barycentric_coordinate(h_curr,halfedge(adj,mesh),curr_alpha);
        prev_p = curr_p;
        curr_p.first=adj;
        curr_p.second= new_bary;
        accumulated += sqrt(squared_distance(construct_point(curr_p,mesh), construct_point(prev_p,mesh)));
        curr_tid = adj;
        h_ref=halfedge(curr_tid,mesh);
        //TODO curr_dir should be normalized everytime (to try to avoid numerical errors)
        curr_dir = compute_new_dir(h_ref,h_curr,prev_p.second,curr_p.second,vpm,mesh);
        curr_flat_tid=init_flat_triangle(h_ref,vpm,mesh);
        flat_p= curr_p.second[0]*curr_flat_tid[0]+curr_p.second[1]*curr_flat_tid[1]+curr_p.second[2]*curr_flat_tid[2];
#ifdef CGAL_DEBUG_BSURF
        std::cout << "  h_curr " << edge(h_curr, mesh)<<"\n";
        std::cout << "  adj " << adj<<"\n";
        Vector_2 intersection_point=new_bary[0]*curr_flat_tid[0]+new_bary[1]*curr_flat_tid[1]+new_bary[2]*curr_flat_tid[2];
        std::cout << "  New intersection point is "<< intersection_point<<std::endl;
        std::cout << "  barycentric coordinates"<<new_bary[0]<<" "<<new_bary[1]<< " "<<new_bary[2]<<"\n";
        std::cout << "  curr_flat_tid"<<std::endl;
        std::cout << "  " << curr_flat_tid[0] << " " << curr_flat_tid[1] << " " << curr_flat_tid[2] << "\n";
        std::cout << "  New_Flat_p "<<flat_p.x()<<" "<< flat_p.y()<<" )\n";
        std::cout << "  New Dir"<<flat_p.x()<<" "<< flat_p.y()<<" "<<flat_p.x()+curr_dir.x()<<" "<<flat_p.y()+ curr_dir.y()<<"\n";
#endif
      }
      result.push_back(curr_p);
#ifdef CGAL_DEBUG_BSURF
      std::cout << "  inter pt: " << construct_point(curr_p, mesh) << "\n";
#endif
    }

    double excess = accumulated - len;
    Point_3 prev_pos = construct_point(*std::next(result.rbegin()),mesh);
    Point_3 last_pos = construct_point(result.back(),mesh);
    double alpha = excess / sqrt((last_pos - prev_pos).squared_length());
    Point_3 pos = barycenter(prev_pos, alpha, last_pos, 1-alpha);
#ifdef CGAL_DEBUG_BSURF
    std::cout << "excess " << excess << std::endl;
    std::cout << "prev_pos " << prev_pos << std::endl;
    std::cout << "last_pos " << last_pos << std::endl;

    std::cout << "curr_tri "<< prev_p.first << std::endl;
    std::cout << "pos " << pos << std::endl;

#endif

    auto [inside, bary] =
        point_in_triangle(vpm,mesh,prev_p.first,pos);
    if (!inside)
      std::cout << "error!This point should be in the triangle" << std::endl;

    result.pop_back();
    prev_p.second=bary;
    result.push_back(prev_p);

    return result;
  }

  // static
  // std::vector<mesh_point>
  // polthier_straightest_geodesic(const vector<vec3i> &triangles, const vector<vec3f> &positions,
  //                               const vector<vec3i> &adjacencies, const vector<vector<int>> &v2t,
  //                               const vector<vector<float>> &angles, const vector<float> &total_angles,
  //                               const mesh_point &p, const vec2f &dir, const float &len)
  // {
  //   auto result = vector<mesh_point>{};
  //   auto accumulated = 0.f;
  //   auto curr_tid = p.face;
  //   auto curr_p = p;
  //   auto prev_p = mesh_point{};
  //   auto curr_dir = dir;

  //   result.push_back(p);

  //   auto k_start = 0;
  //   auto [is_vert, kv] = point_is_vert(p);
  //   auto [is_edge, ke] = point_is_edge(p);
  //   if (is_vert)
  //     k_start = kv;
  //   else if (is_edge)
  //     k_start = ke;

  //   auto count = 0;
  //   while (accumulated < len)
  //   {
  //     ++count;
  //     auto [k, t1] = straightest_path_in_tri(triangles, positions, curr_p,
  //                                           curr_dir, k_start);

  //     auto new_bary = zero3f;
  //     auto point_on_edge = mesh_point{};
  //     if (k != -1) {
  //       new_bary[k] = 1 - t1;
  //       new_bary[(k + 1) % 3] = t1;
  //       point_on_edge = mesh_point{curr_tid, vec2f{new_bary.y, new_bary.z}};
  //     } else {
  //       std::tie(is_edge, ke) = point_is_edge(curr_p, 5e-3);
  //       std::tie(is_vert, kv) = point_is_vert(curr_p, 5e-3);
  //       auto bary = get_bary(curr_p.uv);
  //       if (is_edge) {
  //         k = ke;
  //         t1 = bary[(k + 1) % 3];
  //         point_on_edge = curr_p;
  //       } else if (is_vert) {
  //         auto bary3 = zero3f;
  //         bary3[kv] = 1;
  //         point_on_edge = {curr_p.face, {bary3.y, bary3.z}};
  //       } else {
  //         std::cout << "Error!This should not happen" << std::endl;
  //         return result;
  //       }
  //     }

  //     std::tie(is_vert, kv) = point_is_vert(point_on_edge);

  //     if (is_vert) {
  //       auto vid = triangles[curr_tid][kv];
  //       if (angles[vid].size() == 0)
  //         return result;

  //       accumulated +=
  //           length(eval_position(triangles, positions, curr_p) - positions[vid]);
  //       auto dir3d = normalize(eval_position(triangles, positions, curr_p) -
  //                             positions[vid]);
  //       prev_p = curr_p;
  //       std::tie(curr_dir, curr_p, k_start) =
  //           polthier_condition_at_vert(triangles, positions, adjacencies,
  //                                     total_angles, vid, curr_tid, dir3d);
  //       curr_tid = curr_p.face;
  //       if (curr_tid == -1)
  //         return result;

  //     } else {
  //       auto adj = adjacencies[curr_tid][k];
  //       if (adj == -1)
  //         return result;
  //       auto h = find(adjacencies[adj], curr_tid);

  //       new_bary = zero3f;
  //       new_bary[h] = t1;
  //       new_bary[(h + 1) % 3] = 1 - t1;

  //       prev_p = curr_p;
  //       curr_p = mesh_point{adj, vec2f{new_bary.y, new_bary.z}};
  //       accumulated += length(eval_position(triangles, positions, curr_p) -
  //                             eval_position(triangles, positions, prev_p));

  //       auto T = switch_reference_frame(triangles, positions, adj, curr_tid);
  //       curr_dir = switch_reference_frame_vector(T, curr_dir);

  //       curr_tid = adj;
  //       k_start = h;
  //     }

  //     result.push_back(curr_p);
  //   }

  //   auto excess = accumulated - len;
  //   auto prev_pos = eval_position(triangles, positions, result.rbegin()[1]);
  //   auto last_pos = eval_position(triangles, positions, result.back());
  //   auto alpha = excess / length(last_pos - prev_pos);
  //   auto pos = alpha * prev_pos + (1 - alpha) * last_pos;

  //   auto [inside, bary] =
  //       point_in_triangle(triangles, positions, prev_p.face, pos);
  //   if (!inside)
  //     std::cout << "error!This point should be in the triangle" << std::endl;

  //   result.pop_back();
  //   result.push_back(mesh_point{prev_p.face, bary});

  //   return result;
  // }
#endif

  static
  void
  strip_path(const TriangleMesh& tmesh,
             Face_location<TriangleMesh, FT>& src,
             Face_location<TriangleMesh, FT>& tgt,
             std::vector<halfedge_descriptor>& initial_path)
  {
    // retrieve the vertex id of a face location describing a vertex
    auto is_vertex = [](const Face_location<TriangleMesh, FT>& fl)
    {
      if (fl.second[0]==0 && fl.second[1]==0)
        return 2;
      if (fl.second[1]==0 && fl.second[2]==0)
        return 0;
      if (fl.second[0]==0 && fl.second[2]==0)
        return 1;
      return -1;
    };

    // retrieve the source vertex id of a face location describing a point on a halfedge
    auto is_edge = [](const Face_location<TriangleMesh, FT>& fl)
    {
      if (fl.second[0]==0 && fl.second[1]!=0 && fl.second[2]!=0)
        return 2;
      if (fl.second[1]==0 && fl.second[2]!=0 && fl.second[0]!=0)
        return 0;
      if (fl.second[2]==0 && fl.second[0]!=0 && fl.second[1]!=0)
        return 1;
      return -1;
    };

    // strip initial_path from extra halfedges after src and also update src
    int vi = is_vertex(src);
    if (vi!=-1)
    {
      halfedge_descriptor hsrc=prev(halfedge(src.first, tmesh), tmesh);
      for (int i=0; i<vi; ++i)
        hsrc=next(hsrc, tmesh);
      vertex_descriptor vsrc = target(hsrc, tmesh);
      int last_hi=0;
      for (halfedge_descriptor hip : initial_path)
      {
        if (source(hip, tmesh)!=vsrc && target(hip, tmesh)!=vsrc)
          break;
        ++last_hi;
      }
      initial_path.erase(initial_path.begin(), std::next(initial_path.begin(),last_hi));
      if (initial_path.empty()) return;

      if (last_hi!=0)
      {
        hsrc=halfedge(face(opposite(initial_path[0], tmesh), tmesh), tmesh);
        int vid=0;
        while (source(hsrc, tmesh)!=vsrc)
        {
          hsrc=next(hsrc, tmesh);
          ++vid;
          CGAL_assertion(vid!=3);
        }
        std::array<FT, 3> fl_src = {FT(0),FT(0),FT(0)};
        fl_src[vid]=FT(1);
        src=std::make_pair(face(hsrc, tmesh), fl_src);
      }
    }
    else
    {
      int ei = is_edge(src);
      if (ei!=-1)
      {
        halfedge_descriptor hsrc=prev(halfedge(src.first, tmesh), tmesh);
        for (int i=0; i<ei; ++i)
          hsrc=next(hsrc, tmesh);
        if (opposite(initial_path[0], tmesh)==hsrc)
        {
          initial_path.erase(initial_path.begin(), std::next(initial_path.begin()));
          if (initial_path.empty()) return;

          halfedge_descriptor new_hsrc=prev(halfedge(face(opposite(initial_path[0], tmesh), tmesh), tmesh), tmesh);
          int hid=0;
          while (opposite(new_hsrc, tmesh)!=hsrc)
          {
            new_hsrc=next(new_hsrc, tmesh);
            ++hid;
            CGAL_assertion(hid!=3);
          }
          std::array<FT, 3> fl_src = {FT(0),FT(0),FT(0)};
          fl_src[hid]=src.second[(ei+2)%3];
          fl_src[(hid+2)%3]=src.second[ei];
          src=std::make_pair(face(new_hsrc, tmesh), fl_src);
        }
      }
    }

    // strip initial_path from extra halfedges before tgt and also update tgt
    vi = is_vertex(tgt);
    if (vi!=-1)
    {
      halfedge_descriptor htgt=prev(halfedge(tgt.first, tmesh), tmesh);
      for (int i=0; i<vi; ++i)
        htgt=next(htgt, tmesh);
      vertex_descriptor vtgt = target(htgt, tmesh);
      std::size_t first_hi=initial_path.size();
      for (auto rit=initial_path.rbegin(); rit!=initial_path.rend(); ++rit)
      {
        if (source(*rit, tmesh)!=vtgt && target(*rit, tmesh)!=vtgt)
          break;
        --first_hi;
      }
      bool update_tgt = (first_hi!=initial_path.size());
      initial_path.erase(std::next(initial_path.begin(), first_hi), initial_path.end());
      if (initial_path.empty()) return;

      if (update_tgt)
      {
        halfedge_descriptor new_htgt=halfedge(face(initial_path.back(), tmesh), tmesh);
        int vid=0;
        while (source(new_htgt, tmesh)!=vtgt)
        {
          new_htgt=next(new_htgt, tmesh);
          ++vid;
          CGAL_assertion(vid!=3);
        }
        std::array<FT, 3> fl_tgt = {FT(0),FT(0),FT(0)};
        fl_tgt[vid]=FT(1);
        tgt=std::make_pair(face(new_htgt, tmesh), fl_tgt);
      }
    }
    else
    {
      int ei = is_edge(tgt);
      if (ei!=-1)
      {
        halfedge_descriptor htgt=prev(halfedge(tgt.first, tmesh), tmesh);
        for (int i=0; i<ei; ++i)
          htgt=next(htgt, tmesh);
        if (htgt==initial_path.back())
        {
          initial_path.pop_back();
          if (initial_path.empty()) return;

          halfedge_descriptor new_htgt=prev(halfedge(face(initial_path.back(), tmesh), tmesh), tmesh);
          int hid=0;
          while (opposite(new_htgt, tmesh)!=htgt)
          {
            new_htgt=next(new_htgt, tmesh);
            ++hid;
            CGAL_assertion(hid!=3);
          }
          std::array<FT, 3> fl_tgt = {FT(0),FT(0),FT(0)};
          fl_tgt[hid]=tgt.second[(ei+2)%3];
          fl_tgt[(hid+2)%3]=tgt.second[ei];
          tgt=std::make_pair(face(new_htgt, tmesh), fl_tgt);
        }
      }
    }
  }
}; // end of Locally_shortest_path_imp

template <class K, class TriangleMesh, class VertexPointMap>
struct Bezier_tracing_impl
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
                const Face_location<TriangleMesh, FT>& tgt,const FT& t
#ifndef CGAL_BSURF_USE_DIJKSTRA_SP
               , const Dual_geodesic_solver<FT>& solver
#endif
  )
  {
    std::vector<Edge_location<TriangleMesh, FT>> edge_locations;
    locally_shortest_path<FT>(src,tgt,mesh, edge_locations, solver);
    std::vector<FT> parameters=path_parameters(edge_locations,mesh,src,tgt);
    Face_location<TriangleMesh, FT> point =
      eval_point_on_geodesic(edge_locations,mesh,src,tgt,parameters,t);
    return point;
  }


  static
  std::pair<Bezier_segment<TriangleMesh, FT>,Bezier_segment<TriangleMesh,FT>>
  subdivide_bezier_polygon(const TriangleMesh& mesh,
                           const Bezier_segment<TriangleMesh,FT>& polygon,
                           const FT& t
#ifndef CGAL_BSURF_USE_DIJKSTRA_SP
                           , const Dual_geodesic_solver<FT>& solver
#endif
  )
  {
#ifndef CGAL_BSURF_USE_DIJKSTRA_SP
     Face_location<TriangleMesh, FT> Q0 = geodesic_lerp(mesh, polygon[0], polygon[1], t, solver);
     Face_location<TriangleMesh, FT> Q1 = geodesic_lerp(mesh, polygon[1], polygon[2], t, solver);
     Face_location<TriangleMesh, FT> Q2 = geodesic_lerp(mesh, polygon[2], polygon[3], t, solver);
     Face_location<TriangleMesh, FT> R0 = geodesic_lerp(mesh, Q0, Q1, t, solver);
     Face_location<TriangleMesh, FT> R1 = geodesic_lerp(mesh, Q1, Q2, t, solver);
     Face_location<TriangleMesh, FT> S  = geodesic_lerp(mesh, R0, R1, t, solver);
#else
     Face_location<TriangleMesh, FT> Q0 = geodesic_lerp(mesh, polygon[0], polygon[1], t);
     Face_location<TriangleMesh, FT> Q1 = geodesic_lerp(mesh, polygon[1], polygon[2], t);
     Face_location<TriangleMesh, FT> Q2 = geodesic_lerp(mesh, polygon[2], polygon[3], t);
     Face_location<TriangleMesh, FT> R0 = geodesic_lerp(mesh, Q0, Q1, t);
     Face_location<TriangleMesh, FT> R1 = geodesic_lerp(mesh, Q1, Q2, t);
     Face_location<TriangleMesh, FT> S  = geodesic_lerp(mesh, R0, R1, t);
#endif
    return {{polygon[0], Q0, R0, S}, {S, R1, Q2, polygon[3]}};
  }
};


template <class K, class TriangleMesh, class VertexPointMap, class VertexIndexMap, class FaceIndexMap>
struct Geodesic_circle_impl
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

  struct geodesic_solver {
    struct graph_edge {
      int node = -1;
      double len=DBL_MAX;
    };
    std::vector<std::vector<graph_edge>> graph;
  };

  static
  void connect_nodes(geodesic_solver &solver,
                    const vertex_descriptor& a,
                    const vertex_descriptor& b,
                    const VertexIndexMap& vidmap,
                    const FT& len)
  {
    // TODO: avoid cast
    uint vida=get(vidmap,a);
    uint vidb=get(vidmap,b);
    solver.graph[vida].push_back({static_cast<int>(vidb), len});
    solver.graph[vidb].push_back({static_cast<int>(vida), len});

  }

  static
  double opposite_nodes_arc_length(const VertexPointMap &vpm,
                                   const vertex_descriptor& a,
                                   const vertex_descriptor& c,
                                   const vertex_descriptor& b,
                                   const vertex_descriptor& d)
  {
    // Triangles (a, b, d) and (b, d, c) are connected by (b, d) edge
    // Nodes a and c must be connected.

    Vector_3 ba = get(vpm,a) - get(vpm,b);
    Vector_3 bc = get(vpm,c) - get(vpm,b);
    Vector_3 bd = get(vpm,d) - get(vpm,b);

    Vector_3 ba_norm=ba/sqrt(ba.squared_length());
    Vector_3 bd_norm=bd/sqrt(bd.squared_length());
    Vector_3 bc_norm=bc/sqrt(bc.squared_length());

    double cos_alpha = ba_norm * bd_norm;
    double cos_beta  = bc_norm * bd_norm;
    double sin_alpha = sqrt(std::max(0.0, 1 - cos_alpha * cos_alpha));
    double sin_beta  = sqrt(std::max(0.0, 1 - cos_beta * cos_beta));


    // cos(alpha + beta)
    double cos_alpha_beta = cos_alpha * cos_beta - sin_alpha * sin_beta;
    CGAL_assertion(cos_alpha_beta>-1);
    if (cos_alpha_beta <= -1) return DBL_MAX;

    // law of cosines (generalized Pythagorean theorem)
    double len = ba.squared_length() + bc.squared_length() -
                 sqrt(ba.squared_length()) * sqrt(bc.squared_length()) * 2 * cos_alpha_beta;

     CGAL_assertion(len>0);
    if (len <= 0)
      return DBL_MAX;
    else
      return sqrt(len);
  }

  static
  void connect_opposite_nodes(geodesic_solver& solver,
                              const VertexPointMap &vpm,
                              const TriangleMesh &mesh,
                              const VertexIndexMap& vidmap,
                              const vertex_descriptor& a,
                              const vertex_descriptor& b,
                              const halfedge_descriptor& h)
  {
    vertex_descriptor v0 =target(next(h,mesh),mesh);
    vertex_descriptor v1= target(next(opposite(h,mesh),mesh),mesh);
    auto len = opposite_nodes_arc_length(vpm, v0,v1,a,b);
     //std::cout<<len<<std::endl;
    connect_nodes(solver, v0, v1, vidmap, len);
  }

  static
  geodesic_solver
  make_geodesic_solver(const VertexPointMap &vpm,
                       const VertexIndexMap& vidmap,
                       const TriangleMesh &mesh)
  {
    geodesic_solver solver;
    solver.graph.resize(vertices(mesh).size());
    for (face_descriptor f : faces(mesh)) {
      halfedge_descriptor h=halfedge(f,mesh);
      for(std::size_t i=0;i<3;++i)
    {
      vertex_descriptor a=source(h,mesh);
      vertex_descriptor b=target(h,mesh);

      if(a<b)
        {
          FT len=sqrt(squared_distance(get(vpm,a),get(vpm,b)));
          //std::cout<<len<<std::endl;
          CGAL_assertion(len<DBL_MAX && len>0);
          connect_nodes(solver,a,b,vidmap,len);
        }
        face_descriptor nei=face(opposite(h,mesh),mesh);
        if(f<nei)
            connect_opposite_nodes(solver,vpm,mesh,vidmap,a,b,h);

        h=next(h,mesh);


      }


    }
    return solver;
  }
  // geodesic_solver
  // extended_geodesic_solver(const VertexPointMap& vpm,
  //                          const VertexIndexMap& vidmap,
  //                          const TriangleMesh& mesh,
  //                          const dual_geodesic_solver& dual_solver,
  //                          const int k)
  // {
  //   geodesic_solver solver;
  //   solver.graph.resize(vertices(mesh).size());
  //   for(auto& v:vertices(mesh))
  //   {

  //   }
  // }
  static
  Dual_geodesic_solver<FT>
  make_dual_geodesic_solver(const VertexPointMap &vpm,
                            const FaceIndexMap& tidmap,
                            const TriangleMesh &mesh)
  {
    auto compute_dual_weights=[&mesh,&vpm](const halfedge_descriptor& h)
    {
      using Impl = Locally_shortest_path_imp<K, TriangleMesh, VertexPointMap>;
      std::array<Vector_2,3> flat_tid= Impl::init_flat_triangle(h,vpm,mesh);
      std::array<Vector_2,3> flat_nei= Impl::unfold_face(h,vpm,mesh,flat_tid);

      Vector_2 c0=0.33*(flat_tid[0]+flat_tid[1]+flat_tid[2]);
      Vector_2 c1=0.33*(flat_nei[0]+flat_nei[1]+flat_nei[2]);

      return sqrt((c1 - c0).squared_length());
    };

    Dual_geodesic_solver<FT> solver;
    solver.graph.resize(faces(mesh).size());
    for (auto f : faces(mesh)) {
      halfedge_descriptor h=halfedge(f,mesh);
      int entry=get(tidmap,f);
      for (auto i = 0; i < 3; ++i) {
        solver.graph[entry][i].node = get(tidmap,face(opposite(h,mesh),mesh));
        solver.graph[entry][i].len = compute_dual_weights(h);
        h=next(h,mesh);
      }
    }
    return solver;
  }

  // `update` is a function that is executed during expansion, every time a node
  // is put into queue. `exit` is a function that tells whether to expand the
  // current node or perform early exit.
  template <typename Update, typename Stop, typename Exit>
  static
  void visit_geodesic_graph(std::vector<double> &field, const geodesic_solver &solver,
                            const std::vector<int> &sources, Update &&update,
                            Stop &&stop, Exit &&exit)
  {
    /*
      This algortithm uses the heuristic Small Label Fisrt and Large Label Last
      https://en.wikipedia.org/wiki/Shortest_Path_Faster_Algorithm

      Large Label Last (LLL): When extracting nodes from the queue, pick the
      front one. If it weights more than the average weight of the queue, put
      on the back and check the next node. Continue this way.
      Sometimes average_weight is less than every value due to floating point
      errors (doesn't happen with double precision).

      Small Label First (SLF): When adding a new node to queue, instead of
      always pushing it to the end of the queue, if it weights less than the
      front node of the queue, it is put on front. Otherwise the node is put at
      the end of the queue.
    */

    auto in_queue = std::vector<bool>(solver.graph.size(), false);

    // Cumulative weights of elements in queue. Used to keep track of the
    // average weight of the queue.
    auto cumulative_weight = 0.0;

    // setup queue
    auto queue = std::deque<int>();
    for (auto source : sources) {
      in_queue[source] = true;
      cumulative_weight += field[source];
      queue.push_back(source);
    }

    while (!queue.empty()) {
      auto node = queue.front();
      auto average_weight = (cumulative_weight / queue.size());

      // Large Label Last (see comment at the beginning)
      for (std::size_t tries = 0; tries < queue.size() + 1; tries++) {
        if (field[node] <= average_weight)
          break;
        queue.pop_front();
        queue.push_back(node);
        node = queue.front();
      }

      // Remove node from queue.
      queue.pop_front();
      in_queue[node] = false;
      cumulative_weight -= field[node];

      // Check early exit condition.
      if (exit(node))
        break;
      if (stop(node))
        continue;

      for (auto i = 0; i < (int)solver.graph[node].size(); i++) {
        // Distance of neighbor through this node
        double new_distance = field[node] + solver.graph[node][i].len;

        auto neighbor = solver.graph[node][i].node;

        auto old_distance = field[neighbor];
        if (new_distance >= old_distance)
          continue;

        if (in_queue[neighbor]) {
          // If neighbor already in queue, don't add it.
          // Just update cumulative weight.
          cumulative_weight += new_distance - old_distance;
        } else {
          // If neighbor not in queue, add node to queue using Small Label
          // First (see comment at the beginning).
          if (queue.empty() || (new_distance < field[queue.front()]))
            queue.push_front(neighbor);
          else
            queue.push_back(neighbor);

          // Update queue information.
          in_queue[neighbor] = true;
          cumulative_weight += new_distance;
        }

        // Update distance of neighbor.
        field[neighbor] = new_distance;
        update(node, neighbor, new_distance);
      }
    }
  }

  template <typename Update, typename Stop, typename Exit, typename FT>
  static
  void visit_dual_geodesic_graph(std::vector<double> &field,
                                const Dual_geodesic_solver<FT> &solver,
                                const std::vector<int> &sources,
                                Update &&update,
                                Stop &&stop,
                                Exit &&exit)
  {
    /*
      This algortithm uses the heuristic Small Label Fisrt and Large Label Last
      https://en.wikipedia.org/wiki/Shortest_Path_Faster_Algorithm

      Large Label Last (LLL): When extracting nodes from the queue, pick the
      front one. If it weights more than the average weight of the queue, put
      on the back and check the next node. Continue this way.
      Sometimes average_weight is less than every value due to floating point
      errors (doesn't happen with double precision).

      Small Label First (SLF): When adding a new node to queue, instead of
      always pushing it to the end of the queue, if it weights less than the
      front node of the queue, it is put on front. Otherwise the node is put at
      the end of the queue.
    */

    auto in_queue = std::vector<bool>(solver.graph.size(), false);

    // Cumulative weights of elements in queue. Used to keep track of the
    // average weight of the queue.
    auto cumulative_weight = 0.0;

    // setup queue
    auto queue = std::deque<int>();
    for (auto source : sources) {
      in_queue[source] = true;
      cumulative_weight += field[source];
      queue.push_back(source);
    }

    while (!queue.empty()) {
      auto node = queue.front();
      auto average_weight = (cumulative_weight / queue.size());

      // Large Label Last (see comment at the beginning)
      for (std::size_t tries = 0; tries < queue.size() + 1; tries++) {
        if (field[node] <= average_weight)
          break;
        queue.pop_front();
        queue.push_back(node);
        node = queue.front();
      }

      // Remove node from queue.
      queue.pop_front();
      in_queue[node] = false;
      cumulative_weight -= field[node];

      // Check early exit condition.
      if (exit(node))
        break;
      if (stop(node))
        continue;

      for (auto i = 0; i < (int)solver.graph[node].size(); i++) {
        // Distance of neighbor through this node
        auto new_distance = field[node] + solver.graph[node][i].len;
        auto neighbor = solver.graph[node][i].node;

        auto old_distance = field[neighbor];
        if (new_distance >= old_distance)
          continue;

        if (in_queue[neighbor]) {
          // If neighbor already in queue, don't add it.
          // Just update cumulative weight.
          cumulative_weight += new_distance - old_distance;
        } else {
          // If neighbor not in queue, add node to queue using Small Label
          // First (see comment at the beginning).
          if (queue.empty() || (new_distance < field[queue.front()]))
            queue.push_front(neighbor);
          else
            queue.push_back(neighbor);

          // Update queue information.
          in_queue[neighbor] = true;
          cumulative_weight += new_distance;
        }

        // Update distance of neighbor.
        field[neighbor] = new_distance;
        update(node, neighbor, new_distance);
      }
    }
  }

  static
  std::vector<double>
  compute_geodesic_distances(const geodesic_solver &solver,
                             const VertexIndexMap& vidmap,
                             const std::vector<std::pair<vertex_descriptor, double>> &sources_and_dist)
  {
    auto update = [](int /* node */, int /* neighbor */, double /* new_distance */) {};
    auto stop = [](int /* node */) { return false; };
    auto exit = [](int /* node */) { return false; };

    auto distances = std::vector<double>{};
    distances.assign(solver.graph.size(), DBL_MAX);
    std::vector<int>sources_id((sources_and_dist.size()));
    for (std::size_t i = 0; i < sources_and_dist.size(); ++i) {
      sources_id[i] = get(vidmap,sources_and_dist[i].first);
      distances[sources_id[i]] = sources_and_dist[i].second;
    }
    visit_geodesic_graph(distances, solver, sources_id, update, stop, exit);

    return distances;
  }

  static
  std::vector<double> solve_with_targets(const geodesic_solver& solver,
                                        const VertexIndexMap& vidmap,
                                        const std::vector<std::pair<vertex_descriptor, double>> &sources_and_dist,
                                        const std::vector<std::pair<vertex_descriptor, double>> &targets_and_dist)
  {
    auto update       = [](int node, int neighbor, float new_distance) {};
    auto stop         = [](int node) { return false; };
    double max_distance = DBL_MAX;
    std::vector<int> exit_verts(targets_and_dist.size());
    for (auto i = 0; i < targets_and_dist.size(); ++i) {
        exit_verts[i]=get(vidmap,targets_and_dist[i].first);
    }
    auto exit = [&exit_verts](int node) {
      auto it = find(exit_verts.begin(), exit_verts.end(), node);
      if (it != exit_verts.end())
        exit_verts.erase(it);

    if (exit_verts.empty())
        return true;

      return false;
    };

    auto distances  = std::vector<double>(solver.graph.size(), DBL_MAX);
    std::vector<int>sources_id((sources_and_dist.size()));
    for (auto i = 0; i < sources_and_dist.size(); ++i) {
      sources_id[i] = get(vidmap,sources_and_dist[i].first);
      distances[sources_id[i]] = sources_and_dist[i].second;
    }

    visit_geodesic_graph(distances, solver, sources_id, update, stop, exit);
    return distances;
  }

  //compute the length between two opposite vertices by flattening
  //TODO: handle concave configurations
  static
  double length_by_flattening(const VertexPointMap &vpm,
                            const TriangleMesh &mesh,
                            const halfedge_descriptor& h)
  {
    using Impl = Locally_shortest_path_imp<K, TriangleMesh, VertexPointMap>;
    std::array<Vector_2,3> flat_tid=Impl::init_flat_triangle(h,vpm,mesh);
    std::array<Vector_2,3> flat_nei=Impl::unfold_face(h,vpm,mesh,flat_tid);
    return sqrt((flat_tid[2]-flat_nei[2]).squared_length());
  }

  static
  Point_3 eval_position(const VertexPointMap &vpm,
                        const TriangleMesh &mesh,
                        const Face_location<TriangleMesh, FT>& p)
  {
    halfedge_descriptor h=halfedge(p.first,mesh);
    return CGAL::barycenter(get(vpm, source(h,mesh)), p.second[0], get(vpm, target(h,mesh)), p.second[1], get(vpm, target(next(h,mesh),mesh)), p.second[2]);
  }

  // compute the distance between a point p and some vertices around him
  // TODO: consider to take more vertices (increase accuracy)
  // TODO: handle concave configurations
  static
  std::vector<std::pair<vertex_descriptor, double>>
  nodes_around_point(const VertexPointMap &vpm,
                    const TriangleMesh &mesh,
                    const Face_location<TriangleMesh, FT>& p)
  {
    auto get_vid=[&mesh](const int k,const face_descriptor& tid)
    {
      halfedge_descriptor h=halfedge(tid,mesh);
      switch(k)
      {
        case 0:
          return source(h,mesh);
        case 1:
          return target(h,mesh);
        default:
          return target(next(h,mesh),mesh);
      }
    };

    std::vector<std::pair<vertex_descriptor, double>> nodes;
    nodes.reserve(6);
    auto [is_vert,offset]=Locally_shortest_path_imp<K,TriangleMesh, VertexPointMap>::point_is_vert(p);
    if (is_vert) {
      vertex_descriptor vid = get_vid(offset,p.first);
      nodes.push_back({vid, 0});
    } else {
      face_descriptor tid = p.first;
      Point_3 pos = eval_position(vpm, mesh, p);
      halfedge_descriptor h=halfedge(tid,mesh);
      for (auto i = 0; i < 3; ++i) {
        vertex_descriptor p0 = source(h,mesh);
        //connect to current vertex
        double d = sqrt(squared_distance(get(vpm,p0),pos));
        nodes.push_back(std::make_pair(p0, d));

        //connecting to opposite vertex w.r.t to current halfedge
        vertex_descriptor opp = target(next(opposite(h,mesh),mesh),mesh);
        double l =
            length_by_flattening(vpm,mesh,h);
        nodes.push_back(std::make_pair(opp, l));
        h=next(h,mesh);
      }
    }

    return nodes;
  }

  //compute geodesic distance field from p
  //TODO: can be easiliy extended to more than one source
  static
  std::vector<double>
  compute_geodesic_distances(const geodesic_solver& solver,
                             const VertexPointMap& vpm,
                             const VertexIndexMap& vim,
                             const TriangleMesh &mesh,
                             const Face_location<TriangleMesh, FT>& p)
  {
    std::vector<std::pair<vertex_descriptor,double>> source_nodes=nodes_around_point(vpm,mesh,p);

    return compute_geodesic_distances(solver, vim, source_nodes);
  }

  //compute the geodesic distance field from src, and stop the propagation
  //once tgt is reached
  static
  std::vector<double>
  compute_pruned_geodesic_distances(const geodesic_solver &solver,
                                    const VertexPointMap &vpm,
                                    const VertexIndexMap& vidmap,
                                    const TriangleMesh &mesh,
                                    const Face_location<TriangleMesh, FT>& src,
                                    const Face_location<TriangleMesh, FT>& tgt)
  {
    std::vector<std::pair<vertex_descriptor, double>>
    source_nodes = nodes_around_point(vpm,mesh,src);

    std::vector<std::pair<vertex_descriptor, double>>
    target_nodes = nodes_around_point(vpm,mesh,tgt);

    return solve_with_targets(solver, source_nodes, target_nodes);
  }

  template <class FT>
  static
  std::vector<halfedge_descriptor>
  strip_on_dual_graph(const Dual_geodesic_solver<FT> &solver,
                      const TriangleMesh &mesh,
                      const int src,
                      const int tgt)
  {
    if (src == tgt)
      return {};

    auto common_halfedge = [&mesh](face_descriptor f1, face_descriptor f2) {
      halfedge_descriptor h = halfedge(f1, mesh);
      for (int i = 0; i < 3; ++i) {
        if (face(opposite(h, mesh), mesh) == f2)
          return h;
        h = next(h, mesh);
      }
      CGAL_assertion(!"faces do no share a common edge");
      return halfedge_descriptor();
    };
    // initialize once for all and sparsely cleanup at the end of every solve
    std::vector<int> parents(solver.graph.size(), -1);
    std::vector<double> field(solver.graph.size(), DBL_MAX);
    std::vector<face_descriptor> id_to_face_map(faces(mesh).begin(), faces(mesh).end());

    field[src]=0.0;
    std::vector<int> sources = {src};
    auto update = [&parents](int node, int neighbor, double) {
      parents[neighbor] = node;
    };
    auto stop = [](int) { return false; };
    auto exit = [tgt](int node) { return node==tgt; };

    visit_dual_geodesic_graph(field,solver, sources, update, stop, exit);

    // extract_strip
    std::vector<halfedge_descriptor> strip;
    int node = tgt;
    CGAL_assertion(parents[tgt] != -1);
    //update the result using id_to_face_map
    strip.reserve((int)std::sqrt(parents.size()));
    while (parents[node] != -1) {
      strip.push_back(common_halfedge(id_to_face_map[node],id_to_face_map[parents[node]]));
      node = parents[node];
    }
    std::reverse(strip.begin(),strip.end());
    return strip;
  }
};

#ifdef CGAL_BSURF_USE_DIJKSTRA_SP
  class Dijkstra_end_exception : public std::exception
  {
    const char* what() const throw ()
    {
      return "Dijkstra shortest path: reached the target vertex.";
    }
  };

  template <class Graph>
  class Stop_at_target_Dijkstra_visitor : boost::default_dijkstra_visitor
  {
    typedef typename boost::graph_traits<Graph>::vertex_descriptor vertex_descriptor;
    typedef typename boost::graph_traits<Graph>::edge_descriptor edge_descriptor;

    vertex_descriptor destination_vd;

  public:
    Stop_at_target_Dijkstra_visitor(vertex_descriptor destination_vd)
      : destination_vd(destination_vd)
    { }

    void initialize_vertex(const vertex_descriptor& /*s*/, const Graph& /*mesh*/) const { }
    void examine_vertex(const vertex_descriptor& /*s*/, const Graph& /*mesh*/) const { }
    void examine_edge(const edge_descriptor& /*e*/, const Graph& /*mesh*/) const { }
    void edge_relaxed(const edge_descriptor& /*e*/, const Graph& /*mesh*/) const { }
    void discover_vertex(const vertex_descriptor& /*s*/, const Graph& /*mesh*/) const { }
    void edge_not_relaxed(const edge_descriptor& /*e*/, const Graph& /*mesh*/) const { }
    void finish_vertex(const vertex_descriptor &vd, const Graph& /* mesh*/) const
    {
      if(vd == destination_vd)
        throw Dijkstra_end_exception();
    }
  };
#endif


} // namespace internal
#ifndef CGAL_BSURF_USE_DIJKSTRA_SP
template <class FT>
struct Dual_geodesic_solver
{
  struct Edge {
    int node = -1;
    FT len = DBL_MAX;
  };
  std::vector<std::array<Edge, 3>> graph = {};
};

template <class FT, class TriangleMesh>
void init_geodesic_dual_solver(Dual_geodesic_solver<FT>& solver, const TriangleMesh& tmesh)
{
  //TODO replace with named parameter
  using VPM = typename boost::property_map<TriangleMesh, CGAL::vertex_point_t>::const_type;
  using K =  typename Kernel_traits<typename boost::property_traits<VPM>::value_type>::type;
  VPM vpm = get(CGAL::vertex_point, tmesh);
  typedef typename GetInitializedFaceIndexMap<TriangleMesh, parameters::Default_named_parameters>::const_type FIM;
  typedef typename GetInitializedVertexIndexMap<TriangleMesh, parameters::Default_named_parameters>::const_type VIM;
  const FIM fim = get_initialized_face_index_map(tmesh, parameters::default_values());

  using Impl2 = typename internal::Geodesic_circle_impl<K, TriangleMesh, VPM, VIM, FIM>;
  solver=Impl2::make_dual_geodesic_solver(vpm, fim, tmesh);
}
#endif

// points can be duplicated in the path in case you hit a vertex of the mesh,
// in order to get continuity of edge locations
template <class FT, class TriangleMesh, class EdgeLocationRange>
void locally_shortest_path(Face_location<TriangleMesh, FT> src,
                           Face_location<TriangleMesh, FT> tgt,
                           const TriangleMesh &tmesh,
                           EdgeLocationRange &edge_locations
#ifndef CGAL_BSURF_USE_DIJKSTRA_SP
                           , const Dual_geodesic_solver<FT>& solver
#endif
)
{
  typedef boost::graph_traits<TriangleMesh> BGT;
  typedef typename BGT::halfedge_descriptor halfedge_descriptor;
  typedef typename BGT::vertex_descriptor vertex_descriptor;
  typedef typename BGT::face_descriptor face_descriptor;

  // start by checking if it is not a trivial path
  if (src.first == tgt.first) return;

  auto variant_src = get_descriptor_from_location(src,tmesh);
  auto variant_tgt = get_descriptor_from_location(tgt,tmesh);

  std::vector<face_descriptor> src_visible_face;
  if (const face_descriptor* f_ptr = std::get_if<face_descriptor>(&variant_src))
  {
    src_visible_face.push_back(*f_ptr);
  }
  else
  {
    if (const halfedge_descriptor* h_ptr = std::get_if<halfedge_descriptor>(&variant_src))
    {
      if (!is_border(*h_ptr, tmesh))
        src_visible_face.push_back(face(*h_ptr, tmesh));
      if (!is_border(opposite(*h_ptr, tmesh), tmesh))
        src_visible_face.push_back(face(opposite(*h_ptr, tmesh), tmesh));
    }
    else
    {
      vertex_descriptor v = std::get<vertex_descriptor>(variant_src);
      for(halfedge_descriptor h : halfedges_around_target(v, tmesh))
      {
        if (!is_border(h, tmesh)) src_visible_face.push_back(face(h, tmesh));
      }
    }
  }
  std::sort(src_visible_face.begin(), src_visible_face.end());
  if (const face_descriptor* f_ptr = std::get_if<face_descriptor>(&variant_tgt))
  {
    if (std::find(src_visible_face.begin(), src_visible_face.end(), *f_ptr)!=src_visible_face.end())
      return;
  }
  else
  {
    auto is_visible = [&src_visible_face, &tmesh](halfedge_descriptor h)
    {
      return !is_border(h, tmesh) &&
              std::find(src_visible_face.begin(),
                        src_visible_face.end(), face(h, tmesh))!=src_visible_face.end();
    };
    if (const halfedge_descriptor* h_ptr = std::get_if<halfedge_descriptor>(&variant_tgt))
    {
      if (is_visible(*h_ptr) || is_visible(opposite(*h_ptr, tmesh))) return;
    }
    else
    {
      vertex_descriptor v = std::get<vertex_descriptor>(variant_tgt);
      for(halfedge_descriptor h : halfedges_around_target(v, tmesh))
      {
        if (is_visible(h)) return;
      }
    }
  }



  //TODO replace with named parameter
  using VPM = typename boost::property_map<TriangleMesh, CGAL::vertex_point_t>::const_type;
  using K =  typename Kernel_traits<typename boost::property_traits<VPM>::value_type>::type;
  using Impl = internal::Locally_shortest_path_imp<K, TriangleMesh, VPM>;
  VPM vpm = get(CGAL::vertex_point, tmesh);



// TODO : if (edge(vsrc, vtgt, mesh) || tgt && src on the same edge ) return;


//TODO: handle cases of src and tgt not in the same connected component (assert?)

#ifdef CGAL_BSURF_USE_DIJKSTRA_SP
  typedef typename BGT::edge_descriptor edge_descriptor;
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

 auto compute_dual_weights=[&tmesh,&vpm](const halfedge_descriptor& h)
 {
      auto flat_tid= Impl::init_flat_triangle(h,vpm,tmesh);
      auto flat_nei= Impl::unfold_face(h,vpm,tmesh,flat_tid);

      Vector_2 c0=0.33*(flat_tid[0]+flat_tid[1]+flat_tid[2]);
      Vector_2 c1=0.33*(flat_nei[0]+flat_nei[1]+flat_nei[2]);

      return sqrt((c1 - c0).squared_length());
 };
//TODO: handle boundary edges
  Dual dual(tmesh);

  // TODO: fill the weight map using something better than euclidean distance
  // TODO: the edge map could be precomputed if we know that several queries will be done
  // TODO: construct the dual graph once at the beginning and then use Dijkstra to
  //       to navigate it at every query
  for (edge_descriptor ed : edges(tmesh))
  {
    halfedge_descriptor h=halfedge(ed, tmesh), hopp=opposite(h, tmesh);

    //  distance between centroids of unfolded triangles
    put(weight_map, ed, compute_dual_weights(h));
    // Euclidean distance between centroids
    // put(weight_map, ed,
    //      sqrt(squared_distance(
    //            centroid(get(vpm, source(h, tmesh)), get(vpm, target(h, tmesh)), get(vpm, target(next(h, tmesh), tmesh))),
    //            centroid(get(vpm, source(hopp, tmesh)), get(vpm, target(hopp, tmesh)), get(vpm, target(next(hopp, tmesh), tmesh)))
    //            )));
  }

    internal::Stop_at_target_Dijkstra_visitor<Dual<TriangleMesh>> vis(tgt.first);

    try{
      boost::dijkstra_shortest_paths(dual, src.first,
                                     boost::distance_map(distance_map)
                                         .predecessor_map(predecessor_map)
                                         .weight_map(weight_map)
                                         .visitor(vis));
    }
    catch(const internal::Dijkstra_end_exception&)
    {}

  // TODO try stopping dijkstra as soon tgt is out of the queue.


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
#else
  //TODO: VIM is not needed here
  typedef typename GetInitializedVertexIndexMap<TriangleMesh, parameters::Default_named_parameters>::const_type VIM;
  typedef typename GetInitializedFaceIndexMap<TriangleMesh, parameters::Default_named_parameters>::const_type FIM;
  const FIM fim = get_initialized_face_index_map(tmesh, parameters::default_values());

  using Impl2 = typename internal::Geodesic_circle_impl<K, TriangleMesh, VPM, VIM, FIM>;

  std::vector<halfedge_descriptor> initial_path
    = (solver.graph.empty())
    ? Impl2::strip_on_dual_graph(Impl2::make_dual_geodesic_solver(vpm, fim, tmesh), tmesh, get(fim, src.first), get(fim,tgt.first))
    : Impl2::strip_on_dual_graph(solver, tmesh, get(fim, src.first), get(fim,tgt.first));
#endif

  CGAL_assertion(face(opposite(initial_path.front(), tmesh), tmesh)==src.first);
  CGAL_assertion(face(initial_path.back(), tmesh)==tgt.first);

  // remove extra halfedges from inital_path and update src/tgt to get a minimal sequence
  // in case src and/or tgt are on a vertex or an edge
#ifdef CGAL_DEBUG_BSURF
  std::size_t initial_path_size_before = initial_path.size();
#endif
  Impl::strip_path(tmesh, src, tgt, initial_path);
  if (initial_path.empty()) return;

  CGAL_assertion(face(opposite(initial_path.front(), tmesh), tmesh)==src.first);
  CGAL_assertion(face(initial_path.back(), tmesh)==tgt.first);
#ifdef CGAL_DEBUG_BSURF
  if (initial_path_size_before != initial_path.size())
  {
    std::cout << "initial_path cured: " << initial_path_size_before << " ---> " << initial_path.size() << "\n";
    for (halfedge_descriptor h : initial_path)
    {
      std::cout << "  - " << edge(h,tmesh) << "\n";
    }
    std::cout << "  src=" << src.first << " ("<< src.second[0] << "," << src.second[1] << "," << src.second[2] << ")\n";
    std::cout << "  tgt=" << tgt.first << " ("<< tgt.second[0] << "," << tgt.second[1] << "," << tgt.second[2] << ")\n";
    std::cout <<  "  Updated src/tgt " << construct_point(src, tmesh) << " | " << construct_point(tgt, tmesh) << "\n";
  }
#endif

  // here portals contains 2D coordinates of endpoints of edges in initial_path
  std::vector< std::array<typename K::Vector_2, 2>> portals=Impl::unfold_strip(initial_path,src,tgt,vpm,tmesh);
  std::size_t max_index=0;

  // lerps are barycentric coordinates of the shortest path restricted to the unfolded strip
  std::vector<FT> lerps=Impl::funnel(portals,max_index);
  // TODO: if you comment this if you don't want to shorten the path (option?).
  //       but this part is really fast so maybe does not make sense.
  Impl::straighten_path(portals,lerps,initial_path,src,tgt,vpm,tmesh,max_index);
  CGAL_assertion(lerps.size()==initial_path.size());

  edge_locations.reserve(initial_path.size());
  for(std::size_t i=0; i<initial_path.size(); ++i)
  {
    edge_locations.emplace_back(initial_path[i], make_array(lerps[i], 1.-lerps[i]));
  }
}


//TODO: do we want to also have a way to return Bezier segments?
//      the output is actually bezier segments subdivided.
template <class TriangleMesh, class FT>
std::vector<Face_location<TriangleMesh, FT>>
recursive_de_Casteljau(const TriangleMesh& mesh,
                       const Bezier_segment<TriangleMesh, FT>& control_points,
                       const int num_subdiv
#ifndef CGAL_BSURF_USE_DIJKSTRA_SP
                       , const Dual_geodesic_solver<FT>& solver = Dual_geodesic_solver<FT>()
#endif
                       )
{
  //TODO replace with named parameter
  using VPM = typename boost::property_map<TriangleMesh, CGAL::vertex_point_t>::const_type;
  using K =  typename Kernel_traits<typename boost::property_traits<VPM>::value_type>::type;
  using Impl = internal::Bezier_tracing_impl<K, TriangleMesh, VPM>;

  // init solver if empty
  const Dual_geodesic_solver<FT>* solver_ptr=&solver;
  Dual_geodesic_solver<FT> local_solver;
  if (solver.graph.empty())
  {
    solver_ptr = &local_solver;
    init_geodesic_dual_solver(local_solver, mesh);
  }

  std::vector<Bezier_segment<TriangleMesh, FT>> segments(1,control_points);
  std::vector<Bezier_segment<TriangleMesh, FT>> result;
  for (auto subdivision = 0; subdivision < num_subdiv; ++subdivision)
  {
    result.clear();
    result.reserve(segments.size() * 2);
    for (std::size_t i = 0; i < segments.size(); ++i)
    {
      auto [split0, split1] = Impl::subdivide_bezier_polygon(mesh, segments[i], 0.5, *solver_ptr);
      result.push_back(split0);
      result.push_back(split1);
    }

    CGAL_assertion( 2*segments.size()==result.size() );

    std::swap(segments, result);
  }

  std::vector<Face_location<TriangleMesh, FT>> final_result;
  final_result.reserve(result.size() * 3 + 1);
  final_result.push_back(std::move(result.front()[0]));
  for (Bezier_segment<TriangleMesh, FT>& bs : result)
  {
    for (int i=1;i<4;++i)
      final_result.push_back(std::move(bs[i]));
  }

#if 0
  // nasty trick to build the vector from a pair of iterators
  // using the fact that data in array and vector are contiguous
  return {(Face_location<TriangleMesh, FT>*)segments.data(),
          (Face_location<TriangleMesh, FT>*)segments.data() + segments.size() * 4};
#endif
  return final_result;
}

template <class FT, class TriangleMesh, class VertexDistanceMap>
void approximate_geodesic_distance_field(const Face_location<TriangleMesh, FT>& center,
                                         VertexDistanceMap distance_map,
                                         const TriangleMesh& tmesh)
{
  // TODO: the solver could be init once and used several times for different centers
  //       in particular, it can be tweaked to compute the Voronoi diagram of the initial centers
  //       or geodesic furthest point sampling.
  // TODO: add a parameter for the link size to increase to precision of the approximation of the distance
  //       that is you increase the size of the neighborhood of each vertex and you connect in the graph each vertex to its neighbors
  //       with weight being the geodesic shortest path.
  //       the more neighbors you have, the better is the approximation
  // TODO: graph construction can be done in parallel
  // TODO: concave flattening should be handled to improve the approximation of the distance
  //       (shortest path is not a line in that case)

  //TODO replace with named parameter
  using VPM = typename boost::property_map<TriangleMesh, CGAL::vertex_point_t>::const_type;
  VPM vpm = get(CGAL::vertex_point, tmesh);
  using K =  typename Kernel_traits<typename boost::property_traits<VPM>::value_type>::type;

  typedef typename GetInitializedVertexIndexMap<TriangleMesh, parameters::Default_named_parameters>::const_type VIM;
  const VIM vim = get_initialized_vertex_index_map(tmesh, parameters::default_values());
  typedef typename GetInitializedFaceIndexMap<TriangleMesh, parameters::Default_named_parameters>::const_type FIM;

  using Impl = typename internal::Geodesic_circle_impl<K, TriangleMesh, VPM, VIM, FIM>;

  typename Impl::geodesic_solver solver = Impl::make_geodesic_solver(vpm, vim,tmesh);
  std::vector<double> distances = Impl::compute_geodesic_distances(solver, vpm, vim, tmesh, center);

  for (typename Impl::vertex_descriptor v : vertices(tmesh))
  {
    put(distance_map, v, distances[get(vim,v)]);
  }
}

// TODO: center is put is output, shall we keep it like that?
template <class K, class TriangleMesh>
std::vector<Face_location<TriangleMesh, typename K::FT>>
straightest_geodesic(const Face_location<TriangleMesh, typename K::FT> &src,
                     const typename K::Vector_2& dir,
                     const typename K::FT len,
                     const TriangleMesh &tmesh)
{
  //TODO replace with named parameter
  using VPM = typename boost::property_map<TriangleMesh, CGAL::vertex_point_t>::const_type;
  using Impl = internal::Locally_shortest_path_imp<K, TriangleMesh, VPM>;
  VPM vpm = get(CGAL::vertex_point, tmesh);


  return Impl::straightest_geodesic(src, tmesh, vpm, dir, len);
}

// we don't expect the last point to be duplicated here
// TODO: rename the function, it's not a polygon necessarily
template <class K, class PointRange_2>
std::vector<std::pair<typename K::FT, typename K::FT>>
convert_polygon_to_polar_coordinates(const PointRange_2& polygon,
                                     std::optional<typename K::Point_2> center = std::nullopt)
{
  std::vector<std::pair<typename K::FT, typename K::FT>> result;
  // compute naively the center
  if (center==std::nullopt)
  {
    Bbox_2 bb = bbox_2(polygon.begin(), polygon.end());
    center = typename K::Point_2((bb.xmax()+bb.xmin())/2., (bb.ymax()+bb.ymin())/2);
  }

  bool is_closed = polygon.front()==polygon.back();
  std::size_t nbp=polygon.size();
  if (is_closed) --nbp;
  for (std::size_t i=0; i<nbp; ++i)
  {
    const typename K::Point_2& pt = polygon[i];
    typename K::FT d = approximate_sqrt(squared_distance(*center, pt));
    typename K::FT polar = std::atan2((pt.y()-center->y())/* /d */, (pt.x()-center->x())/* /d */);
    result.emplace_back(d, polar);
    // std::cout << center->x()+d*std::cos(polar) << " " << center->x()+d*std::sin(polar) << " 0\n";
  }
  if (is_closed) result.push_back(result.front());

  return result;
}

//TODO: do we want to handle polylines?
// href of center.first is the y axis in 2D
template <class K, class TriangleMesh>
std::vector<Face_location<TriangleMesh,typename K::FT>>
trace_geodesic_polygon(const Face_location<TriangleMesh, typename K::FT> &center,
                       const std::vector<typename K::Vector_2>& directions,
                       const std::vector<typename K::FT>& lengths,
                       const TriangleMesh &tmesh,
                       const Dual_geodesic_solver<typename K::FT>& solver = {})
{
  std::size_t n=directions.size();
  std::vector<Face_location<TriangleMesh,typename K::FT>> result;
  std::vector<Face_location<TriangleMesh, typename K::FT>> vertices(n);
#ifdef CGAL_DEBUG_BSURF
  std::cout << "trace_geodesic_polygon\n";
#endif
  for(std::size_t i=0;i<n;++i)
  {
    vertices[i]= straightest_geodesic<K>(center,directions[i],lengths[i],tmesh).back();
#ifdef CGAL_DEBUG_BSURF
    std::cout << construct_point(vertices[i], tmesh) << "\n";
#endif
  }
  std::vector<Edge_location<TriangleMesh,typename K::FT>> edge_locations;

  const Dual_geodesic_solver<typename K::FT>* solver_ptr=&solver;
  Dual_geodesic_solver<typename K::FT> local_solver;
  if (solver.graph.empty())
  {
    solver_ptr = &local_solver;
    init_geodesic_dual_solver(local_solver, tmesh);
  }

  for(std::size_t i=0;i<n;++i)
  {
    edge_locations.clear();
    locally_shortest_path<typename K::FT>(vertices[i],vertices[(i+1)%n],tmesh, edge_locations, *solver_ptr);
    result.push_back(vertices[i]);
    for(auto& el : edge_locations)
    {
      result.push_back(to_face_location(el, tmesh));
    }
  }
  result.push_back(vertices[0]);

  return result;
}

//TODO add groups of polygons for better rendering
template <class K, class TriangleMesh>
std::vector<std::vector<Face_location<TriangleMesh,typename K::FT>>>
trace_geodesic_polygons(const Face_location<TriangleMesh, typename K::FT> &center,
                        const std::vector<std::vector<typename K::Point_2>>& polygons,
                        const typename K::FT scaling,
                        const TriangleMesh &tmesh,
                        const Dual_geodesic_solver<typename K::FT>& solver = {})
{
  using VPM = typename boost::property_map<TriangleMesh, CGAL::vertex_point_t>::const_type;
  using Impl = internal::Locally_shortest_path_imp<K, TriangleMesh, VPM>;
  VPM vpm = get(CGAL::vertex_point, tmesh);


  using halfedge_descriptor = typename boost::graph_traits<TriangleMesh>::halfedge_descriptor;
  std::vector<Bbox_2> polygon_bboxes;
  polygon_bboxes.reserve(polygons.size());
  Bbox_2 gbox;
  for (const std::vector<typename K::Point_2>& polygon : polygons)
  {
    polygon_bboxes.push_back( bbox_2(polygon.begin(), polygon.end()) );
    gbox+=polygon_bboxes.back();
  }

  const Dual_geodesic_solver<typename K::FT>* solver_ptr=&solver;
  Dual_geodesic_solver<typename K::FT> local_solver;
  if (solver.graph.empty())
  {
    solver_ptr = &local_solver;
    init_geodesic_dual_solver(local_solver, tmesh);
  }

  std::vector<std::vector<Face_location<TriangleMesh,typename K::FT>>> result(polygons.size());
  typename K::Vector_2 initial_dir(1,0);// TODO: input parameter or 2D PCA of the centers?

  for(std::size_t i=0;i<polygons.size();++i)
  {
    typename K::Vector_2 dir( (-gbox.xmin()-gbox.xmax()+polygon_bboxes[i].xmin()+polygon_bboxes[i].xmax())/2.,
                              (-gbox.ymin()-gbox.ymax()+polygon_bboxes[i].ymin()+polygon_bboxes[i].ymax())/2. );
    std::vector< Face_location<TriangleMesh, typename K::FT> > spath =
      straightest_geodesic<K>(center, dir, scaling * std::sqrt(dir.squared_length()),tmesh);
    Face_location<TriangleMesh, typename K::FT> polygon_center = spath.back();

    std::vector<Edge_location<TriangleMesh, typename K::FT>> shortest_path;
      locally_shortest_path<typename K::FT>(center, polygon_center, tmesh, shortest_path, *solver_ptr);

    // update direction
    typename K::Vector_2 v = initial_dir;
    for(std::size_t i=0;i<shortest_path.size();++i)
    {
      halfedge_descriptor h_curr = halfedge(shortest_path[i].first, tmesh);
      v = Impl::parallel_transport_f2f(h_curr, v, vpm, tmesh);
    }

    std::vector<std::pair<typename K::FT, typename K::FT>> polar_coords =
      convert_polygon_to_polar_coordinates<K>(polygons[i],
                                              typename K::Point_2((polygon_bboxes[i].xmin()+polygon_bboxes[i].xmax())/2.,
                                                                  (polygon_bboxes[i].ymin()+polygon_bboxes[i].ymax())/2.));
    if (polygons[i].front()==polygons[i].back())
      polar_coords.pop_back();

    double theta = atan2(v.y(),v.x());

    std::vector<typename K::Vector_2> directions;
    std::vector<typename K::FT> lens;
    lens.reserve(polar_coords.size());
    directions.reserve(polar_coords.size());

    for (const std::pair<double, double>& coord : polar_coords)
    {
      lens.push_back(scaling * coord.first);
      directions.emplace_back(std::cos(coord.second+theta), std::sin(coord.second+theta));
    }
    result[i] = trace_geodesic_polygon<K>(polygon_center, directions, lens, tmesh, *solver_ptr);
  }

  return result;
}

//TODO add groups of polygons for better rendering
//TODO add a function dedicated to groups of polygons (like a sentence)
//     that takes groups of polyons
template <class K, class TriangleMesh>
std::vector<std::vector<Face_location<TriangleMesh,typename K::FT>>>
trace_geodesic_label(const Face_location<TriangleMesh, typename K::FT> &center,
                     const std::vector<std::vector<typename K::Point_2>>& polygons,
                     const typename K::FT scaling,
                     const TriangleMesh &tmesh,
                     const Dual_geodesic_solver<typename K::FT>& solver = {})
{
  using VPM = typename boost::property_map<TriangleMesh, CGAL::vertex_point_t>::const_type;
  using Impl = internal::Locally_shortest_path_imp<K, TriangleMesh, VPM>;
  VPM vpm = get(CGAL::vertex_point, tmesh);
  using Point_2 = typename K::Point_2;
  using Vector_2 = typename K::Vector_2;


  std::vector<Bbox_2> polygon_bboxes;
  polygon_bboxes.reserve(polygons.size());
  Bbox_2 gbox;
  for (const std::vector<typename K::Point_2>& polygon : polygons)
  {
    polygon_bboxes.push_back( bbox_2(polygon.begin(), polygon.end()) );
    gbox+=polygon_bboxes.back();
  }

  const Dual_geodesic_solver<typename K::FT>* solver_ptr=&solver;
  Dual_geodesic_solver<typename K::FT> local_solver;
  if (solver.graph.empty())
  {
    solver_ptr = &local_solver;
    init_geodesic_dual_solver(local_solver, tmesh);
  }

  std::vector<std::vector<Face_location<TriangleMesh,typename K::FT>>> result(polygons.size());
  // Vector_2 initial_dir(1,0);// TODO: input parameter or 2D PCA of the centers?

  // 1D partition of the letters
  Point_2 c2( (gbox.xmin()+gbox.xmax())/2.,
              (gbox.ymin()+gbox.ymax())/2. );
  Point_2 left_most(gbox.xmin(), c2.y());
  Point_2 right_most(gbox.xmax(), c2.y());
  double len = (gbox.xmax()-gbox.xmin())/2.;

  std::vector< Face_location<TriangleMesh, typename K::FT> > left_path =
    straightest_geodesic<K>(center, left_most-c2, scaling * len, tmesh);
  std::vector< Face_location<TriangleMesh, typename K::FT> > right_path =
    straightest_geodesic<K>(center, right_most-c2, scaling * len, tmesh);

  CGAL_assertion(left_path.size() >=2);
  CGAL_assertion(right_path.size() >=2);

  // TODO: precompute distances along supporting curve + stop when exceeding max distance needed

  for(std::size_t i=0;i<polygons.size();++i)
  {
    Face_location<TriangleMesh, typename K::FT> polygon_center;
    double xc = (polygon_bboxes[i].xmin()+polygon_bboxes[i].xmax())/2.;

    auto get_polygon_center = [&tmesh, &vpm](const std::vector<Face_location<TriangleMesh, typename K::FT>>& path,
                                             double targetd)
    {
      // use left
      double acc=0.;
      std::size_t k=0;
      while(true)
      {
        double len = std::sqrt(squared_distance(construct_point(path[k], tmesh),
                                                construct_point(path[k+1], tmesh)));
        acc+=len;
        if (acc == targetd)
        {
          // TODO: if you land here and path[k] is on an input vertex, it might be that loc_k==loc_k1

          double theta=0;
          Face_location<TriangleMesh, typename K::FT> loc_k=path[k], loc_k1=path[k+1];
          CGAL_assertion_code(bool OK=)
          locate_in_common_face(loc_k, loc_k1, tmesh);
          CGAL_assertion(OK);
          std::array<Vector_2,3> flat_triangle =
            Impl::init_flat_triangle(halfedge(loc_k.first,tmesh),vpm,tmesh);
          Vector_2 src = loc_k.second[0]*flat_triangle[0]+loc_k.second[1]*flat_triangle[1]+loc_k.second[2]*flat_triangle[2];
          Vector_2 tgt = loc_k1.second[0]*flat_triangle[0]+loc_k1.second[1]*flat_triangle[1]+loc_k1.second[2]*flat_triangle[2];

          Vector_2 dir2 = tgt-src;
          theta = atan2(dir2.y(), dir2.x());

          return std::make_pair(path[k+1], theta);
        }

        if (acc > targetd)
        {
          double excess = acc-targetd;

          Face_location<TriangleMesh, typename K::FT> loc_k=path[k], loc_k1=path[k+1];
          CGAL_assertion_code(bool OK=)
          locate_in_common_face(loc_k, loc_k1, tmesh);
          CGAL_assertion(OK);

          Face_location<TriangleMesh, typename K::FT> polygon_center;
          polygon_center.first=loc_k.first;
          double alpha = excess/len;

          for(int ii=0; ii<3;++ii)
            polygon_center.second[ii] = loc_k.second[ii]*alpha+loc_k1.second[ii]*(1.-alpha);

          std::array<Vector_2,3> flat_triangle =
            Impl::init_flat_triangle(halfedge(polygon_center.first,tmesh),vpm,tmesh);
          Vector_2 src = loc_k.second[0]*flat_triangle[0]+loc_k.second[1]*flat_triangle[1]+loc_k.second[2]*flat_triangle[2];
          Vector_2 tgt = loc_k1.second[0]*flat_triangle[0]+loc_k1.second[1]*flat_triangle[1]+loc_k1.second[2]*flat_triangle[2];

          Vector_2 dir2 = tgt-src;
          double theta = atan2(dir2.y(), dir2.x());

          return std::make_pair(polygon_center, theta);
        }

        if (++k==path.size()-1)
        {
          Face_location<TriangleMesh, typename K::FT> loc_k=path[k-1], loc_k1=path[k];
          CGAL_assertion_code(bool OK=)
          locate_in_common_face(loc_k, loc_k1, tmesh);
          CGAL_assertion(OK);
          std::array<Vector_2,3> flat_triangle =
            Impl::init_flat_triangle(halfedge(loc_k.first,tmesh),vpm,tmesh);
          Vector_2 src = loc_k.second[0]*flat_triangle[0]+loc_k.second[1]*flat_triangle[1]+loc_k.second[2]*flat_triangle[2];
          Vector_2 tgt = loc_k1.second[0]*flat_triangle[0]+loc_k1.second[1]*flat_triangle[1]+loc_k1.second[2]*flat_triangle[2];

          Vector_2 dir2 = tgt-src;
          double theta = atan2(dir2.y(), dir2.x());

          return std::make_pair(path.back(), theta);
        }
      }
    };
    double theta;
    if (xc<c2.x())
    {
      std::tie(polygon_center, theta)=get_polygon_center(left_path, scaling * (c2.x()-xc));
      theta=CGAL_PI+theta;
    }
    else
      std::tie(polygon_center, theta)=get_polygon_center(right_path, scaling * (xc-c2.x()));

    std::vector<Edge_location<TriangleMesh, typename K::FT>> shortest_path;
      locally_shortest_path<typename K::FT>(center, polygon_center, tmesh, shortest_path, *solver_ptr);


    std::vector<std::pair<typename K::FT, typename K::FT>> polar_coords =
      convert_polygon_to_polar_coordinates<K>(polygons[i],
                                              typename K::Point_2((polygon_bboxes[i].xmin()+polygon_bboxes[i].xmax())/2.,
                                                                  (gbox.ymin()+gbox.ymax())/2.));

    //already duplicated in trace_geodesic_polygon
    if (polar_coords.front()==polar_coords.back()) polar_coords.pop_back();

    std::vector<Vector_2> directions;
    std::vector<typename K::FT> lens;
    lens.reserve(polar_coords.size());
    directions.reserve(polar_coords.size());

    for (const std::pair<double, double>& coord : polar_coords)
    {
      lens.push_back(scaling * coord.first);
      directions.emplace_back(std::cos(coord.second+theta), std::sin(coord.second+theta));
    }
    result[i] = trace_geodesic_polygon<K>(polygon_center, directions, lens, tmesh, *solver_ptr);
  }

  return result;
}

template <class K, class TriangleMesh>
typename K::FT path_length(const std::vector<Face_location<TriangleMesh,typename K::FT>>& path,
                           const TriangleMesh &tmesh)
{
  std::size_t lpath = path.size();
  if(lpath<2)
    return 0;

  using VPM = typename boost::property_map<TriangleMesh, CGAL::vertex_point_t>::const_type;
  VPM vpm = get(CGAL::vertex_point, tmesh);

  typename K::FT len(0);

  for (std::size_t i=0; i<lpath-1; ++i)
    len += sqrt(squared_distance(construct_point(path[i],tmesh),construct_point(path[i+1],tmesh)));

  return len;
}

//TODO: factorise code
template <class K, class TriangleMesh>
std::vector<std::vector<Face_location<TriangleMesh, typename K::FT>>>
trace_geodesic_label_along_curve(const std::vector<Face_location<TriangleMesh, typename K::FT>>& supporting_curve,
                                 const std::vector<std::vector<typename K::Point_2>>& polygons,
                                 const typename K::FT scaling,
                                 const typename K::FT padding,
                                 const bool is_centered,
                                 const TriangleMesh &tmesh,
                                 const Dual_geodesic_solver<typename K::FT>& solver = {})
{
  using VPM = typename boost::property_map<TriangleMesh, CGAL::vertex_point_t>::const_type;
  using Impl = internal::Locally_shortest_path_imp<K, TriangleMesh, VPM>;
  VPM vpm = get(CGAL::vertex_point, tmesh);
  using Point_2 = typename K::Point_2;
  using Vector_2 = typename K::Vector_2;
  using Face_loc = Face_location<TriangleMesh, typename K::FT>;

  std::vector<Bbox_2> polygon_bboxes;
  polygon_bboxes.reserve(polygons.size());
  Bbox_2 gbox;
  for (const std::vector<typename K::Point_2>& polygon : polygons)
  {
    polygon_bboxes.push_back( bbox_2(polygon.begin(), polygon.end()) );
    gbox+=polygon_bboxes.back();
  }

  const Dual_geodesic_solver<typename K::FT>* solver_ptr=&solver;
  Dual_geodesic_solver<typename K::FT> local_solver;
  if (solver.graph.empty())
  {
    solver_ptr = &local_solver;
    init_geodesic_dual_solver(local_solver, tmesh);
  }

  std::vector<std::vector<Face_loc>> result(polygons.size());
//  Vector_2 initial_dir(1,0);// TODO: input parameter or 2D PCA of the centers?

  // 1D partition of the letters
  Point_2 c2( (gbox.xmin()+gbox.xmax())/2.,
              (gbox.ymin()+gbox.ymax())/2. );
  Point_2 left_most(gbox.xmin(), c2.y());
  Point_2 right_most(gbox.xmax(), c2.y());

  std::vector<typename K::FT> support_len(supporting_curve.size(),0);
  for (std::size_t i=0; i<supporting_curve.size()-1; ++i)
    support_len[i+1]=support_len[i]+std::sqrt(squared_distance(construct_point(supporting_curve[i], tmesh),
                                                               construct_point(supporting_curve[i+1], tmesh)));

  CGAL_assertion( (support_len.back()>padding + scaling * (gbox.xmax()-gbox.xmin())) ||
                  (is_centered && (support_len.back()> scaling * (gbox.xmax()-gbox.xmin()))) );

  typename K::FT pad = is_centered ? (support_len.back()-(scaling*(gbox.xmax()-gbox.xmin())))/2
                                   : padding;

  auto get_polygon_center = [&tmesh, &vpm, &supporting_curve, &support_len](double targetd)
  {
    // use left
    std::size_t k=0;
    while(true) // TODO get rid of the while and std::lower_bound
    {
      double acc = support_len[k];
      if (acc == targetd)
      {
        // note: if the supporting curve goes though supporting_curve[k], that vertex might be duplicated (if it's a path).
        //       But it is not an issue for the code as the 0 lenght segments do not contribute to the distance and
        //       points of loc_k and loc_k1 are expected to be always different.

        double theta=0;
        //TODO: should pick k-1 if k==supporting_curve.size()-1
        Face_location<TriangleMesh, typename K::FT> loc_k=supporting_curve[k], loc_k1=supporting_curve[k+1];
        CGAL_assertion_code(bool OK=)
        locate_in_common_face(loc_k, loc_k1, tmesh);
        CGAL_assertion(OK);
        std::array<Vector_2,3> flat_triangle =
          Impl::init_flat_triangle(halfedge(loc_k.first,tmesh),vpm,tmesh);
        Vector_2 src = loc_k.second[0]*flat_triangle[0]+loc_k.second[1]*flat_triangle[1]+loc_k.second[2]*flat_triangle[2];
        Vector_2 tgt = loc_k1.second[0]*flat_triangle[0]+loc_k1.second[1]*flat_triangle[1]+loc_k1.second[2]*flat_triangle[2];

        Vector_2 dir2 = tgt-src;
        theta = atan2(dir2.y(), dir2.x());

        return std::make_pair(supporting_curve[k], theta);
      }

      if (acc > targetd)
      {
        double excess = acc-targetd;

        Face_location<TriangleMesh, typename K::FT> loc_k=supporting_curve[k-1], loc_k1=supporting_curve[k];
        CGAL_assertion_code(bool OK=)
        locate_in_common_face(loc_k, loc_k1, tmesh);
        CGAL_assertion(OK);

        Face_location<TriangleMesh, typename K::FT> polygon_center;
        polygon_center.first=loc_k.first;
        double alpha = excess/(support_len[k]-support_len[k-1]);

        for(int ii=0; ii<3;++ii)
          polygon_center.second[ii] = loc_k.second[ii]*alpha+loc_k1.second[ii]*(1.-alpha);

        std::array<Vector_2,3> flat_triangle =
          Impl::init_flat_triangle(halfedge(polygon_center.first,tmesh),vpm,tmesh);
        Vector_2 src = loc_k.second[0]*flat_triangle[0]+loc_k.second[1]*flat_triangle[1]+loc_k.second[2]*flat_triangle[2];
        Vector_2 tgt = loc_k1.second[0]*flat_triangle[0]+loc_k1.second[1]*flat_triangle[1]+loc_k1.second[2]*flat_triangle[2];

        Vector_2 dir2 = tgt-src;
        double theta = atan2(dir2.y(), dir2.x());

        return std::make_pair(polygon_center, theta);
      }

      if (++k==supporting_curve.size())
      {
        //TODO: shall we throw an exception instead or return false?

        Face_location<TriangleMesh, typename K::FT> loc_k=supporting_curve[k-2], loc_k1=supporting_curve[k-1];
        CGAL_assertion_code(bool OK=)
        locate_in_common_face(loc_k, loc_k1, tmesh);
        CGAL_assertion(OK);
        std::array<Vector_2,3> flat_triangle =
          Impl::init_flat_triangle(halfedge(loc_k.first,tmesh),vpm,tmesh);
        Vector_2 src = loc_k.second[0]*flat_triangle[0]+loc_k.second[1]*flat_triangle[1]+loc_k.second[2]*flat_triangle[2];
        Vector_2 tgt = loc_k1.second[0]*flat_triangle[0]+loc_k1.second[1]*flat_triangle[1]+loc_k1.second[2]*flat_triangle[2];

        Vector_2 dir2 = tgt-src;
        double theta = atan2(dir2.y(), dir2.x());

        return std::make_pair(supporting_curve.back(), theta);
      }
    }
  };

  for(std::size_t i=0;i<polygons.size();++i)
  {
    Face_location<TriangleMesh, typename K::FT> polygon_center;
    double xc = (polygon_bboxes[i].xmin()+polygon_bboxes[i].xmax())/2.;

    double theta;
    std::tie(polygon_center, theta)=get_polygon_center(scaling * (xc-gbox.xmin())+pad);

    std::vector<std::pair<typename K::FT, typename K::FT>> polar_coords =
      convert_polygon_to_polar_coordinates<K>(polygons[i],
                                              typename K::Point_2((polygon_bboxes[i].xmin()+polygon_bboxes[i].xmax())/2.,
                                                                  (gbox.ymin()+gbox.ymax())/2.));
    //already duplicated in trace_geodesic_polygon
    if (polar_coords.front()==polar_coords.back()) polar_coords.pop_back();

    std::vector<Vector_2> directions;
    std::vector<typename K::FT> lens;
    lens.reserve(polar_coords.size());
    directions.reserve(polar_coords.size());

    for (const std::pair<double, double>& coord : polar_coords)
    {
      lens.push_back(scaling * coord.first);
      directions.emplace_back(std::cos(coord.second+theta), std::sin(coord.second+theta));
    }
    result[i] = trace_geodesic_polygon<K>(polygon_center, directions, lens, tmesh, *solver_ptr);
  }

  return result;
}

template <class K, class TriangleMesh>
std::vector< std::vector<Face_location<TriangleMesh, typename K::FT>> >
trace_bezier_curves(const Face_location<TriangleMesh, typename K::FT> &center,
                    const std::vector<std::array<typename K::Vector_2, 4>>& directions,
                    const std::vector<std::array<typename K::FT, 4>>& lengths,
                    const int num_subdiv,
                    const TriangleMesh &tmesh
#ifndef CGAL_BSURF_USE_DIJKSTRA_SP
                    , const Dual_geodesic_solver<typename K::FT>& solver = {}
#endif
)
{
  using FT = typename K::FT;

  std::size_t n=directions.size();
  std::vector< std::vector<Face_location<TriangleMesh, typename K::FT>> > result(n);
  std::vector<Face_location<TriangleMesh, typename K::FT>> vertices(n);

#ifndef CGAL_BSURF_USE_DIJKSTRA_SP
  const Dual_geodesic_solver<typename K::FT>* solver_ptr=&solver;
  Dual_geodesic_solver<typename K::FT> local_solver;
  if (solver.graph.empty())
  {
    solver_ptr = &local_solver;
    init_geodesic_dual_solver(local_solver, tmesh);
  }
#endif

#ifdef CGAL_DEBUG_BSURF
  std::ofstream debug_cp("/tmp/control_points.xyz");
  std::ofstream debug_ep("/tmp/end_points.xyz");
  debug_cp << std::setprecision(17);
  debug_ep << std::setprecision(17);
#endif
  for (std::size_t i=0; i<n; ++i)
  {
    Bezier_segment<TriangleMesh, FT> control_loc;
    for (int k=0;k<4; ++k)
    {
      control_loc[k] = straightest_geodesic<K>(center,directions[i][k],lengths[i][k],tmesh).back();
    }

    #ifdef CGAL_DEBUG_BSURF
      debug_ep << construct_point(control_loc[0], tmesh) << "\n";
      debug_ep << construct_point(control_loc[3], tmesh) << "\n";
      debug_cp << construct_point(control_loc[1], tmesh) << "\n";
      debug_cp << construct_point(control_loc[2], tmesh) << "\n";
    #endif

    std::vector<Face_location<TriangleMesh, FT>> bezier =
      recursive_de_Casteljau(tmesh, control_loc, num_subdiv
#ifndef CGAL_BSURF_USE_DIJKSTRA_SP
                            , *solver_ptr
#endif
                            );

    result[i].reserve(bezier.size());
    result[i].push_back(bezier[0]);
    for(std::size_t b=0; b<bezier.size()-1; ++b)
    {
      const Face_location<TriangleMesh, FT>& loc = bezier[b];
      const Face_location<TriangleMesh, FT>& loc1 = bezier[b+1];

      // connect the two face locations with shortest path is they are in different faces
      if (loc.first!=loc1.first)
      {
        std::vector<Edge_location<TriangleMesh, FT>> edge_locations;
        locally_shortest_path<FT>(loc, loc1, tmesh, edge_locations, solver);
        result[i].reserve(result[i].size()+edge_locations.size());
        for (const Edge_location<TriangleMesh, FT>& e : edge_locations)
          result[i].push_back(to_face_location(e, tmesh));
      }
      result[i].push_back(loc1);
    }
  }

  return result;
}


template <class K, class TriangleMesh>
typename K::FT path_length(const std::vector<Edge_location<TriangleMesh,typename K::FT>>& path,
                          const Face_location<TriangleMesh, typename K::FT>& src,
                          const Face_location<TriangleMesh, typename K::FT>& tgt,
                          const TriangleMesh &tmesh)
{
  std::size_t lpath = path.size();
  if(lpath==0)
    return sqrt(squared_distance(construct_point(src,tmesh),construct_point(tgt,tmesh)));

  using VPM = typename boost::property_map<TriangleMesh, CGAL::vertex_point_t>::const_type;
  VPM vpm = get(CGAL::vertex_point, tmesh);

  typename K::FT len=sqrt(squared_distance(construct_point(src,tmesh),construct_point(path[0],tmesh)));

  for (std::size_t i=0; i<lpath-1; ++i)
    len += sqrt(squared_distance(construct_point(path[i],tmesh),construct_point(path[i+1],tmesh)));

  len+=sqrt(squared_distance(construct_point(path.back(),tmesh),construct_point(tgt,tmesh)));

  return len;
}

//TODO VNM should be optional as well as FNM
template <class K, class TriangleMesh, class VNM, class FNM, class OutputIterator>
// std::vector<typename boost::graph_traits<TriangleMesh>::vertex_descriptor>
OutputIterator
refine_mesh_along_paths(const std::vector<std::vector<Face_location<TriangleMesh, typename K::FT>>>& paths,
                        TriangleMesh& tmesh,
                        VNM vnm,
                        FNM fnm,
                        OutputIterator out)
{
  // TODO: nothing is done here to identify identical points

  using vertex_descriptor = typename boost::graph_traits<TriangleMesh>::vertex_descriptor;
  using edge_descriptor = typename boost::graph_traits<TriangleMesh>::edge_descriptor;
  using halfedge_descriptor = typename boost::graph_traits<TriangleMesh>::halfedge_descriptor;
  using face_descriptor = typename boost::graph_traits<TriangleMesh>::face_descriptor;
  using Face_loc = Face_location<TriangleMesh, typename K::FT>;
  using dscrptr_vrnt = descriptor_variant<TriangleMesh>;
  using Point_3 = typename K::Point_3;
  using EK = CGAL::Exact_predicates_exact_constructions_kernel;

  // 3 types of points: on edges, on faces, and existing vertices
  typename boost::property_map<TriangleMesh, CGAL::dynamic_face_property_t<std::size_t>>::type
    fid_map = get(CGAL::dynamic_face_property_t<std::size_t>(), tmesh, std::size_t(-1));
  typename boost::property_map<TriangleMesh, CGAL::dynamic_edge_property_t<std::size_t>>::type
    eid_map = get(CGAL::dynamic_edge_property_t<std::size_t>(), tmesh, std::size_t(-1));
  typename boost::property_map<TriangleMesh, CGAL::dynamic_vertex_property_t<std::size_t>>::type
    vid_map = get(CGAL::dynamic_vertex_property_t<std::size_t>(), tmesh, std::size_t(-1));

  std::vector< std::vector<std::pair<std::size_t, std::size_t>> > edges_per_face;
  std::vector< std::vector<std::size_t> > points_per_face;
  std::vector< std::vector<std::size_t> > points_per_edge;
  std::vector< edge_descriptor > edges_to_refine;
  std::vector< face_descriptor > faces_to_refine;
  std::vector< std::pair<std::size_t, vertex_descriptor> > input_vertices;

  std::vector<Face_loc> polyline_locations;

  auto register_point = [&](const dscrptr_vrnt& var, std::size_t id)
  {
    switch(var.index())
    {
      case 1:
      {
        halfedge_descriptor hd=std::get<halfedge_descriptor>(var);
        edge_descriptor ed = edge(hd, tmesh);
        if (get(eid_map, ed)==std::size_t(-1))
        {
          put(eid_map, ed, edges_to_refine.size());
          points_per_edge.emplace_back();
          edges_to_refine.push_back(ed);
        }
        points_per_edge[get(eid_map,ed)].push_back(id);
        return;
      }
      case 2:
      {
        face_descriptor fd=std::get<face_descriptor>(var);
        if (get(fid_map, fd)==std::size_t(-1))
        {
          put(fid_map, fd, faces_to_refine.size());
          points_per_face.emplace_back();
          edges_per_face.emplace_back();
          faces_to_refine.push_back(fd);
        }
        points_per_face[get(fid_map,fd)].push_back(id);
        return;
      }
      default:
        input_vertices.emplace_back(id, std::get<vertex_descriptor>(var));
    }
  };

  auto register_segment = [&](const dscrptr_vrnt& v1, const dscrptr_vrnt& v2,
                              const Face_loc& loc1, const Face_loc& loc2,
                              std::size_t i1, std::size_t i2)
  {
    if (v1.index()==0 && v2.index()==0) return;
    if (v1.index()==1 && v2.index()==1 && v1==v2) return;
    if (v1.index()==0 && v2.index()==1)
    {
      vertex_descriptor vd = std::get<vertex_descriptor>(v1);
      halfedge_descriptor hd = std::get<halfedge_descriptor>(v2);
      if (vd==source(hd, tmesh) || vd==target(hd, tmesh))
        return;
    }
    if (v2.index()==0 && v1.index()==1)
    {
      vertex_descriptor vd = std::get<vertex_descriptor>(v2);
      halfedge_descriptor hd = std::get<halfedge_descriptor>(v1);
      if (vd==source(hd, tmesh) || vd==target(hd, tmesh))
        return;
    }
    Face_loc copy1=loc1, copy2=loc2;
    CGAL_assertion_code(bool OK=)
    locate_in_common_face(copy1, copy2, tmesh);
    CGAL_assertion(OK);
    if (get(fid_map,copy1.first)==std::size_t(-1))
    {
      put(fid_map, copy1.first, faces_to_refine.size());
      points_per_face.emplace_back();
      edges_per_face.emplace_back();
      faces_to_refine.push_back(copy1.first);
    }
    edges_per_face[get(fid_map, copy1.first)].emplace_back(i1,i2);
  };

  for (const std::vector<Face_loc>& path : paths)
  {
    bool closed = are_locations_identical(path.front(), path.back(), tmesh);
    std::size_t nbp = path.size();

    if (closed)  nbp-=1;

    std::size_t prev_id=polyline_locations.size(), first_id=prev_id;
    dscrptr_vrnt prev_var = get_descriptor_from_location(path.front(), tmesh),
                 first_var = prev_var;

    polyline_locations.push_back( path.front() );
    register_point(prev_var, prev_id);

    for (std::size_t i=1; i<nbp; ++i)
    {
      std::size_t id=polyline_locations.size();
      polyline_locations.push_back(path[i]);
      dscrptr_vrnt var = get_descriptor_from_location(path[i], tmesh);
      register_point(var, id);
      register_segment(prev_var, var,
                       path[i-1], path[i],
                       prev_id, id);
      prev_id=id;
      prev_var=var;
    }

    if (closed)
      register_segment(prev_var, first_var,
                       path[path.size()-2], path.front(),
                       prev_id, first_id);
  }

  // pre-compute normals per face
  std::vector<typename K::Vector_3> face_normals(faces_to_refine.size());
  std::vector<std::array<vertex_descriptor, 3>> face_input_vertices(faces_to_refine.size());
  for (std::size_t fid=0; fid<faces_to_refine.size(); ++fid)
  {
    face_descriptor fd = faces_to_refine[fid];
    halfedge_descriptor hd = halfedge(fd, tmesh);
    face_normals[fid]=compute_face_normal(fd, tmesh); // TODO: vpm + compute using EK
    face_input_vertices[fid] = make_array(target(hd, tmesh),
                                          target(next(hd, tmesh), tmesh),
                                          target(prev(hd, tmesh), tmesh));
  }

  //vertices
  std::vector<EK::Point_3> exact_points(polyline_locations.size());
  std::vector<vertex_descriptor> polyline_vertices(
    polyline_locations.size(), boost::graph_traits<TriangleMesh>::null_vertex());

  // collect split point on edges
  std::vector<halfedge_descriptor> hedges_to_refine(edges_to_refine.size());
  std::vector< std::vector<std::pair<Point_3, std::size_t> > > points_on_hedges(edges_to_refine.size());
  for (std::size_t eid=0; eid<edges_to_refine.size(); ++eid)
  {
    edge_descriptor e = edges_to_refine[eid];
    halfedge_descriptor href = halfedge(e, tmesh);
    hedges_to_refine[eid]=href;

    std::vector< std::pair<typename K::FT, typename K::FT> > coordinates;
    coordinates.reserve(points_per_edge[eid].size());
    for (std::size_t pid : points_per_edge[eid])
    {
      halfedge_descriptor hd = halfedge(polyline_locations[pid].first, tmesh);
      int src_id=0;
      if (edge(hd, tmesh)!=e)
      {
        hd=next(hd, tmesh); ++src_id;
        if (edge(hd, tmesh)!=e)
        {
          hd=next(hd, tmesh);
          ++src_id;
        }
      }
      CGAL_assertion(edge(hd, tmesh)==e);
      coordinates.emplace_back(polyline_locations[pid].second[src_id], polyline_locations[pid].second[(src_id+1)%3]);
      if (hd!=href)
      {
        CGAL_assertion( opposite(hd, tmesh)==href );
        std::swap( coordinates.back().first, coordinates.back().second );
      }
    }

    // now sort coordinates
    std::vector<std::size_t> ids(coordinates.size());
    std::iota(ids.begin(), ids.end(), 0);

    std::sort(ids.begin(), ids.end(), [&coordinates](std::size_t i, std::size_t j)
                                      { return coordinates[i].first < coordinates[j].first; });

    points_on_hedges[eid].reserve(ids.size());
    for (std::size_t id : ids)
    {
      points_on_hedges[eid].emplace_back(
        construct_point(polyline_locations[points_per_edge[eid][id]], tmesh),
        points_per_edge[eid][id]);

      exact_points[points_per_edge[eid][id]] =
        construct_point(polyline_locations[points_per_edge[eid][id]], tmesh,
                        parameters::vertex_point_map(make_cartesian_converter_property_map<EK::Point_3>(tmesh.points())));

    }
  }

  CGAL::Cartesian_converter<K, EK> to_exact;

  // add new vertices per face
  for (std::size_t fid=0; fid<faces_to_refine.size(); ++fid)
  {
    typename K::Vector_3 face_normal=get(fnm, faces_to_refine[fid]);
    for(std::size_t vid : points_per_face[fid])
    {
      vertex_descriptor vd = add_vertex(tmesh);
      put(vnm, vd, face_normal);
      put(vid_map, vd, vid);
      tmesh.point(vd)=construct_point(polyline_locations[vid], tmesh);

      CGAL_assertion( polyline_vertices[vid]==boost::graph_traits<TriangleMesh>::null_vertex() );
      polyline_vertices[vid] = vd;
      exact_points[vid] = construct_point(polyline_locations[vid], tmesh,
                            parameters::vertex_point_map(make_cartesian_converter_property_map<EK::Point_3>(tmesh.points())));
    }
  }

  // set existing vertex
  for (const std::pair<std::size_t, vertex_descriptor>& id_and_vd : input_vertices)
  {
    polyline_vertices[id_and_vd.first]=id_and_vd.second;
    put(vid_map, id_and_vd.second, id_and_vd.first);
    exact_points[id_and_vd.first]=to_exact(tmesh.point(id_and_vd.second));
  }

  //Now split edges
  for (std::size_t eid=0; eid<edges_to_refine.size(); ++eid)
  {
    halfedge_descriptor href = hedges_to_refine[eid];
    CGAL_assertion( !is_border_edge(href, tmesh) ); // TODO shall we throw if we reach boundary?
    typename K::Vector_3 edge_normal = 0.5 * (get(fnm, face(href,tmesh))+get(fnm, face(opposite(href, tmesh), tmesh)));

    for (const std::pair<Point_3, std::size_t>& pt_and_id : points_on_hedges[eid])
    {
      CGAL_assertion( polyline_vertices[pt_and_id.second]==boost::graph_traits<TriangleMesh>::null_vertex() );
      href = ::CGAL::Euler::split_edge(href, tmesh);
      polyline_vertices[pt_and_id.second]=target(href, tmesh);
      // TODO: use VPM
      tmesh.point(target(href,tmesh)) = std::move(pt_and_id.first);
      put(vid_map, target(href,tmesh), pt_and_id.second);
      put(vnm, target(href,tmesh), edge_normal);
    }
  }

  // triangulate faces
  using CDT_traits = Projection_traits_3<EK>;
  using Vb = Triangulation_vertex_base_with_info_2<vertex_descriptor,CDT_traits>;
  using Fb = Constrained_triangulation_face_base_2<CDT_traits>;
  using TDS_2 = Triangulation_data_structure_2<Vb,Fb>;
  using CDT = Constrained_Delaunay_triangulation_2<CDT_traits,TDS_2>;
  using CDT_Vertex_handle = typename CDT::Vertex_handle;

  std::size_t nb_verts = polyline_vertices.size();
  polyline_vertices.resize(polyline_vertices.size()+3); // for input triangle vertices
  exact_points.resize(exact_points.size()+3);

  std::vector<CDT_Vertex_handle> id_to_cdt_vhandles(nb_verts);

  for (std::size_t fid=0; fid<faces_to_refine.size(); ++fid)
  {
    CDT_traits traits(to_exact(face_normals[fid]));
    CDT cdt(traits);

    // turn around the face to init the convex hull
    std::array<vertex_descriptor,3> tri_verts =
      make_array(face_input_vertices[fid][0], face_input_vertices[fid][1], face_input_vertices[fid][2]);

    std::array<CDT_Vertex_handle, 3> vhandles;
    vhandles[0]=cdt.insert_outside_affine_hull(to_exact(tmesh.point(tri_verts[0])));
    vhandles[1]=cdt.insert_outside_affine_hull(to_exact(tmesh.point(tri_verts[1])));
    vhandles[2] = cdt.tds().insert_dim_up(cdt.infinite_vertex(), false);
    vhandles[2]->set_point(to_exact(tmesh.point(tri_verts[2])));

    // assign ids to input triangle vertices
    std::size_t offset=nb_verts;
    for (int i=0;i<3;++i)
    {
      vhandles[i]->info()=tri_verts[i];
      if (get(vid_map, tri_verts[i])>=nb_verts) // and not simply -1 as previous triangulation could have set another dummy value
      {
        polyline_vertices[offset]=tri_verts[i];
        put(vid_map, tri_verts[i], offset);
        exact_points[offset]=vhandles[i]->point();
        ++offset;
      }
      else
      {
        id_to_cdt_vhandles[get(vid_map, tri_verts[i])]=vhandles[i];
      }
    }


    face_descriptor fd = faces_to_refine[fid];
    std::array<halfedge_descriptor,3> face_sides;
    for (int i=0; i<3; ++i)
    {
      for (halfedge_descriptor h : halfedges_around_target(face_input_vertices[fid][i], tmesh))
      {
        if (face(h, tmesh) == fd)
        {
          face_sides[i]=h;
          break;
        }
      }
    }

    // insert points on edges
    CDT_Vertex_handle vinf = cdt.infinite_vertex();
    for (int k=0; k<3; ++k)
    {
      if (next(face_sides[k], tmesh) != face_sides[(k+1)%3])
      {
        CDT_Vertex_handle src = vhandles[k], tgt = vhandles[(k+1)%3];
        halfedge_descriptor hcurr=next(face_sides[k], tmesh);
        CGAL_assertion(src->info()==source(hcurr, tmesh));
        CGAL_assertion(tgt->info()==target(face_sides[(k+1)%3], tmesh));

        CDT_Vertex_handle prev = src;
        while(hcurr!=face_sides[(k+1)%3])
        {
          typename CDT::Face_handle fh;
          CGAL_assertion_code(bool ok =)
          cdt.is_face(prev, tgt, vinf, fh);
          CGAL_assertion(ok);
          CGAL_assertion( get(vid_map, target(hcurr,tmesh)) != std::size_t(-1) );
          prev=cdt.insert_in_edge(exact_points[get(vid_map, target(hcurr,tmesh))], fh, fh->index(vinf));
          id_to_cdt_vhandles[get(vid_map, target(hcurr,tmesh))]=prev;
          prev->info() = target(hcurr,tmesh);
          cdt.restore_Delaunay(prev); // TODO maybe not each time but one global?
          CGAL_assertion(cdt.is_valid());
          hcurr=next(hcurr, tmesh);
        }
      }
    }

    // insert points in the face (TODO: we probably don't need spatial sorting)
    for (std::size_t vid : points_per_face[fid])
    {
      id_to_cdt_vhandles[vid] = cdt.insert(exact_points[vid]);
      id_to_cdt_vhandles[vid]->info()=polyline_vertices[vid];
    }



    // insert constrained edges
    for (const std::pair<std::size_t, std::size_t>& pi : edges_per_face[fid])
    {
      cdt.insert_constraint(id_to_cdt_vhandles[pi.first],id_to_cdt_vhandles[pi.second]);
    }

    // register cdt edge -> halfedge
    halfedge_descriptor hd = halfedge(fd, tmesh);
    std::map<std::pair<std::size_t, std::size_t>, halfedge_descriptor> edge_to_hedge;
    // triangle boundary first
    for (halfedge_descriptor h : halfedges_around_face(hd, tmesh))
    {
      edge_to_hedge[std::make_pair(get(vid_map,source(h,tmesh)), get(vid_map,target(h,tmesh)))]=h;
    }
    //grab edges that are not on the convex hull (these have already been created)
    for (typename CDT::Finite_edges_iterator it=cdt.finite_edges_begin();
                                             it!=cdt.finite_edges_end(); ++it)
    {
      typename CDT::Vertex_handle cdt_v0=it->first->vertex( cdt.ccw(it->second) );
      typename CDT::Vertex_handle cdt_v1=it->first->vertex( cdt.cw(it->second) );

      // consider edges not on the convex hull (not on the boundary of the face)
      // and create the corresponding halfedges
      if ( !cdt.is_infinite(it->first->vertex(it->second)) &&
           !cdt.is_infinite(cdt.mirror_vertex(it->first,it->second)) )
      {
        edge_descriptor e=add_edge(tmesh);
        halfedge_descriptor h=halfedge(e,tmesh), h_opp=opposite(h,tmesh);
        std::size_t i0=get(vid_map,cdt_v0->info()), i1=get(vid_map,cdt_v1->info());
        vertex_descriptor v0=polyline_vertices[i0], v1=polyline_vertices[i1];

        set_target(h,v0,tmesh);
        set_target(h_opp,v1,tmesh);
        set_halfedge(v0,h,tmesh);
        set_halfedge(v1,h_opp,tmesh);

        edge_to_hedge[std::make_pair(i0,i1)]=h_opp;
        edge_to_hedge[std::make_pair(i1,i0)]=h;

        if (it->first->is_constrained(it->second))
          *out++=h;
      }
    }

    //grab triangles.
    face_descriptor current_face = fd;
    typename K::Vector_3 face_normal = get(fnm, fd);
    for (typename CDT::Finite_faces_iterator it=cdt.finite_faces_begin(),
                                             it_end=cdt.finite_faces_end();;)
    {
      typename CDT::Vertex_handle cdt_v0=it->vertex(0);
      typename CDT::Vertex_handle cdt_v1=it->vertex(1);
      typename CDT::Vertex_handle cdt_v2=it->vertex(2);

      std::size_t i0=get(vid_map,cdt_v0->info()), i1=get(vid_map,cdt_v1->info()), i2=get(vid_map,cdt_v2->info());

      CGAL_assertion(edge_to_hedge.count(std::make_pair(i0,i1))!= 0);
      CGAL_assertion(edge_to_hedge.count(std::make_pair(i1,i2))!= 0);
      CGAL_assertion(edge_to_hedge.count(std::make_pair(i2,i0))!= 0);

      halfedge_descriptor h01=edge_to_hedge[std::make_pair(i0,i1)];
      halfedge_descriptor h12=edge_to_hedge[std::make_pair(i1,i2)];
      halfedge_descriptor h20=edge_to_hedge[std::make_pair(i2,i0)];

      set_next(h01,h12,tmesh);
      set_next(h12,h20,tmesh);
      set_next(h20,h01,tmesh);

      //update face halfedge
      set_halfedge(current_face,h01,tmesh);

      //update face of halfedges
      set_face(h01,current_face,tmesh);
      set_face(h12,current_face,tmesh);
      set_face(h20,current_face,tmesh);

      if ( ++it!=it_end )
      {
        current_face=add_face(tmesh);
        put(fnm, current_face, face_normal);
      }
      else
        break;
    }
  }

  return out;
}

// TODO: add to doc
// create a polyline from a path on a mesh. It does eliminate duplicated elements in the path
// that are here to make a continuous face transition when the path goes though a mesh vertex
template<class TriangleMesh, class FT, class OutputIterator>
OutputIterator
convert_path_to_polyline(const std::vector<Face_location<TriangleMesh, FT>>& path,
                         const TriangleMesh& tmesh,
                         OutputIterator poly_out)
{
  using vertex_descriptor = typename boost::graph_traits<TriangleMesh>::vertex_descriptor;

  vertex_descriptor vd = boost::graph_traits<TriangleMesh>::null_vertex();
  for (const Face_location<TriangleMesh, FT>& loc : path)
  {
    std::optional<vertex_descriptor> ov = vertex_descriptor_from_location(loc, tmesh);
    if (ov.has_value())
    {
      if (vd==*ov) continue;
      vd=*ov;
    }
    *poly_out++=construct_point(loc, tmesh);
  }
  return poly_out;
}

// same as above but for edge locations
template<class TriangleMesh, class FT, class OutputIterator>
OutputIterator
convert_path_to_polyline(const Face_location<TriangleMesh, FT>& src,
                         const std::vector<Edge_location<TriangleMesh, FT>>& path,
                         const Face_location<TriangleMesh, FT>& tgt,
                         const TriangleMesh& tmesh,
                         OutputIterator poly_out)
{
  using vertex_descriptor = typename boost::graph_traits<TriangleMesh>::vertex_descriptor;

  *poly_out++=construct_point(src, tmesh);
  vertex_descriptor vd = boost::graph_traits<TriangleMesh>::null_vertex();
  for (const Edge_location<TriangleMesh, FT>& loc : path)
  {
    std::optional<vertex_descriptor> ov = vertex_descriptor_from_location(loc, tmesh);
    if (ov.has_value())
    {
      if (vd==*ov) continue;
      vd=*ov;
    }
    *poly_out++=construct_point(loc, tmesh);
  }
  *poly_out++=construct_point(tgt, tmesh);
  return poly_out;
}

// template <class K, class TriangleMesh>
// std::vector<typename K::Vector_2>
// tangent_path_direction(const std::vector<Edge_location<TriangleMesh,typename K::FT>>& path,
//                        const Face_location<TriangleMesh, typename K::FT>& src,
//                        const Face_location<TriangleMesh, typename K::FT>& tgt,
//                        const TriangleMesh &tmesh,const bool initial=true)
// {
//   auto find = [](const std::array<typename K::FT,3> &vec, int x) {
//     for (int i = 0; i < vec.size(); i++)
//       if (vec[i] == x)
//         return i;
//     return -1;
//   };
//  typename K::Vector_2 direction;
//   using VPM = typename boost::property_map<TriangleMesh, CGAL::vertex_point_t>::const_type;
//   VPM vpm = get(CGAL::vertex_point, tmesh);


//  if(initial)
//  {

//     halfedge_descriptor h_ref=halfedge(src.first,mesh);
//     std::array<Vector_2, 3> flat_tid=init_flat_triangle(h_ref,vpm,tmesh);
//     if(path.size()==0)
//     {
//       assert(src.first==tgt.first);//TODO:src and tgt may have different faces because we do not update them when cleaning the strip
//       typename K::Vector_2 flat_src=src.second[0]*flat_tid[0]+src.second[1]*flat_tid[1]+src.second[2]*flat_tid[2];
//       typename K::Vector_2 flat_tgt=tgt.second[0]*flat_tid[0]+tgt.second[1]*flat_tid[1]+tgt.second[2]*flat_tid[2];
//       direction=normalize(flat_tgt-flat_src);
//     }else{
//     halfedge_descriptor h_edge=halfedge(path[0].first,tmesh);
//     int k=0;

//     if(h_edge==prev(h_ref,mesh))
//       k=2;
//     else if(h_edge==next(h_ref,mesh))
//       k=1;
//     else
//       assert(h_edge==h_ref);
//     }










//  }






// }

template <class K, class TriangleMesh>
std::vector<typename K::Point_3>
trace_agap_polygon(const Face_location<TriangleMesh, typename K::FT> &center,
                     const std::vector<typename K::Vector_2>& directions,
                     const std::vector<typename K::FT>& lengths,
                     const TriangleMesh &tmesh,
                     const Dual_geodesic_solver<typename K::FT>& solver = {})
{
  size_t n=directions.size();
  std::vector<typename K::Point_3> result;
  std::vector<Face_location<TriangleMesh, typename K::FT>> vertices(n);
  for(std::size_t i=0;i<n;++i)
    vertices[i]= straightest_geodesic<K>(center,directions[i],lengths[i],tmesh).back();

  std::vector<Edge_location<TriangleMesh,typename K::FT>> edge_locations;

  const Dual_geodesic_solver<typename K::FT>* solver_ptr=&solver;
  Dual_geodesic_solver<typename K::FT> local_solver;
  if (solver.graph.empty())
  {
    solver_ptr = &local_solver;
    init_geodesic_dual_solver(local_solver, tmesh);
  }

  for(std::size_t i=0;i<n;++i)
  {
    edge_locations.clear();
    locally_shortest_path<typename K::FT>(vertices[i],vertices[(i+1)%n],tmesh, edge_locations, *solver_ptr);
    result.push_back(construct_point(vertices[i],tmesh));
    for(auto& el : edge_locations)
    {
      result.push_back(construct_point(el, tmesh));
    }
  }
  result.push_back(construct_point(vertices[0],tmesh));

  return result;

}

} // namespace Polygon_mesh_processing
} // namespace CGAL

#endif
