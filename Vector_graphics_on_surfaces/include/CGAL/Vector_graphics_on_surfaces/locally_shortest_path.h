// Copyright (c) 2023-2026 GeometryFactory and Claudio Mancinelli.
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Claudio Mancinelli and Sébastien Loriot

#ifndef CGAL_POLYGON_MESH_PROCESSING_BSURF_LOCALLY_SHORTEST_PATH_H
#define CGAL_POLYGON_MESH_PROCESSING_BSURF_LOCALLY_SHORTEST_PATH_H

#include <CGAL/license/Vector_graphics_on_surfaces.h>

#include <CGAL/Named_function_parameters.h>
#include <CGAL/Vector_graphics_on_surfaces/internal/utils.h>
#include <CGAL/Vector_graphics_on_surfaces/utils.h>
#include <CGAL/Vector_graphics_on_surfaces/Dual_geodesic_solver.h>
#include <CGAL/boost/graph/named_params_helper.h>
#include <Eigen/Dense>
#include <boost/graph/graph_traits.hpp>

namespace CGAL {
namespace Vector_graphics_on_surfaces {

#ifndef DOXYGEN_RUNNING
template <class FT, class TriangleMesh, class EdgeLocationRange>
void locally_shortest_path(CGAL::Polygon_mesh_processing::Face_location<TriangleMesh, FT> src,
                           CGAL::Polygon_mesh_processing::Face_location<TriangleMesh, FT> tgt,
                           const TriangleMesh &tmesh,
                           EdgeLocationRange &edge_locations
#ifndef CGAL_BSURF_USE_DIJKSTRA_SP
                           , const Dual_geodesic_solver<FT>& solver = Dual_geodesic_solver<FT>()
#endif
);
#endif

/*!
 * \ingroup VGSMiscellaneous
 * array containing the locations of the four control points of a Bézier segment on a triangle mesh
 * `init_geodesic_dual_solver()`.
 * \tparam FT floating point number type (float or double)
 * \tparam TriangleMesh a model of `FaceListGraph` and `EdgeListGraph`
 */
template <class TriangleMesh, class FT>
using Bezier_segment = std::array<CGAL::Polygon_mesh_processing::Face_location<TriangleMesh, FT>, 4>;

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
                 const CGAL::Polygon_mesh_processing::Face_location<TriangleMesh,FT>& src,
                 const CGAL::Polygon_mesh_processing::Face_location<TriangleMesh,FT>& tgt,
                 const TriangleMesh& mesh)
  {
    namespace PMP = CGAL::Polygon_mesh_processing;
    static int i = -1;
    std::cout << "dump current path in path_"+std::to_string(i+1)+".polylines.txt\n";
    std::ofstream out("path_"+std::to_string(++i)+".polylines.txt");
    out << path.size()+2 << " " << PMP::construct_point(src, mesh);
    for(std::size_t i=0; i<path.size(); ++i)
      out << " " << PMP::construct_point(PMP::Edge_location<TriangleMesh, FT>(path[i], make_array(lerps[i], 1.-lerps[i])), mesh);
    out << " " << PMP::construct_point(tgt, mesh) << "\n";
  }
#endif

  static std::array<Vector_2, 2>
  init_source_triangle(halfedge_descriptor hopp,
                       const VertexPointMap &vpm,
                       const TriangleMesh &mesh,
                       CGAL::Polygon_mesh_processing::Face_location<TriangleMesh, FT> src)
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
    tr2d[2] = intersect_circles<K>(tr2d[0], rx, tr2d[1], ry);


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
                       CGAL::Polygon_mesh_processing::Face_location<TriangleMesh, FT> tgt)
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
    tr2d[2] = intersect_circles<K>(tr2d[0], rx, tr2d[1], ry);


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

  static
  std::vector< std::array<Vector_2, 2>>
  unfold_strip(const std::vector<halfedge_descriptor>& initial_path,
              const CGAL::Polygon_mesh_processing::Face_location<TriangleMesh, FT>& src,
              const CGAL::Polygon_mesh_processing::Face_location<TriangleMesh, FT>& tgt,
              const VertexPointMap &vpm, const TriangleMesh &mesh)
  {
    std::size_t s=initial_path.size();
    std::vector<std::array<Vector_2, 2>> result(s+1); // Vector_2 should be Point_2 as they are funnel endpoints in 2D (after flattening)
    result[0]=init_source_triangle(initial_path[0], vpm, mesh, src);
#ifdef CGAL_DEBUG_BSURF
    std::cout << "unfolding faces\n";
#endif
    for(std::size_t i=1;i<s;++i)
      result[i]=unfold_face<K>(initial_path[i-1],initial_path[i], vpm, mesh, result[i-1]);

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
                       CGAL::Polygon_mesh_processing::Face_location<TriangleMesh, FT>& src,
                       CGAL::Polygon_mesh_processing::Face_location<TriangleMesh, FT>& tgt,
                       const VertexPointMap &vpm, const TriangleMesh &mesh, std::size_t index)
  {
    namespace PMP = CGAL::Polygon_mesh_processing;

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
      std::cout << "src = " << PMP::construct_point(src, mesh) << "\n";
      std::cout << "tgt = " << PMP::construct_point(tgt, mesh) << "\n";
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

    std::array<Vector_2,3> flat_from = init_flat_triangle<K>(h,vpm,mesh);
     std::array<Vector_2,3> flat_to = unfold_face<K>(h,vpm,mesh,flat_from);
    Vector_2 bary = Vector_2{0.333, 0.333};
    Vector_2 c0 = 0.33*(flat_from[0]+flat_from[1]+flat_from[2]);
    Vector_2 c1 = 0.33*(flat_to[0]+flat_to[1]+flat_to[2]);
    Vector_2 e0 = flat_from[0] - c0;
    Vector_2 e1 = flat_to[0] - c1;

    Vector_2 w = c1 - c0;
    FT phi_ij = angle(e0, w);
    if (e0.x()*w.y()-e0.y()*w.x() < 0)
      phi_ij = 2 * CGAL_PI - phi_ij;
    w *= -1;
    FT phi_ji = angle(e1, w);

    if (e1.x()*w.y()-e1.y()*w.x() < 0)
      phi_ji = 2 * CGAL_PI - phi_ji;

    Vector_3 n_from= face_normal(vpm, mesh, from);
    std::vector<Vector_3> e_from = polar_basis(vpm, mesh, from);
    double teta = angle(e_from, v);
    if (cross_product(e_from, v)*n_from < 0)
      teta = 2 * CGAL_PI - teta;

    std::vector<Vector_3> e_to = polar_basis(vpm, mesh, from ,to);
    Vector_3 n_to = face_normal(vpm, mesh, to);
    double rot = teta + phi_ji + CGAL_PI - phi_ij;
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

    std::array<Vector_2,3> flat_curr = init_flat_triangle<K>(h_ref,vpm,mesh);
    std::array<Vector_2,3> flat_prev = unfold_face<K>(h_edge,vpm,mesh,flat_curr,k);

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

  Vector_3 parallel_transport_along_path(const std::vector<CGAL::Polygon_mesh_processing::Edge_location<TriangleMesh, FT>>& edge_locations,
                                         const VertexPointMap &vpm,
                                         const TriangleMesh& mesh,
                                         const CGAL::Polygon_mesh_processing::Face_location<TriangleMesh, FT>& src,
                                         const CGAL::Polygon_mesh_processing::Face_location<TriangleMesh, FT>& tgt,
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
#if 0
  static
  std::vector<Vector_3>
  get_3D_basis_at_point(const CGAL::Polygon_mesh_processing::Face_location<TriangleMesh, FT>& p,
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
#endif

/////////////////////////////////////////
/////////////////////////////////////////
/////////////////////////////////////////
/////////////////////////////////////////

  static
  void
  strip_path(const TriangleMesh& tmesh,
             CGAL::Polygon_mesh_processing::Face_location<TriangleMesh, FT>& src,
             CGAL::Polygon_mesh_processing::Face_location<TriangleMesh, FT>& tgt,
             std::vector<halfedge_descriptor>& initial_path)
  {
    namespace PMP = CGAL::Polygon_mesh_processing;

    // retrieve the vertex id of a face location describing a vertex
    auto is_vertex = [](const PMP::Face_location<TriangleMesh, FT>& fl)
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
    auto is_edge = [](const PMP::Face_location<TriangleMesh, FT>& fl)
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

} // namespace internal

/*!
 * \ingroup VGSFunctions
 * computes an approximated geodesic shortest path between two locations on
 * `tmesh`. The points `src` and `tgt` must be on the same connected component.
 * \tparam TriangleMesh a model of `FaceListGraph` and `EdgeListGraph`
 * \tparam FT floating point number type (float or double)
 * \tparam EdgeLocationRange a model of `BackInsertionSequence` whose value type `CGAL::Polygon_mesh_processing::Edge_location<FT>`.
 * \param src source of the path
 * \param tgt target of the path
 * \param tmesh input triangle mesh to compute the path on
 * \param edge_locations contains the path as a sequence of edge locations.
 *                       Two consecutive edges `e_k` and `e_kp1` stored in `edge_locations` are such that
 *                       `face(halfedge(e_k, tmesh), tmesh) == face(opposite(halfedge(e_kp1, tmesh), tmesh), tmesh))`.
 *                       In particular, it means that if the path goes through a vertex of `tmesh`, several
 *                       edge locations will be reported to maintain this property.
 *                       Additionally, if `src` is in the interior of a face `f`, then the first edge location `e_0`
 *                       of `edge_locations` is such that `f == face(opposite(halfedge(e_0, tmesh), tmesh), tmesh))`.
 *                       Similarly, if `tgt` is in the interior of a face `f`, then the last edge location `e_n`
 *                       of `edge_locations` is such that `f == face(halfedge(e_n, tmesh), tmesh)`.
 * \param solver container for the precomputed information. If not initialized, it will be initialized internally.
 * \todo add named parameters
 * \todo should we have halfedge location instead?
 */
template <class FT, class TriangleMesh, class EdgeLocationRange>
void locally_shortest_path(CGAL::Polygon_mesh_processing::Face_location<TriangleMesh, FT> src,
                           CGAL::Polygon_mesh_processing::Face_location<TriangleMesh, FT> tgt,
                           const TriangleMesh &tmesh,
                           EdgeLocationRange &edge_locations
#ifndef CGAL_BSURF_USE_DIJKSTRA_SP
                           , const Dual_geodesic_solver<FT>& solver
#endif
)
{
  namespace PMP = CGAL::Polygon_mesh_processing;

  typedef boost::graph_traits<TriangleMesh> BGT;
  typedef typename BGT::halfedge_descriptor halfedge_descriptor;
  typedef typename BGT::vertex_descriptor vertex_descriptor;
  typedef typename BGT::face_descriptor face_descriptor;

  // start by checking if it is not a trivial path
  if (src.first == tgt.first) return;

  auto variant_src = PMP::get_descriptor_from_location(src,tmesh);
  auto variant_tgt = PMP::get_descriptor_from_location(tgt,tmesh);

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
      auto flat_tid= internal::init_flat_triangle<K>(h,vpm,tmesh);
      auto flat_nei= internal::unfold_face<K>(h,vpm,tmesh,flat_tid);

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

  // TODO try stopping Dijkstra as soon tgt is out of the queue.


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

  // remove extra halfedges from initial_path and update src/tgt to get a minimal sequence
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
    std::cout <<  "  Updated src/tgt " << PMP::construct_point(src, tmesh) << " | " << PMP::construct_point(tgt, tmesh) << "\n";
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

} // namespace Vector_graphics_on_surfaces
} // namespace CGAL

#endif // CGAL_POLYGON_MESH_PROCESSING_BSURF_LOCALLY_SHORTEST_PATH_H
