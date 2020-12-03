// Copyright (c) 2019 GeometryFactory SARL (France).
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
// Author(s)     : Simon Giraudot

#ifndef CGAL_KSR_3_DATA_STRUCTURE_H
#define CGAL_KSR_3_DATA_STRUCTURE_H

// #include <CGAL/license/Kinetic_shape_reconstruction.h>

// CGAL includes.
#include <CGAL/Delaunay_triangulation_2.h>

// Internal includes.
#include <CGAL/KSR/utils.h>
#include <CGAL/KSR/debug.h>

#include <CGAL/KSR_3/Support_plane.h>
#include <CGAL/KSR_3/Intersection_graph.h>

namespace CGAL {
namespace KSR_3 {

template<typename GeomTraits>
class Data_structure {

public:
  using Kernel = GeomTraits;

private:
  using FT          = typename Kernel::FT;
  using Point_2     = typename Kernel::Point_2;
  using Point_3     = typename Kernel::Point_3;
  using Segment_2   = typename Kernel::Segment_2;
  using Segment_3   = typename Kernel::Segment_3;
  using Vector_2    = typename Kernel::Vector_2;
  using Vector_3    = typename Kernel::Vector_3;
  using Direction_2 = typename Kernel::Direction_2;
  using Triangle_2  = typename Kernel::Triangle_2;
  using Line_2      = typename Kernel::Line_2;
  using Line_3      = typename Kernel::Line_3;
  using Plane_3     = typename Kernel::Plane_3;

  using Support_plane      = KSR_3::Support_plane<Kernel>;
  using Intersection_graph = KSR_3::Intersection_graph<Kernel>;

  using Mesh           = typename Support_plane::Mesh;
  using Vertex_index   = typename Mesh::Vertex_index;
  using Face_index     = typename Mesh::Face_index;
  using Edge_index     = typename Mesh::Edge_index;
  using Halfedge_index = typename Mesh::Halfedge_index;

  using Polygon_2 = CGAL::Polygon_2<Kernel>;

public:
  using PVertex = std::pair<KSR::size_t, Vertex_index>;
  using PFace   = std::pair<KSR::size_t, Face_index>;
  using PEdge   = std::pair<KSR::size_t, Edge_index>;

  template<typename PSimplex>
  struct Make_PSimplex {
    using argument_type = typename PSimplex::second_type;
    using result_type   = PSimplex;

    const KSR::size_t support_plane_idx;
    Make_PSimplex(const KSR::size_t sp_idx) :
    support_plane_idx(sp_idx)
    { }

    const result_type operator()(const argument_type& arg) const {
      return result_type(support_plane_idx, arg);
    }
  };

  using PVertex_iterator =
    boost::transform_iterator<Make_PSimplex<PVertex>, typename Mesh::Vertex_range::iterator>;
  using PVertices = CGAL::Iterator_range<PVertex_iterator>;

  using PFace_iterator =
    boost::transform_iterator<Make_PSimplex<PFace>, typename Mesh::Face_range::iterator>;
  using PFaces = CGAL::Iterator_range<PFace_iterator>;

  using PEdge_iterator =
    boost::transform_iterator<Make_PSimplex<PEdge>, typename Mesh::Edge_range::iterator>;
  using PEdges = CGAL::Iterator_range<PEdge_iterator>;

  struct Halfedge_to_pvertex {
    using argument_type = Halfedge_index;
    using result_type   = PVertex;

    const KSR::size_t support_plane_idx;
    const Mesh& mesh;

    Halfedge_to_pvertex(const KSR::size_t sp_idx, const Mesh& m) :
    support_plane_idx(sp_idx),
    mesh(m)
    { }

    const result_type operator()(const argument_type& arg) const {
      return result_type(support_plane_idx, mesh.target(arg));
    }
  };

  using PVertex_of_pface_iterator =
    boost::transform_iterator<Halfedge_to_pvertex, CGAL::Halfedge_around_face_iterator<Mesh> >;
  using PVertices_of_pface = CGAL::Iterator_range<PVertex_of_pface_iterator>;

  struct Halfedge_to_pedge {
    using argument_type = Halfedge_index;
    using result_type = PEdge;

    const KSR::size_t support_plane_idx;
    const Mesh& mesh;

    Halfedge_to_pedge(const KSR::size_t sp_idx, const Mesh& m) :
    support_plane_idx(sp_idx),
    mesh(m)
    { }

    const result_type operator()(const argument_type& arg) const {
      return result_type(support_plane_idx, mesh.edge(arg));
    }
  };

  using PEdge_around_pvertex_iterator =
    boost::transform_iterator<Halfedge_to_pedge, CGAL::Halfedge_around_target_iterator<Mesh> >;
  using PEdges_around_pvertex = CGAL::Iterator_range<PEdge_around_pvertex_iterator>;

  using PEdge_of_pface_iterator =
    boost::transform_iterator<Halfedge_to_pedge, CGAL::Halfedge_around_face_iterator<Mesh> >;
  using PEdges_of_pface = CGAL::Iterator_range<PEdge_of_pface_iterator>;

  using IVertex = typename Intersection_graph::Vertex_descriptor;
  using IEdge   = typename Intersection_graph::Edge_descriptor;

  struct Volume_cell {
    std::vector<PFace> pfaces;
    std::vector<int> neighbors;
    std::set<PVertex> pvertices;
    std::size_t index = std::size_t(-1);
    Point_3 centroid;

    void add_pface(const PFace& pface, const int neighbor) {
      pfaces.push_back(pface);
      neighbors.push_back(neighbor);
    }
    void set_index(const std::size_t idx) {
      index = idx;
    }
    void set_centroid(const Point_3& point) {
      centroid = point;
    }
  };

private:
  KSR::vector<Support_plane> m_support_planes;
  Intersection_graph m_intersection_graph;
  std::vector<Volume_cell> m_volumes;
  FT m_current_time;
  FT m_previous_time;

public:
  Data_structure() :
  m_current_time(FT(0)),
  m_previous_time(FT(0))
  { }

  void clear() {
    m_support_planes.clear();
    m_intersection_graph.clear();
    m_volumes.clear();
    m_current_time = FT(0);
  }

  template<typename DS>
  void convert(DS& ds) {

    ds.clear();
    ds.resize(number_of_support_planes());
    CGAL_assertion(ds.number_of_support_planes() == number_of_support_planes());

    m_intersection_graph.convert(ds.igraph());
    for (KSR::size_t i = 0; i < number_of_support_planes(); ++i) {
      m_support_planes[i].convert(m_intersection_graph, ds.support_planes()[i]);
    }
  }

  KSR::vector<Support_plane>& support_planes() {
    return m_support_planes;
  }

  Intersection_graph& igraph() {
    return m_intersection_graph;
  }

  void resize(const KSR::size_t number_of_items) {
    m_support_planes.resize(number_of_items);
  }

  // TODO: It looks like here we lose precision during the conversion because
  // KSR::size_t is usually smaller than std::size_t!
  void reserve(const std::size_t number_of_polygons) {
    m_support_planes.reserve(static_cast<KSR::size_t>(number_of_polygons) + 6);
  }

  const FT& current_time() const { return m_current_time; }

  void set_last_event_time(const PVertex& pvertex, const FT time) {
    support_plane(pvertex).set_last_event_time(pvertex.second, time);
  }

  const FT last_event_time(const PVertex& pvertex) {
    return support_plane(pvertex).last_event_time(pvertex.second);
  }

  const std::vector<Volume_cell>& polyhedrons() const {
    return m_volumes;
  }

  /*******************************
  **      SUPPORT PLANES        **
  ********************************/

  const KSR::size_t number_of_support_planes() const {
    return m_support_planes.size();
  }

  const bool is_bbox_support_plane(const KSR::size_t support_plane_idx) const {
    return (support_plane_idx < 6);
  }

  const bool is_mesh_valid(
    const bool check_simplicity,
    const bool check_convexity,
    const KSR::size_t support_plane_idx) const {

    const bool is_valid = mesh(support_plane_idx).is_valid();
    if (!is_valid) {
      return false;
    }

    // bbox faces may have multiple equal points after converting from exact to inexact!
    if (support_plane_idx < 6) {
      return true;
    }

    for (const auto pface : pfaces(support_plane_idx)) {
      std::function<Point_2(PVertex)> unary_f =
      [&](const PVertex& pvertex) -> Point_2 {
        return point_2(pvertex);
      };

      const Polygon_2 polygon(
        boost::make_transform_iterator(pvertices_of_pface(pface).begin(), unary_f),
        boost::make_transform_iterator(pvertices_of_pface(pface).end(), unary_f));

      // Use only with an exact kernel!
      if (check_simplicity && !polygon.is_simple()) {
        const std::string msg = "ERROR: pface " + str(pface) + " is not simple!";
        CGAL_assertion_msg(false, msg.c_str());
        return false;
      }

      // Use only with an exact kernel!
      if (check_convexity && !polygon.is_convex()) {
        const std::string msg = "ERROR: pface " + str(pface) + " is not convex!";
        CGAL_assertion_msg(false, msg.c_str());
        return false;
      }

      auto prev = null_pvertex();
      for (const auto pvertex : pvertices_of_pface(pface)) {
        if (prev == null_pvertex()) {
          prev = pvertex;
          continue;
        }

        if (point_2(prev) == point_2(pvertex) &&
          direction(prev) == direction(pvertex)) {

          const std::string msg = "ERROR: pface " + str(pface) +
          " has two consequent identical vertices "
          + str(prev) + " and " + str(pvertex) + "!";
          CGAL_assertion_msg(false, msg.c_str());
          return false;
        }
        prev = pvertex;
      }
    }
    return true;
  }

  void check_integrity(
    const bool check_simplicity = false,
    const bool check_convexity  = false) const {

    for (KSR::size_t i = 0; i < number_of_support_planes(); ++i) {
      if (!is_mesh_valid(check_simplicity, check_convexity, i)) {
        const std::string msg = "ERROR: mesh " + std::to_string(i) + " is not valid!";
        CGAL_assertion_msg(false, msg.c_str());
      }

      for (const auto& iedge : this->iedges(i)) {
        const auto& iplanes = this->intersected_planes(iedge);
        if (iplanes.find(i) == iplanes.end()) {

          const std::string msg = "ERROR: support_plane " + std::to_string(i) +
          " is intersected by " + str(iedge) +
          " but it claims it does not intersect it!";
          CGAL_assertion_msg(false, msg.c_str());
        }
      }
    }

    for (const auto iedge : this->iedges()) {
      const auto& iplanes = this->intersected_planes(iedge);
      for (const auto support_plane_idx : iplanes) {

        const auto& sp_iedges = this->iedges(support_plane_idx);
        if (sp_iedges.find(iedge) == sp_iedges.end()) {

          const std::string msg = "ERROR: iedge " + str(iedge) +
          " intersects support plane " + std::to_string(support_plane_idx) +
          " but it claims it is not intersected by it!";
          CGAL_assertion_msg(false, msg.c_str());
        }
      }
    }
  }

  template<typename PointRange>
  const KSR::size_t add_support_plane(const PointRange& polygon) {

    const Support_plane new_support_plane(polygon);
    KSR::size_t support_plane_idx = KSR::no_element();
    for (KSR::size_t i = 0; i < number_of_support_planes(); ++i) {
      if (new_support_plane == support_plane(i)) {
        support_plane_idx = i;
        break;
      }
    }

    if (support_plane_idx == KSR::no_element()) {
      support_plane_idx = number_of_support_planes();
      m_support_planes.push_back(new_support_plane);
    }

    // Intersect planes with the bounding box.
    if (support_plane_idx >= 6) {

      Point_3 point;
      Point_3 centroid_3 = CGAL::ORIGIN;
      std::vector< std::pair<IEdge, Point_3> > intersections;

      for (const IEdge iedge : m_intersection_graph.edges()) {
        if (!KSR::intersection(
          support_plane(support_plane_idx).plane(), segment_3(iedge), point)) {
          continue;
        }

        centroid_3 = CGAL::barycenter(
          centroid_3, static_cast<FT>(intersections.size()), point, FT(1));
        intersections.push_back(std::make_pair(iedge, point));
      }

      Point_2 centroid_2 = support_plane(support_plane_idx).to_2d(centroid_3);
      std::sort(intersections.begin(), intersections.end(),
      [&] (const std::pair<IEdge, Point_3>& a, const std::pair<IEdge, Point_3>& b) -> bool {
        const auto a2 = support_plane(support_plane_idx).to_2d(a.second);
        const auto b2 = support_plane(support_plane_idx).to_2d(b.second);
        const Segment_2 sega(centroid_2, a2);
        const Segment_2 segb(centroid_2, b2);
        return ( Direction_2(sega) < Direction_2(segb) );
      });

      KSR::vector<KSR::size_t> common_planes_idx;
      std::map<KSR::size_t, KSR::size_t> map_lines_idx;
      KSR::vector<IVertex> vertices;

      const std::size_t n = intersections.size();
      vertices.reserve(n);

      for (std::size_t i = 0; i < n; ++i) {
        const auto& iedge0 = intersections[i].first;
        const auto& iedge1 = intersections[(i + 1) % n].first;

        KSR::size_t common_plane_idx = KSR::no_element();
        std::set_intersection(
          m_intersection_graph.intersected_planes(iedge0).begin(),
          m_intersection_graph.intersected_planes(iedge0).end(),
          m_intersection_graph.intersected_planes(iedge1).begin(),
          m_intersection_graph.intersected_planes(iedge1).end(),
          boost::make_function_output_iterator(
            [&](const KSR::size_t& idx) -> void {
              if (idx < 6) {
                CGAL_assertion(common_plane_idx == KSR::no_element());
                common_plane_idx = idx;
              }
            }
          )
        );
        CGAL_assertion(common_plane_idx != KSR::no_element());
        common_planes_idx.push_back(common_plane_idx);

        typename std::map<KSR::size_t, KSR::size_t>::iterator iter;
        const auto pair = map_lines_idx.insert(std::make_pair(common_plane_idx, KSR::no_element()));
        const bool is_inserted = pair.second;
        if (is_inserted) {
          pair.first->second = m_intersection_graph.add_line();
        }
        vertices.push_back(m_intersection_graph.add_vertex(
          intersections[i].second).first);
      }
      CGAL_assertion(vertices.size() == n);

      for (std::size_t i = 0; i < n; ++i) {
        const auto& iplanes = m_intersection_graph.intersected_planes(intersections[i].first);
        for (const KSR::size_t sp_idx : iplanes) {
          support_plane(sp_idx).iedges().erase(intersections[i].first);
        }
        const auto edges = m_intersection_graph.split_edge(
          intersections[i].first, vertices[i]);

        const auto& iplanes_1 = m_intersection_graph.intersected_planes(edges.first);
        for (const KSR::size_t sp_idx : iplanes_1) {
          support_plane(sp_idx).iedges().insert(edges.first);
        }

        const auto& iplanes_2 = m_intersection_graph.intersected_planes(edges.second);
        for (const KSR::size_t sp_idx : iplanes_2) {
          support_plane(sp_idx).iedges().insert(edges.second);
        }

        const auto new_edge = m_intersection_graph.add_edge(
          vertices[i], vertices[(i + 1) % n], support_plane_idx).first;
        m_intersection_graph.intersected_planes(new_edge).insert(common_planes_idx[i]);
        m_intersection_graph.set_line(new_edge, map_lines_idx[common_planes_idx[i]]);

        support_plane(support_plane_idx).iedges().insert(new_edge);
        support_plane(common_planes_idx[i]).iedges().insert(new_edge);
      }
    }
    return support_plane_idx;
  }

  template<typename PointRange>
  void add_bbox_polygon(const PointRange& polygon) {

    const KSR::size_t support_plane_idx = add_support_plane(polygon);

    std::array<IVertex, 4> ivertices;
    std::array<Point_2, 4> points;
    for (std::size_t i = 0; i < 4; ++i) {
      points[i] = support_plane(support_plane_idx).to_2d(polygon[i]);
      ivertices[i] = m_intersection_graph.add_vertex(polygon[i]).first;
    }

    const auto vertices =
      support_plane(support_plane_idx).add_bbox_polygon(points, ivertices);

    for (std::size_t i = 0; i < 4; ++i) {
      const auto pair = m_intersection_graph.add_edge (ivertices[i], ivertices[(i+1)%4], support_plane_idx);
      const auto& iedge = pair.first;
      const bool is_inserted = pair.second;
      if (is_inserted) {
        m_intersection_graph.set_line(iedge, m_intersection_graph.add_line());
      }

      support_plane(support_plane_idx).set_iedge(vertices[i], vertices[(i + 1) % 4], iedge);
      support_plane(support_plane_idx).iedges().insert(iedge);
    }
  }

  template<typename PointRange>
  void add_input_polygon(const PointRange& polygon, const KSR::size_t input_idx) {

    const KSR::size_t support_plane_idx = add_support_plane(polygon);
    std::vector<Point_2> points;
    points.reserve(polygon.size());
    for (const auto& point : polygon) {
      const Point_3 converted(
        static_cast<FT>(point.x()),
        static_cast<FT>(point.y()),
        static_cast<FT>(point.z()));
      points.push_back(support_plane(support_plane_idx).to_2d(converted));
    }

    // const auto centroid = CGAL::centroid(points.begin(), points.end());

    using TRI = CGAL::Delaunay_triangulation_2<Kernel>;
    TRI tri(points.begin(), points.end());
    std::vector<Triangle_2> triangles;
    triangles.reserve(tri.number_of_faces());
    for (auto fit = tri.finite_faces_begin(); fit != tri.finite_faces_end(); ++fit) {
      triangles.push_back(Triangle_2(
          fit->vertex(0)->point(), fit->vertex(1)->point(), fit->vertex(2)->point()));
    }
    const auto centroid = CGAL::centroid(triangles.begin(), triangles.end());

    std::sort(points.begin(), points.end(),
    [&](const Point_2& a, const Point_2& b) -> bool {
      const Segment_2 sega(centroid, a);
      const Segment_2 segb(centroid, b);
      return ( Direction_2(sega) < Direction_2(segb) );
    });
    support_plane(support_plane_idx).add_input_polygon(points, centroid, input_idx);
  }

  /*******************************
  **        PSimplices          **
  ********************************/

  static PVertex null_pvertex() { return PVertex(KSR::no_element(), Vertex_index()); }
  static PEdge   null_pedge()   { return   PEdge(KSR::no_element(),   Edge_index()); }
  static PFace   null_pface()   { return   PFace(KSR::no_element(),   Face_index()); }

  const PVertices pvertices(const KSR::size_t support_plane_idx) const {
    return PVertices(
      boost::make_transform_iterator(
        mesh(support_plane_idx).vertices().begin(),
        Make_PSimplex<PVertex>(support_plane_idx)),
      boost::make_transform_iterator(
        mesh(support_plane_idx).vertices().end(),
        Make_PSimplex<PVertex>(support_plane_idx)));
  }

  const PEdges pedges(const KSR::size_t support_plane_idx) const {
    return PEdges(
      boost::make_transform_iterator(
        mesh(support_plane_idx).edges().begin(),
        Make_PSimplex<PEdge>(support_plane_idx)),
      boost::make_transform_iterator(
        mesh(support_plane_idx).edges().end(),
        Make_PSimplex<PEdge>(support_plane_idx)));
  }

  const PFaces pfaces(const KSR::size_t support_plane_idx) const {
    return PFaces(
      boost::make_transform_iterator(
        mesh(support_plane_idx).faces().begin(),
        Make_PSimplex<PFace>(support_plane_idx)),
      boost::make_transform_iterator(
        mesh(support_plane_idx).faces().end(),
        Make_PSimplex<PFace>(support_plane_idx)));
  }

  // Get prev and next of free vertex
  const PVertex prev(const PVertex& pvertex) const {
    return PVertex(pvertex.first, support_plane(pvertex).prev(pvertex.second));
  }
  const PVertex next(const PVertex& pvertex) const {
    return PVertex(pvertex.first, support_plane(pvertex).next(pvertex.second));
  }

  // Get prev and next of constrained vertex.
  const std::pair<PVertex, PVertex> prev_and_next(const PVertex& pvertex) const {

    std::pair<PVertex, PVertex> out(null_pvertex(), null_pvertex());
    for (const auto& he : halfedges_around_target(
      halfedge(pvertex.second, mesh(pvertex)), mesh(pvertex))) {

      const auto iedge = support_plane(pvertex).iedge(mesh(pvertex).edge(he));
      if (iedge == this->iedge(pvertex)) {
        continue;
      }
      if (out.first == null_pvertex()) {
        out.first  = PVertex(pvertex.first, mesh(pvertex).source(he));
      } else {
        out.second = PVertex(pvertex.first, mesh(pvertex).source(he));
        return out;
      }
    }
    return out;
  }

  const std::pair<PVertex, PVertex> border_prev_and_next(const PVertex& pvertex) const {

    // std::cout << point_3(pvertex) << std::endl;
    auto he = mesh(pvertex).halfedge(pvertex.second);
    const auto end = he;

    // std::cout << point_3(PVertex(pvertex.first, mesh(pvertex).source(he))) << std::endl;
    // std::cout << point_3(PVertex(pvertex.first, mesh(pvertex).target(he))) << std::endl;

    // If the assertion below fails, it probably means that we need to circulate
    // longer until we hit the border edge!

    std::size_t count = 0;
    while (true) {
      if (mesh(pvertex).face(he) != Face_index()) {
        he = mesh(pvertex).prev(mesh(pvertex).opposite(he));

        // std::cout << point_3(PVertex(pvertex.first, mesh(pvertex).source(he))) << std::endl;
        // std::cout << point_3(PVertex(pvertex.first, mesh(pvertex).target(he))) << std::endl;

        ++count;
      } else { break; }

      // std::cout << "count: " << count << std::endl;
      CGAL_assertion(count <= 2);
      if (he == end) {
        CGAL_assertion_msg(false, "ERROR: BORDER HALFEDGE IS NOT FOUND, FULL CIRCLE!");
        break;
      }
      if (count == 100) {
        CGAL_assertion_msg(false, "ERROR: BORDER HALFEDGE IS NOT FOUND, LIMIT ITERATIONS!");
        break;
      }
    }

    CGAL_assertion(mesh(pvertex).face(he) == Face_index());
    return std::make_pair(
      PVertex(pvertex.first, mesh(pvertex).source(he)),
      PVertex(pvertex.first, mesh(pvertex).target(mesh(pvertex).next(he))));
  }

  const PVertex add_pvertex(const KSR::size_t support_plane_idx, const Point_2& point) {

    CGAL_assertion(support_plane_idx != KSR::uninitialized());
    CGAL_assertion(support_plane_idx != KSR::no_element());

    auto& m = mesh(support_plane_idx);
    const auto vi = m.add_vertex(point);
    CGAL_assertion(vi != typename Support_plane::Mesh::Vertex_index());
    return PVertex(support_plane_idx, vi);
  }

  template<typename VertexRange>
  const PFace add_pface(const VertexRange& pvertices) {

    const auto support_plane_idx = pvertices.front().first;
    CGAL_assertion(support_plane_idx != KSR::uninitialized());
    CGAL_assertion(support_plane_idx != KSR::no_element());

    auto& m = mesh(support_plane_idx);
    const auto range = CGAL::make_range(
      boost::make_transform_iterator(pvertices.begin(),
      CGAL::Property_map_to_unary_function<CGAL::Second_of_pair_property_map<PVertex> >()),
      boost::make_transform_iterator(pvertices.end(),
      CGAL::Property_map_to_unary_function<CGAL::Second_of_pair_property_map<PVertex> >()));
    const auto fi = m.add_face(range);
    CGAL_assertion(fi != Support_plane::Mesh::null_face());
    return PFace(support_plane_idx, fi);
  }

  void clear_polygon_faces(const KSR::size_t support_plane_idx) {
    Mesh& m = mesh(support_plane_idx);
    for (const auto& fi : m.faces()) {
      m.remove_face(fi);
    }
    for (const auto& ei : m.edges()) {
      m.remove_edge(ei);
    }
    for (const auto& vi : m.vertices()) {
      m.set_halfedge(vi, Halfedge_index());
    }
  }

  const PVertex source(const PEdge& pedge) const {
    return PVertex(pedge.first, mesh(pedge).source(mesh(pedge).halfedge(pedge.second)));
  }
  const PVertex target(const PEdge& pedge) const {
    return PVertex(pedge.first, mesh(pedge).target(mesh(pedge).halfedge(pedge.second)));
  }
  const PVertex opposite(const PEdge& pedge, const PVertex& pvertex) const {

    if (mesh(pedge).target(mesh(pedge).halfedge(pedge.second)) == pvertex.second) {
      return PVertex(pedge.first, mesh(pedge).source(mesh(pedge).halfedge(pedge.second)));
    }
    CGAL_assertion(mesh(pedge).source(mesh(pedge).halfedge(pedge.second)) == pvertex.second);
    return PVertex(pedge.first, mesh(pedge).target(mesh(pedge).halfedge(pedge.second)));
  }

  const PFace pface_of_pvertex(const PVertex& pvertex) const {
    return PFace(pvertex.first, support_plane(pvertex).face(pvertex.second));
  }

  const std::pair<PFace, PFace> pfaces_of_pvertex(const PVertex& pvertex) const {

    std::pair<PFace, PFace> out(null_pface(), null_pface());
    std::tie(out.first.second, out.second.second) =
      support_plane(pvertex).faces(pvertex.second);
    if (out.first.second != Face_index()) {
      out.first.first = pvertex.first;
    }
    if (out.second.second != Face_index()) {
      out.second.first = pvertex.first;
    }
    return out;
  }

  const PVertices_of_pface pvertices_of_pface(const PFace& pface) const {

    return PVertices_of_pface(
      boost::make_transform_iterator(
        halfedges_around_face(halfedge(pface.second, mesh(pface)), mesh(pface)).begin(),
        Halfedge_to_pvertex(pface.first, mesh(pface))),
      boost::make_transform_iterator(
        halfedges_around_face(halfedge(pface.second, mesh(pface)), mesh(pface)).end(),
        Halfedge_to_pvertex(pface.first, mesh(pface))));
  }

  const PEdges_of_pface pedges_of_pface(const PFace& pface) const {

    return PEdges_of_pface(
      boost::make_transform_iterator(
        halfedges_around_face(halfedge(pface.second, mesh(pface)), mesh(pface)).begin(),
        Halfedge_to_pedge(pface.first, mesh(pface))),
      boost::make_transform_iterator(
        halfedges_around_face(halfedge(pface.second, mesh(pface)), mesh(pface)).end(),
        Halfedge_to_pedge(pface.first, mesh(pface))));
  }

  const PEdges_around_pvertex pedges_around_pvertex(const PVertex& pvertex) const {
    return PEdges_around_pvertex(
      boost::make_transform_iterator(
        halfedges_around_target(halfedge(pvertex.second, mesh(pvertex)), mesh(pvertex)).begin(),
        Halfedge_to_pedge(pvertex.first, mesh(pvertex))),
      boost::make_transform_iterator(
        halfedges_around_target(halfedge(pvertex.second, mesh(pvertex)), mesh(pvertex)).end(),
        Halfedge_to_pedge(pvertex.first, mesh(pvertex))));
  }

  const KSR::size_t& input(const PFace& pface) const{ return support_plane(pface).input(pface.second); }
  KSR::size_t& input(const PFace& pface) { return support_plane(pface).input(pface.second); }

  const unsigned int& k(const PFace& pface) const { return support_plane(pface).k(pface.second); }
  unsigned int& k(const PFace& pface) { return support_plane(pface).k(pface.second); }

  const bool is_frozen(const PVertex& pvertex) const { return support_plane(pvertex).is_frozen(pvertex.second); }

  const Vector_2& direction(const PVertex& pvertex) const { return support_plane(pvertex).direction(pvertex.second); }
  Vector_2& direction(const PVertex& pvertex) { return support_plane(pvertex).direction(pvertex.second); }

  const FT speed(const PVertex& pvertex) { return support_plane(pvertex).speed(pvertex.second); }

  const bool is_active(const PVertex& pvertex) const { return support_plane(pvertex).is_active(pvertex.second); }

  void deactivate(const PVertex& pvertex) {

    support_plane(pvertex).set_active(pvertex.second, false);
    if (iedge(pvertex) != null_iedge()) {
      m_intersection_graph.is_active(iedge(pvertex)) = false;
    }
    // std::cout << str(pvertex) << " ";
    if (ivertex(pvertex) != null_ivertex()) {
      // std::cout << " ivertex: " << point_3(ivertex(pvertex));
      m_intersection_graph.is_active(ivertex(pvertex)) = false;
    }
    // std::cout << std::endl;
  }

  void activate(const PVertex& pvertex) {

    support_plane(pvertex).set_active(pvertex.second, true);
    if (iedge(pvertex) != null_iedge()) {
      m_intersection_graph.is_active(iedge(pvertex)) = true;
    }
    if (ivertex(pvertex) != null_ivertex()) {
      m_intersection_graph.is_active(ivertex(pvertex)) = true;
    }
  }

  /*******************************
  **          ISimplices        **
  ********************************/

  static IVertex null_ivertex() { return Intersection_graph::null_ivertex(); }
  static IEdge null_iedge() { return Intersection_graph::null_iedge(); }

  decltype(auto) ivertices() const { return m_intersection_graph.vertices(); }
  decltype(auto) iedges() const { return m_intersection_graph.edges(); }

  const KSR::size_t nb_intersection_lines() const { return m_intersection_graph.nb_lines(); }
  const KSR::size_t line_idx(const IEdge& iedge) const { return m_intersection_graph.line(iedge); }
  const KSR::size_t line_idx(const PVertex& pvertex) const { return line_idx(iedge(pvertex)); }

  const IVertex add_ivertex(const Point_3& point, const KSR::Idx_set& support_planes_idx) {

    KSR::Idx_vector vec_planes;
    std::copy(
      support_planes_idx.begin(),
      support_planes_idx.end(),
      std::back_inserter(vec_planes));
    const auto pair = m_intersection_graph.add_vertex(point, vec_planes);
    const auto ivertex = pair.first;
    return ivertex;
  }

  void add_iedge(const KSR::Idx_set& support_planes_idx, KSR::vector<IVertex>& vertices) {

    const auto source = m_intersection_graph.point_3(vertices.front());
    std::sort(vertices.begin(), vertices.end(),
      [&](const IVertex& a, const IVertex& b) -> bool {
        const auto ap = m_intersection_graph.point_3(a);
        const auto bp = m_intersection_graph.point_3(b);
        const auto sq_dist_a = CGAL::squared_distance(source, ap);
        const auto sq_dist_b = CGAL::squared_distance(source, bp);
        return (sq_dist_a < sq_dist_b);
      }
    );

    KSR::size_t line_idx = m_intersection_graph.add_line();
    for (KSR::size_t i = 0; i < vertices.size() - 1; ++i) {

      const auto pair = m_intersection_graph.add_edge(
        vertices[i], vertices[i + 1], support_planes_idx);
      const auto iedge = pair.first;
      const auto is_inserted = pair.second;
      CGAL_assertion(is_inserted);
      m_intersection_graph.set_line(iedge, line_idx);

      for (const auto support_plane_idx : support_planes_idx) {
        support_plane(support_plane_idx).iedges().insert(iedge);
      }
    }
  }

  const IVertex source(const IEdge& edge) const { return m_intersection_graph.source(edge); }
  const IVertex target(const IEdge& edge) const { return m_intersection_graph.target(edge); }

  const IVertex opposite(const IEdge& edge, const IVertex& ivertex) const {
    const auto out = source(edge);
    if (out == ivertex) {
      return target(edge);
    }
    CGAL_assertion(target(edge) == ivertex);
    return out;
  }

  decltype(auto) incident_iedges(const IVertex& ivertex) const {
    return m_intersection_graph.incident_edges(ivertex);
  }

  const std::set<IEdge>& iedges(const KSR::size_t support_plane_idx) const {
    return support_plane(support_plane_idx).iedges();
  }

  const KSR::Idx_set& intersected_planes(const IEdge& iedge) const {
    return m_intersection_graph.intersected_planes(iedge);
  }

  const KSR::Idx_set intersected_planes(
    const IVertex& ivertex, const bool keep_bbox = true) const {

    KSR::Idx_set out;
    for (const auto incident_iedge : incident_iedges(ivertex)) {
      for (const auto support_plane_idx : intersected_planes(incident_iedge)) {
        if (!keep_bbox && support_plane_idx < 6) {
          continue;
        }
        out.insert(support_plane_idx);
      }
    }
    return out;
  }

  const bool is_iedge(const IVertex& source, const IVertex& target) const {
    return m_intersection_graph.is_edge(source, target);
  }

  const bool is_active(const IEdge& iedge) const {
    return m_intersection_graph.is_active(iedge);
  }
  const bool is_active(const IVertex& ivertex) const {
    return m_intersection_graph.is_active(ivertex);
  }

  const bool is_bbox_iedge(const IEdge& edge) const {

    for (const auto support_plane_idx : m_intersection_graph.intersected_planes(edge)) {
      if (support_plane_idx < 6) {
        return true;
      }
    }
    return false;
  }

  /*******************************
  **        Connectivity        **
  ********************************/

  const bool has_ivertex(const PVertex& pvertex) const { return support_plane(pvertex).has_ivertex(pvertex.second); }
  const IVertex ivertex(const PVertex& pvertex) const { return support_plane(pvertex).ivertex(pvertex.second); }

  const bool has_iedge(const PVertex& pvertex) const { return support_plane(pvertex).has_iedge(pvertex.second); }
  const IEdge iedge(const PVertex& pvertex) const { return support_plane(pvertex).iedge(pvertex.second); }

  const bool has_iedge(const PEdge& pedge) const { return support_plane(pedge).has_iedge(pedge.second); }
  const IEdge iedge(const PEdge& pedge) const { return support_plane(pedge).iedge(pedge.second); }

  void connect(const PVertex& pvertex, const IVertex& ivertex) { support_plane(pvertex).set_ivertex(pvertex.second, ivertex); }
  void connect(const PVertex& pvertex, const IEdge& iedge) { support_plane(pvertex).set_iedge(pvertex.second, iedge); }
  void connect(const PVertex& a, const PVertex& b, const IEdge& iedge) { support_plane(a).set_iedge(a.second, b.second, iedge); }
  void connect(const PEdge& pedge, const IEdge& iedge) { support_plane(pedge).set_iedge(pedge.second, iedge); }

  const IVertex disconnect_ivertex(const PVertex& pvertex) {
    const auto out = ivertex(pvertex);
    support_plane(pvertex).set_ivertex(pvertex.second, null_ivertex());
    return out;
  }

  const IEdge disconnect_iedge(const PVertex& pvertex) {
    const auto out = iedge(pvertex);
    support_plane(pvertex).set_iedge(pvertex.second, null_iedge());
    return out;
  }

  struct Queue_element {
    PVertex previous;
    PVertex pvertex;
    bool front;
    bool previous_was_free;

    Queue_element(
      const PVertex& previous, const PVertex& pvertex,
      const bool front, const bool previous_was_free) :
    previous(previous), pvertex(pvertex),
    front(front), previous_was_free(previous_was_free)
    { }
  };

  const std::vector<PVertex> pvertices_around_ivertex(
    const PVertex& pvertex, const IVertex& ivertex) const {

    // std::cout.precision(20);
    std::deque<PVertex> vertices;
    vertices.push_back(pvertex);

    std::cout << "came from: " <<
    str(iedge(pvertex)) << " " << segment_3(iedge(pvertex)) << std::endl;

    std::queue<Queue_element> todo;
    PVertex prev, next;
    std::tie(prev, next) = border_prev_and_next(pvertex);
    // std::cout << "prev in: " << str(prev) << " " << point_3(prev) << std::endl;
    // std::cout << "next in: " << str(next) << " " << point_3(next) << std::endl;
    // std::cout << "curr in: " << str(pvertex) << " " << point_3(pvertex) << std::endl;

    todo.push(Queue_element(pvertex, prev, true, false));
    todo.push(Queue_element(pvertex, next, false, false));

    while (!todo.empty()) {

      // std::cout << std::endl;
      auto previous = todo.front().previous;
      auto current = todo.front().pvertex;
      bool front = todo.front().front;
      bool previous_was_free = todo.front().previous_was_free;
      todo.pop();

      auto iedge = this->iedge(current);
      bool is_free = (iedge == null_iedge());
      // std::cout << "is free 1: " << is_free << std::endl;

      // std::cout << "iedge: " << segment_3(iedge) << std::endl;
      if (!is_free && source(iedge) != ivertex && target(iedge) != ivertex) {
        // std::cout << "is free 2: " << is_free << std::endl;
        is_free = true;
      }

      if (!is_free) {

        auto other = source(iedge);
        if (other == ivertex) {
          other = target(iedge);
        } else {
          CGAL_assertion(target(iedge) == ivertex);
        }

        // Filter backwards vertex.
        const Vector_2 dir1 = direction(current);
        // std::cout << "dir1: " << dir1 << std::endl;
        const Vector_2 dir2(
          point_2(current.first, other), point_2(current.first, ivertex));
        // std::cout << "dir2: " << dir2 << std::endl;
        const FT dot_product = dir1 * dir2;
        // std::cout << "dot: " << dot_product << std::endl;

        if (dot_product < FT(0)) {
          std::cerr << str(current) << " is backwards" << std::endl;
          // std::cout << point_3(current) << std::endl;
          is_free = true;
        }

        if (is_frozen(current)) {
          std::cerr << str(current) << " is frozen" << std::endl;
          // std::cout << point_3(current) << std::endl;
          is_free = true;
        }
        // std::cout << "is free 3: " << is_free << std::endl;
      }

      if (previous_was_free && is_free) {
        std::cerr << str(current) << " has no iedge, stopping there" << std::endl;
        // std::cout << point_3(current) << std::endl;
        continue;
      }

      if (is_free) {
        std::cerr << str(current) << " has no iedge" << std::endl;
        // std::cout << point_3(current) << std::endl;
      }
      else {
        std::cerr << str(current) << " has iedge " << str(iedge)
                  << " from " << str(source(iedge)) << " to " << str(target(iedge)) << std::endl;
        // std::cout << segment_3(iedge) << std::endl;
        // std::cout << point_3(current) << std::endl;
      }

      if (front) {
        vertices.push_front (current);
        // std::cout << "pushed front" << std::endl;
      }
      else {
        vertices.push_back (current);
        // std::cout << "pushed back" << std::endl;
      }

      std::tie(prev, next) = border_prev_and_next(current);
      if (prev == previous) {
        CGAL_assertion(next != previous);
        todo.push(Queue_element(current, next, front, is_free));
        // std::cout << "pushed next" << std::endl;
      }
      else {
        todo.push(Queue_element(current, prev, front, is_free));
        // std::cout << "pushed prev" << std::endl;
      }
    }

    std::vector<PVertex> out;
    out.reserve(vertices.size());
    std::copy(vertices.begin(), vertices.end(), std::back_inserter(out));

    std::cout << "*** Found vertices:";
    for (const auto& pv : out)
      std::cout << " " << str(pv);
    std::cout << std::endl;
    return out;
  }

  /*******************************
  **        Conversions         **
  ********************************/

  const Point_2 to_2d(const KSR::size_t support_plane_idx, const IVertex& ivertex) const {
    return support_plane(support_plane_idx).to_2d(point_3(ivertex));
  }
  const Segment_2 to_2d(const KSR::size_t support_plane_idx, const Segment_3& segment_3) const {
    return support_plane(support_plane_idx).to_2d(segment_3);
  }

  const Point_2 point_2(const PVertex& pvertex, const FT time) const {
    return support_plane(pvertex).point_2(pvertex.second, time);
  }
  const Point_2 point_2(const PVertex& pvertex) const {
    return point_2(pvertex, m_current_time);
  }
  const Point_2 point_2(const KSR::size_t support_plane_idx, const IVertex& ivertex) const {
    return support_plane(support_plane_idx).to_2d(point_3(ivertex));
  }

  const Segment_2 segment_2(const KSR::size_t support_plane_idx, const IEdge& iedge) const {
    return support_plane(support_plane_idx).to_2d(segment_3(iedge));
  }

  const Point_3 to_3d(const KSR::size_t support_plane_idx, const Point_2& point_2) const {
    return support_plane(support_plane_idx).to_3d(point_2);
  }

  const Point_3 point_3(const PVertex& pvertex, const FT time) const {
    return support_plane(pvertex).point_3(pvertex.second, time);
  }
  const Point_3 point_3(const PVertex& pvertex) const {
    return point_3(pvertex, m_current_time);
  }
  const Point_3 point_3(const IVertex& vertex) const {
    return m_intersection_graph.point_3(vertex);
  }

  const Segment_3 segment_3(const PEdge& pedge, const FT time) const {
    return support_plane(pedge).segment_3(pedge.second, time);
  }
  const Segment_3 segment_3(const PEdge& pedge) const {
    return segment_3 (pedge, m_current_time);
  }
  const Segment_3 segment_3(const IEdge& edge) const {
    return m_intersection_graph.segment_3(edge);
  }

  /*******************************
  **          Predicates        **
  ********************************/

  // Check if there is a collision with another polygon.
  const std::pair<bool, bool> collision_occured (
    const PVertex& pvertex, const IEdge& iedge) const {

    // const FT tol = FT(1) / FT(100000);
    bool collision = false;
    for (const auto support_plane_idx : intersected_planes(iedge)) {
      if (support_plane_idx < 6) {
        return std::make_pair(true, true);
      }

      for (const auto pedge : pedges(support_plane_idx)) {
        if (this->iedge(pedge) == iedge) {
          const auto pedge_segment = Segment_3(point_3(source(pedge)), point_3(target(pedge)));

          const Segment_3 source_to_pvertex(pedge_segment.source(), point_3(pvertex));
          // if (CGAL::sqrt(source_to_pvertex.squared_length()) < tol) {
          //   // std::cout << "WARNING: POINTS ARE ALMOST EQUAL!" << std::endl;
          //   collision = true;
          //   break;
          // }

          // std::cout << point_3(source(pedge)) << std::endl;
          // std::cout << point_3(target(pedge)) << std::endl;
          // std::cout << point_3(pvertex) << std::endl;

          const FT dot_product = pedge_segment.to_vector() * source_to_pvertex.to_vector();
          if (dot_product < FT(0)) {
            continue;
          }
          // std::cout << source_to_pvertex.squared_length() << std::endl;
          // std::cout << pedge_segment.squared_length() << std::endl;

          if (pedge_segment.squared_length() == FT(0)) {
            std::cout << "ERROR: SOURCE_TO_PVERTEX/PEDGE SEGMENT SQ LENGTH = "
            << source_to_pvertex.squared_length() << std::endl;
          }
          CGAL_assertion(pedge_segment.squared_length() != FT(0));

          // if (pedge_segment.squared_length() == FT(0)) {
          //   if (pedge_segment.source() == point_3(pvertex)) {
          //     collision = false;
          //     break;
          //   }
          // }
          // CGAL_assertion(pedge_segment.squared_length() != FT(0));

          if (source_to_pvertex.squared_length() <= pedge_segment.squared_length()) {
            collision = true;
            break;
          }
        }
      }
    }
    return std::make_pair(collision, false);
  }

  const std::pair<bool, bool> is_occupied(
    const PVertex& pvertex,
    const IVertex& ivertex,
    const IEdge& query_iedge) {

    const auto pair = is_occupied(pvertex, query_iedge);
    const bool has_polygon = pair.first;
    const bool is_bbox_reached = pair.second;

    if (is_bbox_reached) return std::make_pair(true, true);
    CGAL_assertion(!is_bbox_reached);
    if (!has_polygon) {
      // std::cout << "NO POLYGON DETECTED" << std::endl;
      return std::make_pair(false, false);
    }
    CGAL_assertion(has_polygon);

    CGAL_assertion(ivertex != null_ivertex());
    std::set<PEdge> pedges;
    get_occupied_pedges(pvertex, query_iedge, pedges);
    for (const auto& pedge : pedges) {
      CGAL_assertion(pedge != null_pedge());
      // std::cout << "PEDGE: " << segment_3(pedge) << std::endl;

      const auto source = this->source(pedge);
      const auto target = this->target(pedge);
      if (this->ivertex(source) == ivertex || this->ivertex(target) == ivertex) {
        return std::make_pair(true, false);
      }
    }
    return std::make_pair(false, false);
  }

  void get_occupied_pedges(
    const PVertex& pvertex,
    const IEdge& query_iedge,
    std::set<PEdge>& pedges) {

    pedges.clear();
    for (const auto plane_idx : intersected_planes(query_iedge)) {
      if (plane_idx == pvertex.first) continue; // current plane
      if (plane_idx < 6) continue; // bbox plane

      for (const auto pedge : this->pedges(plane_idx)) {
        if (iedge(pedge) == query_iedge) {
          pedges.insert(pedge);
        }
      }
    }
    CGAL_assertion(pedges.size() > 0);
  }

  const std::pair<bool, bool> is_occupied(
    const PVertex& pvertex,
    const IEdge& query_iedge) {

    CGAL_assertion(query_iedge != null_iedge());
    // std::cout << str(query_iedge) << " " << segment_3(query_iedge) << std::endl;
    KSR::size_t num_adjacent_faces = 0;
    for (const auto plane_idx : intersected_planes(query_iedge)) {
      if (plane_idx == pvertex.first) continue; // current plane
      if (plane_idx < 6) return std::make_pair(true, true); // bbox plane

      for (const auto pedge : pedges(plane_idx)) {
        if (!has_iedge(pedge)) continue;

        // std::cout << str(iedge(pedge)) << std::endl;
        if (iedge(pedge) == query_iedge) {
          const auto& m = mesh(plane_idx);
          const auto he = m.halfedge(pedge.second);
          const auto op = m.opposite(he);
          const auto face1 = m.face(he);
          const auto face2 = m.face(op);
          if (face1 != Support_plane::Mesh::null_face()) ++num_adjacent_faces;
          if (face2 != Support_plane::Mesh::null_face()) ++num_adjacent_faces;
        }
      }
    }

    std::cout << "num adjacent faces: " << num_adjacent_faces << std::endl;
    if (num_adjacent_faces <= 1)
      return std::make_pair(false, false);
    return std::make_pair(true, false);
  }

  /*******************************
  **    Operations on polygons  **
  ********************************/

  const PVertex crop_polygon(const PVertex& pvertex, const IEdge& iedge) {

    std::cout << "*** Cropping " << str(pvertex) << " along " << str(iedge) << std::endl;

    Point_2 future_point_a, future_point_b;
    Vector_2 direction_a, direction_b;

    compute_future_points_and_directions (pvertex, iedge,
                                          future_point_a, future_point_b,
                                          direction_a, direction_b);

    PEdge pedge (pvertex.first, support_plane(pvertex).split_vertex(pvertex.second));
    CGAL_assertion (source(pedge) == pvertex || target(pedge) == pvertex);

    PVertex other = opposite(pedge, pvertex);

    std::cout << "*** New edge " << str(pedge) << " between " << str(pvertex)
                     << " and " << str(other) << std::endl;

    connect (pedge, iedge);
    connect (pvertex, iedge);
    connect (other, iedge);

    support_plane(pvertex).set_point (pvertex.second, future_point_a);
    support_plane(other).set_point (other.second, future_point_b);

    direction(pvertex) = direction_a;
    direction(other) = direction_b;

    // std::cout << "pvertex: " << point_3(pvertex) << std::endl;
    // std::cout << "pvertex dir: " << direction_a << std::endl;
    // std::cout << "other: " << point_3(other) << std::endl;
    // std::cout << "other dir: " << direction_b << std::endl;

    std::cout << "New vertices: " << str(other) << std::endl;
    return other;
  }

  std::array<PVertex, 3> propagate_polygon (
    const unsigned int last_k,
    const PVertex& pvertex, const IEdge& iedge)
  {
    std::cout << "*** Propagating " << str(pvertex) << " along " << str(iedge) << std::endl;

    Point_2 original_point = point_2 (pvertex, 0);
    Vector_2 original_direction = direction(pvertex);

    PVertex other = crop_polygon (pvertex, iedge);

    PVertex propagated = add_pvertex (pvertex.first, original_point);
    direction(propagated) = original_direction;

    std::array<PVertex, 3> pvertices;

    pvertices[0] = pvertex;
    pvertices[1] = other;
    pvertices[2] = propagated;

    PFace new_pface = add_pface (pvertices);
    this->k(new_pface) = last_k;
    CGAL_assertion (new_pface.second != Face_index());

    std::cout << "*** New face " << lstr(new_pface) << std::endl;

    return pvertices;
  }

  void crop_polygon (const PVertex& pv0, const PVertex& pv1, const IEdge& iedge)
  {
    std::cout << "*** Cropping " << str(pv0) << "/" << str(pv1) << " along " << str(iedge) << std::endl;

    std::cout.precision(20);
    // std::cout << "pv0: " << point_3(pv0) << std::endl;
    // std::cout << "pv1: " << point_3(pv1) << std::endl;

    Point_2 future_point;
    Vector_2 future_direction;
    // const Line_2 iedge_line = segment_2(pv0.first, iedge).supporting_line();
    CGAL_assertion(pv0.first == pv1.first);

    {
      // const Point_2 pinit = iedge_line.projection(point_2(pv0, m_current_time));
      // const Point_2 future_point = iedge_line.projection(point_2(pv0, m_current_time + FT(1)));

      const PVertex prev(pv0.first, support_plane(pv0).prev(pv0.second));
      const PVertex next(pv0.first, support_plane(pv0).next(pv0.second));

      if (prev == pv1)
        compute_future_point_and_direction(0, pv0, next, iedge, future_point, future_direction);
      else {
        CGAL_assertion(next == pv1);
        compute_future_point_and_direction(0, pv0, prev, iedge, future_point, future_direction);
      }

      direction(pv0) = future_direction;
      std::cout << "pv0 dir: " << direction(pv0) << std::endl;
      support_plane(pv0).set_point(pv0.second, future_point);
      connect(pv0, iedge);
    }

    {
      // const Point_2 pinit = iedge_line.projection(point_2(pv1, m_current_time));
      // const Point_2 future_point = iedge_line.projection(point_2(pv1, m_current_time + FT(1)));

      const PVertex prev(pv1.first, support_plane(pv1).prev(pv1.second));
      const PVertex next(pv1.first, support_plane(pv1).next(pv1.second));

      if (prev == pv0)
        compute_future_point_and_direction(0, pv1, next, iedge, future_point, future_direction);
      else {
        CGAL_assertion(next == pv0);
        compute_future_point_and_direction(0, pv1, prev, iedge, future_point, future_direction);
      }

      direction(pv1) = future_direction;
      std::cout << "pv1 dir: " << direction(pv1) << std::endl;
      support_plane(pv1).set_point(pv1.second, future_point);
      connect(pv1, iedge);
    }

    const PEdge pedge(pv0.first, support_plane(pv0).edge(pv0.second, pv1.second));
    connect(pedge, iedge);
  }

  std::pair<PVertex, PVertex> propagate_polygon(
    const unsigned int, // last_k,
    const PVertex& pvertex, const PVertex& pother, const IEdge& iedge)
  {
    std::cout << "*** Propagating " << str(pvertex) << "/" << str(pother) << " along " << str(iedge) << std::endl;
    CGAL_assertion_msg(false, "TODO: PROPAGATE POLYGON VIA THE EDGE!");
    return std::make_pair(null_pvertex(), null_pvertex());
  }

  bool transfer_vertex (const PVertex& pvertex, const PVertex& pother)
  {
    std::cout << "*** Transfering " << str(pother) << " through " << str(pvertex) << std::endl;

    // If pvertex is adjacent to one or two.
    PFace source_face, target_face;
    std::tie (source_face, target_face) = pfaces_of_pvertex (pvertex);

    PFace common_pface = pface_of_pvertex (pother);

    if (common_pface == target_face)
      std::swap (source_face, target_face);
    CGAL_assertion (common_pface == source_face);

    std::cout << "*** Initial faces: " << lstr(source_face)
                     << " and " << lstr(target_face) << std::endl;

    // std::cout << "pv: " << point_3(pvertex) << std::endl;
    // std::cout << "p0: " << point_3(pother) << std::endl;
    // if (source_face != null_pface())
    //   std::cout << "source pface: " << centroid_of_pface(source_face) << std::endl;
    // if (target_face != null_pface())
    //   std::cout << "target pface: " << centroid_of_pface(target_face) << std::endl;

    PVertex pthird = next(pother);
    if (pthird == pvertex)
      pthird = prev(pother);

    if (target_face == null_pface())
    {
      Vector_2 new_direction;

      Line_2 iedge_line = segment_2(pother.first, iedge(pvertex)).supporting_line();
      Point_2 pinit = iedge_line.projection(point_2 (pother, m_current_time));

      Line_2 future_line (point_2 (pother, m_current_time + 1),
                          point_2 (pthird, m_current_time + 1));

      Point_2 future_point = KSR::intersection<Point_2>(future_line, iedge_line);

      direction(pvertex) = Vector_2 (pinit, future_point);
      support_plane(pvertex).set_point (pvertex.second,
                                        pinit - direction(pvertex) * m_current_time);

      const auto he = mesh(pvertex).halfedge(pother.second, pvertex.second);
      CGAL::Euler::join_vertex(he, mesh(pvertex));
    }
    else
    {
      IEdge iedge = disconnect_iedge (pvertex);
      // std::cerr << "disconnect " << str(pvertex) << " from " << str(iedge) << std::endl;

      PEdge pedge = null_pedge();
      for (PEdge pe : pedges_around_pvertex (pvertex))
        if (this->iedge(pe) == iedge)
        {
          pedge = pe;
          break;
        }
      CGAL_assertion (pedge != null_pedge());

      auto he = mesh(pedge).halfedge(pedge.second);
      if (mesh(pedge).face(he) != common_pface.second)
        he = mesh(pedge).opposite(he);
      CGAL_assertion (mesh(pedge).face(he) == common_pface.second);

      if (mesh(pedge).target(he) == pvertex.second)
      {
        // std::cerr << "shift target." << std::endl;
        CGAL::Euler::shift_target (he, mesh(pedge));
      }
      else
      {
        CGAL_assertion(mesh(pedge).source(he) == pvertex.second);
        // std::cerr << "shift source." << std::endl;
        CGAL::Euler::shift_source (he, mesh(pedge));
      }

      Vector_2 new_direction;

      Line_2 iedge_line = segment_2(pother.first, iedge).supporting_line();
      Point_2 pinit = iedge_line.projection(point_2 (pother, m_current_time));

      direction(pvertex) = direction(pother);
      support_plane(pother).set_point (pvertex.second,
                                       pinit - direction(pvertex) * m_current_time);

      Line_2 future_line (point_2 (pvertex, m_current_time + 1),
                          point_2 (pthird, m_current_time + 1));

      Point_2 future_point = KSR::intersection<Point_2>(future_line, iedge_line);

      direction(pother) = Vector_2 (pinit, future_point);
      support_plane(pother).set_point (pother.second,
                                       pinit - direction(pother) * m_current_time);

      // std::cerr << "connect " << str(pother) << " to " << str(iedge) << std::endl;
      connect (pother, iedge);
    }

    std::cout << "*** New faces: " << lstr(source_face)
                     << " and " << lstr(target_face) << std::endl;

    return (target_face != null_pface());
  }

  void merge_pvertices (const PVertex& pvertex, const PVertex& pother)
  {
    std::cout << "*** Merging " << str(pvertex) << " with " << str(pother) << std::endl;

    const auto he = mesh(pvertex).halfedge(pother.second, pvertex.second);
    disconnect_ivertex (pother);
    CGAL::Euler::join_vertex(he, mesh(pvertex));
  }

  std::vector<PVertex> merge_pvertices_on_ivertex (const FT min_time,
                                                   const FT max_time,
                                                   const PVertex& event_pvertex,
                                                   std::vector<PVertex>& pvertices,
                                                   const IVertex& ivertex,
                                                   std::vector<IEdge>& crossed)
  {
    crossed.clear();
    KSR::size_t support_plane_idx = pvertices.front().first;

    PVertex prev = pvertices.front();
    PVertex next = pvertices.back();

    IEdge prev_iedge = null_iedge(), next_iedge = null_iedge();

    // std::cout << "start from: " << segment_3(iedge(pvertices[1])) << std::endl;
    // std::ofstream("start_from.polylines.txt")
    // << "2 " << segment_3(iedge(pvertices[1])) << std::endl;

    std::cout << "came from: " <<
    str(iedge(pvertices[1])) << " " << segment_3(iedge(pvertices[1])) << std::endl;

    // Copy front/back to remember position/direction.
    PVertex front, back;
    if (pvertices.size() < 3) {
      CGAL_assertion_msg(false, "ERROR: INVALID CASE!");
    } else if (pvertices.size() == 3 || pvertices.size() == 4) {

      // BUG: In this case, the point that is duplicated twice is not always copied.
      // To fix it, we copy the second point not from the original vertex but from the first
      // copy of that vertex.

      const auto& initial = pvertices[1];
      front = PVertex(support_plane_idx, support_plane(support_plane_idx).duplicate_vertex(initial.second));
      support_plane(support_plane_idx).set_point(
        front.second, support_plane(support_plane_idx).get_point(initial.second));
      back  = PVertex(support_plane_idx, support_plane(support_plane_idx).duplicate_vertex(front.second));
      support_plane(support_plane_idx).set_point(
        back.second, support_plane(support_plane_idx).get_point(front.second));

    } else if (pvertices.size() >= 5) {

      const auto& initial1 = pvertices[1];
      front = PVertex(support_plane_idx, support_plane(support_plane_idx).duplicate_vertex(initial1.second));
      support_plane(support_plane_idx).set_point(
        front.second, support_plane(support_plane_idx).get_point(initial1.second));

      const auto& initial2 = pvertices[pvertices.size() - 2];
      back  = PVertex(support_plane_idx, support_plane(support_plane_idx).duplicate_vertex(initial2.second));
      support_plane(support_plane_idx).set_point(
        back.second, support_plane(support_plane_idx).get_point(initial2.second));

    } else {
      CGAL_assertion_msg(false, "ERROR: INVALID CASE!");
    }

    // auto pvertex_to_point =
    //   [&](const PVertex& a) -> Point_2 {
    //     return point_2(a);
    //   };

    // PFace fprev = pface_of_pvertex(prev);
    // Point_2 pprev = CGAL::centroid
    //   (boost::make_transform_iterator (pvertices_of_pface(fprev).begin(), pvertex_to_point),
    //    boost::make_transform_iterator (pvertices_of_pface(fprev).end(), pvertex_to_point));
    // PFace fnext = pface_of_pvertex(next);
    // Point_2 pnext = CGAL::centroid
    //   (boost::make_transform_iterator (pvertices_of_pface(fnext).begin(), pvertex_to_point),
    //    boost::make_transform_iterator (pvertices_of_pface(fnext).end(), pvertex_to_point));

    bool was_swapped = false;
    // if (CGAL::orientation(pprev, point_2(support_plane_idx, ivertex), pnext) == CGAL::LEFT_TURN) {
    //   std::cout << "Swapped!" << std::endl;
    //   was_swapped = true;
    //   std::swap(prev, next);
    //   std::swap(front, back);
    // }

    // Freeze vertices.
    for (std::size_t i = 1; i < pvertices.size() - 1; ++i) {
      PVertex& pvertex = pvertices[i];
      Point_2 point = point_2(support_plane_idx, ivertex);
      support_plane(pvertex).direction(pvertex.second) = CGAL::NULL_VECTOR;
      support_plane(pvertex).set_point(pvertex.second, point);
    }

    PVertex pvertex = pvertices[1];
    connect (pvertex, ivertex);

    std::cout << "*** Frozen vertex: " << str(pvertex) << std::endl;
    // std::cout << point_3(pvertex) << std::endl;
    // std::cout << "*** Removed vertices:";

    // Join vertices
    for (std::size_t i = 2; i < pvertices.size() - 1; ++ i)
    {
      // std::cout << " " << str(pvertices[i]) << std::endl;
      // std::cout << point_3(pvertices[i]) << std::endl;

      const auto he = mesh(support_plane_idx).halfedge(pvertices[i].second, pvertex.second);
      disconnect_ivertex (pvertices[i]);
      CGAL::Euler::join_vertex(he, mesh(support_plane_idx));
    }
    // std::cout << std::endl;

    auto i_iedges = incident_iedges(ivertex);
    std::vector<std::pair<IEdge, Direction_2> > iedges;
    std::copy (i_iedges.begin(), i_iedges.end(),
               boost::make_function_output_iterator
               ([&](const IEdge& ie) -> void
                {
                  if (intersected_planes(ie).find (support_plane_idx)
                      == intersected_planes(ie).end())
                    return;

                  Direction_2 dir (point_2 (support_plane_idx, opposite (ie, ivertex))
                                   - point_2 (support_plane_idx, ivertex));
                  iedges.push_back (std::make_pair (ie, dir));
                }));

    std::sort (iedges.begin(), iedges.end(),
               [&](const std::pair<IEdge, Direction_2>& a,
                   const std::pair<IEdge, Direction_2>& b) -> bool
               {
                 return a.second < b.second;
               });
    CGAL_assertion(iedges.size() != 0);

    std::cout.precision(20);
    std::cout << "Prev = "  << point_3 (prev)  << " / " << direction (prev)  << std::endl
              << "Front = " << point_3 (front) << " / " << direction (front) << std::endl
              << "Back = "  << point_3 (back)  << " / " << direction (back)  << std::endl
              << "Next = "  << point_3 (next)  << " / " << direction (next)  << std::endl;

    // std::cout << (iedge(next) != null_iedge()) << std::endl;
    // std::cout << "source: " << point_3(source(iedge(next))) << std::endl;
    // std::cout << "target: " << point_3(target(iedge(next))) << std::endl;
    // std::cout << "ivertex: " << point_3(ivertex) << std::endl;

    bool back_constrained = false;
    if ((iedge(next) != null_iedge()
         && (source(iedge(next)) == ivertex || target(iedge(next)) == ivertex))
        || (this->ivertex(next) != null_ivertex()
            && is_iedge (this->ivertex(next), ivertex)))
      back_constrained = true;

    bool front_constrained = false;
    if ((iedge(prev) != null_iedge()
         && (source(iedge(prev)) == ivertex || target(iedge(prev)) == ivertex))
        || (this->ivertex(prev) != null_ivertex()
            && is_iedge (this->ivertex(prev), ivertex)))
      front_constrained = true;

    std::vector<PVertex> new_vertices;
    if (back_constrained && front_constrained) // Closing case
    {
      std::cout << "*** Closing case" << std::endl;
    }
    else if (back_constrained) // Border case
    {
      std::cout << "*** Back border case" << std::endl;

      CGAL_assertion(has_iedge(pvertex));
      // std::ofstream("limit.polylines.txt")
      // << "2 " << segment_3(iedge(pvertex)) << std::endl;
      const KSR::size_t other_side_limit = line_idx(pvertex);

      // const Direction_2 dir(point_2(prev) - point_2(pvertex));

      const FT prev_time = last_event_time(prev);
      CGAL_assertion(prev_time < m_current_time);
      CGAL_assertion(prev_time >= FT(0));

      const auto pp_last = point_2(prev, prev_time);
      const auto pp_curr = point_2(prev, m_current_time);
      const auto dirp = Vector_2(pp_last, pp_curr);
      const auto tmp_prev = pp_curr - dirp / FT(10);

      const Direction_2 tmp_dir(tmp_prev - point_2(pvertex.first, ivertex));
      // std::cout << "tmp_dir: " << to_3d(prev.first, tmp_prev) << std::endl;

      std::reverse(iedges.begin(), iedges.end());

      // std::cout << "initial iedges: " << std::endl;
      // for (const auto& iedge : iedges) {
      //   std::cout << segment_3(iedge.first) << std::endl;
      // }

      KSR::size_t first_idx = KSR::no_element();
      for (std::size_t i = 0; i < iedges.size(); ++i) {
        if (tmp_dir.counterclockwise_in_between(
          iedges[(i + 1) % iedges.size()].second, iedges[i].second)) {

          first_idx = (i + 1) % iedges.size();
          break;
        }
      }

      // std::cout << "first: " << segment_3(iedges[first_idx].first) << std::endl;
      // std::ofstream("first.polylines.txt")
      // << "2 " << segment_3(iedges[first_idx].first) << std::endl;

      CGAL_assertion(first_idx != KSR::no_element());
      crossed.clear();

      KSR::size_t iedge_idx = first_idx;
      std::size_t iter = 0;
      while (true) {
        const IEdge& iedge = iedges[iedge_idx].first;
        bool collision, bbox_reached;
        std::tie(collision, bbox_reached) = collision_occured(pvertex, iedge);
        const bool limit_reached = (line_idx(iedge) == other_side_limit);
        std::cout << "limit/bbox: " << limit_reached << "/" << bbox_reached << std::endl;

        // std::cout << "next: " << segment_3(iedge) << std::endl;
        // std::ofstream("next" + std::to_string(iter) + ".polylines.txt")
        // << "2 " << segment_3(iedge) << std::endl;
        crossed.push_back(iedge);
        if (limit_reached || bbox_reached) {
          break;
        }
        iedge_idx = (iedge_idx + 1) % iedges.size();
        if (iter == 100) {
          CGAL_assertion_msg(false, "ERROR: BACK STD WHY SO MANY ITERATIONS?");
        } ++iter;
      }

      std::vector<IEdge> all_crossed;
      iedge_idx = first_idx; iter = 0;
      while (true) {
        const IEdge& iedge = iedges[iedge_idx].first;
        bool limit_reached, bbox_reached;
        std::tie(limit_reached, bbox_reached) = is_occupied(pvertex, iedge);
        all_crossed.push_back(iedge);
        if (limit_reached || bbox_reached) {
          break;
        }
        iedge_idx = (iedge_idx + 1) % iedges.size();
        if (iter == 100) {
          CGAL_assertion_msg(false, "ERROR: BACK LIMIT WHY SO MANY ITERATIONS?");
        } ++iter;
      }

      CGAL_assertion(crossed.size() != 0);
      std::cerr << "IEdges crossed = " << crossed.size() << std::endl;
      for (const auto& iedge : crossed)
        std::cout << segment_3(iedge) << std::endl;
      CGAL_assertion(crossed[0] == all_crossed[0]);

      std::vector<Point_2> future_points;
      std::vector<Vector_2> future_directions;

      if (crossed.size() > all_crossed.size()) {
        future_points.resize(crossed.size());
        future_directions.resize(crossed.size());
        for (std::size_t i = 0; i < crossed.size(); ++i) {
          const bool is_parallel = compute_future_point_and_direction(
            i, back, prev, crossed[i], future_points[i], future_directions[i]);
          if (is_parallel) {
            if (is_intersecting_iedge(min_time, max_time, prev, crossed[i])) {
              prev_iedge = crossed[i];
            }
          }
        }
      } else {
        future_points.resize(all_crossed.size());
        future_directions.resize(all_crossed.size());
        for (std::size_t i = 0; i < all_crossed.size(); ++i) {
          const bool is_parallel = compute_future_point_and_direction(
            i, back, prev, all_crossed[i], future_points[i], future_directions[i]);
          if (is_parallel) {
            if (is_intersecting_iedge(min_time, max_time, prev, all_crossed[i])) {
              prev_iedge = all_crossed[i];
            }
          }
        }
      }

      for (std::size_t i = 0; i < iedges.size(); ++i) {
        // std::cout << "back saved " << str(iedges[i].first) << std::endl;
        Point_2 future_point;
        Vector_2 future_direction;
        compute_future_point_and_direction(
            i, back, prev, iedges[i].first, future_point, future_direction);
        m_points[std::make_pair(pvertex.first, iedges[i].first)] = future_point;
        m_directions[std::make_pair(pvertex.first, iedges[i].first)] = future_direction;
      }

      PVertex previous = null_pvertex();
      for (std::size_t i = 0; i < crossed.size(); ++i) {
        if (i == 0) // crop
        {
          std::cout << "Cropping" << std::endl;
          PVertex cropped;
          if (prev_iedge != null_iedge() && prev_iedge == crossed[i]) {
            std::cout << "prev parallel case" << std::endl;

            cropped = prev;
            Point_2 future_point; Vector_2 future_direction;
            const auto pair = this->border_prev_and_next(prev);
            const auto pprev = pair.first;
            compute_future_point_and_direction(
              i, prev, pprev, prev_iedge, future_point, future_direction);
            future_points[i] = future_point;
            future_directions[i] = future_direction;

          } else {
            std::cout << "standard case" << std::endl;
            cropped = PVertex(support_plane_idx, support_plane(pvertex).split_edge(pvertex.second, prev.second));
            // future_point = future_points[i];
            // future_direction = future_directions[i];
          }

          const PEdge pedge(support_plane_idx, support_plane(pvertex).edge(pvertex.second, cropped.second));
          new_vertices.push_back(cropped);

          connect(pedge, crossed[i]);
          connect(cropped, crossed[i]);

          support_plane(cropped).set_point(cropped.second, future_points[i]);
          direction(cropped) = future_directions[i];
          previous = cropped;
          // std::cerr << "cropped point -> direction: " << point_2 (cropped) << " -> " << direction(cropped) << std::endl;
          std::cout << "cropped: " << point_3(cropped) << std::endl;
        }
        else // create triangle face
        {
          CGAL_assertion_msg(i == 1,
          "TODO: BACK, CAN WE HAVE MORE THAN 1 NEW FACE? IF YES, I SHOULD CHECK K FOR EACH!");
          bool is_occupied_edge, bbox_reached;
          std::tie(is_occupied_edge, bbox_reached) = is_occupied(pvertex, ivertex, crossed[i - 1]);
          // std::tie(is_occupied_edge, bbox_reached) = collision_occured(pvertex, crossed[0]);
          std::cout << "is already occupied / bbox: " << is_occupied_edge << "/" << bbox_reached << std::endl;

          // Stop.
          const auto pface = pface_of_pvertex(pvertex);
          std::cout << "k intersections: " << this->k(pface) << std::endl;
          if (bbox_reached) {
            std::cout << "stop bbox" << std::endl;
            CGAL_assertion_msg(false, "ERROR: THIS CASE CANNOT HAPPEN!");
            break;
          } else if (is_occupied_edge && this->k(pface) == 1) {
            std::cout << "stop k" << std::endl;
            break;
          }

          // Create a new face.
          std::cout << "adding new face!" << std::endl;
          if (is_occupied_edge && this->k(pface) > 1) {
            std::cout << "continue k > 1" << std::endl;
            this->k(pface)--;
          } else {
            std::cout << "continue k = 1" << std::endl;
          }
          CGAL_assertion(this->k(pface) >= 1);

          const PVertex propagated = add_pvertex(pvertex.first, future_points[i]);
          direction(propagated) = future_directions[i];
          new_vertices.push_back(propagated);

          std::cout << "propagated: " << point_3(propagated) << std::endl;
          const PFace new_pface = add_pface(std::array<PVertex, 3>{pvertex, propagated, previous});
          this->k(new_pface) = this->k(pface);
          previous = propagated;

          const PEdge pedge(support_plane_idx, support_plane(pvertex).edge(pvertex.second, propagated.second));
          connect(pedge, crossed[i]);
          connect(propagated, crossed[i]);
        }
      }

      if (crossed.size() < all_crossed.size()) {
        const std::size_t csize = crossed.size();
        const std::size_t asize = all_crossed.size();

        std::cout << "crossed size: " << csize << std::endl;
        std::cout << "all_crossed size: " << asize << std::endl;

        const std::size_t num_extra_faces = asize - csize;
        CGAL_assertion(num_extra_faces > 0);
        if (num_extra_faces == 1) {

          std::cout << "adding extra face!" << std::endl;
          PVertex propagated = find_opposite_pvertex(pvertex, ivertex, all_crossed.back());
          if (propagated == null_pvertex()) {
            CGAL_assertion_msg(false, "TODO: BACK, NULL PROPAGATED CASE!");
          } else {

            std::cout << "propagated: " << point_3(propagated) << std::endl;
            CGAL_assertion(num_extra_faces == 1);

            // Old code.
            // const PFace pface = pface_of_pvertex(pvertex);
            // const PFace new_pface = add_pface(std::array<PVertex, 3>{pvertex, propagated, previous});
            // this->k(new_pface) = this->k(pface);
            // previous = propagated;
            // continue;

            // New code!
            const auto opposite = this->opposite(all_crossed.back(), ivertex);
            const auto& mesh = this->mesh(pvertex);
            auto he = mesh.halfedge(pvertex.second, propagated.second);
            const PEdge qpedge(pvertex.first, mesh.edge(he));
            std::cout << "qpedge: " << segment_3(qpedge) << std::endl;

            PFace target_pface = null_pface();
            for (const auto pface : pfaces(pvertex.first)) {
              for (const auto pedge : pedges_of_pface(pface)) {
                if (pedge == qpedge) {
                  target_pface = pface;
                  break;
                }
              }
            }
            CGAL_assertion(target_pface != null_pface());

            const auto tt = pedges_of_pface(target_pface);
            std::vector<PEdge> pedges;
            pedges.reserve(tt.size());
            for (const auto t : tt) pedges.push_back(t);

            PEdge other_pedge = null_pedge();
            for (std::size_t j = 0; j < pedges.size(); ++j) {
              if (pedges[j] == qpedge) {
                const std::size_t jp = (j + 1) % pedges.size();
                const std::size_t jm = (j + pedges.size() - 1) % pedges.size();
                const auto& pedge1 = pedges[jm];
                const auto& pedge2 = pedges[jp];
                const auto iv1 = this->ivertex(this->target(pedge1));
                const auto iv2 = this->ivertex(this->source(pedge2));
                if (iv1 == opposite) {
                  CGAL_assertion(iv2 != opposite);
                  other_pedge = pedge1;
                  break;
                } else if (iv2 == opposite) {
                  CGAL_assertion(iv1 != opposite);
                  other_pedge = pedge2;
                  break;
                } else {
                  CGAL_assertion_msg(false, "ERROR: WRONG CASE!");
                }
              }
            }
            CGAL_assertion(other_pedge != null_pedge());
            std::cout << "other pedge: " << segment_3(other_pedge) << std::endl;

            IEdge other_iedge;
            const auto& iedges = support_plane(pvertex).iedges();
            CGAL_assertion(has_iedge(other_pedge));
            const auto query_iedge = this->iedge(other_pedge);
            for (const auto& iedge : iedges) {
              if (iedge == query_iedge) continue;
              if (this->source(iedge) == opposite || this->target(iedge) == opposite) {
                if (line_idx(query_iedge) == line_idx(iedge)) {
                  other_iedge = iedge;
                }
              }
            }
            CGAL_assertion(other_iedge != null_iedge());
            std::cout << "other iedge: " << segment_3(other_iedge) << std::endl;

            CGAL_assertion(m_points.find(std::make_pair(pvertex.first, other_iedge)) != m_points.end());
            CGAL_assertion(m_directions.find(std::make_pair(pvertex.first, other_iedge)) != m_directions.end());
            const Point_2 future_point = m_points.at(std::make_pair(pvertex.first, other_iedge));
            const Vector_2 future_direction = m_directions.at(std::make_pair(pvertex.first, other_iedge));

            auto tmp = future_direction;
            tmp = KSR::normalize(tmp);
            std::cout << "future: " << to_3d(pvertex.first, point_2(propagated) + m_current_time * tmp) << std::endl;

            const auto before = propagated;
            propagated = add_pvertex(propagated.first, future_point);
            direction(propagated) = future_direction;
            new_vertices.push_back(propagated);

            std::cout << "before: " << point_3(before) << std::endl;
            std::cout << "propagated: " << point_3(propagated) << std::endl;

            const PFace pface = pface_of_pvertex(pvertex);
            const PFace new_pface = add_pface(std::array<PVertex, 4>{pvertex, before, propagated, previous});
            this->k(new_pface) = this->k(pface);
            previous = propagated;

            // CGAL_assertion_msg(false, "DEBUG THIS CASE!");

            const PEdge pedge(support_plane_idx, support_plane(pvertex).edge(before.second, propagated.second));
            connect(pedge, other_iedge);
            connect(propagated, other_iedge);
            crossed.push_back(other_iedge);
          }
        } else {
          CGAL_assertion_msg(false, "TODO: BACK, CROSSED < LIMIT, MULTIPLE FACES!");
        }
      }

      if (crossed.size() == all_crossed.size()) {
        // continue...
      }

      if (crossed.size() > all_crossed.size()) {
        // continue ..
        // std::cout << "crossed size: " << crossed.size() << std::endl;
        // std::cout << "all crossed size: " << all_crossed.size() << std::endl;
        // CGAL_assertion_msg(false, "TODO: BACK CROSSED > LIMIT!");
      }
    }
    else if (front_constrained) // Border case
    {
      std::cout << "*** Front border case" << std::endl;

      CGAL_assertion(has_iedge(pvertex));
      // std::ofstream("limit.polylines.txt")
      // << "2 " << segment_3(iedge(pvertex)) << std::endl;
      const KSR::size_t other_side_limit = line_idx(pvertex);

      // const Direction_2 dir(point_2(next) - point_2(pvertex));

      const FT next_time = last_event_time(next);
      std::cout << next_time << std::endl;
      std::cout << m_current_time << std::endl;
      CGAL_assertion(next_time < m_current_time);
      CGAL_assertion(next_time >= FT(0));

      const auto pn_last = point_2(next, next_time);
      const auto pn_curr = point_2(next, m_current_time);
      const auto dirn = Vector_2(pn_last, pn_curr);
      const auto tmp_next = pn_curr - dirn / FT(10);

      const Direction_2 tmp_dir(tmp_next - point_2(pvertex.first, ivertex));
      // std::cout << "tmp_dir: " << to_3d(next.first, tmp_next) << std::endl;

      if (was_swapped) {
        std::reverse(iedges.begin(), iedges.end());
      }

      // std::cout << "initial iedges: " << std::endl;
      // for (const auto& iedge : iedges) {
      //   std::cout << segment_3(iedge.first) << std::endl;
      // }

      KSR::size_t first_idx = KSR::no_element();
      for (std::size_t i = 0; i < iedges.size(); ++ i)
      {
        if (!was_swapped) {
          if (tmp_dir.counterclockwise_in_between(
            iedges[i].second, iedges[(i + 1) % iedges.size()].second)) {
            first_idx = (i + 1) % iedges.size();
            break;
          }
        } else {
          if (tmp_dir.counterclockwise_in_between(
            iedges[(i + 1) % iedges.size()].second, iedges[i].second)) {
            first_idx = (i + 1) % iedges.size();
            break;
          }
        }
      }

      // std::cout << "first: " << segment_3(iedges[first_idx].first) << std::endl;
      // std::ofstream("first.polylines.txt")
      // << "2 " << segment_3(iedges[first_idx].first) << std::endl;

      CGAL_assertion(first_idx != KSR::no_element());
      crossed.clear();

      KSR::size_t iedge_idx = first_idx;
      std::size_t iter = 0;
      while (true) {
        const IEdge& iedge = iedges[iedge_idx].first;
        bool collision, bbox_reached;
        std::tie(collision, bbox_reached) = collision_occured(pvertex, iedge);
        const bool limit_reached = (line_idx(iedge) == other_side_limit);
        std::cout << "limit/bbox: " << limit_reached << "/" << bbox_reached << std::endl;

        // std::cout << "next: " << segment_3(iedge) << std::endl;
        // std::ofstream("next" + std::to_string(iter) + ".polylines.txt")
        // << "2 " << segment_3(iedge) << std::endl;
        crossed.push_back(iedge);
        if (limit_reached || bbox_reached) {
          break;
        }
        iedge_idx = (iedge_idx + 1) % iedges.size();
        if (iter == 100) {
          CGAL_assertion_msg(false, "ERROR: FRONT STD WHY SO MANY ITERATIONS?");
        } ++iter;
      }

      std::vector<IEdge> all_crossed;
      iedge_idx = first_idx; iter = 0;
      while (true) {
        const IEdge& iedge = iedges[iedge_idx].first;
        bool limit_reached, bbox_reached;
        std::tie(limit_reached, bbox_reached) = is_occupied(pvertex, iedge);
        all_crossed.push_back(iedge);
        if (limit_reached || bbox_reached) {
          break;
        }
        iedge_idx = (iedge_idx + 1) % iedges.size();
        if (iter == 100) {
          CGAL_assertion_msg(false, "ERROR: FRONT LIMIT WHY SO MANY ITERATIONS?");
        } ++iter;
      }

      CGAL_assertion(crossed.size() != 0);
      std::cerr << "IEdges crossed = " << crossed.size() << std::endl;
      for (const auto& iedge : crossed)
        std::cout << segment_3(iedge) << std::endl;
      CGAL_assertion(crossed[0] == all_crossed[0]);

      std::vector<Point_2> future_points;
      std::vector<Vector_2> future_directions;

      if (crossed.size() > all_crossed.size()) {
        future_points.resize(crossed.size());
        future_directions.resize(crossed.size());
        for (std::size_t i = 0; i < crossed.size(); ++i) {
          const bool is_parallel = compute_future_point_and_direction(
            i, front, next, crossed[i], future_points[i], future_directions[i]);

          if (is_parallel) {
            if (is_intersecting_iedge(min_time, max_time, next, crossed[i])) {
              next_iedge = crossed[i];
            }
          }
        }
      } else {
        future_points.resize(all_crossed.size());
        future_directions.resize(all_crossed.size());
        for (std::size_t i = 0; i < all_crossed.size(); ++i) {
          const bool is_parallel = compute_future_point_and_direction(
            i, front, next, all_crossed[i], future_points[i], future_directions[i]);

          if (is_parallel) {
            if (is_intersecting_iedge(min_time, max_time, next, all_crossed[i])) {
              next_iedge = all_crossed[i];
            }
          }
        }
      }

      for (std::size_t i = 0; i < iedges.size(); ++i) {
        // std::cout << "front saved " << str(iedges[i].first) << std::endl;
        Point_2 future_point;
        Vector_2 future_direction;
        compute_future_point_and_direction(
            i, front, next, iedges[i].first, future_point, future_direction);
        m_points[std::make_pair(pvertex.first, iedges[i].first)] = future_point;
        m_directions[std::make_pair(pvertex.first, iedges[i].first)] = future_direction;
      }

      PVertex previous = null_pvertex();
      for (std::size_t i = 0; i < crossed.size(); ++i) {
        if (i == 0) // crop
        {
          std::cout << "Cropping" << std::endl;
          PVertex cropped;
          if (next_iedge != null_iedge() && next_iedge == crossed[i]) {
            std::cout << "next parallel case" << std::endl;

            cropped = next;
            Point_2 future_point; Vector_2 future_direction;
            const auto pair = this->border_prev_and_next(next);
            const auto nnext = pair.second;
            compute_future_point_and_direction(
              i, next, nnext, next_iedge, future_point, future_direction);
            future_points[i] = future_point;
            future_directions[i] = future_direction;

          } else {
            std::cout << "standard case" << std::endl;
            cropped = PVertex(support_plane_idx, support_plane(pvertex).split_edge(pvertex.second, next.second));
            // future_point = future_points[i];
            // future_direction = future_directions[i];
          }

          const PEdge pedge(support_plane_idx, support_plane(pvertex).edge(pvertex.second, cropped.second));
          CGAL_assertion(cropped != pvertex);
          new_vertices.push_back(cropped);

          connect(pedge, crossed[i]);
          connect(cropped, crossed[i]);

          support_plane(cropped).set_point(cropped.second, future_points[i]);
          direction(cropped) = future_directions[i];
          previous = cropped;
          // std::cerr << point_2 (cropped) << " -> " << direction(cropped) << std::endl;
        }
        else // create triangle face
        {
          CGAL_assertion_msg(i == 1,
          "TODO: FRONT, CAN WE HAVE MORE THAN 1 NEW FACE? IF YES, I SHOULD CHECK K FOR EACH!");
          bool is_occupied_edge, bbox_reached;
          std::tie(is_occupied_edge, bbox_reached) = is_occupied(pvertex, ivertex, crossed[i - 1]);
          // std::tie(is_occupied_edge, bbox_reached) = collision_occured(pvertex, crossed[0]);
          std::cout << "is already occupied / bbox: " << is_occupied_edge << "/" << bbox_reached << std::endl;

          // Stop.
          const auto pface = pface_of_pvertex(pvertex);
          std::cout << "k intersections: " << this->k(pface) << std::endl;
          if (bbox_reached) {
            std::cout << "stop bbox" << std::endl;
            CGAL_assertion_msg(false, "ERROR: THIS CASE CANNOT HAPPEN!");
            break;
          } else if (is_occupied_edge && this->k(pface) == 1) {
            std::cout << "stop k" << std::endl;
            break;
          }

          // Create a new face.
          std::cout << "adding new face!" << std::endl;
          if (is_occupied_edge && this->k(pface) > 1) {
            std::cout << "continue k > 1" << std::endl;
            this->k(pface)--;
          } else {
            std::cout << "continue k = 1" << std::endl;
          }
          CGAL_assertion(this->k(pface) >= 1);

          const PVertex propagated = add_pvertex(pvertex.first, future_points[i]);
          direction(propagated) = future_directions[i];
          new_vertices.push_back(propagated);

          std::cout << "propagated: " << point_3(propagated) << std::endl;
          const PFace new_pface = add_pface(std::array<PVertex, 3>{pvertex, previous, propagated});
          this->k(new_pface) = this->k(pface);
          previous = propagated;

          const PEdge pedge(support_plane_idx, support_plane(pvertex).edge(pvertex.second, propagated.second));
          connect(pedge, crossed[i]);
          connect(propagated, crossed[i]);
        }
      }

      if (crossed.size() < all_crossed.size()) {
        const std::size_t csize = crossed.size();
        const std::size_t asize = all_crossed.size();

        std::cout << "crossed size: " << csize << std::endl;
        std::cout << "all_crossed size: " << asize << std::endl;

        const std::size_t num_extra_faces = asize - csize;
        CGAL_assertion(num_extra_faces != 0);

        bool is_ok = true;
        if (num_extra_faces == 2) {

          PVertex propagated = find_opposite_pvertex(pvertex, ivertex, all_crossed.back());
          const auto opposite = this->opposite(all_crossed.back(), ivertex);
          const auto& mesh = this->mesh(pvertex);
          auto he = mesh.halfedge(pvertex.second, propagated.second);
          const PEdge qpedge(pvertex.first, mesh.edge(he));
          std::cout << "qpedge: " << segment_3(qpedge) << std::endl;

          PFace target_pface = null_pface();
          for (const auto pface : pfaces(pvertex.first)) {
            for (const auto pedge : pedges_of_pface(pface)) {
              if (pedge == qpedge) {
                target_pface = pface;
                break;
              }
            }
          }
          CGAL_assertion(target_pface != null_pface());

          const auto tt = pedges_of_pface(target_pface);
          std::vector<PEdge> pedges;
          pedges.reserve(tt.size());
          for (const auto t : tt) pedges.push_back(t);

          PEdge other_pedge = null_pedge();
          for (std::size_t j = 0; j < pedges.size(); ++j) {
            if (pedges[j] == qpedge) {
              const std::size_t jp = (j + 1) % pedges.size();
              const std::size_t jm = (j + pedges.size() - 1) % pedges.size();
              const auto& pedge1 = pedges[jm];
              const auto& pedge2 = pedges[jp];
              const auto iv1 = this->ivertex(this->target(pedge1));
              const auto iv2 = this->ivertex(this->source(pedge2));
              if (iv1 == opposite) {
                CGAL_assertion(iv2 != opposite);
                other_pedge = pedge1;
                break;
              } else if (iv2 == opposite) {
                CGAL_assertion(iv1 != opposite);
                other_pedge = pedge2;
                break;
              } else {
                CGAL_assertion_msg(false, "ERROR: WRONG CASE!");
              }
            }
          }
          CGAL_assertion(other_pedge != null_pedge());
          std::cout << "other pedge: " << segment_3(other_pedge) << std::endl;

          IEdge other_iedge;
          const auto& iedges = support_plane(pvertex).iedges();
          CGAL_assertion(has_iedge(other_pedge));
          const auto query_iedge = this->iedge(other_pedge);
          for (const auto& iedge : iedges) {
            if (iedge == query_iedge) continue;
            if (this->source(iedge) == opposite || this->target(iedge) == opposite) {
              if (line_idx(query_iedge) == line_idx(iedge)) {
                other_iedge = iedge;
              }
            }
          }
          CGAL_assertion(other_iedge != null_iedge());
          std::cout << "other iedge: " << segment_3(other_iedge) << std::endl;

          if (!is_occupied(propagated, other_iedge).first) {
            is_ok = false;
          }
        }

        if (is_ok) {
        if (num_extra_faces < 3) {

          // CGAL_assertion_msg(false, "TODO: FRONT, CROSSED < LIMIT, 1 or 2 FACES!");
          CGAL_assertion(future_points.size() == asize);
          CGAL_assertion(future_directions.size() == asize);

          for (std::size_t i = csize; i < asize; ++i) {
            std::cout << "adding extra face!" << std::endl;

            PVertex propagated = find_opposite_pvertex(pvertex, ivertex, all_crossed[i]);
            if (propagated == null_pvertex()) {

              const Line_2 iedge_line = segment_2(pvertex.first, all_crossed[i]).supporting_line();
              const Point_2 pinit = iedge_line.projection(point_2(pvertex));

              const IVertex opposite = this->opposite(all_crossed[i], ivertex);
              Point_2 future_point = to_2d(pvertex.first, opposite);

              Vector_2 future_direction = Vector_2(pinit, future_point);
              future_point = pinit - m_current_time * future_direction;

              // auto tmp = future_direction;
              // tmp = KSR::normalize(tmp);
              // std::cout << "future: " << to_3d(pvertex.first, pinit + m_current_time * tmp) << std::endl;

              const FT dot_product = future_direction * future_directions[i];
              if (dot_product < FT(0)) {
                future_direction = -future_directions[i];
                future_point = pinit - m_current_time * future_direction;
              } else {
                future_direction = future_directions[i];
                future_point = future_points[i];
              }

              // future_point = future_points[i]; // old, does not work
              // future_direction = future_directions[i]; // old, does not work

              propagated = add_pvertex(pvertex.first, future_point);
              direction(propagated) = future_direction;
              new_vertices.push_back(propagated);

              std::cout << "propagated null: " << point_3(propagated) << std::endl;
              const PFace pface = pface_of_pvertex(pvertex);
              const PFace new_pface = add_pface(std::array<PVertex, 3>{pvertex, previous, propagated});
              this->k(new_pface) = this->k(pface);
              previous = propagated;

              const PEdge pedge(support_plane_idx, support_plane(pvertex).edge(pvertex.second, propagated.second));
              connect(pedge, all_crossed[i]);
              connect(propagated, all_crossed[i]);
              crossed.push_back(all_crossed[i]); // remove events from this one

              // CGAL_assertion_msg(false, "TODO: FRONT, NULL PROPAGATED CASE!");

            } else {

              std::cout << "propagated std: " << point_3(propagated) << std::endl;
              CGAL_assertion(i == asize - 1);

              // Old code!
              // const PFace pface = pface_of_pvertex(pvertex);
              // const PFace new_pface = add_pface(std::array<PVertex, 3>{pvertex, previous, propagated});
              // this->k(new_pface) = this->k(pface);
              // previous = propagated;
              // continue;

              // New code!
              const auto opposite = this->opposite(all_crossed[i], ivertex);
              const auto& mesh = this->mesh(pvertex);
              auto he = mesh.halfedge(pvertex.second, propagated.second);
              const PEdge qpedge(pvertex.first, mesh.edge(he));
              std::cout << "qpedge: " << segment_3(qpedge) << std::endl;

              PFace target_pface = null_pface();
              for (const auto pface : pfaces(pvertex.first)) {
                for (const auto pedge : pedges_of_pface(pface)) {
                  if (pedge == qpedge) {
                    target_pface = pface;
                    break;
                  }
                }
              }
              CGAL_assertion(target_pface != null_pface());

              const auto tt = pedges_of_pface(target_pface);
              std::vector<PEdge> pedges;
              pedges.reserve(tt.size());
              for (const auto t : tt) pedges.push_back(t);

              PEdge other_pedge = null_pedge();
              for (std::size_t j = 0; j < pedges.size(); ++j) {
                if (pedges[j] == qpedge) {
                  const std::size_t jp = (j + 1) % pedges.size();
                  const std::size_t jm = (j + pedges.size() - 1) % pedges.size();
                  const auto& pedge1 = pedges[jm];
                  const auto& pedge2 = pedges[jp];
                  const auto iv1 = this->ivertex(this->target(pedge1));
                  const auto iv2 = this->ivertex(this->source(pedge2));
                  if (iv1 == opposite) {
                    CGAL_assertion(iv2 != opposite);
                    other_pedge = pedge1;
                    break;
                  } else if (iv2 == opposite) {
                    CGAL_assertion(iv1 != opposite);
                    other_pedge = pedge2;
                    break;
                  } else {
                    CGAL_assertion_msg(false, "ERROR: WRONG CASE!");
                  }
                }
              }
              CGAL_assertion(other_pedge != null_pedge());
              std::cout << "other pedge: " << segment_3(other_pedge) << std::endl;

              IEdge other_iedge;
              const auto& iedges = support_plane(pvertex).iedges();
              CGAL_assertion(has_iedge(other_pedge));
              const auto query_iedge = this->iedge(other_pedge);
              for (const auto& iedge : iedges) {
                if (iedge == query_iedge) continue;
                if (this->source(iedge) == opposite || this->target(iedge) == opposite) {
                  if (line_idx(query_iedge) == line_idx(iedge)) {
                    other_iedge = iedge;
                  }
                }
              }
              CGAL_assertion(other_iedge != null_iedge());
              std::cout << "other iedge: " << segment_3(other_iedge) << std::endl;

              // if (!is_occupied(propagated, other_iedge).first) {
              //   break;
              // }

              CGAL_assertion(m_points.find(std::make_pair(pvertex.first, other_iedge)) != m_points.end());
              CGAL_assertion(m_directions.find(std::make_pair(pvertex.first, other_iedge)) != m_directions.end());
              Point_2 future_point = m_points.at(std::make_pair(pvertex.first, other_iedge));
              Vector_2 future_direction = m_directions.at(std::make_pair(pvertex.first, other_iedge));

              // const Point_2 pinit = point_2(propagated);
              // Point_2 future_point = point_2(pvertex.first, this->opposite(other_iedge, opposite));
              // Vector_2 future_direction = Vector_2(pinit, future_point);
              // future_point = pinit - m_current_time * future_direction;

              auto tmp = future_direction;
              tmp = KSR::normalize(tmp);
              std::cout << "future: " << to_3d(pvertex.first, point_2(propagated) + m_current_time * tmp) << std::endl;

              const auto before = propagated;
              propagated = add_pvertex(propagated.first, future_point);
              direction(propagated) = future_direction;


              std::cout << "before: " << point_3(before) << std::endl;
              std::cout << "propagated: " << point_3(propagated) << std::endl;

              std::size_t count = 0;
              for (const IVertex& iver : { this->source(other_iedge), this->target(other_iedge) }) {
                const Point_2 pi = to_2d(propagated.first, iver);
                const Segment_2 sv(
                  point_2(propagated, m_current_time),
                  point_2(propagated, max_time));
                if (sv.to_vector() * Vector_2(sv.source(), pi) < FT(0)) {
                  ++count;
                  continue;
                }
              }
              if (count == 2) {

                const PFace pface = pface_of_pvertex(pvertex);
                const PFace new_pface = add_pface(std::array<PVertex, 3>{pvertex, previous, before});
                this->k(new_pface) = this->k(pface);
                previous = before;
                support_plane(pvertex.first).remove_vertex(propagated.second);
                // break;

                // const Point_2 pinit = point_2(before);
                // future_point = point_2(pvertex.first, this->opposite(other_iedge, opposite));
                // future_direction = Vector_2(pinit, future_point);
                // future_point = pinit - m_current_time * future_direction;

                // support_plane(propagated).set_point(propagated.second, future_point);
                // direction(propagated) = future_direction;

                CGAL_assertion_msg(false, "TODO!");
              } else {
                new_vertices.push_back(propagated);
              }

              const PFace pface = pface_of_pvertex(pvertex);
              const PFace new_pface = add_pface(std::array<PVertex, 4>{pvertex, previous, propagated, before});
              this->k(new_pface) = this->k(pface);
              previous = propagated;

              // CGAL_assertion_msg(false, "DEBUG THIS CASE!");

              const PEdge pedge(support_plane_idx, support_plane(pvertex).edge(before.second, propagated.second));
              connect(pedge, other_iedge);
              connect(propagated, other_iedge);
              crossed.push_back(other_iedge);
            }
          }
          // CGAL_assertion_msg(false, "TODO: TEST THIS LOOP!");
        } else {

          std::cout << "crossed size: " << crossed.size() << std::endl;
          std::cout << "all crossed size: " << all_crossed.size() << std::endl;
          for (const auto& iedge : all_crossed) {
            std::cout << segment_3(iedge) << std::endl;
          }

          CGAL_assertion_msg(false, "TODO: FRONT, CROSSED < LIMIT, MULTIPLE FACES!");
        }
        } else {

        }
      }

      if (crossed.size() == all_crossed.size()) {
        // continue...
      }

      if (crossed.size() > all_crossed.size()) {
        // continue ..
        // std::cout << "crossed size: " << crossed.size() << std::endl;
        // std::cout << "all crossed size: " << all_crossed.size() << std::endl;
        // CGAL_assertion_msg(false, "TODO: FRONT, CROSSED > LIMIT!");
      }
    }
    else // Open case
    {
      std::cout << "*** Open case" << std::endl;

      // const Direction_2 dir_prev(point_2(prev) - point_2(pvertex));
      // const Direction_2 dir_next(point_2(next) - point_2(pvertex));

      const FT prev_time = last_event_time(prev);
      const FT next_time = last_event_time(next);
      CGAL_assertion(prev_time < m_current_time);
      CGAL_assertion(next_time < m_current_time);
      CGAL_assertion(prev_time >= FT(0));
      CGAL_assertion(next_time >= FT(0));

      const auto pp_last = point_2(prev, prev_time);
      const auto pp_curr = point_2(prev, m_current_time);
      const auto dirp = Vector_2(pp_last, pp_curr);
      const auto tmp_prev = pp_curr - dirp / FT(10);

      const auto pn_last = point_2(next, next_time);
      const auto pn_curr = point_2(next, m_current_time);
      const auto dirn = Vector_2(pn_last, pn_curr);
      const auto tmp_next = pn_curr - dirn / FT(10);

      const Direction_2 dir_prev(tmp_prev - point_2(pvertex.first, ivertex));
      const Direction_2 dir_next(tmp_next - point_2(pvertex.first, ivertex));

      // std::cout << to_3d(prev.first, tmp_prev) << std::endl;
      // std::cout << to_3d(next.first, tmp_next) << std::endl;

      // std::cout << "initial iedges: " << std::endl;
      // for (const auto& iedge : iedges) {
      //   std::cout << segment_3(iedge.first) << std::endl;
      // }

      KSR::size_t first_idx = KSR::no_element();
      for (std::size_t i = 0; i < iedges.size(); ++i) {
        if (dir_next.counterclockwise_in_between(
          iedges[i].second, iedges[(i + 1) % iedges.size()].second)) {

          first_idx = (i + 1) % iedges.size();
          break;
        }
      }

      CGAL_assertion(first_idx != KSR::no_element());
      crossed.clear();

      // std::cout << "first: " << segment_3(iedges[first_idx].first) << std::endl;
      // std::ofstream("first.polylines.txt")
      // << "2 " << segment_3(iedges[first_idx].first) << std::endl;

      KSR::size_t iedge_idx = first_idx;
      std::size_t iter = 0;
      while (true)
      {
        const IEdge& iedge = iedges[iedge_idx].first;
        const Direction_2& dir = iedges[iedge_idx].second;

        if (!dir.counterclockwise_in_between (dir_next, dir_prev))
          break;

        // std::cout << "next: " << segment_3(iedge) << std::endl;
        // std::ofstream("next" + std::to_string(iter) + ".polylines.txt")
        // << "2 " << segment_3(iedge) << std::endl;
        crossed.push_back(iedge);

        iedge_idx = (iedge_idx + 1) % iedges.size();

        if (iter == 100) {
          CGAL_assertion_msg(false, "ERROR: OPEN WHY SO MANY ITERATIONS?");
        }
        ++iter;
      }

      CGAL_assertion(crossed.size() != 0);
      std::cerr << "IEdges crossed = " << crossed.size() << std::endl;
      for (const auto& iedge : crossed)
        std::cout << segment_3(iedge) << std::endl;

      std::vector<Point_2> future_points(crossed.size());
      std::vector<Vector_2> future_directions(crossed.size());
      for (std::size_t i = 0; i < crossed.size(); ++i) {
        const bool is_parallel = compute_future_point_and_direction(
          pvertex, prev, next, crossed[i], future_points[i], future_directions[i]);

        if (is_parallel) {
          if (is_intersecting_iedge(min_time, max_time, prev, crossed[i])) {
            prev_iedge = crossed[i];
          }
          if (is_intersecting_iedge(min_time, max_time, next, crossed[i])) {
            next_iedge = crossed[i];
          }
        }
      }
      CGAL_assertion(future_points.size() == crossed.size());
      CGAL_assertion(future_directions.size() == crossed.size());

      for (std::size_t i = 0; i < iedges.size(); ++i) {
        // std::cout << "open saved " << str(iedges[i].first) << std::endl;
        Point_2 future_point;
        Vector_2 future_direction;
        compute_future_point_and_direction(
            pvertex, prev, next, iedges[i].first, future_point, future_direction);
        m_points[std::make_pair(pvertex.first, iedges[i].first)] = future_point;
        m_directions[std::make_pair(pvertex.first, iedges[i].first)] = future_direction;
      }

      {
        PVertex cropped;
        if (next_iedge != null_iedge() && next_iedge == crossed.front()) {
          std::cout << "next parallel case" << std::endl;

          cropped = next;
          Point_2 future_point; Vector_2 future_direction;
          const auto pair = this->border_prev_and_next(next);
          const auto nnext = pair.second;
          compute_future_point_and_direction(
            0, next, nnext, next_iedge, future_point, future_direction);
          future_points[0] = future_point;
          future_directions[0] = future_direction;

        } else {
          std::cout << "standard case" << std::endl;
          cropped = PVertex(support_plane_idx, support_plane(pvertex).split_edge(pvertex.second, next.second));
          // future_point = future_points.front();
          // future_direction = future_directions.front();
        }

        const PEdge pedge(support_plane_idx, support_plane(pvertex).edge(pvertex.second, cropped.second));
        new_vertices.push_back(cropped);

        connect(pedge, crossed.front());
        connect(cropped, crossed.front());

        support_plane(cropped).set_point(cropped.second, future_points.front());
        direction(cropped) = future_directions.front();
        std::cout << direction(cropped) << std::endl;
        std::cout << "cropped 1: " << point_3(cropped) << std::endl;
      }

      for (std::size_t i = 1; i < crossed.size() - 1; ++i)
      {
        const PVertex propagated = add_pvertex(pvertex.first, future_points[i]);
        direction(propagated) = future_directions[i];
        connect(propagated, crossed[i]);
        new_vertices.push_back(propagated);
        std::cout << "propagated " << std::to_string(i) << ": " << point_3(propagated) << std::endl;
      }

      {
        PVertex cropped;
        if (prev_iedge != null_iedge() && prev_iedge == crossed.back()) {
          std::cout << "prev parallel case" << std::endl;

          cropped = prev;
          Point_2 future_point; Vector_2 future_direction;
          const auto pair = this->border_prev_and_next(prev);
          const auto pprev = pair.first;
          compute_future_point_and_direction(
            0, prev, pprev, prev_iedge, future_point, future_direction);
          future_points[future_points.size() - 1] = future_point;
          future_directions[future_directions.size() - 1] = future_direction;

        } else {
          std::cout << "standard case" << std::endl;
          cropped = PVertex(support_plane_idx, support_plane(pvertex).split_edge(pvertex.second, prev.second));
          // future_point = future_points.back();
          // future_direction = future_directions.back();
        }

        const PEdge pedge(support_plane_idx, support_plane(pvertex).edge(pvertex.second, cropped.second));
        new_vertices.push_back(cropped);

        connect(pedge, crossed.back());
        connect(cropped, crossed.back());

        support_plane(cropped).set_point(cropped.second, future_points.back());
        direction(cropped) = future_directions.back();
        std::cout << direction(cropped) << std::endl;
        std::cout << "cropped 2: " << point_3(cropped) << std::endl;
      }

      std::cerr << new_vertices.size() << " new vertice(s)" << std::endl;
      CGAL_assertion(new_vertices.size() == crossed.size());

      bool is_occupied_edge_back, bbox_reached_back;
      std::tie(is_occupied_edge_back, bbox_reached_back) = is_occupied(pvertex, ivertex, crossed.back());
      // std::tie(is_occupied_edge_back, bbox_reached_back) = collision_occured(pvertex, crossed.back());
      std::cout << "is already occupied back / bbox: " << is_occupied_edge_back << "/" << bbox_reached_back << std::endl;

      bool is_occupied_edge_front, bbox_reached_front;
      std::tie(is_occupied_edge_front, bbox_reached_front) = is_occupied(pvertex, ivertex, crossed.front());
      // std::tie(is_occupied_edge_front, bbox_reached_front) = collision_occured(pvertex, crossed.front());
      std::cout << "is already occupied front / bbox: " << is_occupied_edge_front << "/" << bbox_reached_front << std::endl;

      const auto pface = pface_of_pvertex(pvertex);
      std::cout << "k intersections: " << this->k(pface) << std::endl;
      if (bbox_reached_back) {

        CGAL_assertion(bbox_reached_front);
        std::cout << "stop bbox back" << std::endl;

      } else if (bbox_reached_front) {

        CGAL_assertion(bbox_reached_back);
        std::cout << "stop bbox front" << std::endl;

      } else if ((is_occupied_edge_back && is_occupied_edge_front) && this->k(pface) == 1) {

        std::cout << "stop back && front k = 1" << std::endl;

      } else if ((is_occupied_edge_back && is_occupied_edge_front) && this->k(pface) > 1) {

        this->k(pface)--;
        CGAL_assertion(this->k(pface) >= 1);
        add_new_faces(this->k(pface), pvertex, ivertex, new_vertices, pface, crossed);
        std::cout << "continue back && front k > 1" << std::endl;

      } else if ((!is_occupied_edge_back && !is_occupied_edge_front)) {

        add_new_faces(this->k(pface), pvertex, ivertex, new_vertices, pface, crossed);
        std::cout << "continue !back && !front" << std::endl;

      } else if (is_occupied_edge_back || is_occupied_edge_front) {

        add_new_faces(this->k(pface), pvertex, ivertex, new_vertices, pface, crossed);
        std::cout << "continue back || front" << std::endl;

        // std::cout << "centroid pface: " << centroid_of_pface(pface) << std::endl;
        // CGAL_assertion_msg(false, "TEST THIS CASE: BACK || FRONT!");

      } else {
        CGAL_assertion_msg(false, "TODO: ADD NEW OPEN CASE! DO NOT FORGET TO UPDATE K!");
      }

      // for (std::size_t i = 1; i < crossed.size() - 1; ++i) {
      //   PEdge pedge(support_plane_idx, support_plane(pvertex).edge(pvertex.second, new_vertices[i].second));
      //   connect(pedge, crossed[i]);
      //   connect(new_vertices[i], crossed[i]);
      // }
    }

    support_plane(support_plane_idx).remove_vertex(front.second);
    support_plane(support_plane_idx).remove_vertex(back.second);

    // push also remaining vertex so that its events are recomputed
    // std::cout << "pushing new pv: " << str(pvertex) << std::endl;
    // std::cout << "pv direction: " << direction(pvertex) << std::endl;
    new_vertices.push_back(pvertex);
    crossed.push_back(iedge(pvertex));

    std::cout << "*** New vertices:";
    for (const PVertex& pv : new_vertices)
      std::cout << " " << str(pv);
    std::cout << std::endl;

    // for (const PVertex& pv : new_vertices)
    //   std::cout << point_3(pv) << std::endl;

    // if (has_iedge(prev) && !is_frozen(prev)) {
    // // if (iedge(prev) != iedge(pvertex)) {
    //   std::cout << "pushing new prev: " << str(prev) << std::endl;
    //   new_vertices.push_back (prev);
    // }

    // if (has_iedge(next) && !is_frozen(next)) {
    // // if (back_constrained) {
    //   std::cout << "pushing new next: " << str(next) << std::endl;
    //   new_vertices.push_back (next);
    // }

    return new_vertices;
  }

  void add_new_faces(
    const unsigned int k,
    const PVertex& pvertex,
    const IVertex& ivertex,
    const std::vector<PVertex>& new_vertices,
    const PFace& pface,
    const std::vector<IEdge>& crossed) {

    CGAL_assertion(new_vertices.size() >= 2);
    CGAL_assertion(crossed.size() == new_vertices.size());

    std::size_t num_added_faces = 0;
    for (std::size_t i = 0; i < new_vertices.size() - 1; ++i) {

      if (i >= 1) {
        // bool is_occupied_edge, bbox_reached;
        // std::tie(is_occupied_edge, bbox_reached) = is_occupied(pvertex, ivertex, crossed[i]);

        // if (bbox_reached) return;
        // if (is_occupied_edge && this->k(pface) == 1) return;
        // if (is_occupied_edge && this->k(pface) > 1) {
        //   this->k(pface)--;
        //   CGAL_assertion(this->k(pface) >= 1);
        // }

        const PEdge pedge(pvertex.first, support_plane(pvertex).edge(pvertex.second, new_vertices[i].second));
        connect(pedge, crossed[i]);
        connect(new_vertices[i], crossed[i]);
      }

      std::cout << "adding new face" << std::endl;
      const PFace new_pface = add_pface(std::array<PVertex, 3>{new_vertices[i], new_vertices[i + 1], pvertex});
      this->k(new_pface) = k;
      ++num_added_faces;
    }
    CGAL_assertion(num_added_faces > 0);
    CGAL_assertion_msg(num_added_faces == 1,
    "TODO: OPEN, CAN WE HAVE MORE THAN 1 NEW FACE? IF YES, I SHOULD CHECK K FOR EACH!");
  }

  const PVertex find_opposite_pvertex(
    const PVertex& pvertex,
    const IVertex& ivertex,
    const IEdge& iedge) const {

    std::set<PEdge> pedges;
    const KSR::size_t support_plane_idx = pvertex.first;
    // std::cout << "query: " << segment_3(iedge) << " : " << str(iedge) << std::endl;
    for (const auto pedge : this->pedges(support_plane_idx)) {
      if (!has_iedge(pedge)) continue;
      // std::cout << "other: " << segment_3(pedge) << " : " << str(this->iedge(pedge)) << std::endl;
      if (this->iedge(pedge) == iedge) {
        pedges.insert(pedge);
      }
    }
    // CGAL_assertion(pedges.size() == 1);

    for (const auto& pedge : pedges) {
      CGAL_assertion(pedge != null_pedge());

      const PVertex source = this->source(pedge);
      const PVertex target = this->target(pedge);

      if (source == pvertex) {
        // std::cout << "here1" << std::endl;
        return target;
      }
      if (target == pvertex) {
        // std::cout << "here2" << std::endl;
        return source;
      }
      CGAL_assertion(source != pvertex && target != pvertex);

      if (this->ivertex(source) == ivertex) {
        // std::cout << "here3" << std::endl;
        return target;
      }
      if (this->ivertex(target) == ivertex) {
        // std::cout << "here4" << std::endl;
        return source;
      }
      CGAL_assertion(this->ivertex(source) != ivertex);
      CGAL_assertion(this->ivertex(target) != ivertex);

      CGAL_assertion_msg(false, "TODO: CAN WE HAVE THIS CASE?");

      // const auto s = point_2(source);
      // const auto t = point_2(target);
      // const auto ipoint = point_2(support_plane_idx, ivertex);
      // Vector_2 vec1(s, t);
      // Vector_2 vec2(s, ipoint);
      // vec1 = KSR::normalize(vec1);
      // vec2 = KSR::normalize(vec2);

      // const FT dot_product = vec1 * vec2;
      // if (dot_product < FT(0)) {
      //   return target;
      // } else {
      //   return source;
      // }
    }
    return null_pvertex();
  }

  void check_bbox() {

    for (KSR::size_t i = 0; i < 6; ++i) {
      const auto pfaces = this->pfaces(i);
      for (const auto pface : pfaces) {
        for (const auto pedge : pedges_of_pface(pface)) {
          CGAL_assertion_msg(has_iedge(pedge), "ERROR: BBOX EDGE IS MISSING AN IEDGE!");
        }
        for (const auto pvertex : pvertices_of_pface(pface)) {
          CGAL_assertion_msg(has_ivertex(pvertex), "ERROR: BBOX VERTEX IS MISSING AN IVERTEX!");
        }
      }
    }
  }

  void check_interior() {

    for (KSR::size_t i = 6; i < number_of_support_planes(); ++i) {
      const auto pfaces = this->pfaces(i);
      for (const auto pface : pfaces) {
        for (const auto pedge : pedges_of_pface(pface)) {
          if (!has_iedge(pedge)) {
            dump_pedge(pedge, "debug-pedge");
          }
          CGAL_assertion_msg(has_iedge(pedge), "ERROR: INTERIOR EDGE IS MISSING AN IEDGE!");
        }
        for (const auto pvertex : pvertices_of_pface(pface)) {
          CGAL_assertion_msg(has_ivertex(pvertex), "ERROR: INTERIOR VERTEX IS MISSING AN IVERTEX!");
        }
      }
    }
  }

  void check_vertices() {

    for (const auto vertex : m_intersection_graph.vertices()) {
      const auto nedges = m_intersection_graph.incident_edges(vertex);
      if (nedges.size() <= 2)
        std::cerr << "ERROR: current num edges = " << nedges.size() << std::endl;
      CGAL_assertion_msg(nedges.size() > 2,
      "ERROR: VERTEX MUST HAVE AT LEAST 3 NEIGHBORS!");
    }
  }

  void check_edges() {

    std::vector<PFace> nfaces;
    for (const auto edge : m_intersection_graph.edges()) {
      incident_faces(edge, nfaces);
      if (nfaces.size() == 1) {
        std::cout << segment_3(edge) << std::endl;
        std::cerr << "ERROR: current num faces = " << nfaces.size() << std::endl;
      }
      CGAL_assertion_msg(nfaces.size() != 1,
      "ERROR: EDGE MUST HAVE 0 OR AT LEAST 2 NEIGHBORS!");
    }
  }

  void check_faces() {

    for (std::size_t i = 0; i < number_of_support_planes(); ++i) {
      const auto pfaces = this->pfaces(i);
      for (const auto pface : pfaces) {
        const auto nvolumes = incident_volumes(pface);
        if (nvolumes.size() == 0 || nvolumes.size() > 2)
          std::cerr << "current num volumes = " << nvolumes.size() << std::endl;
        CGAL_assertion_msg(nvolumes.size() == 1 || nvolumes.size() == 2,
        "ERROR: FACE MUST HAVE 1 OR 2 NEIGHBORS!");
      }
    }
  }

  void create_polyhedrons() {

    std::cout.precision(20);
    // for (std::size_t i = 0; i < number_of_support_planes(); ++i)
    //   std::cout << "num faces sp " << i << ": " << pfaces(i).size() << std::endl;

    check_bbox();
    check_interior();
    check_vertices();
    check_edges();
    create_volumes();
    check_faces();
  }

  void create_volumes() {

    // Initialize an empty volume map.
    m_volumes.clear();
    std::map<int, Point_3> centroids;
    std::map<PFace, std::pair<int, int> > map_volumes;
    for (std::size_t i = 0; i < number_of_support_planes(); ++i) {
      const auto pfaces = this->pfaces(i);
      for (const auto pface : pfaces)
        map_volumes[pface] = std::make_pair(-1, -1);
    }

    // First, traverse only boundary volumes.
    bool is_found_new_volume = false;
    std::size_t volume_size = 0;
    int num_volumes  = 0;
    int volume_index = 0;
    for (std::size_t i = 0; i < 6; ++i) {
      const auto pfaces = this->pfaces(i);
      for (const auto pface : pfaces) {
        CGAL_assertion(pface.first < 6);
        std::tie(is_found_new_volume, volume_size) = traverse_boundary_volume(
          pface, volume_index, num_volumes, map_volumes, centroids);
        if (is_found_new_volume) {
          check_volume(volume_index, volume_size, map_volumes);
          ++volume_index;
        }
      }
    }
    std::cout << "* found boundary volumes: "<< volume_index << std::endl;
    num_volumes = volume_index;
    CGAL_assertion(num_volumes > 0);

    // Then traverse all other volumes if any.
    bool quit = true;
    do {
      quit = true;
      const int before = volume_index;
      for (std::size_t i = 6; i < number_of_support_planes(); ++i) {
        const auto pfaces = this->pfaces(i);
        for (const auto pface : pfaces) {
          CGAL_assertion(pface.first >= 6);
          std::tie(is_found_new_volume, volume_size) = traverse_interior_volume(
            pface, volume_index, num_volumes, map_volumes, centroids);
          if (is_found_new_volume) {
            quit = false;
            check_volume(volume_index, volume_size, map_volumes);
            ++volume_index;
          }
        }
      }
      const int after = volume_index;
      std::cout << "* found interior volumes: "<< after - before << std::endl;
      CGAL_assertion(after >= before);
      num_volumes = volume_index;

    } while (!quit);

    // Now, set final polyhedrons and their neighbors.
    for (const auto& item : map_volumes) {
      const auto& pface = item.first;
      const auto& pair  = item.second;

      if (pair.first == -1) {
        dump_pface(pface, "face-debug");
        std::cout << "DEBUG face: " << str(pface) << " "   << std::endl;
        std::cout << "DEBUG  map: " << pair.first << " : " << pair.second << std::endl;
      }

      CGAL_assertion(pair.first != -1);
      if (m_volumes.size() <= pair.first)
        m_volumes.resize(pair.first + 1);
      m_volumes[pair.first].add_pface(pface, pair.second);

      if (pface.first < 6 && pair.second == -1) continue;
      CGAL_assertion(pair.second != -1);
      if (m_volumes.size() <= pair.second)
        m_volumes.resize(pair.second + 1);
      m_volumes[pair.second].add_pface(pface, pair.first);
    }
    for (auto& volume : m_volumes)
      create_cell_pvertices(volume);
    std::cout << "* created polyhedrons: " << m_volumes.size() << std::endl;

    dump_polyhedrons(*this, "polyhedrons/final");
    for (std::size_t i = 0; i < m_volumes.size(); ++i) {
      const auto& volume = m_volumes[i];
      CGAL_assertion(volume.pfaces.size() > 3);
      std::cout <<
      " POLYHEDRON " << std::to_string(i) << ": "
      " pvertices: " << volume.pvertices.size() <<
      " pfaces: "    << volume.pfaces.size()    << std::endl;
    }

    CGAL_assertion(m_volumes.size() == centroids.size());
    for (std::size_t i = 0; i < m_volumes.size(); ++i) {
      auto& volume = m_volumes[i];
      volume.set_index(i);
      volume.set_centroid(centroids.at(i));
    }
  }

  void check_volume(
    const int volume_index,
    const std::size_t volume_size,
    const std::map<PFace, std::pair<int, int> >& map_volumes) const {

    std::vector<PFace> pfaces;
    for (const auto& item : map_volumes) {
      const auto& pface = item.first;
      const auto& pair  = item.second;
      if (pair.first == volume_index || pair.second == volume_index) {
        pfaces.push_back(pface);
      }
    }

    const bool is_broken_volume = is_volume_degenerate(pfaces);
    if (is_broken_volume) {
      dump_polyhedron(*this, pfaces, "polyhedrons/degenerate");
    }
    CGAL_assertion(!is_broken_volume);
    CGAL_assertion(pfaces.size() == volume_size);
  }

  const bool is_volume_degenerate(
    const std::vector<PFace>& pfaces) const {

    for (const auto& pface : pfaces) {
      const auto pedges = pedges_of_pface(pface);
      const std::size_t n = pedges.size();

      std::size_t count = 0;
      for (const auto pedge : pedges) {
        CGAL_assertion(has_iedge(pedge));
        const auto iedge = this->iedge(pedge);
        const std::size_t num_found = find_adjacent_pfaces(pface, iedge, pfaces);
        if (num_found == 1) ++count;
      }
      if (count != n) return true;
    }
    return false;
  }

  const std::size_t find_adjacent_pfaces(
    const PFace& current,
    const IEdge& query,
    const std::vector<PFace>& pfaces) const {

    std::size_t num_found = 0;
    for (const auto& pface : pfaces) {
      if (pface == current) continue;
      const auto pedges = pedges_of_pface(pface);
      for (const auto pedge : pedges) {
        CGAL_assertion(has_iedge(pedge));
        const auto iedge = this->iedge(pedge);
        if (iedge == query) ++num_found;
      }
    }
    return num_found;
  }

  const std::pair<bool, std::size_t> traverse_boundary_volume(
    const PFace& pface,
    const int volume_index,
    const int num_volumes,
    std::map<PFace, std::pair<int, int> >& map_volumes,
    std::map<int, Point_3>& centroids) const {

    CGAL_assertion(num_volumes  == 0);
    CGAL_assertion(volume_index >= 0);
    if (pface.first >= 6) return std::make_pair(false, 0);
    CGAL_assertion(pface.first < 6);
    const auto& pair = map_volumes.at(pface);
    CGAL_assertion(pair.second == -1);
    if (pair.first != -1) return std::make_pair(false, 0);
    CGAL_assertion(pair.first == -1);

    std::deque<PFace> queue;
    queue.push_front(pface);

    Point_3 volume_centroid;
    std::size_t volume_size = 0;

    while (!queue.empty()) {
      // print_queue(volume_index, queue);
      const auto query = queue.front();
      queue.pop_front();
      propagate_pface(
        false, query, volume_index, num_volumes, centroids,
        volume_size, volume_centroid, map_volumes, queue);
    }
    std::cout << "- FOUND VOLUME " << volume_index << ", (SIZE/BARYCENTER): "
    << volume_size << " / " << volume_centroid << std::endl;
    centroids[volume_index] = volume_centroid;
    return std::make_pair(true, volume_size);
  }

  const std::pair<bool, std::size_t> traverse_interior_volume(
    const PFace& pface,
    const int volume_index,
    const int num_volumes,
    std::map<PFace, std::pair<int, int> >& map_volumes,
    std::map<int, Point_3>& centroids) const {

    CGAL_assertion(volume_index > 0);
    CGAL_assertion(volume_index >= num_volumes);

    if (pface.first < 6) return std::make_pair(false, 0);
    CGAL_assertion(pface.first >= 6);
    const auto& pair = map_volumes.at(pface);
    if (pair.second != -1) {
      CGAL_assertion(pair.first != -1);
      return std::make_pair(false, 0);
    }
    CGAL_assertion(pair.second == -1);
    if (pair.first == -1) {
      CGAL_assertion(pair.second == -1);
      return std::make_pair(false, 0);
    }
    CGAL_assertion(pair.first != -1);
    if (pair.first >= num_volumes) return std::make_pair(false, 0);
    CGAL_assertion(pair.first < num_volumes);

    std::deque<PFace> queue;
    queue.push_front(pface);

    Point_3 volume_centroid;
    std::size_t volume_size = 0;

    while (!queue.empty()) {
      // print_queue(volume_index, queue);
      const auto query = queue.front();
      queue.pop_front();
      propagate_pface(
        false, query, volume_index, num_volumes, centroids,
        volume_size, volume_centroid, map_volumes, queue);
    }
    std::cout << "- FOUND VOLUME " << volume_index << ", (SIZE/BARYCENTER): "
    << volume_size << " / " << volume_centroid << std::endl;
    centroids[volume_index] = volume_centroid;
    return std::make_pair(true, volume_size);
  }

  void print_queue(
    const int volume_index,
    const std::deque<PFace>& queue) const {

    // if (volume_index != -1) return;
    std::cout << "QUEUE: " << std::endl;
    for (const auto& pface : queue) {
      std::cout << volume_index << " "
      << pface.first << " " << pface.second << std::endl;
    }
  }

  void propagate_pface(
    const bool verbose,
    const PFace& pface,
    const int volume_index,
    const int num_volumes,
    const std::map<int, Point_3>& centroids,
    std::size_t& volume_size,
    Point_3& volume_centroid,
    std::map<PFace, std::pair<int, int> >& map_volumes,
    std::deque<PFace>& queue) const {

    const bool is_boundary = is_boundary_pface(
      pface, volume_index, num_volumes, map_volumes);
    if (is_boundary) {
      propagate_boundary_pface(
        verbose, pface, volume_index, num_volumes, centroids,
        volume_size, volume_centroid, map_volumes, queue);
    } else {
      propagate_interior_pface(
        verbose, pface, volume_index, num_volumes, centroids,
        volume_size, volume_centroid, map_volumes, queue);
    }
  }

  const bool is_boundary_pface(
    const PFace& pface,
    const int volume_index,
    const int num_volumes,
    const std::map<PFace, std::pair<int, int> >& map_volumes) const {

    CGAL_assertion(volume_index >= 0);
    if (pface.first < 6) return true;
    CGAL_assertion(pface.first >= 6);
    if (num_volumes == 0) return false;
    CGAL_assertion(num_volumes  > 0);
    CGAL_assertion(volume_index > 0);
    CGAL_assertion(volume_index >= num_volumes);

    const auto& pair = map_volumes.at(pface);
    if (pair.first == -1) {
      CGAL_assertion(pair.second == -1);
      return false;
    }
    CGAL_assertion(pair.first != -1);
    if (pair.first < num_volumes) return true;
    CGAL_assertion(pair.first >= num_volumes);
    return false;
  }

  void propagate_boundary_pface(
    const bool verbose,
    const PFace& pface,
    const int volume_index,
    const int num_volumes,
    const std::map<int, Point_3>& centroids,
    std::size_t& volume_size,
    Point_3& volume_centroid,
    std::map<PFace, std::pair<int, int> >& map_volumes,
    std::deque<PFace>& queue) const {

    auto& pair = map_volumes.at(pface);
    if (pair.first >= num_volumes) return;
    CGAL_assertion(pair.first < num_volumes);
    if (pair.second != -1) {
      CGAL_assertion(pair.first != -1);
      return;
    }
    CGAL_assertion(pair.second == -1);

    if (pair.first == -1) {
      pair.first  = volume_index;
    } else {
      pair.second = volume_index;
    }

    const Point_3 center = centroid_of_pface(pface);
    volume_centroid = CGAL::barycenter(
      volume_centroid, static_cast<FT>(volume_size), center, FT(1));

    Point_3 saved_centroid;
    bool is_saved_centroid = false;
    if (num_volumes > 0 && volume_size == 0) {

      is_saved_centroid = true;
      saved_centroid = volume_centroid;

      CGAL_assertion(pair.first < num_volumes);
      CGAL_assertion(centroids.find(pair.first) != centroids.end());
      const auto& other_centroid = centroids.at(pair.first);
      const auto plane = plane_of_pface(pface);
      auto vec1 = plane.orthogonal_vector();
      vec1 = KSR::normalize(vec1);
      auto vec2 = Vector_3(center, other_centroid);
      vec2 = KSR::normalize(vec2);

      const FT d = FT(1) / FT(100000); // TODO: CAN WE AVOID THIS VALUE?
      const FT dot_product = vec1 * vec2;

      if (dot_product < FT(0)) {
        volume_centroid += d * vec1;
      } else {
        volume_centroid -= d * vec1;
      }
    }

    // std::cout << "volume centroid: " << volume_centroid << std::endl;
    ++volume_size;

    if (verbose) {
      // std::cout << "BND PFACE MAP: (" <<
      // pair.first << ", " << pair.second << ")" << std::endl;
      std::cout << "DUMPING BND PFACE: " <<
        std::to_string(volume_index) + "-" +
        std::to_string(pface.first) + "-" +
        std::to_string(pface.second) << std::endl;
      dump_pface(pface, "bnd-pface-" +
        std::to_string(volume_index) + "-" +
        std::to_string(pface.first) + "-" +
        std::to_string(pface.second));
    }

    std::vector<PFace> nfaces, bnd_nfaces, int_nfaces, all_nfaces;
    const auto pedges = pedges_of_pface(pface);
    for (const auto pedge : pedges) {
      CGAL_assertion(has_iedge(pedge));
      incident_faces(this->iedge(pedge), nfaces);
      split_pfaces(
        pface, volume_index, num_volumes, map_volumes, nfaces,
        bnd_nfaces, int_nfaces, all_nfaces);

      if (num_volumes == 0) {
        CGAL_assertion(bnd_nfaces.size() == 1);
        CGAL_assertion(int_nfaces.size() == 0 || int_nfaces.size() == 1);
      }

      if (int_nfaces.size() == 1) {
        queue.push_back(int_nfaces[0]);
        continue;
      }

      if (int_nfaces.size() == 0 && bnd_nfaces.size() == 1) {
        queue.push_front(bnd_nfaces[0]);
        continue;
      }

      if (all_nfaces.size() == 0) {
        dump_info(pface, pedge, nfaces);
        std::cout << "DEBUG: num nfaces: " << nfaces.size() << std::endl;
      }
      CGAL_assertion(all_nfaces.size() > 0);

      const auto found_nface = find_using_2d_directions(
      volume_index, volume_centroid, pface, pedge, all_nfaces);
      if (found_nface == null_pface()) continue;

      if (is_boundary_pface(
        found_nface, volume_index, num_volumes, map_volumes)) {
        queue.push_front(found_nface);
      } else {
        queue.push_back(found_nface);
      }
    }

    if (is_saved_centroid) {
      volume_centroid = saved_centroid;
      // CGAL_assertion_msg(false, "TODO: test this volume centroid!");
    }
  }

  void propagate_interior_pface(
    const bool verbose,
    const PFace& pface,
    const int volume_index,
    const int num_volumes,
    const std::map<int, Point_3>& centroids,
    std::size_t& volume_size,
    Point_3& volume_centroid,
    std::map<PFace, std::pair<int, int> >& map_volumes,
    std::deque<PFace>& queue) const {

    CGAL_assertion(num_volumes >= 0);
    auto& pair = map_volumes.at(pface);
    if (pair.first != -1 && pair.second != -1) return;
    CGAL_assertion(pair.second == -1);
    if (pair.first == volume_index) return;
    CGAL_assertion(pair.first != volume_index);
    if (pair.first != -1) {
      pair.second = volume_index;
    } else {
      pair.first  = volume_index;
    }

    const Point_3 center = centroid_of_pface(pface);
    volume_centroid = CGAL::barycenter(
      volume_centroid, static_cast<FT>(volume_size), center, FT(1));
    // std::cout << "volume centroid: " << volume_centroid << std::endl;
    ++volume_size;

    if (verbose) {
      // std::cout << "INT PFACE MAP: (" <<
      // pair.first << ", " << pair.second << ")" << std::endl;
      std::cout << "DUMPING INT PFACE: " <<
        std::to_string(volume_index) + "-" +
        std::to_string(pface.first) + "-" +
        std::to_string(pface.second) << std::endl;
      dump_pface(pface, "int-pface-" +
        std::to_string(volume_index) + "-" +
        std::to_string(pface.first) + "-" +
        std::to_string(pface.second));
    }

    std::vector<PFace> nfaces, bnd_nfaces, int_nfaces, all_nfaces;
    const auto pedges = pedges_of_pface(pface);
    for (const auto pedge : pedges) {
      CGAL_assertion(has_iedge(pedge));
      incident_faces(this->iedge(pedge), nfaces);
      split_pfaces(
        pface, volume_index, num_volumes, map_volumes, nfaces,
        bnd_nfaces, int_nfaces, all_nfaces);

      if (all_nfaces.size() == 0) {
        dump_info(pface, pedge, nfaces);
        std::cout << "DEBUG: num nfaces: " << nfaces.size() << std::endl;
      }
      CGAL_assertion(all_nfaces.size() > 0);

      const auto found_nface = find_using_2d_directions(
      volume_index, volume_centroid, pface, pedge, all_nfaces);
      if (found_nface == null_pface()) continue;

      if (is_boundary_pface(
        found_nface, volume_index, num_volumes, map_volumes)) {
        queue.push_front(found_nface);
      } else {
        queue.push_back(found_nface);
      }
    }
  }

  void split_pfaces(
    const PFace& current,
    const int volume_index,
    const int num_volumes,
    const std::map<PFace, std::pair<int, int> >& map_volumes,
    const std::vector<PFace>& pfaces,
    std::vector<PFace>& bnd_pfaces,
    std::vector<PFace>& int_pfaces,
    std::vector<PFace>& all_pfaces) const {

    bnd_pfaces.clear();
    int_pfaces.clear();
    all_pfaces.clear();
    for (const auto& pface : pfaces) {
      if (pface == current) continue;
      CGAL_assertion(pface != current);
      all_pfaces.push_back(pface);

      const auto& pair = map_volumes.at(pface);
      if (num_volumes > 0 && pair.first != -1) {
        if (pair.first < num_volumes && pair.second != -1) {
          if (pair.second < num_volumes) {
            continue;
          }
          CGAL_assertion(pair.second >= num_volumes);
        }
      }
      if (is_boundary_pface(
        pface, volume_index, num_volumes, map_volumes))  {
        bnd_pfaces.push_back(pface);
      } else {
        int_pfaces.push_back(pface);
      }
    }
  }

  const PFace find_using_2d_directions(
    const int volume_index,
    const Point_3& volume_centroid,
    const PFace& pface,
    const PEdge& pedge,
    const std::vector<PFace>& nfaces) const {

    CGAL_assertion(nfaces.size() > 0);
    if (nfaces.size() == 1) return nfaces[0];
    const bool is_debug = false;
      // ( volume_index == 31 &&
      //   pface.first == 8 &&
      //   static_cast<std::size_t>(pface.second) == 7);

    if (is_debug) {
      dump_info(pface, pedge, nfaces);
    }
    CGAL_assertion(nfaces.size() > 1);

    Point_3 center = centroid_of_pface(pface);
    const Segment_3 segment = segment_3(pedge);
    const Line_3 line(segment.source(), segment.target());
    Point_3 midp = CGAL::midpoint(segment.source(), segment.target());
    Vector_3 norm(segment.source(), segment.target());
    norm = KSR::normalize(norm);
    const Plane_3 plane(midp, norm);

    std::vector<Point_3> points;
    points.reserve(nfaces.size() + 2);

    points.push_back(midp);
    points.push_back(center);
    for (const auto& nface : nfaces) {
      center = centroid_of_pface(nface);
      points.push_back(center);
    }
    CGAL_assertion(points.size() >= 3);

    for (auto& point : points) {
      point = plane.projection(point);
    }

    if (is_debug) {
      dump_frame(points, "polyhedrons/directions-init");
    }

    const FT cx = volume_centroid.x();
    const FT cy = volume_centroid.y();
    const FT cz = volume_centroid.z();
    FT d = (norm.x() * cx + norm.y() * cy + norm.z() * cz + plane.d());

    const Plane_3 tr_plane(midp + norm * d, norm);
    Point_3 inter;
    const bool is_intersection_found = KSR::intersection(line, tr_plane, inter);
    if (!is_intersection_found) {
      std::cout << "d = " << d << std::endl;
    }
    CGAL_assertion(is_intersection_found);
    d = KSR::distance(midp, inter);
    norm = Vector_3(midp, inter);
    norm = KSR::normalize(norm);
    for (auto& point : points) {
      point += norm * d;
    }

    if (is_debug) {
      auto extended = points;
      extended.push_back(volume_centroid);
      dump_frame(extended, "polyhedrons/directions");
    }

    std::vector< std::pair<Direction_2, PFace> > dir_edges;
    dir_edges.reserve(nfaces.size() + 1);

    const Point_2 proj_0 = plane.to_2d(points[0]);
    for (std::size_t i = 1; i < points.size(); ++i) {
      const Point_2 proj_i = plane.to_2d(points[i]);
      const Vector_2 vec(proj_0, proj_i);
      if (i == 1) {
        dir_edges.push_back(std::make_pair(Direction_2(vec), pface));
      } else {
        dir_edges.push_back(std::make_pair(Direction_2(vec), nfaces[i - 2]));
      }
    }
    CGAL_assertion(dir_edges.size() == nfaces.size() + 1);

    const Point_2 proj_vc = plane.to_2d(volume_centroid);
    const Vector_2 vec(proj_0, proj_vc);
    const Direction_2 ref_dir(vec);

    std::sort(dir_edges.begin(), dir_edges.end(), [&](
      const std::pair<Direction_2, PFace>& p,
      const std::pair<Direction_2, PFace>& q) -> bool {
        return p.first < q.first;
      }
    );

    const std::size_t n = dir_edges.size();
    for (std::size_t i = 0; i < n; ++i) {
      if (dir_edges[i].second == pface) {

        const std::size_t im = (i + n - 1) % n;
        const std::size_t ip = (i + 1) % n;

        const auto& dir_prev = dir_edges[im].first;
        const auto& dir_curr = dir_edges[i].first;
        const auto& dir_next = dir_edges[ip].first;

        if (is_debug) {
          dump_pface(dir_edges[im].second, "prev");
          dump_pface(dir_edges[ip].second, "next");
        }

        if (ref_dir.counterclockwise_in_between(dir_prev, dir_curr)) {
          if (is_debug) {
            std::cout << "found prev" << std::endl;
            exit(EXIT_SUCCESS);
          }
          return dir_edges[im].second;
        } else if (ref_dir.counterclockwise_in_between(dir_curr, dir_next)) {
          if (is_debug) {
            std::cout << "found next" << std::endl;
            exit(EXIT_SUCCESS);
          }
          return dir_edges[ip].second;
        } else {
          // return null_pface();
          dump_info(pface, pedge, nfaces);
          dump_frame(points, "polyhedrons/directions-init");
          auto extended = points;
          extended.push_back(volume_centroid);
          dump_frame(extended, "polyhedrons/directions");
          CGAL_assertion_msg(false, "ERROR: WRONG ORIENTATION!");
        }
      }
    }

    CGAL_assertion_msg(false, "ERROR: NEXT PFACE IS NOT FOUND!");
    return null_pface();
  }

  const Point_3 centroid_of_pface(const PFace& pface) const {

    const std::function<Point_3(PVertex)> unary_f =
    [&](const PVertex& pvertex) -> Point_3 {
      return point_3(pvertex);
    };
    const std::vector<Point_3> polygon(
      boost::make_transform_iterator(pvertices_of_pface(pface).begin(), unary_f),
      boost::make_transform_iterator(pvertices_of_pface(pface).end()  , unary_f));
    CGAL_assertion(polygon.size() >= 3);
    return CGAL::centroid(polygon.begin(), polygon.end());
  }

  const Plane_3 plane_of_pface(const PFace& pface) const {

    const std::function<Point_3(PVertex)> unary_f =
    [&](const PVertex& pvertex) -> Point_3 {
      return point_3(pvertex);
    };
    const std::vector<Point_3> polygon(
      boost::make_transform_iterator(pvertices_of_pface(pface).begin(), unary_f),
      boost::make_transform_iterator(pvertices_of_pface(pface).end()  , unary_f));
    CGAL_assertion(polygon.size() >= 3);
    return Plane_3(polygon[0], polygon[1], polygon[2]);
  }

  void create_cell_pvertices(Volume_cell& cell) {
    cell.pvertices.clear();
    for (const auto& pface : cell.pfaces) {
      for (const auto pvertex : pvertices_of_pface(pface)) {
        cell.pvertices.insert(pvertex);
      }
    }
  }

  void dump_pface(
    const PFace& pface,
    const std::string name) const {

    std::vector<Point_3> polygon;
    std::vector< std::vector<Point_3> > polygons;
    for (const auto pvertex : pvertices_of_pface(pface)) {
      polygon.push_back(point_3(pvertex));
    }
    polygons.push_back(polygon);
    KSR_3::Saver<Kernel> saver;
    saver.export_polygon_soup_3(polygons, "polyhedrons/" + name);
  }

  void dump_pedge(
    const PEdge& pedge,
    const std::string name) const {

    const std::vector<Segment_3> segments = { segment_3(pedge) };
    KSR_3::Saver<Kernel> saver;
    saver.export_segments_3(segments, "polyhedrons/" + name);
  }

  void dump_info(
    const PFace& pface,
    const PEdge& pedge,
    const std::vector<PFace>& nfaces) const {

    std::cout << "DEBUG: number of found nfaces: " << nfaces.size() << std::endl;
    dump_pface(pface, "face-curr");
    dump_pedge(pedge, "face-edge");
    for (std::size_t i = 0; i < nfaces.size(); ++i) {
      dump_pface(nfaces[i], "nface-" + std::to_string(i));
    }
  }

  void dump_frame(
    const std::vector<Point_3>& points,
    const std::string name) const {

    std::vector<Segment_3> segments;
    segments.reserve(points.size() - 1);
    for (std::size_t i = 1; i < points.size(); ++i)
      segments.push_back(Segment_3(points[0], points[i]));
    KSR_3::Saver<Kernel> saver;
    saver.export_segments_3(segments, name);
  }

  const std::vector<Volume_cell> incident_volumes(const PFace& query_pface) const {
    std::vector<Volume_cell> nvolumes;
    for (const auto& volume : m_volumes) {
      for (const auto& pface : volume.pfaces) {
        if (pface == query_pface) nvolumes.push_back(volume);
      }
    }
    return nvolumes;
  }

  void incident_faces(const IEdge& query_iedge, std::vector<PFace>& nfaces) const {

    nfaces.clear();
    for (const auto plane_idx : intersected_planes(query_iedge)) {
      for (const auto pedge : pedges(plane_idx)) {
        if (iedge(pedge) == query_iedge) {
          const auto& m = mesh(plane_idx);
          const auto he = m.halfedge(pedge.second);
          const auto op = m.opposite(he);
          const auto face1 = m.face(he);
          const auto face2 = m.face(op);
          if (face1 != Support_plane::Mesh::null_face()) {
            nfaces.push_back(PFace(plane_idx, face1));
          }
          if (face2 != Support_plane::Mesh::null_face()) {
            nfaces.push_back(PFace(plane_idx, face2));
          }
        }
      }
    }
  }

  void update_positions(const FT time) {
    m_previous_time = m_current_time;
    m_current_time = time;
  }

  const FT previous_time() const {
    return m_previous_time;
  }

  // Strings.
  inline const std::string str(const PVertex& pvertex) const {
    return "PVertex(" + std::to_string(pvertex.first) + ":v" + std::to_string(pvertex.second) + ")";
  }
  inline const std::string str(const PEdge& pedge) const {
    return "PEdge(" + std::to_string(pedge.first) + ":e" + std::to_string(pedge.second) + ")";
  }
  inline const std::string str(const PFace& pface) const {
    return "PFace(" + std::to_string(pface.first) + ":f" + std::to_string(pface.second) + ")";
  }
  inline const std::string str(const IVertex& ivertex) const {
    return "IVertex(" + std::to_string(ivertex) + ")";
  }
  inline const std::string str(const IEdge& iedge) const {
    std::ostringstream oss; oss << "IEdge" << iedge; return oss.str();
  }

  inline const std::string lstr(const PFace& pface) const {

    if (pface == null_pface()) {
      return "PFace(null)";
    }
    std::string out = "PFace(" + std::to_string(pface.first) + ":f" + std::to_string(pface.second) + ")[";
    for (const auto pvertex : pvertices_of_pface(pface)) {
      out += "v" + std::to_string(pvertex.second);
    }
    out += "]";
    return out;
  }

  inline const std::string lstr(const PEdge& pedge) const {
    return "PEdge(" + std::to_string(pedge.first) + ":e" + std::to_string(pedge.second)
      + ")[v" + std::to_string(source(pedge).second) + "->v" + std::to_string(target(pedge).second) + "]";
  }

  template<typename PSimplex>
  const Support_plane& support_plane(const PSimplex& psimplex) const { return support_plane(psimplex.first); }
  const Support_plane& support_plane(const KSR::size_t idx) const { return m_support_planes[idx]; }

  template<typename PSimplex>
  Support_plane& support_plane(const PSimplex& psimplex) { return support_plane(psimplex.first); }
  Support_plane& support_plane(const KSR::size_t idx) { return m_support_planes[idx]; }

private:
  std::map< std::pair<KSR::size_t, IEdge>, Point_2>  m_points;
  std::map< std::pair<KSR::size_t, IEdge>, Vector_2> m_directions;

  template<typename PSimplex>
  const Mesh& mesh(const PSimplex& psimplex) const { return mesh(psimplex.first); }
  const Mesh& mesh(const KSR::size_t support_plane_idx) const { return support_plane(support_plane_idx).mesh(); }

  template<typename PSimplex>
  Mesh& mesh(const PSimplex& psimplex) { return mesh(psimplex.first); }
  Mesh& mesh(const KSR::size_t support_plane_idx) { return support_plane(support_plane_idx).mesh(); }

  void compute_future_points_and_directions(
    const PVertex& pvertex, const IEdge& iedge,
    Point_2& future_point_a, Point_2& future_point_b,
    Vector_2& direction_a, Vector_2& direction_b) const {

    const PVertex prev(pvertex.first, support_plane(pvertex).prev(pvertex.second));
    const PVertex next(pvertex.first, support_plane(pvertex).next(pvertex.second));
    const auto& curr = pvertex;

    const Line_2 iedge_line = segment_2(pvertex.first, iedge).supporting_line();
    const Point_2 pinit = iedge_line.projection(point_2(pvertex));

    // std::cout << "iedge segment: " << segment_3(iedge) << std::endl;

    const auto prev_p = point_2(prev);
    const auto next_p = point_2(next);
    const auto curr_p = point_2(curr);

    // std::cout << "prev: " << point_3(prev) << std::endl;
    // std::cout << "next: " << point_3(next) << std::endl;
    // std::cout << "curr: " << point_3(curr) << std::endl;

    const Line_2 future_line_prev(
      point_2(prev, m_current_time + FT(1)),
      point_2(curr, m_current_time + FT(1)));
    const Line_2 future_line_next(
      point_2(next, m_current_time + FT(1)),
      point_2(curr, m_current_time + FT(1)));

    // std::cout << "future line prev: " <<
    // Segment_3(
    //   to_3d(pvertex.first, point_2(prev, m_current_time + FT(1))),
    //   to_3d(pvertex.first, point_2(curr, m_current_time + FT(1)))) << std::endl;
    // std::cout << "future line next: " <<
    // Segment_3(
    //   to_3d(pvertex.first, point_2(next, m_current_time + FT(1))),
    //   to_3d(pvertex.first, point_2(curr, m_current_time + FT(1)))) << std::endl;

    const Vector_2 current_vec_prev(prev_p, curr_p);
    const Vector_2 current_vec_next(next_p, curr_p);

    const auto source_p = point_2(pvertex.first, source(iedge));
    const auto target_p = point_2(pvertex.first, target(iedge));
    const Vector_2 iedge_vec(source_p, target_p);

    const FT tol = FT(1) / FT(100000);
    FT m1 = FT(100000), m2 = FT(100000), m3 = FT(100000);

    const FT prev_d = (curr_p.x() - prev_p.x());
    const FT next_d = (curr_p.x() - next_p.x());
    const FT edge_d = (target_p.x() - source_p.x());

    if (CGAL::abs(prev_d) > tol)
      m1 = (curr_p.y() - prev_p.y()) / prev_d;
    if (CGAL::abs(next_d) > tol)
      m2 = (curr_p.y() - next_p.y()) / next_d;
    if (CGAL::abs(edge_d) > tol)
      m3 = (target_p.y() - source_p.y()) / edge_d;

    // std::cout << "prev slope: "  << m1 << std::endl;
    // std::cout << "next slope: "  << m2 << std::endl;
    // std::cout << "iedge slope: " << m3 << std::endl;

    if (CGAL::abs(m1 - m3) < tol) {
      std::cout << "prev parallel lines" << std::endl;
      const FT prev_dot = current_vec_prev * iedge_vec;
      if (prev_dot < FT(0)) {
        std::cout << "prev moves backwards" << std::endl;
        future_point_a = target_p;
      } else {
        std::cout << "prev moves forwards" << std::endl;
        future_point_a = source_p;
      }
    } else {

      std::cout << "prev intersected lines" << std::endl;
      const bool a_found = KSR::intersection(future_line_prev, iedge_line, future_point_a);
      if (!a_found)
      {
        std::cerr << "Warning: a not found" << std::endl;
        future_point_b = pinit + (pinit - future_point_a);
      }
    }

    direction_a = Vector_2(pinit, future_point_a);
    future_point_a = pinit - m_current_time * direction_a;
    std::cout << "future point a: " << to_3d(pvertex.first, future_point_a + m_current_time * direction_a) << std::endl;
    std::cout << "dir a: " << direction_a << std::endl;

    if (CGAL::abs(m2 - m3) < tol) {
      std::cout << "next parallel lines" << std::endl;
      const FT next_dot = current_vec_next * iedge_vec;
      if (next_dot < FT(0)) {
        std::cout << "next moves backwards" << std::endl;
        future_point_b = target_p;
      } else {
        std::cout << "next moves forwards" << std::endl;
        future_point_b = source_p;
      }

    } else {

      std::cout << "next intersected lines" << std::endl;
      const bool b_found = KSR::intersection(future_line_next, iedge_line, future_point_b);
      if (!b_found)
      {
        std::cerr << "Warning: b not found" << std::endl;
        future_point_a = pinit + (pinit - future_point_b);
      }
    }

    direction_b = Vector_2(pinit, future_point_b);
    future_point_b = pinit - m_current_time * direction_b;
    std::cout << "future point b: " << to_3d(pvertex.first, future_point_b + m_current_time * direction_b) << std::endl;
    std::cout << "dir b: " << direction_b << std::endl;
  }

  const bool compute_future_point_and_direction(
    const std::size_t idx,
    const PVertex& pvertex, const PVertex& next, // back prev
    const IEdge& iedge,
    Point_2& future_point, Vector_2& future_direction) const {

    bool is_parallel = false;
    // if (this->iedge(pvertex) != null_iedge()
    //     && line_idx(pvertex) == line_idx(iedge))
    // {
    //   std::cerr << "found limit" << std::endl;
    //   future_point = point_2(pvertex, FT(0));
    //   future_direction = this->direction(pvertex);
    //   return is_parallel;
    // }

    const Line_2 iedge_line = segment_2(pvertex.first, iedge).supporting_line();
    const Point_2 pinit = iedge_line.projection(point_2(pvertex));

    const auto& curr = pvertex;
    const auto next_p = point_2(next);
    const auto curr_p = point_2(curr);

    const Line_2 future_line_next(
      point_2(next, m_current_time + FT(1)),
      point_2(curr, m_current_time + FT(1)));
    const Vector_2 current_vec_next(next_p, curr_p);

    const auto source_p = point_2(pvertex.first, source(iedge));
    const auto target_p = point_2(pvertex.first, target(iedge));
    const Vector_2 iedge_vec(source_p, target_p);

    const FT tol = FT(1) / FT(100000);
    FT m2 = FT(100000), m3 = FT(100000);

    const FT next_d = (curr_p.x() - next_p.x());
    const FT edge_d = (target_p.x() - source_p.x());

    if (CGAL::abs(next_d) > tol)
      m2 = (curr_p.y() - next_p.y()) / next_d;
    if (CGAL::abs(edge_d) > tol)
      m3 = (target_p.y() - source_p.y()) / edge_d;

    // std::cout << "m2: " << m2 << std::endl;
    // std::cout << "m3: " << m3 << std::endl;

    if (CGAL::abs(m2 - m3) < tol) {
      std::cout << "back/front parallel lines" << std::endl;

      is_parallel = true;
      const FT next_dot = current_vec_next * iedge_vec;
      if (next_dot < FT(0)) {
        std::cout << "back/front moves backwards" << std::endl;
        future_point = target_p;
        // std::cout << point_3(target(iedge)) << std::endl;
      } else {
        std::cout << "back/front moves forwards" << std::endl;
        future_point = source_p;
        // std::cout << point_3(source(iedge)) << std::endl;
      }

    } else {
      std::cout << "back/front intersected lines" << std::endl;
      future_point = KSR::intersection<Point_2>(future_line_next, iedge_line);
    }

    future_direction = Vector_2(pinit, future_point);

    // std::cout << "prev: " << point_3(next, m_current_time + FT(1)) << std::endl;
    // std::cout << "back: " << point_3(curr, m_current_time + FT(1)) << std::endl;

    // auto tmp = future_direction;
    // tmp = KSR::normalize(tmp);
    // std::cout << "future: " << to_3d(pvertex.first, pinit + m_current_time * tmp) << std::endl;

    future_point = pinit - m_current_time * future_direction;
    return is_parallel;
  }

  const bool compute_future_point_and_direction(
    const PVertex& pvertex,
    const PVertex& prev, const PVertex& next,
    const IEdge& iedge,
    Point_2& future_point, Vector_2& future_direction) const {

    const Line_2 iedge_line = segment_2(pvertex.first, iedge).supporting_line();
    const auto pv_point = point_2(pvertex);
    const Point_2 pinit = iedge_line.projection(pv_point);

    const auto& curr = prev;
    const auto next_p = point_2(next);
    const auto curr_p = point_2(curr);

    const Line_2 future_line_next(
      point_2(next, m_current_time + FT(1)),
      point_2(curr, m_current_time + FT(1)));

    const auto source_p = point_2(pvertex.first, source(iedge));
    const auto target_p = point_2(pvertex.first, target(iedge));
    const Vector_2 iedge_vec(source_p, target_p);

    const FT tol = FT(1) / FT(100000);
    FT m2 = FT(100000), m3 = FT(100000);

    const FT next_d = (curr_p.x() - next_p.x());
    const FT edge_d = (target_p.x() - source_p.x());

    if (CGAL::abs(next_d) > tol)
      m2 = (curr_p.y() - next_p.y()) / next_d;
    if (CGAL::abs(edge_d) > tol)
      m3 = (target_p.y() - source_p.y()) / edge_d;

    // std::cout << "m2: " << m2 << std::endl;
    // std::cout << "m3: " << m3 << std::endl;
    // std::cout << "mm: " << m2 - m3 << std::endl;

    bool is_parallel = false;
    if (CGAL::abs(m2 - m3) < tol) {
      std::cout << "open parallel lines" << std::endl;

      is_parallel = true;
      if (source_p == pv_point)
        future_point = target_p;
      else
        future_point = source_p;

    } else {
      std::cout << "open intersected lines" << std::endl;
      future_point = KSR::intersection<Point_2>(future_line_next, iedge_line);
    }

    future_direction = Vector_2(pinit, future_point);
    future_point = pinit - m_current_time * future_direction;

    // auto tmp = future_direction;
    // tmp = KSR::normalize(tmp);
    // std::cout << "future: " << to_3d(pvertex.first, pinit + m_current_time * tmp) << std::endl;

    return is_parallel;
  }

  const bool is_intersecting_iedge(
    const FT min_time, const FT max_time,
    const PVertex& pvertex, const IEdge& iedge) {

    const FT time_step = (max_time - min_time) / FT(100);
    const FT time_1 = m_current_time - time_step;
    const FT time_2 = m_current_time + time_step;
    CGAL_assertion(time_1 != time_2);

    const Segment_2 pv_seg(
      point_2(pvertex, time_1), point_2(pvertex, time_2));
    const auto pv_bbox = pv_seg.bbox();

    const auto iedge_seg  = segment_2(pvertex.first, iedge);
    const auto iedge_bbox = iedge_seg.bbox();

    if (has_iedge(pvertex)) {
      std::cout << "constrained pvertex case" << std::endl;
      return false;
    }

    if (!is_active(pvertex)) {
      std::cout << "pvertex no active case" << std::endl;
      return false;
    }

    if (!is_active(iedge)) {
      std::cout << "iedge no active case" << std::endl;
      return false;
    }

    if (!CGAL::do_overlap(pv_bbox, iedge_bbox)) {
      std::cout << "no overlap case" << std::endl;
      return false;
    }

    Point_2 point;
    if (!KSR::intersection(pv_seg, iedge_seg, point)) {
      std::cout << "no intersection case" << std::endl;
      return false;
    }

    std::cout << "found intersection" << std::endl;
    return true;
  }
};

}} // namespace CGAL::KSR_3

#endif // CGAL_KSR_3_DATA_STRUCTURE_H
