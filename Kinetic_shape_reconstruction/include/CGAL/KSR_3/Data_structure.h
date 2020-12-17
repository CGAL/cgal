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

  struct Halfedge_to_pface {
    using argument_type = Halfedge_index;
    using result_type = PFace;

    const KSR::size_t support_plane_idx;
    const Mesh& mesh;

    Halfedge_to_pface(const KSR::size_t sp_idx, const Mesh& m) :
    support_plane_idx(sp_idx),
    mesh(m)
    { }

    const result_type operator()(const argument_type& arg) const {
      return result_type(support_plane_idx, mesh.face(arg));
    }
  };

  using PEdge_around_pvertex_iterator =
    boost::transform_iterator<Halfedge_to_pedge, CGAL::Halfedge_around_target_iterator<Mesh> >;
  using PEdges_around_pvertex = CGAL::Iterator_range<PEdge_around_pvertex_iterator>;

  using PEdge_of_pface_iterator =
    boost::transform_iterator<Halfedge_to_pedge, CGAL::Halfedge_around_face_iterator<Mesh> >;
  using PEdges_of_pface = CGAL::Iterator_range<PEdge_of_pface_iterator>;

  using PFace_around_pvertex_iterator =
    boost::transform_iterator<Halfedge_to_pface, CGAL::Halfedge_around_target_iterator<Mesh> >;
  using PFaces_around_pvertex = CGAL::Iterator_range<PFace_around_pvertex_iterator>;

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
  std::map< std::pair<KSR::size_t, IEdge>, Point_2>  m_points;
  std::map< std::pair<KSR::size_t, IEdge>, Vector_2> m_directions;
  KSR::vector<Support_plane> m_support_planes;
  Intersection_graph m_intersection_graph;

  FT m_current_time;
  FT m_previous_time;
  bool m_verbose;

  std::vector<Volume_cell> m_volumes;
  std::size_t m_num_volume_levels;
  std::map<int, std::size_t> m_volume_level_map;

public:
  Data_structure(const bool verbose) :
  m_current_time(FT(0)),
  m_previous_time(FT(0)),
  m_verbose(verbose),
  m_num_volume_levels(0)
  { }

  void clear() {
    m_points.clear();
    m_directions.clear();
    m_support_planes.clear();
    m_intersection_graph.clear();

    m_current_time  = FT(0);
    m_previous_time = FT(0);

    m_volumes.clear();
    m_num_volume_levels = 0;
    m_volume_level_map.clear();
  }

  const std::size_t number_of_volume_levels() const {
    return m_num_volume_levels;
  }

  const std::size_t number_of_volumes(const int volume_level) const {
    if (volume_level < 0) {
      return m_volumes.size();
    }
    CGAL_assertion(volume_level >= 0);
    return m_volume_level_map.at(volume_level);
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

  /*******************************
  **          GENERAL           **
  ********************************/

  const KSR::vector<Support_plane>& support_planes() const { return m_support_planes; }
  KSR::vector<Support_plane>& support_planes() { return m_support_planes; }

  const Intersection_graph& igraph() const { return m_intersection_graph; }
  Intersection_graph& igraph() { return m_intersection_graph; }

  void resize(const KSR::size_t number_of_items) {
    m_support_planes.resize(number_of_items);
  }

  // TODO: It looks like here we lose precision during the conversion because
  // KSR::size_t is usually smaller than std::size_t!
  void reserve(const std::size_t number_of_polygons) {
    m_support_planes.reserve(static_cast<KSR::size_t>(number_of_polygons) + 6);
  }

  const FT current_time() const { return m_current_time; }
  const FT previous_time() const { return m_previous_time; }

  void update_positions(const FT time) {
    m_previous_time = m_current_time;
    m_current_time = time;
  }

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

  template<typename PSimplex>
  const Support_plane& support_plane(const PSimplex& psimplex) const { return support_plane(psimplex.first); }
  const Support_plane& support_plane(const KSR::size_t idx) const { return m_support_planes[idx]; }

  template<typename PSimplex>
  Support_plane& support_plane(const PSimplex& psimplex) { return support_plane(psimplex.first); }
  Support_plane& support_plane(const KSR::size_t idx) { return m_support_planes[idx]; }

  template<typename PSimplex>
  const Mesh& mesh(const PSimplex& psimplex) const { return mesh(psimplex.first); }
  const Mesh& mesh(const KSR::size_t support_plane_idx) const { return support_plane(support_plane_idx).mesh(); }

  template<typename PSimplex>
  Mesh& mesh(const PSimplex& psimplex) { return mesh(psimplex.first); }
  Mesh& mesh(const KSR::size_t support_plane_idx) { return support_plane(support_plane_idx).mesh(); }

  const KSR::size_t number_of_support_planes() const {
    return m_support_planes.size();
  }

  const bool is_bbox_support_plane(const KSR::size_t support_plane_idx) const {
    return (support_plane_idx < 6);
  }

  template<typename PointRange>
  const KSR::size_t add_support_plane(const PointRange& polygon) {

    const Support_plane new_support_plane(polygon);
    KSR::size_t support_plane_idx = KSR::no_element();
    bool found_coplanar_polygons = false;
    for (KSR::size_t i = 0; i < number_of_support_planes(); ++i) {
      if (new_support_plane == support_plane(i)) {
        found_coplanar_polygons = true;
        support_plane_idx = i;
        return support_plane_idx;
      }
    }
    CGAL_assertion_msg(!found_coplanar_polygons,
    "ERROR: NO COPLANAR POLYGONS HERE!");

    if (support_plane_idx == KSR::no_element()) {
      support_plane_idx = number_of_support_planes();
      m_support_planes.push_back(new_support_plane);
    }

    intersect_with_bbox(support_plane_idx);
    return support_plane_idx;
  }

  void intersect_with_bbox(const KSR::size_t support_plane_idx) {
    if (support_plane_idx < 6) return;

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
      const auto pair = m_intersection_graph.add_edge(ivertices[i], ivertices[(i+1)%4], support_plane_idx);
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
  void add_input_polygon(
    const PointRange& polygon, const KSR::size_t input_index) {

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
    const auto centroid = sort_points_by_direction(points);
    std::vector<KSR::size_t> input_indices;
    input_indices.push_back(input_index);
    support_plane(support_plane_idx).
      add_input_polygon(points, centroid, input_indices);
  }

  const Point_2 sort_points_by_direction(
    std::vector<Point_2>& points) const {

    // Naive version.
    // const auto centroid = CGAL::centroid(points.begin(), points.end());

    // Better version.
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
    return centroid;
  }

  void add_input_polygon(
    const KSR::size_t support_plane_idx,
    const std::vector<KSR::size_t>& input_indices,
    std::vector<Point_2>& points) {

    const auto centroid = sort_points_by_direction(points);
    support_plane(support_plane_idx).
      add_input_polygon(points, centroid, input_indices);
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

  // Get prev and next pvertices of the free pvertex.
  const PVertex prev(const PVertex& pvertex) const {
    return PVertex(pvertex.first, support_plane(pvertex).prev(pvertex.second));
  }
  const PVertex next(const PVertex& pvertex) const {
    return PVertex(pvertex.first, support_plane(pvertex).next(pvertex.second));
  }

  // Get prev and next pvertices of the constrained pvertex.
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

  const PFaces_around_pvertex pfaces_around_pvertex(const PVertex& pvertex) const {

    return PFaces_around_pvertex(
      boost::make_transform_iterator(
        halfedges_around_target(halfedge(pvertex.second, mesh(pvertex)), mesh(pvertex)).begin(),
        Halfedge_to_pface(pvertex.first, mesh(pvertex))),
      boost::make_transform_iterator(
        halfedges_around_target(halfedge(pvertex.second, mesh(pvertex)), mesh(pvertex)).end(),
        Halfedge_to_pface(pvertex.first, mesh(pvertex))));
  }

  void non_null_pfaces_around_pvertex(
    const PVertex& pvertex, std::vector<PFace>& pfaces) const {

    pfaces.clear();
    const auto nfaces = pfaces_around_pvertex(pvertex);
    for (const auto pface : nfaces) {
      if (pface.second == Support_plane::Mesh::null_face()) continue;
      pfaces.push_back(pface);
    }
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

  const std::vector<KSR::size_t>& input(const PFace& pface) const{ return support_plane(pface).input(pface.second); }
  std::vector<KSR::size_t>& input(const PFace& pface) { return support_plane(pface).input(pface.second); }

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
  **          STRINGS           **
  ********************************/

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

  /*******************************
  **        CONNECTIVITY        **
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

    if (m_verbose) {
      std::cout << "** searching pvertices around "
      << str(pvertex) << " wrt " << str(ivertex) << std::endl;
    }

    // std::cout.precision(20);
    std::deque<PVertex> vertices;
    vertices.push_back(pvertex);

    if (m_verbose) {
      std::cout << "- came from: " <<
      str(iedge(pvertex)) << " " << segment_3(iedge(pvertex)) << std::endl;
    }

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
          if (m_verbose) {
            std::cout << "- " << str(current) << " is backwards" << std::endl;
            // std::cout << point_3(current) << std::endl;
          }
          is_free = true;
        }

        if (is_frozen(current)) {
          if (m_verbose) {
            std::cout << "- " << str(current) << " is frozen" << std::endl;
            // std::cout << point_3(current) << std::endl;
          }
          is_free = true;
        }
        // std::cout << "is free 3: " << is_free << std::endl;
      }

      if (previous_was_free && is_free) {
        if (m_verbose) {
          std::cout << "- " << str(current) << " has no iedge, stopping there" << std::endl;
          // std::cout << point_3(current) << std::endl;
        }
        continue;
      }

      if (is_free) {
        if (m_verbose) {
          std::cout << "- " << str(current) << " has no iedge" << std::endl;
          // std::cout << point_3(current) << std::endl;
        }
      }
      else {
        if (m_verbose) {
          std::cout << "- " << str(current) << " has iedge " << str(iedge)
          << " from " << str(source(iedge)) << " to " << str(target(iedge)) << std::endl;
          // std::cout << segment_3(iedge) << std::endl;
          // std::cout << point_3(current) << std::endl;
        }
      }

      if (front) {
        vertices.push_front(current);
        // std::cout << "pushed front" << std::endl;
      }
      else {
        vertices.push_back(current);
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

    if (m_verbose) {
      std::cout << "- found pvertices:";
      for (const auto& pv : out)
        std::cout << " " << str(pv);
      std::cout << std::endl;
    }
    return out;
  }

  /*******************************
  **        CONVERSIONS         **
  ********************************/

  const Point_2 to_2d(const KSR::size_t support_plane_idx, const IVertex& ivertex) const {
    return support_plane(support_plane_idx).to_2d(point_3(ivertex));
  }
  const Segment_2 to_2d(const KSR::size_t support_plane_idx, const Segment_3& segment_3) const {
    return support_plane(support_plane_idx).to_2d(segment_3);
  }
  const Point_2 to_2d(const KSR::size_t support_plane_idx, const Point_3& point_3) const {
    return support_plane(support_plane_idx).to_2d(point_3);
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
  **          PREDICATES        **
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
          //   std::cerr << "WARNING: POINTS ARE ALMOST EQUAL!" << std::endl;
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
            std::cerr << "ERROR: SOURCE_TO_PVERTEX/PEDGE SEGMENT SQ LENGTH = "
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

    // std::cout << "num adjacent faces: " << num_adjacent_faces << std::endl;
    if (num_adjacent_faces <= 1)
      return std::make_pair(false, false);
    return std::make_pair(true, false);
  }

  /*******************************
  **    OPERATIONS ON POLYGONS  **
  ********************************/

  const PVertex crop_pvertex_along_iedge(
    const PVertex& pvertex, const IEdge& iedge) {

    std::cout.precision(20);
    if (m_verbose) {
      std::cout << "** cropping " <<
      str(pvertex) << " along " << str(iedge) << std::endl;
    }

    Point_2 future_point_a, future_point_b;
    Vector_2 future_direction_a, future_direction_b;
    compute_future_points_and_directions(
      pvertex, iedge,
      future_point_a, future_point_b,
      future_direction_a, future_direction_b);

    const PEdge pedge(pvertex.first, support_plane(pvertex).split_vertex(pvertex.second));
    CGAL_assertion(source(pedge) == pvertex || target(pedge) == pvertex);
    const PVertex pother = opposite(pedge, pvertex);

    if (m_verbose) {
      std::cout << "- new pedge: " << str(pedge) << " between "
      << str(pvertex) << " and " << str(pother) << std::endl;
    }

    connect(pedge, iedge);
    connect(pvertex, iedge);
    connect(pother, iedge);

    support_plane(pvertex).set_point(pvertex.second, future_point_a);
    support_plane(pother).set_point(pother.second, future_point_b);
    direction(pvertex) = future_direction_a;
    direction(pother) = future_direction_b;

    // std::cout << "pvertex: "     << point_3(pvertex) << std::endl;
    // std::cout << "pvertex dir: " << future_direction_a << std::endl;
    // std::cout << "pother: "      << point_3(pother) << std::endl;
    // std::cout << "pother dir: "  << future_direction_b << std::endl;

    if (m_verbose) {
      std::cout << "- new pvertices: " << str(pother) << std::endl;
    }
    return pother;
  }

  const std::array<PVertex, 3> propagate_pvertex_beyond_iedge(
    const PVertex& pvertex, const IEdge& iedge,
    const unsigned int k) {

    std::cout.precision(20);
    if (m_verbose) {
      std::cout << "** propagating " <<
      str(pvertex) << " beyond " << str(iedge) << std::endl;
    }

    const Point_2 original_point = point_2(pvertex, FT(0));
    const Vector_2 original_direction = direction(pvertex);
    const PVertex pother = crop_pvertex_along_iedge(pvertex, iedge);

    const PVertex propagated = add_pvertex(pvertex.first, original_point);
    direction(propagated) = original_direction;

    std::array<PVertex, 3> pvertices;
    pvertices[0] = pvertex;
    pvertices[1] = pother;
    pvertices[2] = propagated;

    const PFace new_pface = add_pface(pvertices);
    if (m_verbose) std::cout << "- new pface: " << lstr(new_pface) << std::endl;
    this->k(new_pface) = k;
    CGAL_assertion(new_pface.second != Face_index());
    return pvertices;
  }

  void crop_pedge_along_iedge(
    const PVertex& pvertex, const PVertex& pother, const IEdge& iedge) {

    std::cout.precision(20);
    if (m_verbose) {
      std::cout << "** cropping pedge [" << str(pvertex) << "-" << str(pother)
      << "] along " << str(iedge) << std::endl;
    }

    // std::cout << "pvertex: " << point_3(pvertex) << std::endl;
    // std::cout << "pother: "  << point_3(pother) << std::endl;

    Point_2 future_point;
    Vector_2 future_direction;
    CGAL_assertion(pvertex.first == pother.first);

    // Crop first pvertex.
    {
      const PVertex prev(pvertex.first, support_plane(pvertex).prev(pvertex.second));
      const PVertex next(pvertex.first, support_plane(pvertex).next(pvertex.second));

      if (prev == pother) {
        compute_future_point_and_direction(0, pvertex, next, iedge, future_point, future_direction);
      } else {
        CGAL_assertion(next == pother);
        compute_future_point_and_direction(0, pvertex, prev, iedge, future_point, future_direction);
      }

      direction(pvertex) = future_direction;
      if (m_verbose) {
        std::cout << "- pvertex direction: " << direction(pvertex) << std::endl;
      }
      support_plane(pvertex).set_point(pvertex.second, future_point);
      connect(pvertex, iedge);
    }

    // Crop second pvertex.
    {
      const PVertex prev(pother.first, support_plane(pother).prev(pother.second));
      const PVertex next(pother.first, support_plane(pother).next(pother.second));

      if (prev == pvertex) {
        compute_future_point_and_direction(0, pother, next, iedge, future_point, future_direction);
      } else {
        CGAL_assertion(next == pvertex);
        compute_future_point_and_direction(0, pother, prev, iedge, future_point, future_direction);
      }

      direction(pother) = future_direction;
      if (m_verbose) {
        std::cout << "- pother direction: " << direction(pother) << std::endl;
      }
      support_plane(pother).set_point(pother.second, future_point);
      connect(pother, iedge);
    }

    const PEdge pedge(pvertex.first, support_plane(pvertex).edge(pvertex.second, pother.second));
    connect(pedge, iedge);
  }

  const std::pair<PVertex, PVertex> propagate_pedge_beyond_iedge(
    const PVertex& pvertex, const PVertex& pother, const IEdge& iedge,
    const unsigned int /* k */) {

    std::cout.precision(20);
    if (m_verbose) {
      std::cout << "** propagating pedge [" << str(pvertex) << "-" << str(pother)
      << "] along " << str(iedge) << std::endl;
    }
    CGAL_assertion_msg(false, "TODO: PROPAGATE PEDGE BEYOND IEDGE!");
    return std::make_pair(null_pvertex(), null_pvertex());
  }

  const bool transfer_pvertex_via_iedge(
    const PVertex& pvertex, const PVertex& pother) {

    std::cout.precision(20);
    if (m_verbose) {

      std::cout << "** transfering " <<
      str(pother) << " through " << str(pvertex) << " via "
      << str(this->iedge(pvertex)) << std::endl;
    }

    // If pvertex is adjacent to one or two.
    PFace source_pface, target_pface;
    std::tie(source_pface, target_pface) = pfaces_of_pvertex(pvertex);
    const PFace common_pface = pface_of_pvertex(pother);

    if (common_pface == target_pface) {
      std::swap(source_pface, target_pface);
    }
    CGAL_assertion(common_pface == source_pface);

    if (m_verbose) {
      std::cout << "- initial pfaces: " <<
      lstr(source_pface) << " and " << lstr(target_pface) << std::endl;
    }

    // std::cout << "pvertex: " << point_3(pvertex) << std::endl;
    // std::cout << "pother: "  << point_3(pother) << std::endl;

    // if (source_pface != null_pface()) {
    //   std::cout << "source pface center: " << centroid_of_pface(source_pface) << std::endl;
    // }
    // if (target_pface != null_pface()) {
    //   std::cout << "target pface center: " << centroid_of_pface(target_pface) << std::endl;
    // }

    PVertex pthird = next(pother);
    if (pthird == pvertex) {
      pthird = prev(pother);
    }

    CGAL_assertion(has_iedge(pvertex));
    if (target_pface == null_pface()) {

      const Line_2 iedge_line = segment_2(pother.first, this->iedge(pvertex)).supporting_line();
      const Point_2 pinit = iedge_line.projection(point_2(pother, m_current_time));

      const Line_2 future_line(
        point_2(pother, m_current_time + FT(1)),
        point_2(pthird, m_current_time + FT(1)));
      CGAL_assertion_msg(!CGAL::parallel(future_line, iedge_line),
      "TODO: TRANSFER PVERTEX, HANDLE CASE WITH PARALLEL LINES!");
      Point_2 future_point = KSR::intersection<Point_2>(future_line, iedge_line);

      const Vector_2 future_direction(pinit, future_point);
      direction(pvertex) = future_direction;
      future_point = pinit - future_direction * m_current_time;
      support_plane(pvertex).set_point(pvertex.second, future_point);

      const auto he = mesh(pvertex).halfedge(pother.second, pvertex.second);
      CGAL::Euler::join_vertex(he, mesh(pvertex));

    } else {

      // std::cout << "disconnecting " <<
      // str(pvertex) << " from " << str(this->iedge(pvertex)) << std::endl;
      const IEdge iedge = disconnect_iedge(pvertex);
      PEdge pedge = null_pedge();
      for (const auto edge : pedges_around_pvertex(pvertex)) {
        if (this->iedge(edge) == iedge) {
          pedge = edge;
          break;
        }
      }
      CGAL_assertion(pedge != null_pedge());

      auto he = mesh(pedge).halfedge(pedge.second);
      if (mesh(pedge).face(he) != common_pface.second) {
        he = mesh(pedge).opposite(he);
      }
      CGAL_assertion(mesh(pedge).face(he) == common_pface.second);

      if (mesh(pedge).target(he) == pvertex.second) {
        // std::cout << "shifting target" << std::endl;
        CGAL::Euler::shift_target(he, mesh(pedge));
      } else {
        CGAL_assertion(mesh(pedge).source(he) == pvertex.second);
        // std::cout << "shifting source" << std::endl;
        CGAL::Euler::shift_source(he, mesh(pedge));
      }

      const Line_2 iedge_line = segment_2(pother.first, iedge).supporting_line();
      const Point_2 pinit = iedge_line.projection(point_2(pother, m_current_time));

      direction(pvertex) = direction(pother);
      support_plane(pother).set_point(
        pvertex.second, pinit - direction(pvertex) * m_current_time);

      const Line_2 future_line(
        point_2(pvertex, m_current_time + FT(1)),
        point_2(pthird , m_current_time + FT(1)));
      CGAL_assertion_msg(!CGAL::parallel(future_line, iedge_line),
      "TODO: TRANSFER PVERTEX, HANDLE CASE WITH PARALLEL LINES!");
      Point_2 future_point = KSR::intersection<Point_2>(future_line, iedge_line);

      const Vector_2 future_direction(pinit, future_point);
      direction(pother) = future_direction;
      future_point = pinit - future_direction * m_current_time;
      support_plane(pother).set_point(pother.second, future_point);

      // std::cout << "connecting " << str(pother) << " to " << str(iedge) << std::endl;
      connect(pother, iedge);
    }

    if (m_verbose) {
      std::cout << "- new pfaces: " <<
      lstr(source_pface) << " and " << lstr(target_pface) << std::endl;
    }
    return (target_pface != null_pface());
  }

  const std::vector<PVertex> merge_pvertices_on_ivertex(
    const FT min_time, const FT max_time,
    const PVertex& event_pvertex,
    const IVertex& ivertex,
    const std::vector<PVertex>& pvertices,
    std::vector<IEdge>& crossed) {

    std::cout.precision(20);
    if (m_verbose) {
      std::cout << "** merging " << str(event_pvertex) << " on " << str(ivertex) << std::endl;
    }

    // std::cout << "event pvertex: " << point_3(event_pvertex) << std::endl;
    // std::cout << "ivertex: " << point_3(ivertex) << std::endl;

    CGAL_assertion(pvertices.size() >= 3);
    const KSR::size_t support_plane_idx = pvertices.front().first;
    const PVertex prev = pvertices.front();
    const PVertex next = pvertices.back();
    const PVertex pvertex = pvertices[1];

    if (m_verbose) {
      std::cout << "- starting from: " <<
      str(iedge(pvertex)) << " " << segment_3(iedge(pvertex)) << std::endl;
    }

    // Copy front/back to remember position/direction.
    PVertex front, back;
    if (pvertices.size() < 3) {
      CGAL_assertion_msg(false, "ERROR: INVALID CASE!");
    } else if (pvertices.size() == 3 || pvertices.size() == 4) {

      // BUG: In this case, the point that is duplicated twice is not always copied.
      // To fix it, we copy the second point not from the original vertex but from the first
      // copy of that vertex.

      const auto& initial = pvertex;
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

    if (m_verbose) {
      std::cout << "- neighbors: " << std::endl <<
      "prev = " << point_3(prev)  << " / " << direction(prev)  << std::endl <<
      "fron = " << point_3(front) << " / " << direction(front) << std::endl <<
      "back = " << point_3(back)  << " / " << direction(back)  << std::endl <<
      "next = " << point_3(next)  << " / " << direction(next)  << std::endl;
    }

    // Freeze pvertices.
    const Point_2 ipoint = point_2(support_plane_idx, ivertex);
    for (std::size_t i = 1; i < pvertices.size() - 1; ++i) {
      const PVertex& curr = pvertices[i];
      support_plane(curr).direction(curr.second) = CGAL::NULL_VECTOR;
      support_plane(curr).set_point(curr.second, ipoint);
    }
    connect(pvertex, ivertex);
    if (m_verbose) {
      std::cout << "- frozen pvertex: "
      << str(pvertex) << " : " << point_3(pvertex) << std::endl;
    }

    // Join pvertices.
    // std::cout << "removed pvertices:";
    for (std::size_t i = 2; i < pvertices.size() - 1; ++i) {
      // std::cout << " " << str(pvertices[i]) << std::endl;
      // std::cout << point_3(pvertices[i]) << std::endl;
      const auto he = mesh(support_plane_idx).halfedge(pvertices[i].second, pvertex.second);
      disconnect_ivertex(pvertices[i]);
      CGAL::Euler::join_vertex(he, mesh(support_plane_idx));
    }
    // std::cout << std::endl;

    // Get all connected iedges.
    auto inc_iedges = this->incident_iedges(ivertex);
    std::vector< std::pair<IEdge, Direction_2> > iedges;
    std::copy(inc_iedges.begin(), inc_iedges.end(),
      boost::make_function_output_iterator(
        [&](const IEdge& inc_iedge) -> void {
          const auto iplanes = this->intersected_planes(inc_iedge);
          if (iplanes.find(support_plane_idx) == iplanes.end()) {
            return;
          }
          const Direction_2 direction(
            point_2(support_plane_idx, opposite(inc_iedge, ivertex)) -
            point_2(support_plane_idx, ivertex));
          iedges.push_back(std::make_pair(inc_iedge, direction));
        }
      )
    );

    std::sort(iedges.begin(), iedges.end(),
      [&](const std::pair<IEdge, Direction_2>& a,
          const std::pair<IEdge, Direction_2>& b) -> bool {
        return a.second < b.second;
      }
    );
    CGAL_assertion(iedges.size() != 0);

    // Get sub-event type.
    bool back_constrained = false;
    if (
      (iedge(next) != null_iedge() && (source(iedge(next)) == ivertex || target(iedge(next)) == ivertex)) ||
      (this->ivertex(next) != null_ivertex() && is_iedge(this->ivertex(next), ivertex))) {
      back_constrained = true;
    }

    bool front_constrained = false;
    if (
      (iedge(prev) != null_iedge() && (source(iedge(prev)) == ivertex || target(iedge(prev)) == ivertex)) ||
      (this->ivertex(prev) != null_ivertex() && is_iedge(this->ivertex(prev), ivertex))) {
      front_constrained = true;
    }

    if (back_constrained && !front_constrained) {
      std::reverse(iedges.begin(), iedges.end());
    }

    if (m_verbose) {
      std::cout << "- initial iedges: " << std::endl;
      for (const auto& iedge : iedges) {
        std::cout << segment_3(iedge.first) << std::endl;
      }
    }

    // Handle sub-events.
    crossed.clear();
    std::vector<PVertex> new_pvertices;

    if (back_constrained && front_constrained) {
      apply_closing_case();
    } else if (back_constrained) {
      apply_back_border_case(
        min_time, max_time,
        pvertex, ivertex,
        prev, back,
        iedges, crossed, new_pvertices);
    } else if (front_constrained) {
      apply_front_border_case(
        min_time, max_time,
        pvertex, ivertex,
        next, front,
        iedges, crossed, new_pvertices);
    } else {
      apply_open_case(
        min_time, max_time,
        pvertex, ivertex,
        prev, next,
        iedges, crossed, new_pvertices);
    }

    support_plane(support_plane_idx).remove_vertex(front.second);
    support_plane(support_plane_idx).remove_vertex(back.second);

    // Push also remaining vertex so that its events are recomputed.
    new_pvertices.push_back(pvertex);
    crossed.push_back(iedge(pvertex));

    if (m_verbose) {
      std::cout << "- new pvertices:";
      for (const PVertex& pv : new_pvertices)
        std::cout << " " << str(pv);
      std::cout << std::endl;
    }
    return new_pvertices;
  }

  void apply_closing_case() {

    std::cout.precision(20);
    if (m_verbose) {
      std::cout << "*** CLOSING CASE" << std::endl;
    }
  }

  void apply_back_border_case(
    const FT min_time, const FT max_time,
    const PVertex& pvertex,
    const IVertex& ivertex,
    const PVertex& prev,
    const PVertex& back,
    const std::vector< std::pair<IEdge, Direction_2> >& iedges,
    std::vector<IEdge>& crossed,
    std::vector<PVertex>& new_pvertices) {

    std::cout.precision(20);
    if (m_verbose) {
      std::cout << "*** BACK BORDER CASE" << std::endl;
    }

    // We use this modification in order to avoid collinear directions.
    CGAL_assertion(has_iedge(pvertex));
    const KSR::size_t other_side_limit = line_idx(pvertex);
    const FT prev_time = last_event_time(prev);
    CGAL_assertion(prev_time < m_current_time);
    CGAL_assertion(prev_time >= FT(0));

    const auto pp_last = point_2(prev, prev_time);
    const auto pp_curr = point_2(prev, m_current_time);
    const auto dirp = Vector_2(pp_last, pp_curr);
    const auto shifted_prev = pp_curr - dirp / FT(10);
    // std::cout << "shifted prev: " << to_3d(pvertex.first, shifted_prev) << std::endl;

    const auto ipoint = point_2(pvertex.first, ivertex);
    const Direction_2 ref_direction_prev(shifted_prev - ipoint);

    // Find the first iedge.
    KSR::size_t first_idx = KSR::no_element();
    const std::size_t n = iedges.size();
    for (std::size_t i = 0; i < n; ++i) {
      const std::size_t ip = (i + 1) % n;

      const auto& i_dir  = iedges[i].second;
      const auto& ip_dir = iedges[ip].second;
      if (ref_direction_prev.counterclockwise_in_between(ip_dir, i_dir)) {
        first_idx = ip;
        break;
      }
    }
    CGAL_assertion(first_idx != KSR::no_element());
    // std::cout << "first: " << segment_3(iedges[first_idx].first) << std::endl;

    // Find all crossed iedges.
    CGAL_assertion(crossed.size() == 0);
    KSR::size_t iedge_idx = first_idx;
    std::size_t iteration = 0;
    while (true) {
      const auto& iedge = iedges[iedge_idx].first;
      // std::cout << "next: " << segment_3(iedge) << std::endl;

      const bool bbox_reached  = ( collision_occured(pvertex, iedge)   ).second;
      const bool limit_reached = ( line_idx(iedge) == other_side_limit );
      if (m_verbose) {
        std::cout << "- limit/bbox: " << limit_reached << "/" << bbox_reached << std::endl;
      }

      crossed.push_back(iedge);
      if (limit_reached || bbox_reached) {
        break;
      }
      iedge_idx = (iedge_idx + 1) % n;
      if (iteration == 100) {
        CGAL_assertion_msg(false, "ERROR: BACK, WHY SO MANY ITERATIONS?");
      } ++iteration;
    }

    CGAL_assertion(crossed.size() != 0);
    if (m_verbose) {
      std::cout << "- crossed " << crossed.size() << " iedges:" << std::endl;
      for (const auto& iedge : crossed) {
        std::cout << segment_3(iedge) << std::endl;
      }
    }

    // Compute future points and directions.
    std::vector<Point_2> future_points(crossed.size());
    std::vector<Vector_2> future_directions(crossed.size());

    IEdge prev_iedge = null_iedge();
    for (std::size_t i = 0; i < crossed.size(); ++i) {
      const bool is_parallel = compute_future_point_and_direction(
        i, back, prev, crossed[i], future_points[i], future_directions[i]);
      if (is_parallel) {
        if (is_intersecting_iedge(min_time, max_time, prev, crossed[i])) {
          CGAL_assertion_msg(i == 0, "TODO: BACK, CAN WE HAVE NON-ZERO I HERE?");
          prev_iedge = crossed[i];
        }
      }
    }

    // Crop/propagate the pvertex.
    PVertex previous = null_pvertex();
    for (std::size_t i = 0; i < crossed.size(); ++i) {
      if (i == 0) {

        if (m_verbose) std::cout << "- cropping" << std::endl;
        PVertex cropped;
        if (prev_iedge != null_iedge() && prev_iedge == crossed[i]) {
          if (m_verbose) std::cout << "- prev, parallel case" << std::endl;

          // In case, we are parallel, we update the future point and direction.
          cropped = prev;
          Point_2 future_point; Vector_2 future_direction;
          const auto pprev = ( border_prev_and_next(prev) ).first;
          compute_future_point_and_direction(
            i, prev, pprev, prev_iedge, future_point, future_direction);
          future_points[i] = future_point;
          future_directions[i] = future_direction;

        } else {
          if (m_verbose) std::cout << "- prev, standard case" << std::endl;
          cropped = PVertex(pvertex.first, support_plane(pvertex).split_edge(pvertex.second, prev.second));
        }

        const PEdge pedge(pvertex.first, support_plane(pvertex).edge(pvertex.second, cropped.second));
        CGAL_assertion(cropped != pvertex);
        new_pvertices.push_back(cropped);

        connect(pedge, crossed[i]);
        connect(cropped, crossed[i]);

        support_plane(cropped).set_point(cropped.second, future_points[i]);
        direction(cropped) = future_directions[i];
        previous = cropped;
        if (m_verbose) std::cout << "- cropped: " << point_3(cropped) << std::endl;

      } else {
        if (m_verbose) std::cout << "- propagating" << std::endl;
        CGAL_assertion_msg(i == 1,
        "TODO: BACK, CAN WE HAVE MORE THAN 1 NEW PFACE? IF YES, I SHOULD CHECK K FOR EACH!");

        // Now, we check if we should add a new pface.
        bool is_occupied_edge, bbox_reached;
        std::tie(is_occupied_edge, bbox_reached) = is_occupied(pvertex, ivertex, crossed[i - 1]);
        if (m_verbose) {
          std::cout << "- is already occupied / bbox: "
          << is_occupied_edge << "/" << bbox_reached << std::endl;
        }

        // Stop propagating.
        const auto pface = pface_of_pvertex(pvertex);
        if (m_verbose) std::cout << "- k intersections befor: " << this->k(pface) << std::endl;
        if (bbox_reached) {
          if (m_verbose) std::cout << "- stop bbox" << std::endl;
          CGAL_assertion_msg(false, "ERROR: BACK, THIS CASE CANNOT HAPPEN!");
          break;
        } else if (is_occupied_edge && this->k(pface) == 1) {
          if (m_verbose) std::cout << "- stop k" << std::endl;
          break;
        }

        // Create a new pface.
        if (m_verbose) std::cout << "- adding new pface" << std::endl;
        if (is_occupied_edge && this->k(pface) > 1) {
          if (m_verbose) std::cout << "- continue k > 1" << std::endl;
          this->k(pface)--;
        } else {
          if (m_verbose) std::cout << "- continue k = 1" << std::endl;
        }
        CGAL_assertion(this->k(pface) >= 1);

        if (m_verbose) {
          // std::cout << "PFACE: " << centroid_of_pface(pface) << std::endl;
          std::cout << "- k intersections after: " << this->k(pface) << std::endl;
        }

        const PVertex propagated = add_pvertex(pvertex.first, future_points[i]);
        direction(propagated) = future_directions[i];
        CGAL_assertion(propagated != pvertex);
        new_pvertices.push_back(propagated);

        if (m_verbose) std::cout << "- propagated: " << point_3(propagated) << std::endl;
        const PFace new_pface = add_pface(std::array<PVertex, 3>{pvertex, propagated, previous});
        if (m_verbose) std::cout << "- new pface: " << lstr(new_pface) << std::endl;
        this->k(new_pface) = this->k(pface);
        previous = propagated;

        const PEdge pedge(pvertex.first, support_plane(pvertex).edge(pvertex.second, propagated.second));
        connect(pedge, crossed[i]);
        connect(propagated, crossed[i]);
      }
    }
  }

  void apply_front_border_case(
    const FT min_time, const FT max_time,
    const PVertex& pvertex,
    const IVertex& ivertex,
    const PVertex& next,
    const PVertex& front,
    const std::vector< std::pair<IEdge, Direction_2> >& iedges,
    std::vector<IEdge>& crossed,
    std::vector<PVertex>& new_pvertices) {

    std::cout.precision(20);
    if (m_verbose) {
      std::cout << "*** FRONT BORDER CASE" << std::endl;
    }

    // We use this modification in order to avoid collinear directions.
    CGAL_assertion(has_iedge(pvertex));
    const KSR::size_t other_side_limit = line_idx(pvertex);
    const FT next_time = last_event_time(next);
    CGAL_assertion(next_time < m_current_time);
    CGAL_assertion(next_time >= FT(0));

    const auto pn_last = point_2(next, next_time);
    const auto pn_curr = point_2(next, m_current_time);
    const auto dirn = Vector_2(pn_last, pn_curr);
    const auto shifted_next = pn_curr - dirn / FT(10);
    // std::cout << "shifted next: " << to_3d(pvertex.first, shifted_next) << std::endl;

    const auto ipoint = point_2(pvertex.first, ivertex);
    const Direction_2 ref_direction_next(shifted_next - ipoint);

    // Find the first iedge.
    KSR::size_t first_idx = KSR::no_element();
    const std::size_t n = iedges.size();
    for (std::size_t i = 0; i < n; ++i) {
      const std::size_t ip = (i + 1) % n;

      const auto& i_dir  = iedges[i].second;
      const auto& ip_dir = iedges[ip].second;
      if (ref_direction_next.counterclockwise_in_between(i_dir, ip_dir)) {
        first_idx = ip;
        break;
      }
    }
    CGAL_assertion(first_idx != KSR::no_element());
    // std::cout << "first: " << segment_3(iedges[first_idx].first) << std::endl;

    // Find all crossed iedges.
    CGAL_assertion(crossed.size() == 0);
    KSR::size_t iedge_idx = first_idx;
    std::size_t iteration = 0;
    while (true) {
      const auto& iedge = iedges[iedge_idx].first;
      // std::cout << "next: " << segment_3(iedge) << std::endl;

      const bool bbox_reached  = ( collision_occured(pvertex, iedge)   ).second;
      const bool limit_reached = ( line_idx(iedge) == other_side_limit );
      if (m_verbose) {
        std::cout << "- limit/bbox: " << limit_reached << "/" << bbox_reached << std::endl;
      }

      crossed.push_back(iedge);
      if (limit_reached || bbox_reached) {
        break;
      }
      iedge_idx = (iedge_idx + 1) % n;
      if (iteration == 100) {
        CGAL_assertion_msg(false, "ERROR: FRONT, WHY SO MANY ITERATIONS?");
      } ++iteration;
    }

    CGAL_assertion(crossed.size() != 0);
    if (m_verbose) {
      std::cout << "- crossed " << crossed.size() << " iedges:" << std::endl;
      for (const auto& iedge : crossed) {
        std::cout << segment_3(iedge) << std::endl;
      }
    }

    // Compute future points and directions.
    std::vector<Point_2> future_points(crossed.size());
    std::vector<Vector_2> future_directions(crossed.size());

    IEdge next_iedge = null_iedge();
    for (std::size_t i = 0; i < crossed.size(); ++i) {
      const bool is_parallel = compute_future_point_and_direction(
        i, front, next, crossed[i], future_points[i], future_directions[i]);
      if (is_parallel) {
        if (is_intersecting_iedge(min_time, max_time, next, crossed[i])) {
          CGAL_assertion_msg(i == 0, "TODO: FRONT, CAN WE HAVE NON-ZERO I HERE?");
          next_iedge = crossed[i];
        }
      }
    }

    // Crop/propagate the pvertex.
    PVertex previous = null_pvertex();
    for (std::size_t i = 0; i < crossed.size(); ++i) {
      if (i == 0) {
        if (m_verbose) std::cout << "- cropping" << std::endl;

        PVertex cropped;
        if (next_iedge != null_iedge() && next_iedge == crossed[i]) {
          if (m_verbose) std::cout << "- next, parallel case" << std::endl;

          // In case, we are parallel, we update the future point and direction.
          cropped = next;
          Point_2 future_point; Vector_2 future_direction;
          const auto nnext = ( border_prev_and_next(next) ).second;
          compute_future_point_and_direction(
            i, next, nnext, next_iedge, future_point, future_direction);
          future_points[i] = future_point;
          future_directions[i] = future_direction;

        } else {
          if (m_verbose) std::cout << "- next, standard case" << std::endl;
          cropped = PVertex(pvertex.first, support_plane(pvertex).split_edge(pvertex.second, next.second));
        }

        const PEdge pedge(pvertex.first, support_plane(pvertex).edge(pvertex.second, cropped.second));
        CGAL_assertion(cropped != pvertex);
        new_pvertices.push_back(cropped);

        connect(pedge, crossed[i]);
        connect(cropped, crossed[i]);

        support_plane(cropped).set_point(cropped.second, future_points[i]);
        direction(cropped) = future_directions[i];
        previous = cropped;
        if (m_verbose) std::cout << "- cropped: " << point_3(cropped) << std::endl;

      } else {
        if (m_verbose) std::cout << "- propagating" << std::endl;
        CGAL_assertion_msg(i == 1,
        "TODO: FRONT, CAN WE HAVE MORE THAN 1 NEW PFACE? IF YES, I SHOULD CHECK K FOR EACH!");

        // Now, we check if we should add a new pface.
        bool is_occupied_edge, bbox_reached;
        std::tie(is_occupied_edge, bbox_reached) = is_occupied(pvertex, ivertex, crossed[i - 1]);
        if (m_verbose) {
          std::cout << "- is already occupied / bbox: " << is_occupied_edge << "/" << bbox_reached << std::endl;
        }

        // Stop propagating.
        const auto pface = pface_of_pvertex(pvertex);
        if (m_verbose) std::cout << "- k intersections befor: " << this->k(pface) << std::endl;
        if (bbox_reached) {
          if (m_verbose) std::cout << "- stop bbox" << std::endl;
          CGAL_assertion_msg(false, "ERROR: FRONT, THIS CASE CANNOT HAPPEN!");
          break;
        } else if (is_occupied_edge && this->k(pface) == 1) {
          if (m_verbose) std::cout << "- stop k" << std::endl;
          break;
        }

        // Create a new pface.
        if (m_verbose) std::cout << "- adding new pface" << std::endl;
        if (is_occupied_edge && this->k(pface) > 1) {
          if (m_verbose) std::cout << "- continue k > 1" << std::endl;
          this->k(pface)--;
        } else {
          if (m_verbose) std::cout << "- continue k = 1" << std::endl;
        }
        CGAL_assertion(this->k(pface) >= 1);

        if (m_verbose) {
          // std::cout << "PFACE: " << centroid_of_pface(pface) << std::endl;
          std::cout << "- k intersections after: " << this->k(pface) << std::endl;
        }

        const PVertex propagated = add_pvertex(pvertex.first, future_points[i]);
        direction(propagated) = future_directions[i];
        CGAL_assertion(propagated != pvertex);
        new_pvertices.push_back(propagated);

        if (m_verbose) std::cout << "- propagated: " << point_3(propagated) << std::endl;
        const PFace new_pface = add_pface(std::array<PVertex, 3>{pvertex, previous, propagated});
        if (m_verbose) std::cout << "- new pface: " << lstr(new_pface) << std::endl;
        this->k(new_pface) = this->k(pface);
        previous = propagated;

        const PEdge pedge(pvertex.first, support_plane(pvertex).edge(pvertex.second, propagated.second));
        connect(pedge, crossed[i]);
        connect(propagated, crossed[i]);
      }
    }
  }

  void apply_open_case(
    const FT min_time, const FT max_time,
    const PVertex& pvertex,
    const IVertex& ivertex,
    const PVertex& prev,
    const PVertex& next,
    const std::vector< std::pair<IEdge, Direction_2> >& iedges,
    std::vector<IEdge>& crossed,
    std::vector<PVertex>& new_pvertices) {

    std::cout.precision(20);
    if (m_verbose) {
      std::cout << "*** OPEN CASE" << std::endl;
    }

    // We use this modification in order to avoid collinear directions.
    const FT prev_time = last_event_time(prev);
    const FT next_time = last_event_time(next);
    CGAL_assertion(prev_time < m_current_time);
    CGAL_assertion(next_time < m_current_time);
    CGAL_assertion(prev_time >= FT(0));
    CGAL_assertion(next_time >= FT(0));

    const auto pp_last = point_2(prev, prev_time);
    const auto pp_curr = point_2(prev, m_current_time);
    const auto dirp = Vector_2(pp_last, pp_curr);
    const auto shifted_prev = pp_curr - dirp / FT(10);

    const auto pn_last = point_2(next, next_time);
    const auto pn_curr = point_2(next, m_current_time);
    const auto dirn = Vector_2(pn_last, pn_curr);
    const auto shifted_next = pn_curr - dirn / FT(10);

    // std::cout << "shifted prev: " << to_3d(pvertex.first, shifted_prev) << std::endl;
    // std::cout << "shifted next: " << to_3d(pvertex.first, shifted_next) << std::endl;

    const auto ipoint = point_2(pvertex.first, ivertex);
    const Direction_2 ref_direction_prev(shifted_prev - ipoint);
    const Direction_2 ref_direction_next(shifted_next - ipoint);

    // Find the first iedge.
    KSR::size_t first_idx = KSR::no_element();
    const std::size_t n = iedges.size();
    for (std::size_t i = 0; i < n; ++i) {
      const std::size_t ip = (i + 1) % n;

      const auto& i_dir  = iedges[i].second;
      const auto& ip_dir = iedges[ip].second;
      if (ref_direction_next.counterclockwise_in_between(i_dir, ip_dir)) {
        first_idx = ip;
        break;
      }
    }
    CGAL_assertion(first_idx != KSR::no_element());
    // std::cout << "first: " << segment_3(iedges[first_idx].first) << std::endl;

    // Find all crossed iedges.
    CGAL_assertion(crossed.size() == 0);
    KSR::size_t iedge_idx = first_idx;
    std::size_t iteration = 0;
    while (true) {
      const auto& iedge = iedges[iedge_idx].first;
      // std::cout << "next: " << segment_3(iedge) << std::endl;

      const auto& ref_direction = iedges[iedge_idx].second;
      if (!ref_direction.counterclockwise_in_between(
        ref_direction_next, ref_direction_prev)) {
        break;
      }

      crossed.push_back(iedge);
      iedge_idx = (iedge_idx + 1) % n;
      if (iteration == 100) {
        CGAL_assertion_msg(false, "ERROR: OPEN, WHY SO MANY ITERATIONS?");
      } ++iteration;
    }

    CGAL_assertion(crossed.size() != 0);
    if (m_verbose) {
      std::cout << "- crossed " << crossed.size() << " iedges: " << std::endl;
      for (const auto& iedge : crossed) {
        std::cout << segment_3(iedge) << std::endl;
      }
    }

    // Compute future points and directions.
    std::vector<Point_2> future_points(crossed.size());
    std::vector<Vector_2> future_directions(crossed.size());

    IEdge prev_iedge = null_iedge(), next_iedge = null_iedge();
    for (std::size_t i = 0; i < crossed.size(); ++i) {
      const bool is_parallel = compute_future_point_and_direction(
        pvertex, prev, next, crossed[i], future_points[i], future_directions[i]);
      if (is_parallel) {
        if (is_intersecting_iedge(min_time, max_time, prev, crossed[i])) {
          CGAL_assertion_msg(i == crossed.size() - 1, "TODO: FRONT, CAN WE HAVE OTHER I HERE?");
          prev_iedge = crossed[i];
        }
        if (is_intersecting_iedge(min_time, max_time, next, crossed[i])) {
          CGAL_assertion_msg(i == 0, "TODO: FRONT, CAN WE HAVE OTHER I HERE?");
          next_iedge = crossed[i];
        }
      }
    }

    // Crop/propagate the pvertex.
    { // first crop
      PVertex cropped;
      if (next_iedge != null_iedge() && next_iedge == crossed.front()) {
        if (m_verbose) std::cout << "- next, parallel case" << std::endl;

        // In case, we are parallel, we update the future point and direction.
        cropped = next;
        Point_2 future_point; Vector_2 future_direction;
        const auto nnext = ( border_prev_and_next(next) ).second;
        compute_future_point_and_direction(
          0, next, nnext, next_iedge, future_point, future_direction);
        future_points[0] = future_point;
        future_directions[0] = future_direction;

      } else {
        if (m_verbose) std::cout << "- next, standard case" << std::endl;
        cropped = PVertex(pvertex.first, support_plane(pvertex).split_edge(pvertex.second, next.second));
      }

      const PEdge pedge(pvertex.first, support_plane(pvertex).edge(pvertex.second, cropped.second));
      CGAL_assertion(cropped != pvertex);
      new_pvertices.push_back(cropped);

      connect(pedge, crossed.front());
      connect(cropped, crossed.front());

      support_plane(cropped).set_point(cropped.second, future_points.front());
      direction(cropped) = future_directions.front();
      if (m_verbose) std::cout << "- cropped 1: " << point_3(cropped) << std::endl;
    }

    { // then propagate
      for (std::size_t i = 1; i < crossed.size() - 1; ++i) {
        const PVertex propagated = add_pvertex(pvertex.first, future_points[i]);
        direction(propagated) = future_directions[i];
        connect(propagated, crossed[i]);

        CGAL_assertion(propagated != pvertex);
        new_pvertices.push_back(propagated);
        if (m_verbose) {
          std::cout << "- propagated " << std::to_string(i) << ": " << point_3(propagated) << std::endl;
        }
      }
    }

    { // then crop again
      PVertex cropped;
      if (prev_iedge != null_iedge() && prev_iedge == crossed.back()) {
        if (m_verbose) std::cout << "- prev, parallel case" << std::endl;

        // In case, we are parallel, we update the future point and direction.
        cropped = prev;
        Point_2 future_point; Vector_2 future_direction;
        const auto pprev = ( border_prev_and_next(prev) ).first;
        compute_future_point_and_direction(
          0, prev, pprev, prev_iedge, future_point, future_direction);
        future_points[future_points.size() - 1] = future_point;
        future_directions[future_directions.size() - 1] = future_direction;

      } else {
        if (m_verbose) std::cout << "- prev, standard case" << std::endl;
        cropped = PVertex(pvertex.first, support_plane(pvertex).split_edge(pvertex.second, prev.second));
      }

      const PEdge pedge(pvertex.first, support_plane(pvertex).edge(pvertex.second, cropped.second));
      CGAL_assertion(cropped != pvertex);
      new_pvertices.push_back(cropped);

      connect(pedge, crossed.back());
      connect(cropped, crossed.back());

      support_plane(cropped).set_point(cropped.second, future_points.back());
      direction(cropped) = future_directions.back();
      if (m_verbose) std::cout << "- cropped 2: " << point_3(cropped) << std::endl;
    }

    if (m_verbose) std::cout << "- new pvertices size: " << new_pvertices.size() << std::endl;
    CGAL_assertion(new_pvertices.size() == crossed.size());

    // Now, we check if we should add new pfaces.
    bool is_occupied_edge_back, bbox_reached_back;
    std::tie(is_occupied_edge_back, bbox_reached_back) = is_occupied(pvertex, ivertex, crossed.back());
    if (m_verbose) {
      std::cout << "- is already occupied back / bbox: " << is_occupied_edge_back << "/" << bbox_reached_back << std::endl;
    }

    bool is_occupied_edge_front, bbox_reached_front;
    std::tie(is_occupied_edge_front, bbox_reached_front) = is_occupied(pvertex, ivertex, crossed.front());
    if (m_verbose) {
      std::cout << "- is already occupied fron / bbox: " << is_occupied_edge_front << "/" << bbox_reached_front << std::endl;
    }

    const auto pface = pface_of_pvertex(pvertex);
    if (m_verbose) std::cout << "- k intersections befor: " << this->k(pface) << std::endl;
    if (bbox_reached_back) {

      CGAL_assertion(bbox_reached_front);
      if (m_verbose) std::cout << "- stop bbox back" << std::endl;

    } else if (bbox_reached_front) {

      CGAL_assertion(bbox_reached_back);
      if (m_verbose) std::cout << "- stop bbox front" << std::endl;

    } else if ((is_occupied_edge_back && is_occupied_edge_front) && this->k(pface) == 1) {

      if (m_verbose) std::cout << "- stop back && front k = 1" << std::endl;

    } else if ((is_occupied_edge_back && is_occupied_edge_front) && this->k(pface) > 1) {

      this->k(pface)--;
      CGAL_assertion(this->k(pface) >= 1);
      add_new_pfaces(this->k(pface), pvertex, ivertex, new_pvertices, pface, crossed);
      if (m_verbose) std::cout << "- continue back && front k > 1" << std::endl;

    } else if ((!is_occupied_edge_back && !is_occupied_edge_front)) {

      add_new_pfaces(this->k(pface), pvertex, ivertex, new_pvertices, pface, crossed);
      if (m_verbose) std::cout << "- continue !back && !front" << std::endl;

    } else if (is_occupied_edge_back || is_occupied_edge_front) {

      add_new_pfaces(this->k(pface), pvertex, ivertex, new_pvertices, pface, crossed);
      if (m_verbose) std::cout << "- continue back || front" << std::endl;

      // std::cout << "pver pface: " << str(pface_of_pvertex(pvertex))      << std::endl;
      // std::cout << "back pface: " << str(pface_of_pvertex(pvertices[1])) << std::endl;
      // std::cout << "fron pface: " << str(pface_of_pvertex(pvertices[2])) << std::endl;
      // CGAL_assertion_msg(false, "TEST THIS CASE: BACK || FRONT!");

    } else {
      CGAL_assertion_msg(false, "TODO: ADD NEW OPEN CASE! DO NOT FORGET TO UPDATE K!");
    }

    if (m_verbose) {
      // std::cout << "PFACE: " << centroid_of_pface(pface) << std::endl;
      std::cout << "- k intersections after: " << this->k(pface) << std::endl;
    }
  }

  void add_new_pfaces(
    const unsigned int k,
    const PVertex& pvertex,
    const IVertex& ivertex,
    const std::vector<PVertex>& new_pvertices,
    const PFace& pface,
    const std::vector<IEdge>& crossed) {

    CGAL_assertion(new_pvertices.size() >= 2);
    CGAL_assertion(crossed.size() == new_pvertices.size());

    std::size_t num_added_pfaces = 0;
    for (std::size_t i = 0; i < new_pvertices.size() - 1; ++i) {

      if (i >= 1) {
        // bool is_occupied_edge, bbox_reached;
        // std::tie(is_occupied_edge, bbox_reached) = is_occupied(pvertex, ivertex, crossed[i]);

        // if (bbox_reached) return;
        // if (is_occupied_edge && this->k(pface) == 1) return;
        // if (is_occupied_edge && this->k(pface) > 1) {
        //   this->k(pface)--;
        //   CGAL_assertion(this->k(pface) >= 1);
        // }

        const PEdge pedge(pvertex.first, support_plane(pvertex).edge(pvertex.second, new_pvertices[i].second));
        connect(pedge, crossed[i]);
        connect(new_pvertices[i], crossed[i]);
      }

      if (m_verbose) {
        std::cout << "- adding new pface" << std::endl;
      }
      const PFace new_pface = add_pface(std::array<PVertex, 3>{new_pvertices[i], new_pvertices[i + 1], pvertex});
      if (m_verbose) std::cout << "- new pface: " << lstr(new_pface) << std::endl;
      this->k(new_pface) = k;
      ++num_added_pfaces;
    }
    CGAL_assertion(num_added_pfaces > 0);
    CGAL_assertion_msg(num_added_pfaces == 1,
    "TODO: OPEN, CAN WE HAVE MORE THAN 1 NEW PFACE? IF YES, I SHOULD CHECK K FOR EACH!");
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

      CGAL_assertion_msg(false,
      "TODO: FIND_OPPOSITE_PVERTEX, CAN WE HAVE THIS CASE?");

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

  /*******************************
  **          CLEANING          **
  ********************************/

  void finalize() {

    bool quit = true;
    std::size_t num_removed_faces = 0;
    do {
      quit = true;
      for (const auto iedge : m_intersection_graph.edges()) {
        const std::size_t num_faces = check_edge(iedge);
        if (num_faces != 0) {
          num_removed_faces += num_faces;
          quit = false; break;
        }
      }
    } while (!quit);
    if (m_verbose) {
      std::cout << "* number of removed hanging faces: " << num_removed_faces << std::endl;
    }
    // CGAL_assertion_msg(false, "TODO: DEBUG THIS FUNCTION!");

    // TODO: Should I also implement here the part that removes all
    // identical pfaces within the same support plane? If the k intersection
    // criteria works well, that should not be necessary!
  }

  const std::size_t check_edge(const IEdge& iedge) {

    std::vector<PFace> pfaces;
    std::size_t num_removed_pfaces = 0;
    incident_faces(iedge, pfaces);
    if (pfaces.size() == 1) {
      return remove_pfaces(iedge, pfaces[0], false);
    }
    if (pfaces.size() == 2) {
      const auto& pface0 = pfaces[0];
      const auto& pface1 = pfaces[1];
      if (pface0.first >= 6 && pface1.first >= 6 && pface0.first != pface1.first) {
        return remove_pfaces(iedge, pface0, false);
      }
    }
    return num_removed_pfaces;
  }

  const std::size_t remove_pfaces(
    const IEdge& init_iedge, const PFace& init_pface, const bool stop) {

    std::set<PFace> unique;
    std::vector< std::pair<Halfedge_index, PFace> > nfaces;
    const Halfedge_index init_he = find_crossing_he(init_iedge, init_pface);
    add_pfaces(init_he, init_pface, unique, nfaces);

    if (m_verbose) {
      std::cout << "* found faces to remove: " << nfaces.size() << std::endl;
    }

    std::size_t num_removed_pfaces = 0;
    for (const auto& item : nfaces) {
      const auto& he = item.first;
      const auto& nface = item.second;
      const bool success = remove_pface(he, nface);
      if (success) ++num_removed_pfaces;
    }
    CGAL_assertion(num_removed_pfaces == nfaces.size());
    if (stop) CGAL_assertion_msg(false, "TODO: DEBUG THIS FUNCTION!");
    return num_removed_pfaces;
  }

  const Halfedge_index find_crossing_he(
    const IEdge& iedge, const PFace& pface) {

    const auto& mesh = this->mesh(pface.first);
    const auto pedges = pedges_of_pface(pface);
    bool found_pedge = false;
    for (const auto pedge : pedges) {
      CGAL_assertion(has_iedge(pedge));
      if (this->iedge(pedge) == iedge) {
        found_pedge = true;

        const auto he = mesh.halfedge(pedge.second);
        const auto op = mesh.opposite(he);
        const auto face1 = mesh.face(he);
        const auto face2 = mesh.face(op);
        const bool has_face1 = (face1 != Support_plane::Mesh::null_face());
        const bool has_face2 = (face2 != Support_plane::Mesh::null_face());
        if (!has_face1) {
          return op;
        } else if (!has_face2) {
          return he;
        } else {
          CGAL_assertion_msg(false, "ERROR: CROSSING HE IS NOT FOUND!");
        }
      }
    }
    CGAL_assertion(found_pedge);
    return Halfedge_index();
  }

  void add_pfaces(
    const Halfedge_index crossing_he,
    const PFace& pface,
    std::set<PFace>& unique,
    std::vector< std::pair<Halfedge_index, PFace> >& nfaces) {

    const auto pair = unique.insert(pface);
    if (!pair.second) return;

    CGAL_assertion(crossing_he != Halfedge_index());
    CGAL_assertion(pface != null_pface());
    CGAL_assertion(pface.second != Support_plane::Mesh::null_face());
    nfaces.push_back(std::make_pair(crossing_he, pface));

    const auto& mesh = this->mesh(pface.first);
    const auto pedges = pedges_of_pface(pface);
    for (const auto pedge : pedges) {
      CGAL_assertion(has_iedge(pedge));

      const PVertex pvertex(pface.first, 0);
      bool is_occupied_edge, bbox_reached;
      std::tie(is_occupied_edge, bbox_reached) = is_occupied(pvertex, this->iedge(pedge));
      if (is_occupied_edge || bbox_reached) continue;

      const auto he = mesh.halfedge(pedge.second);
      const auto op = mesh.opposite(he);
      const auto face1 = mesh.face(he);
      const auto face2 = mesh.face(op);

      const auto nface1 = PFace(pface.first, face1);
      const auto nface2 = PFace(pface.first, face2);
      const bool has_nface1 = (face1 != Support_plane::Mesh::null_face());
      const bool has_nface2 = (face2 != Support_plane::Mesh::null_face());

      if (nface1 == pface) {
        if (has_nface2) {
          // std::cout << "adding nface2" << std::endl;
          add_pfaces(op, nface2, unique, nfaces);
        }
        continue;
      }
      if (nface2 == pface) {
        if (has_nface1) {
          // std::cout << "adding nface1" << std::endl;
          add_pfaces(he, nface1, unique, nfaces);
        }
        continue;
      }
      CGAL_assertion_msg(false, "ERROR: NO PFACE FOUND!");
    }
  }

  const bool remove_pface(const Halfedge_index he, const PFace& pface) {

    const std::string plane_idx = std::to_string(pface.first);
    const std::string face_idx  = std::to_string(pface.second);

    // std::cout << "removing " << str(pface) << std::endl;
    // dump_pface(*this, pface, "removed-pface-" + plane_idx + "-" + face_idx);

    auto& mesh = this->mesh(pface.first);
    CGAL::Euler::remove_face(he, mesh);
    return true;
  }

  /*******************************
  **    CHECKING PROPERTIES     **
  ********************************/

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
            dump_pedge(*this, pedge, "debug-pedge");
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
        std::cerr << "ERROR: CURRENT NUMBER OF EDGES = " << nedges.size() << std::endl;
      CGAL_assertion_msg(nedges.size() > 2,
      "ERROR: VERTEX MUST HAVE AT LEAST 3 NEIGHBORS!");
    }
  }

  void check_edges() {

    std::vector<PFace> nfaces;
    for (const auto edge : m_intersection_graph.edges()) {
      incident_faces(edge, nfaces);
      if (nfaces.size() == 1) {
        std::cerr << segment_3(edge) << std::endl;
        std::cerr << "ERROR: CURRENT NUMBER OF FACES = " << nfaces.size() << std::endl;
      }
      CGAL_assertion_msg(nfaces.size() != 1,
      "ERROR: EDGE MUST HAVE 0 OR AT LEAST 2 NEIGHBORS!");
    }
  }

  void check_faces() {

    for (KSR::size_t i = 0; i < number_of_support_planes(); ++i) {
      const auto pfaces = this->pfaces(i);
      for (const auto pface : pfaces) {
        const auto nvolumes = incident_volumes(pface);
        if (nvolumes.size() == 0 || nvolumes.size() > 2)
          std::cout << "ERROR: CURRENT NUMBER OF VOLUMES = " << nvolumes.size() << std::endl;
        CGAL_assertion_msg(nvolumes.size() == 1 || nvolumes.size() == 2,
        "ERROR: FACE MUST HAVE 1 OR 2 NEIGHBORS!");
      }
    }
  }

  const bool is_mesh_valid(
    const bool check_simplicity,
    const bool check_convexity,
    const bool check_equal_faces,
    const KSR::size_t support_plane_idx) const {

    const bool is_valid = mesh(support_plane_idx).is_valid();
    if (!is_valid) {
      return false;
    }

    // Note: bbox faces may have multiple equal points after converting from exact to inexact!
    if (support_plane_idx < 6) {
      return true;
    }

    const auto pfaces = this->pfaces(support_plane_idx);
    for (const auto pface : pfaces) {
      std::function<Point_2(PVertex)> unary_f =
      [&](const PVertex& pvertex) -> Point_2 {
        return point_2(pvertex);
      };

      const auto pvertices = pvertices_of_pface(pface);
      const Polygon_2 polygon(
        boost::make_transform_iterator(pvertices.begin(), unary_f),
        boost::make_transform_iterator(pvertices.end(), unary_f));

      // Very slow!
      if (check_equal_faces) {
        for (const auto oface : pfaces) {
          if (oface == pface) continue;

          const auto overtices = pvertices_of_pface(oface);
          if (overtices.size() != pvertices.size()) continue;
          const Polygon_2 oolygon(
            boost::make_transform_iterator(overtices.begin(), unary_f),
            boost::make_transform_iterator(overtices.end(), unary_f));

          std::set<std::size_t> unique;
          std::size_t num_overtices = 0;
          for (const auto& ppoint : polygon) {
            std::size_t count = 0;
            for (const auto& opoint : oolygon) {
              if (CGAL::squared_distance(ppoint, opoint) < KSR::tolerance<FT>()) {
                const auto res = unique.insert(count);
                const bool is_inserted = res.second;
                if (is_inserted) { ++num_overtices; }
                ++count; break;
              } else {
                ++count;
              }
            }
          }

          if (num_overtices == pvertices.size()) {

            std::cout << "pvertices: " << std::endl;
            for (const auto pvertex : pvertices) {
              std::cout << str(pvertex) << std::endl;
            }

            std::cout << "overtices: " << std::endl;
            for (const auto overtex : overtices) {
              std::cout << str(overtex) << std::endl;
            }

            dump_pface(*this, pface, "pface");
            dump_pface(*this, oface, "oface");
            dump_2d_surface_mesh(*this, support_plane_idx,
            "iter-10000-surface-mesh-" + std::to_string(support_plane_idx));
            const std::string msg = "ERROR: MESH " + std::to_string(support_plane_idx) +
            " HAS TWO EQUAL/OVERLAPPING PFACES: " + str(pface) + " AND " + str(oface) + "!";
            CGAL_assertion_msg(false, msg.c_str());
            return false;
          }
        }
      }

      // Use only with an exact kernel!
      if (check_simplicity && !polygon.is_simple()) {
        const std::string msg = "ERROR: PFACE " + str(pface) + " IS NOT SIMPLE!";
        CGAL_assertion_msg(false, msg.c_str());
        return false;
      }

      // Use only with an exact kernel!
      if (check_convexity && !polygon.is_convex()) {
        const std::string msg = "ERROR: PFACE " + str(pface) + " IS NOT CONVEX!";
        CGAL_assertion_msg(false, msg.c_str());
        return false;
      }

      auto prev = null_pvertex();
      for (const auto pvertex : pvertices) {
        if (prev == null_pvertex()) {
          prev = pvertex;
          continue;
        }

        if (point_2(prev) == point_2(pvertex) &&
          direction(prev) == direction(pvertex)) {

          const std::string msg = "ERROR: PFACE " + str(pface) +
          " HAS TWO CONSEQUENT IDENTICAL VERTICES "
          + str(prev) + " AND " + str(pvertex) + "!";
          CGAL_assertion_msg(false, msg.c_str());
          return false;
        }
        prev = pvertex;
      }
    }
    return true;
  }

  void check_integrity(
    const bool check_simplicity  = false,
    const bool check_convexity   = false,
    const bool check_equal_faces = true) const {

    for (KSR::size_t i = 0; i < number_of_support_planes(); ++i) {
      if (!is_mesh_valid(check_simplicity, check_convexity, check_equal_faces, i)) {
        const std::string msg = "ERROR: MESH " + std::to_string(i) + " IS NOT VALID!";
        CGAL_assertion_msg(false, msg.c_str());
      }

      for (const auto& iedge : this->iedges(i)) {
        const auto& iplanes = this->intersected_planes(iedge);
        if (iplanes.find(i) == iplanes.end()) {

          const std::string msg = "ERROR: SUPPORT PLANE " + std::to_string(i) +
          " IS INTERSECTED BY " + str(iedge) +
          " BUT IT CLAIMS IT DOES NOT INTERSECT IT!";
          CGAL_assertion_msg(false, msg.c_str());
        }
      }
    }

    for (const auto iedge : this->iedges()) {
      const auto& iplanes = this->intersected_planes(iedge);
      for (const auto support_plane_idx : iplanes) {

        const auto& sp_iedges = this->iedges(support_plane_idx);
        if (sp_iedges.find(iedge) == sp_iedges.end()) {

          const std::string msg = "ERROR: IEDGE " + str(iedge) +
          " INTERSECTS SUPPORT PLANE " + std::to_string(support_plane_idx) +
          " BUT IT CLAIMS IT IS NOT INTERSECTED BY IT!";
          CGAL_assertion_msg(false, msg.c_str());
        }
      }
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
      if (count != n) {
        std::cout << "current num neighbors " << count << " != " << n << std::endl;
        dump_info(*this, pface, *pedges.begin(), pfaces);
        return true;
      }
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

  /*******************************
  **     EXTRACTING VOLUMES     **
  ********************************/

  void create_polyhedrons() {

    std::cout.precision(20);
    // for (KSR::size_t i = 0; i < number_of_support_planes(); ++i)
    //   std::cout << "num pfaces sp " << i << ": " << pfaces(i).size() << std::endl;

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
    for (KSR::size_t i = 0; i < number_of_support_planes(); ++i) {
      const auto pfaces = this->pfaces(i);
      for (const auto pface : pfaces)
        map_volumes[pface] = std::make_pair(-1, -1);
    }

    // First, traverse only boundary volumes.
    // TODO: SORT HERE BY PFACE AREA!
    // Actually, we should sort by both number of edges and area!
    bool is_found_new_volume = false;
    std::size_t volume_size = 0;
    int num_volumes  = 0;
    int volume_index = 0;
    int volume_level = 0;
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
    if (m_verbose) {
      std::cout << "* found boundary volumes: "<< volume_index << std::endl;
    }
    num_volumes = volume_index;
    m_volume_level_map[volume_level] =
      static_cast<std::size_t>(num_volumes);
    ++volume_level;
    CGAL_assertion(num_volumes > 0);

    // Then traverse all other volumes if any.
    std::vector<PFace> other_pfaces;
    for (KSR::size_t i = 6; i < number_of_support_planes(); ++i) {
      const auto pfaces = this->pfaces(i);
      for (const auto pface : pfaces) {
        CGAL_assertion(pface.first >= 6);
        other_pfaces.push_back(pface);
      }
    }

    // TODO: SORT HERE BY PFACE AREA!
    // Actually, we should sort by both number of edges and area!
    std::sort(other_pfaces.begin(), other_pfaces.end(),
      [&](const PFace& pface1, const PFace& pface2) -> bool {
        const auto pedges1 = pedges_of_pface(pface1);
        const auto pedges2 = pedges_of_pface(pface2);
        return pedges1.size() > pedges2.size();
      }
    );

    bool quit = true;
    do {
      quit = true;
      const int before = volume_index;
      for (const auto& other_pface : other_pfaces) {
        std::tie(is_found_new_volume, volume_size) = traverse_interior_volume(
          other_pface, volume_index, num_volumes, map_volumes, centroids);
        if (is_found_new_volume) {
          quit = false;
          check_volume(volume_index, volume_size, map_volumes);
          ++volume_index;
        }
      }
      const int after = volume_index;
      if (m_verbose) {
        std::cout << "* found interior volumes: "<< after - before << std::endl;
      }
      CGAL_assertion(after >= before);
      num_volumes = volume_index;
      m_volume_level_map[volume_level] =
        static_cast<std::size_t>(num_volumes);
      ++volume_level;

    } while (!quit);
    m_num_volume_levels = volume_level;

    // Now, set final polyhedrons and their neighbors.
    for (const auto& item : map_volumes) {
      const auto& pface = item.first;
      const auto& pair  = item.second;

      if (pair.first == -1) {
        dump_pface(*this, pface, "face-debug");
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

    if (m_verbose) {
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
    }

    CGAL_assertion(m_volumes.size() == centroids.size());
    for (std::size_t i = 0; i < m_volumes.size(); ++i) {
      auto& volume = m_volumes[i];
      volume.set_index(i);
      volume.set_centroid(centroids.at(i));
    }
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

    if (m_verbose) {
      std::cout << "- FOUND VOLUME " << volume_index << ", (SIZE/BARYCENTER): "
      << volume_size << " / " << volume_centroid << std::endl;
    }
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
    if (m_verbose) {
      std::cout << "- FOUND VOLUME " << volume_index << ", (SIZE/BARYCENTER): "
      << volume_size << " / " << volume_centroid << std::endl;
    }
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

    Point_3 centroid = centroid_of_pface(pface);
    if (num_volumes > 0) {
      // std::cout << "SHIFTING CENTROID" << std::endl;

      CGAL_assertion(pair.first < num_volumes);
      CGAL_assertion(centroids.find(pair.first) != centroids.end());
      const auto& other_centroid = centroids.at(pair.first);
      const auto plane = plane_of_pface(pface);
      auto vec1 = plane.orthogonal_vector();
      vec1 = KSR::normalize(vec1);
      auto vec2 = Vector_3(centroid, other_centroid);
      vec2 = KSR::normalize(vec2);

      const FT d = FT(1) / FT(100000); // TODO: CAN WE AVOID THIS VALUE?
      const FT dot_product = vec1 * vec2;

      if (dot_product < FT(0)) {
        centroid += d * vec1;
      } else {
        centroid -= d * vec1;
      }
      volume_centroid = CGAL::barycenter(
        volume_centroid, static_cast<FT>(volume_size), centroid, FT(1));

    } else {
      volume_centroid = CGAL::barycenter(
        volume_centroid, static_cast<FT>(volume_size), centroid, FT(1));
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
      dump_pface(*this, pface, "bnd-pface-" +
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
        dump_info(*this, pface, pedge, nfaces);
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

    if (verbose) {
      std::cout << "pface: " << str(pface) << std::endl;
      std::cout << "pair: " <<
      std::to_string(pair.first) << "/" << std::to_string(pair.second) << std::endl;
    }

    const Point_3 centroid = centroid_of_pface(pface);
    volume_centroid = CGAL::barycenter(
      volume_centroid, static_cast<FT>(volume_size), centroid, FT(1));
    // std::cout << "volume centroid: " << volume_centroid << std::endl;
    ++volume_size;

    if (verbose) {
      // std::cout << "INT PFACE MAP: (" <<
      // pair.first << ", " << pair.second << ")" << std::endl;
      std::cout << "DUMPING INT PFACE: " <<
        std::to_string(volume_index) + "-" +
        std::to_string(pface.first) + "-" +
        std::to_string(pface.second) << std::endl;
      dump_pface(*this, pface, "int-pface-" +
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
        dump_info(*this, pface, pedge, nfaces);
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
      dump_info(*this, pface, pedge, nfaces);
    }
    CGAL_assertion(nfaces.size() > 1);

    Point_3 center = centroid_of_pface(pface);
    const Segment_3 segment = segment_3(pedge);
    const Line_3 line(segment.source(), segment.target());
    Point_3 midp = CGAL::midpoint(segment.source(), segment.target());
    // std::cout << "midp: " << midp << std::endl;
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

    // std::cout << "1 d: " << d << std::endl;
    // std::cout << "1 norm: " << norm << std::endl;
    const Plane_3 tr_plane(midp + norm * d, norm);
    Point_3 inter;
    const bool is_intersection_found = KSR::intersection(line, tr_plane, inter);
    if (!is_intersection_found) {
      std::cout << "d = " << d << std::endl;
    }
    CGAL_assertion(is_intersection_found);
    // std::cout << "inter: " << inter << std::endl;

    d = KSR::distance(midp, inter);
    norm = Vector_3(midp, inter);
    // std::cout << "2 d: " << d << std::endl;
    // std::cout << "2 norm: " << norm << std::endl;

    if (d != FT(0)) {
      CGAL_assertion(norm != Vector_3(FT(0), FT(0), FT(0)));
      norm = KSR::normalize(norm);
      for (auto& point : points) {
        point += norm * d;
      }
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
          dump_pface(*this, dir_edges[im].second, "prev");
          dump_pface(*this, dir_edges[ip].second, "next");
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
          dump_info(*this, pface, pedge, nfaces);
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

  void create_cell_pvertices(Volume_cell& cell) {
    cell.pvertices.clear();
    for (const auto& pface : cell.pfaces) {
      for (const auto pvertex : pvertices_of_pface(pface)) {
        cell.pvertices.insert(pvertex);
      }
    }
  }

private:

  /*******************************
  **   FUTURE POINTS AND DIRS   **
  ********************************/

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

    // std::cout << "prev slope: " << m1 << std::endl;
    // std::cout << "next slope: " << m2 << std::endl;
    // std::cout << "iedg slope: " << m3 << std::endl;

    if (CGAL::abs(m1 - m3) < tol) {
      if (m_verbose) std::cout << "- prev parallel lines" << std::endl;
      const FT prev_dot = current_vec_prev * iedge_vec;
      if (prev_dot < FT(0)) {
        if (m_verbose) std::cout << "- prev moves backwards" << std::endl;
        future_point_a = target_p;
      } else {
        if (m_verbose) std::cout << "- prev moves forwards" << std::endl;
        future_point_a = source_p;
      }
    } else {

      if (m_verbose) std::cout << "- prev intersected lines" << std::endl;
      const bool a_found = KSR::intersection(future_line_prev, iedge_line, future_point_a);
      if (!a_found) {
        std::cerr << "WARNING: A IS NOT FOUND!" << std::endl;
        future_point_b = pinit + (pinit - future_point_a);
      }
    }

    direction_a = Vector_2(pinit, future_point_a);
    future_point_a = pinit - m_current_time * direction_a;
    if (m_verbose) {
      std::cout << "- prev future point a: " <<
      to_3d(pvertex.first, future_point_a + m_current_time * direction_a) << std::endl;
      std::cout << "- prev future direction a: " << direction_a << std::endl;
    }

    if (CGAL::abs(m2 - m3) < tol) {
      if (m_verbose) std::cout << "- next parallel lines" << std::endl;
      const FT next_dot = current_vec_next * iedge_vec;
      if (next_dot < FT(0)) {
        if (m_verbose) std::cout << "- next moves backwards" << std::endl;
        future_point_b = target_p;
      } else {
        if (m_verbose) std::cout << "- next moves forwards" << std::endl;
        future_point_b = source_p;
      }

    } else {

      if (m_verbose) std::cout << "- next intersected lines" << std::endl;
      const bool b_found = KSR::intersection(future_line_next, iedge_line, future_point_b);
      if (!b_found) {
        std::cerr << "WARNING: B IS NOT FOUND!" << std::endl;
        future_point_a = pinit + (pinit - future_point_b);
      }
    }

    direction_b = Vector_2(pinit, future_point_b);
    future_point_b = pinit - m_current_time * direction_b;
    if (m_verbose) {
      std::cout << "- next future point b: " <<
      to_3d(pvertex.first, future_point_b + m_current_time * direction_b) << std::endl;
      std::cout << "- next furure direction b: " << direction_b << std::endl;
    }
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
    //   std::cout << "found limit" << std::endl;
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
      if (m_verbose) std::cout << "- back/front parallel lines" << std::endl;

      is_parallel = true;
      const FT next_dot = current_vec_next * iedge_vec;
      if (next_dot < FT(0)) {
        if (m_verbose) std::cout << "- back/front moves backwards" << std::endl;
        future_point = target_p;
        // std::cout << point_3(target(iedge)) << std::endl;
      } else {
        if (m_verbose) std::cout << "- back/front moves forwards" << std::endl;
        future_point = source_p;
        // std::cout << point_3(source(iedge)) << std::endl;
      }

    } else {
      if (m_verbose) std::cout << "- back/front intersected lines" << std::endl;
      future_point = KSR::intersection<Point_2>(future_line_next, iedge_line);
    }

    // std::cout << "prev: " << point_3(next, m_current_time + FT(1)) << std::endl;
    // std::cout << "back: " << point_3(curr, m_current_time + FT(1)) << std::endl;

    future_direction = Vector_2(pinit, future_point);
    future_point = pinit - m_current_time * future_direction;

    // auto tmp = future_direction;
    // tmp = KSR::normalize(tmp);
    // std::cout << "future tmp: " << to_3d(pvertex.first, pinit + m_current_time * tmp) << std::endl;

    if (m_verbose) {
      std::cout << "- back/front future point: " <<
      to_3d(pvertex.first, future_point + m_current_time * future_direction) << std::endl;
      std::cout << "- back/front furure direction: " << future_direction << std::endl;
    }
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
      if (m_verbose) std::cout << "- open parallel lines" << std::endl;

      is_parallel = true;
      if (source_p == pv_point)
        future_point = target_p;
      else
        future_point = source_p;

    } else {
      if (m_verbose) std::cout << "- open intersected lines" << std::endl;
      future_point = KSR::intersection<Point_2>(future_line_next, iedge_line);
    }

    future_direction = Vector_2(pinit, future_point);
    future_point = pinit - m_current_time * future_direction;

    // auto tmp = future_direction;
    // tmp = KSR::normalize(tmp);
    // std::cout << "future tmp: " << to_3d(pvertex.first, pinit + m_current_time * tmp) << std::endl;

    if (m_verbose) {
      std::cout << "- open future point: " <<
      to_3d(pvertex.first, future_point + m_current_time * future_direction) << std::endl;
      std::cout << "- open furure direction: " << future_direction << std::endl;
    }
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
      if (m_verbose) std::cout << "* constrained pvertex case" << std::endl;
      return false;
    }

    if (!is_active(pvertex)) {
      if (m_verbose) std::cout << "* pvertex no active case" << std::endl;
      return false;
    }

    if (!is_active(iedge)) {
      if (m_verbose) std::cout << "* iedge no active case" << std::endl;
      return false;
    }

    if (!CGAL::do_overlap(pv_bbox, iedge_bbox)) {
      if (m_verbose) std::cout << "* no overlap case" << std::endl;
      return false;
    }

    Point_2 point;
    if (!KSR::intersection(pv_seg, iedge_seg, point)) {
      if (m_verbose) std::cout << "* no intersection case" << std::endl;
      return false;
    }

    if (m_verbose) std::cout << "* found intersection" << std::endl;
    return true;
  }
};

} // namespace KSR_3
} // namespace CGAL

#endif // CGAL_KSR_3_DATA_STRUCTURE_H
