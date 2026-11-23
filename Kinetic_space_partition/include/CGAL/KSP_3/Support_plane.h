// Copyright (c) 2023 GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Sven Oesau, Florent Lafarge, Dmitry Anisimov, Simon Giraudot

#ifndef CGAL_KSP_3_SUPPORT_PLANE_H
#define CGAL_KSP_3_SUPPORT_PLANE_H

#include <CGAL/license/Kinetic_space_partition.h>

// CGAL includes.
#include <CGAL/Surface_mesh.h>
#include <CGAL/centroid.h>

// Internal includes.
#include <CGAL/KSP/utils.h>
#include <CGAL/KSP_3/Intersection_graph.h>

namespace CGAL {
namespace KSP_3 {
namespace internal {

#ifdef DOXYGEN_RUNNING
#else

template<typename GeomTraits, typename IntersectionKernel>
class Support_plane {

public:
  using Kernel = GeomTraits;
  using Intersection_kernel = IntersectionKernel;
  using To_exact = CGAL::Cartesian_converter<Kernel, Intersection_kernel>;
  using From_exact = CGAL::Cartesian_converter<Intersection_kernel, Kernel>;

  using FT = typename Kernel::FT;
  using Point_2 = typename Kernel::Point_2;
  using Point_3 = typename Kernel::Point_3;
  using Vector_2 = typename Kernel::Vector_2;
  using Vector_3 = typename Kernel::Vector_3;
  using Direction_2 = typename Kernel::Direction_2;
  using Segment_2 = typename Kernel::Segment_2;
  using Segment_3 = typename Kernel::Segment_3;
  using Line_2 = typename Kernel::Line_2;
  using Line_3 = typename Kernel::Line_3;
  using Plane_3 = typename Kernel::Plane_3;
  using Triangle_2 = typename Kernel::Triangle_2;
  using IkFT = typename Intersection_kernel::FT;
  using IkVector_2 = typename Intersection_kernel::Vector_2;
  using IkLine_2 = typename Intersection_kernel::Line_2;
  using IkPoint_2 = typename Intersection_kernel::Point_2;
  using IkSegment_2 = typename Intersection_kernel::Segment_2;

  using Mesh = CGAL::Surface_mesh<Point_2>;
  using Intersection_graph = CGAL::KSP_3::internal::Intersection_graph<Kernel, Intersection_kernel>;
  using Bbox_2 = CGAL::Bbox_2;

  using IVertex = typename Intersection_graph::Vertex_descriptor;
  using IEdge = typename Intersection_graph::Edge_descriptor;
  using IFace = typename Intersection_graph::Face_descriptor;

  using IEdge_set = typename Intersection_graph::IEdge_set;

  using Vertex_index = typename Mesh::Vertex_index;
  using Face_index = typename Mesh::Face_index;
  using Edge_index = typename Mesh::Edge_index;
  using Halfedge_index = typename Mesh::Halfedge_index;

  using V_ivertex_map = typename Mesh::template Property_map<Vertex_index, IVertex>;
  using V_iedge_map = typename Mesh::template Property_map<Vertex_index, IEdge>;
  using E_iedge_map = typename Mesh::template Property_map<Edge_index, IEdge>;
  using F_index_map = typename Mesh::template Property_map<Face_index, std::vector<std::size_t> >;
  using F_bool_map = typename Mesh::template Property_map<Face_index, bool>;

  struct Cached_event {
    IkFT t;
    IkFT u;
    IkPoint_2 p;
    IEdge e = Intersection_graph::null_iedge();// for unconstrained vertices
    IVertex d  = Intersection_graph::null_ivertex();// for constrained vertices
  };

  struct Vertex_event {
    Vertex_event() : vertex(-1), other(-1), time(-1), crossed_edge(Intersection_graph::null_iedge()), destination(Intersection_graph::null_ivertex()) {}
    Vertex_event(Cached_event& ce, std::size_t v) : vertex(v), other(-1), time(ce.t), u(ce.u), p(ce.p), crossed_edge(ce.e), destination(ce.d) {}
    std::size_t vertex, other;
    IkFT time;
    IkFT u;
    IkPoint_2 p;
    IEdge crossed_edge;
    IVertex destination; // Target ivertex for constrained vertices.
  };

  struct Vertex {
    Vertex() : sp_idx(std::size_t(-1)), queued_events(0) {}
    Vertex(std::size_t sp_idx, const IkPoint_2& p0, const IkVector_2& v, IkFT t_init = 0)
      : sp_idx(sp_idx), p0(p0), v(v), t_init(t_init), ivertex(-1), itarget(-1), moving(true), k(0), face(std::size_t(-1)), other(-1), constraints(),
      constraint_edge(IVertex(-1), IVertex(-1), nullptr), other_constraint_edge(IVertex(-1), IVertex(-1), nullptr), queued_events(0) {
    }

    void stop(const Vertex_event &event) {
      moving = false;
      t_init = event.time;
      v = IkVector_2(0, 0);
      p0 = event.p;
      ivertex = itarget;

      known_intersections.clear();
      cached_events.clear();
      queued_events = 0;
    }

    std::size_t sp_idx;
    IkPoint_2 p0;
    IkVector_2 v;
    IkFT t_init;
    IVertex ivertex, itarget;
    bool moving;
    int k;
    std::size_t face;
    std::size_t other; // If the vertex is constrained, it may have a mirrored vertex in the neighbor face.
    std::set<std::size_t> constraints; // line indices
    IEdge constraint_edge, other_constraint_edge;
    // Reference from edge back to moving vertices is required.
    std::unordered_map<std::size_t, Cached_event> known_intersections; // Cache to avoid recalculating intersections.
    std::list<Cached_event> cached_events;
    int queued_events;
  };

  // A polygon expanding with a homothety about its centroid, i.e., the vertex p moves with the velocity p - centroid.
  struct Input_polygon {
    std::vector<typename Intersection_kernel::Point_2> vertices; // Counterclockwise around the centroid.
    std::vector<typename Intersection_kernel::Vector_2> velocities;
    typename Intersection_kernel::Point_2 centroid;
    std::size_t input; // Index of the input polygon. For merged polygons, the index of one of the merged input polygons.
  };

  struct Data {
    Data(std::vector<Vertex> &vertices)
      : mesh(),
      v_ivertex_map(mesh.template add_property_map<Vertex_index, IVertex>("v:ivertex", Intersection_graph::null_ivertex()).first),
      v_iedge_map(mesh.template add_property_map<Vertex_index, IEdge>("v:iedge", Intersection_graph::null_iedge()).first),
      e_iedge_map(mesh.template add_property_map<Edge_index, IEdge>("e:iedge", Intersection_graph::null_iedge()).first),
      input_map(mesh.template add_property_map<Face_index, std::vector<std::size_t> >("f:input", std::vector<std::size_t>()).first),
      f_initial_map(mesh.template add_property_map<Face_index, bool >("f:initial", false).first),
      vertices(vertices) {}

    bool is_bbox;
    Point_2 centroid;
    IkPoint_2 ikcentroid;
    Plane_3 plane;
    typename Intersection_kernel::Plane_3 exact_plane;
    Mesh mesh;

    V_ivertex_map v_ivertex_map;
    V_iedge_map v_iedge_map;
    E_iedge_map e_iedge_map;
    F_index_map input_map;
    F_bool_map f_initial_map;
    std::map<IEdge, std::pair<IFace, IFace> > iedge2ifaces;
    // Ordered containers: the iteration order of lines and polygons determines the order in which polygons are split
    // and events are queued. An unordered container would make the result depend on the standard library.
    std::map<std::size_t, typename Intersection_kernel::Line_2> lines, bbox_lines;
    std::set<IFace> ifaces; // All ifaces in the support plane
    std::map<IVertex, Vertex_index> ivertex2pvertex;
    IEdge_set iedges;

    std::vector<Point_2> original_vertices;
    std::vector<typename Intersection_kernel::Point_2> exact_vertices;
    std::vector<typename Intersection_kernel::Direction_2> exact_directions;
    std::map<IFace, std::list<std::size_t> > polygons;
    // Input polygons of the support plane. Coplanar input polygons are propagated separately. Polygons sharing an iface
    // at t = 0 are merged within the iface into the convex hull of their pieces, which is added as an additional polygon.
    std::vector<Input_polygon> input_polygons;
    std::map<IFace, std::size_t> piece_polygon; // Index into input_polygons of the polygon each piece in polygons originates from.
    // Input polygons lying in a bbox support plane: index of the input polygon of the partition and vertices.
    std::vector<std::pair<std::size_t, std::vector<typename Intersection_kernel::Point_2> > > coinciding_polygons;
    std::unordered_map<IFace, Face_index> iface2pface;
    std::unordered_map<IFace, std::list<IEdge>> borders;
    std::vector<Vertex> &vertices;


    std::size_t actual_input_polygon;

    int k;
  };
public:

private:

  std::shared_ptr<Data> m_data;

public:
  //Support_plane() : m_data(std::make_shared<Data>()) {}

  template<typename PointRange>
  Support_plane(const PointRange& polygon, const bool is_bbox, typename Intersection_kernel::Plane_3 plane, std::vector<Vertex> &vertices)
    : m_data(std::make_shared<Data>(vertices)) {

    From_exact from_exact;

    m_data->k = 0;
    m_data->plane = from_exact(plane);
    m_data->exact_plane = plane;
    m_data->is_bbox = is_bbox;
    m_data->actual_input_polygon = static_cast<std::size_t>(- 1);

    std::vector<IkPoint_2> points;
    points.reserve(polygon.size());
    for (const auto& point : polygon) {
      points.push_back(to_2d(point));
    }
    CGAL_assertion(points.size() == polygon.size());

    m_data->ikcentroid = CGAL::centroid(points.begin(), points.end(), CGAL::Dimension_tag<0>());
    m_data->centroid = from_exact(m_data->ikcentroid);

    add_property_maps();
  }

  Support_plane(const std::vector<typename Intersection_kernel::Point_3>& polygon, const bool is_bbox, std::vector<Vertex>& vertices)
    : m_data(std::make_shared<Data>(vertices)) {

    From_exact from_exact;

    m_data->k = 0;
    m_data->exact_plane = typename Intersection_kernel::Plane_3(polygon[0], polygon[1], polygon[2]);
    m_data->plane = from_exact(m_data->exact_plane);
    m_data->is_bbox = is_bbox;
    m_data->actual_input_polygon = static_cast<std::size_t>(- 1);

    std::vector<IkPoint_2> points;
    points.reserve(polygon.size());
    for (const auto& point : polygon) {
      points.push_back(to_2d(point));
    }
    CGAL_assertion_code(const std::size_t n = points.size();)
    CGAL_assertion(n == polygon.size());
    CGAL_assertion(n != 0);

    m_data->ikcentroid = CGAL::centroid(points.begin(), points.end(), CGAL::Dimension_tag<0>());
    m_data->centroid = from_exact(m_data->ikcentroid);

    add_property_maps();
  }

  void add_property_maps() {
    m_data->v_ivertex_map = m_data->mesh.template add_property_map<Vertex_index, IVertex>("v:ivertex", Intersection_graph::null_ivertex()).first;
    m_data->v_iedge_map = m_data->mesh.template add_property_map<Vertex_index, IEdge>("v:iedge", Intersection_graph::null_iedge()).first;
    m_data->e_iedge_map = m_data->mesh.template add_property_map<Edge_index, IEdge>("e:iedge", Intersection_graph::null_iedge()).first;
    m_data->input_map = m_data->mesh.template add_property_map<Face_index, std::vector<std::size_t> >("f:input", std::vector<std::size_t>()).first;
    m_data->f_initial_map = m_data->mesh.template add_property_map<Face_index, bool >("f:initial", false).first;
  }

  // Returns a copy of the support plane that does not share its data with this one.
  // The copy still refers to the same kinetic vertices and the same intersection graph.
  Support_plane deep_copy() const {
    Support_plane copy(*this);
    copy.m_data = std::make_shared<Data>(*m_data);
    copy.link_property_maps();
    return copy;
  }

  void link_property_maps() {
    m_data->v_ivertex_map = m_data->mesh.template property_map<Vertex_index, IVertex>("v:ivertex").value();
    m_data->v_iedge_map = m_data->mesh.template property_map<Vertex_index, IEdge>("v:iedge").value();
    m_data->e_iedge_map = m_data->mesh.template property_map<Edge_index, IEdge>("e:iedge").value();
    m_data->input_map = m_data->mesh.template property_map<Face_index, std::vector<std::size_t> >("f:input").value();
    m_data->f_initial_map = m_data->mesh.template property_map<Face_index, bool >("f:initial").value();
  }

  void centroid(Point_2& c) {
    if (m_data->original_vertices.size() < 2)
      return;
    std::vector<Triangle_2> tris(m_data->original_vertices.size() - 2);

    for (std::size_t i = 2; i < m_data->original_vertices.size(); i++) {
      tris[i - 2] = Triangle_2(m_data->original_vertices[0], m_data->original_vertices[i - 1], m_data->original_vertices[i]);
    }

    c = CGAL::centroid(tris.begin(), tris.end(), CGAL::Dimension_tag<2>());
  }

  Data& data() { return *m_data; }

  void clear_pfaces() {
    m_data->mesh.clear();
    add_property_maps();
  }

  static void check_edge_is_turning(const Vertex& v1, const Vertex& v2) {
    IkFT t_min = ((v1.t_init < v2.t_init) ? v2.t_init : v1.t_init) + 1;
    IkPoint_2 v1p = v1.p0 + (t_min - v1.t_init) * v1.v;
    IkPoint_2 v2p = v2.p0 + (t_min - v2.t_init) * v2.v;
    IkPoint_2 v1_t = v1p + v1.v;
    IkPoint_2 v2_t = v2p + v2.v;
    IkVector_2 e1 = v2p - v1p;
    IkVector_2 e2 = v2_t - v1_t;
    CGAL_assertion((e1.x() * e2.y() - e1.y() * e2.x()) == 0); // Collinear
  }

  static void check_edge_is_turning(const Vertex& v1, const Vertex& v2, IkVector_2& v) {
    IkFT t_min = ((v1.t_init < v2.t_init) ? v2.t_init : v1.t_init) + 1;
    IkPoint_2 v1p = v1.p0 + (t_min - v1.t_init) * v1.v;
    IkPoint_2 v2p = v2.p0 + (t_min - v2.t_init) * v2.v;
    IkPoint_2 v1_t = v1p + v1.v;
    IkPoint_2 v2_t = v2p + v2.v;
    IkPoint_2 v1_2t = v1p + v1.v + v1.v;
    IkPoint_2 v2_2t = v2p + v2.v + v2.v;
    IkVector_2 e1 = v2p - v1p;
    IkVector_2 e2 = v2_t - v1_t;
    IkVector_2 e3 = v2_2t - v1_2t;
    CGAL_assertion((e1.x() * e2.y() - e1.y() * e2.x()) == 0); // Collinear
    CGAL_assertion((e3.x() * e2.y() - e3.y() * e2.x()) == 0); // Collinear
    v = e1;
  }

  static IkVector_2 calculate_edge_speed(const Vertex& a, const Vertex& b, const IkPoint_2& intersection, const IkLine_2& line) {
    CGAL_assertion_code(check_edge_is_turning(a, b));
    IkPoint_2 a_t = a.p0 + a.v;
    IkPoint_2 b_t = b.p0 + b.v;

    bool half = (a_t == b_t);
    if (half) {
      a_t = a.p0 + a.v / FT(2);
      b_t = b.p0 + b.v / FT(2);
    }

    IkPoint_2 i_t;
    auto res = CGAL::intersection(line, IkLine_2(a_t, b_t));
    CGAL::assign(i_t, res);

    IkVector_2 edge_vec = i_t - intersection;
    if (half)
      edge_vec = edge_vec * FT(2);

    return edge_vec;
  }

  static IkVector_2 calculate_edge_speed(const Vertex& a, const Vertex& b, const IkPoint_2& intersection, const IkLine_2& line, const IkFT &t) {
    CGAL_assertion_code(check_edge_is_turning(a, b));
    IkPoint_2 a_t = a.p0 + (1 + t - a.t_init) * a.v;
    IkPoint_2 b_t = b.p0 + (1 + t - b.t_init) * b.v;

    bool half = (a_t == b_t);
    if (half) {
      a_t = a.p0 + (0.5 + t - a.t_init) * a.v;
      b_t = b.p0 + (0.5 + t - b.t_init) * b.v;
    }

    IkPoint_2 i_t;
    auto res = CGAL::intersection(line, IkLine_2(a_t, b_t));
    CGAL::assign(i_t, res);

    IkVector_2 edge_vec = i_t - intersection;
    if (half)
      edge_vec = edge_vec * FT(2);

    return edge_vec;
  }

  const std::array<Vertex_index, 4>
    add_bbox_polygon(
      const std::array<Point_2, 4>& points,
      const std::array<IVertex, 4>& ivertices) {

    CGAL_assertion(CGAL::is_simple_2(points.begin(), points.end()));
    CGAL_assertion(CGAL::is_convex_2(points.begin(), points.end()));

    std::array<Vertex_index, 4> vertices;
    for (std::size_t i = 0; i < 4; ++i) {
      const auto vi = m_data->mesh.add_vertex(points[i]);
      m_data->v_ivertex_map[vi] = ivertices[i];
      vertices[i] = vi;
    }

    const auto fi = m_data->mesh.add_face(vertices);
    CGAL_assertion(fi != Mesh::null_face());
    auto& input_vec = m_data->input_map[fi];
    CGAL_assertion(input_vec.empty());
    input_vec.push_back(std::size_t(-1));
    return vertices;
  }

  // Adds an input polygon to the support plane and returns its index in input_polygons().
  // The first input polygon also defines the polygon data of the support plane (exact_vertices, centroid, ...).
  template<typename Pair>
  std::size_t add_input_polygon(
    const std::vector<Pair>& points, std::size_t input) {

    std::vector<typename Intersection_kernel::Point_2> pts;
    pts.reserve(points.size());
    for (const auto& pair : points)
      pts.push_back(pair.first);
    const std::size_t idx = add_polygon(pts, input);
    if (idx == 0)
      set_polygon_data(points);
    return idx;
  }

  // Adds a polygon expanding with a homothety about the centroid of its vertices and returns its index in input_polygons().
  // The polygon has to be convex.
  std::size_t add_polygon(const std::vector<typename Intersection_kernel::Point_2>& points, std::size_t input) {
    using IkVector_2 = typename Intersection_kernel::Vector_2;
    CGAL_assertion(points.size() >= 3);
    Input_polygon poly;
    poly.input = input;
    poly.centroid = CGAL::centroid(points.begin(), points.end(), CGAL::Dimension_tag<0>());

    std::vector<std::pair<std::size_t, typename Intersection_kernel::Direction_2> > dir_vec;
    dir_vec.reserve(points.size());
    for (std::size_t i = 0; i < points.size(); i++)
      dir_vec.push_back(std::make_pair(i, typename Intersection_kernel::Direction_2(IkVector_2(poly.centroid, points[i]))));

    std::sort(dir_vec.begin(), dir_vec.end(), [](const auto& a, const auto& b) { return a.second < b.second; });

    for (const auto& d : dir_vec) {
      poly.vertices.push_back(points[d.first]);
      poly.velocities.push_back(IkVector_2(poly.centroid, points[d.first]));
    }

    m_data->input_polygons.push_back(std::move(poly));
    return m_data->input_polygons.size() - 1;
  }

  std::size_t number_of_input_polygons() const { return m_data->input_polygons.size(); }
  const Input_polygon& input_polygon(std::size_t idx) const { return m_data->input_polygons[idx]; }

  // Index of the input polygon the piece in the iface originates from.
  std::size_t piece_polygon(IFace face) const {
    auto it = m_data->piece_polygon.find(face);
    CGAL_assertion(it != m_data->piece_polygon.end());
    return (it == m_data->piece_polygon.end()) ? 0 : it->second;
  }

  void set_piece_polygon(IFace face, std::size_t polygon) { m_data->piece_polygon[face] = polygon; }

  template<typename Pair>
  void set_polygon_data(
    const std::vector<Pair>& points) {

    CGAL_assertion(is_simple_polygon(points));
    CGAL_assertion(is_valid_polygon(points));

    From_exact from_exact;

    CGAL_assertion(points.size() >= 3);

    std::vector<Vertex_index> vertices;
    const std::size_t n = points.size();
    CGAL_assertion(n >= 3);
    vertices.reserve(n);
    m_data->original_vertices.resize(n);
    m_data->exact_vertices.resize(n);
    m_data->exact_directions.resize(n);

    for (std::size_t i = 0; i < n; ++i) {
      m_data->exact_vertices[i] = points[i].first;
    }

    m_data->ikcentroid = CGAL::centroid(m_data->exact_vertices.begin(), m_data->exact_vertices.end(), CGAL::Dimension_tag<0>());
    m_data->centroid = from_exact(m_data->ikcentroid);

    std::vector<typename Intersection_kernel::Vector_2> directions;
    directions.reserve(n);

    std::vector<std::pair<std::size_t, typename Intersection_kernel::Direction_2> > dir_vec;

    for (const auto& pair : points) {
      const auto& point = pair.first;
      directions.push_back(typename Intersection_kernel::Vector_2(m_data->ikcentroid, point));
    }
    CGAL_assertion(directions.size() == n);

    dir_vec.reserve(n);
    for (std::size_t i = 0; i < n; i++)
      dir_vec.push_back(std::pair<std::size_t, typename Intersection_kernel::Direction_2>(i, directions[i]));

    std::sort(dir_vec.begin(), dir_vec.end(),
      [&](const std::pair<std::size_t, typename Intersection_kernel::Direction_2>& a,
        const std::pair<std::size_t, typename Intersection_kernel::Direction_2>& b) -> bool {
          return a.second < b.second;
      });

    for (std::size_t i = 0; i < n; ++i) {
      const auto& point = points[dir_vec[i].first].first;
      m_data->original_vertices[i] = from_exact(point);
      m_data->exact_vertices[i] = point;
      m_data->exact_directions[i] = dir_vec[i].second;
    }

    // The exact directions are sorted. The approximated directions may coincide for small polygons.
    CGAL_assertion_code(
    for (std::size_t i = 0; i < m_data->exact_directions.size(); i++) {
      for (std::size_t j = 0; j < m_data->exact_directions.size(); j++) {
        if (j < i)
          CGAL_assertion(m_data->exact_directions[j] < m_data->exact_directions[i]);
        if (j > i)
          CGAL_assertion(m_data->exact_directions[i] < m_data->exact_directions[j]);
      }
    });
  }

  void set_input_polygon(std::size_t input_polygon_idx) {
    m_data->actual_input_polygon = input_polygon_idx;
  }

  template<typename Pair>
  bool is_valid_polygon(const std::vector<Pair>& polygon) const {
    for (std::size_t i = 0; i < polygon.size(); ++i) {
      const std::size_t ip = (i + 1) % polygon.size();
      const auto& p = polygon[i].first;
      const auto& q = polygon[ip].first;
      const bool is_equal_zero = (p == q);
      CGAL_assertion_msg(!is_equal_zero,
        "ERROR: WE HAVE EQUAL POINTS IN THE INPUT POLYGON!");
      if (is_equal_zero) return false;
    }
    return true;
  }

  template<typename Pair>
  bool is_simple_polygon(const std::vector<Pair>& points) const {
    // Decided with the exact points, as a small polygon may not be simple after rounding.
    std::vector<typename Intersection_kernel::Point_2> polygon;
    polygon.reserve(points.size());
    for (const auto& pair : points)
      polygon.push_back(pair.first);
    return CGAL::is_simple_2(polygon.begin(), polygon.end(), Intersection_kernel());
  }

  const Plane_3& plane() const { return m_data->plane; }
  const typename Intersection_kernel::Plane_3& exact_plane() const { return m_data->exact_plane; }
  const Point_2& centroid() const { return m_data->centroid; }
  bool is_bbox() const { return m_data->is_bbox; }
  std::map<IVertex, Vertex_index>& ivertex2pvertex() { return m_data->ivertex2pvertex; }

  const Mesh& mesh() const { return m_data->mesh; }
  Mesh& mesh() { return m_data->mesh; }

  void set_point(const Vertex_index& vi, const Point_2& point) {
    m_data->mesh.point(vi) = point;
  }

  void add_neighbor(IEdge edge, IFace face) {
    std::pair<IEdge, std::pair<IFace, IFace>> neighbor(edge, std::pair<IFace, IFace>(face, Intersection_graph::null_iface()));
    auto pair = m_data->iedge2ifaces.insert(neighbor);
    m_data->ifaces.insert(face);
    if (!pair.second) {
      CGAL_assertion(pair.first->second.first != Intersection_graph::null_iface());
      CGAL_assertion(pair.first->second.second == Intersection_graph::null_iface());
      pair.first->second.second = face;
    }
  }

  IFace iface(IEdge edge) {
    auto it = m_data->iedge2ifaces.find(edge);
    if (it == m_data->iedge2ifaces.end())
      return Intersection_graph::null_iface();
    else return it->second.first;
  }

  IFace other(IEdge edge, IFace face) {
    auto it = m_data->iedge2ifaces.find(edge);
    if (it == m_data->iedge2ifaces.end())
      return Intersection_graph::null_iface();
    if (it->second.first == face)
      return it->second.second;
    else
      return it->second.first;
  }

  const Vertex_index prev(const Vertex_index& vi) const {
    return m_data->mesh.source(m_data->mesh.halfedge(vi));
  }
  const Vertex_index next(const Vertex_index& vi) const {
    return m_data->mesh.target(m_data->mesh.next(m_data->mesh.halfedge(vi)));
  }

  const Face_index face(const Vertex_index& vi) const {

    auto out = m_data->mesh.face(m_data->mesh.halfedge(vi));
    if (out == Face_index()) {
      out = m_data->mesh.face(m_data->mesh.opposite(m_data->mesh.halfedge(vi)));
    }
    CGAL_assertion(out != Face_index());
    return out;
  }

  void add_face(IFace &face, Face_index fi) {
    m_data->iface2pface[face] = fi;
  }

  const Face_index face(IFace &face) {
    auto it = m_data->iface2pface.find(face);
    if (it == m_data->iface2pface.end())
      return Face_index();
    else
      return it->second;
  }

  const std::pair<Face_index, Face_index> faces(const Vertex_index& vi) const {

    for (const auto& he : halfedges_around_target(halfedge(vi, m_data->mesh), m_data->mesh)) {
      if (has_iedge(m_data->mesh.edge(he))) {
        return std::make_pair(
          m_data->mesh.face(he), m_data->mesh.face(m_data->mesh.opposite(he)));
      }
    }
    CGAL_assertion_msg(false, "ERROR: NO CONSTRAINED EDGE FOUND!");
    return std::make_pair(Face_index(), Face_index());
  }

  const std::pair<Face_index, Face_index> faces(const Halfedge_index& he) const {

    if (has_iedge(m_data->mesh.edge(he))) {
      return std::make_pair(
        m_data->mesh.face(he), m_data->mesh.face(m_data->mesh.opposite(he)));
    }
    CGAL_assertion_msg(false, "ERROR: NO CONSTRAINED EDGE FOUND!");
    return std::make_pair(Face_index(), Face_index());
  }

  const Point_2 point_2(const Vertex_index& vi) const {
    return m_data->mesh.point(vi);
  }

  const Point_3 point_3(const Vertex_index& vi) const {
    return to_3d(m_data->mesh.point(vi));
  }

  const Segment_3 segment_3(const Edge_index& ei) const {
    return Segment_3(
      point_3(m_data->mesh.source(m_data->mesh.halfedge(ei))),
      point_3(m_data->mesh.target(m_data->mesh.halfedge(ei))));
  }

  void set_iedge(
    const Vertex_index& v0, const Vertex_index& v1, const IEdge& iedge) const {

    const auto he = m_data->mesh.halfedge(v0, v1);
    CGAL_assertion(he != Halfedge_index());
    const auto ei = m_data->mesh.edge(he);
    m_data->e_iedge_map[ei] = iedge;
  }

  void set_ivertex(const Vertex_index& vi, const IVertex& ivertex) const {
    m_data->v_ivertex_map[vi] = ivertex;
  }

  void set_iedge(const Vertex_index& vi, const IEdge& iedge) const {
    m_data->v_iedge_map[vi] = iedge;
  }

  void set_iedge(const Edge_index& ei, const IEdge& iedge) const {
    m_data->e_iedge_map[ei] = iedge;
  }

  const IEdge& iedge(const Edge_index& ei) const {
    return m_data->e_iedge_map[ei];
  }

  const IEdge& iedge(const Vertex_index& vi) const {
    return m_data->v_iedge_map[vi];
  }

  const IVertex& ivertex(const Vertex_index& vi) const {
    return m_data->v_ivertex_map[vi];
  }

  bool has_iedge(const Edge_index& ei) const {
    return (m_data->e_iedge_map[ei] != Intersection_graph::null_iedge());
  }
  bool has_iedge(const Vertex_index& vi) const {
    return (m_data->v_iedge_map[vi] != Intersection_graph::null_iedge());
  }
  bool has_ivertex(const Vertex_index& vi) const {
    return (m_data->v_ivertex_map[vi] != Intersection_graph::null_ivertex());
  }

  const std::vector<std::size_t>& input(const Face_index& fi) const { return m_data->input_map[fi]; }
  std::vector<std::size_t>& input(const Face_index& fi) { return m_data->input_map[fi]; }

  bool is_initial(const Face_index& fi) const { return m_data->f_initial_map[fi]; }

  void set_initial(const Face_index& fi) { m_data->f_initial_map[fi] = true; }

  const int& k() const { return m_data->k; }
  int& k() { return m_data->k; }

  const std::set<IFace>& ifaces() const { return m_data->ifaces; }

  const IEdge_set& iedges() const { return m_data->iedges; }
  IEdge_set& iedges() { return m_data->iedges; }

  const Point_2 to_2d(const Point_3& point) const {
    return m_data->plane.to_2d(point);
  }

  const Vector_2 to_2d(const Vector_3& vec) const {
    return Vector_2(
      m_data->plane.to_2d(Point_3(0, 0, 0)),
      m_data->plane.to_2d(Point_3(0, 0, 0) + vec));
  }

  template<class IK = Intersection_kernel>
  auto to_2d(const typename Intersection_kernel::Point_3& point) const
    -> std::enable_if_t<!std::is_same_v<Kernel, IK>, const typename Intersection_kernel::Point_2> {
    return m_data->exact_plane.to_2d(point);
  }

  const Line_2 to_2d(const Line_3& line) const {
    return Line_2(
      m_data->plane.to_2d(line.point()),
      m_data->plane.to_2d(line.point() + line.to_vector()));
  }

  template <typename IK = Intersection_kernel>
  auto to_2d(const typename Intersection_kernel::Line_3& line) const
      -> std::enable_if_t<!std::is_same_v<Kernel, IK>, const typename Intersection_kernel::Line_2> {
    return typename Intersection_kernel::Line_2(
      m_data->exact_plane.to_2d(line.point()),
      m_data->exact_plane.to_2d(line.point() + line.to_vector()));
  }

  const Segment_2 to_2d(const Segment_3& segment) const {
    return Segment_2(
      m_data->plane.to_2d(segment.source()),
      m_data->plane.to_2d(segment.target()));
  }

  template <typename IK = Intersection_kernel>
  auto to_2d(const typename Intersection_kernel::Segment_3& segment) const
      -> std::enable_if_t<!std::is_same_v<Kernel, IK>, const typename Intersection_kernel::Segment_2> {
    return typename Intersection_kernel::Segment_2(
      m_data->exact_plane.to_2d(segment.source()),
      m_data->exact_plane.to_2d(segment.target()));
  }

  const Point_3 to_3d(const Point_2& point) const {
    return m_data->plane.to_3d(point);
  }

  template <typename IK = Intersection_kernel>
  auto to_3d(const typename Intersection_kernel::Point_2& point) const
    -> std::enable_if_t<!std::is_same_v<Kernel, IK>, const typename Intersection_kernel::Point_3> {
    return m_data->exact_plane.to_3d(point);
  }

  const Edge_index edge(const Vertex_index& v0, const Vertex_index& v1) {
    return m_data->mesh.edge(m_data->mesh.halfedge(v0, v1));
  }

  const Edge_index add_edge(
    const Vertex_index& v0, const Vertex_index& v1, const IEdge& iedge) {

    const auto out = m_data->mesh.edge(m_data->mesh.add_edge(v0, v1));
    m_data->e_iedge_map[out] = iedge;
    return out;
  }

  const Vertex_index add_vertex(const Point_2& point) {
    return m_data->mesh.add_vertex(point);
  }
};

template<typename GeomTraits, typename IntersectionKernel>
bool operator==(const Support_plane<GeomTraits, IntersectionKernel>& a, const Support_plane<GeomTraits, IntersectionKernel>& b) {

  if (a.is_bbox() || b.is_bbox()) {
    return false;
  }

  // Coplanar input polygons share a support plane, independent of the orientation of their planes.
  return (a.exact_plane() == b.exact_plane() || a.exact_plane() == b.exact_plane().opposite());
}

#endif //DOXYGEN_RUNNING

} // namespace internal
} // namespace KSP_3
} // namespace CGAL

#endif // CGAL_KSP_3_SUPPORT_LINE_H
