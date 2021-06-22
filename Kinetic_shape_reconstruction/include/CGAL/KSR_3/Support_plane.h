// Copyright (c) 2019 GeometryFactory SARL (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Simon Giraudot

#ifndef CGAL_KSR_3_SUPPORT_PLANE_H
#define CGAL_KSR_3_SUPPORT_PLANE_H

// #include <CGAL/license/Kinetic_shape_reconstruction.h>

// CGAL includes.
#include <CGAL/Surface_mesh.h>

// Internal includes.
#include <CGAL/KSR/utils.h>
#include <CGAL/KSR_3/Intersection_graph.h>

namespace CGAL {
namespace KSR_3 {

template<typename GeomTraits>
class Support_plane {

public:
  using Kernel = GeomTraits;

  using FT        = typename Kernel::FT;
  using Point_2   = typename Kernel::Point_2;
  using Point_3   = typename Kernel::Point_3;
  using Vector_2  = typename Kernel::Vector_2;
  using Vector_3  = typename Kernel::Vector_3;
  using Segment_2 = typename Kernel::Segment_2;
  using Segment_3 = typename Kernel::Segment_3;
  using Line_2    = typename Kernel::Line_2;
  using Line_3    = typename Kernel::Line_3;
  using Plane_3   = typename Kernel::Plane_3;

  using Mesh = CGAL::Surface_mesh<Point_2>;
  using Intersection_graph = KSR_3::Intersection_graph<Kernel>;
  using Bbox_2 = CGAL::Bbox_2;

  using IVertex = typename Intersection_graph::Vertex_descriptor;
  using IEdge   = typename Intersection_graph::Edge_descriptor;

  using Vertex_index   = typename Mesh::Vertex_index;
  using Face_index     = typename Mesh::Face_index;
  using Edge_index     = typename Mesh::Edge_index;
  using Halfedge_index = typename Mesh::Halfedge_index;

  using V_vector_map   = typename Mesh::template Property_map<Vertex_index, Vector_2>;
  using V_ivertex_map  = typename Mesh::template Property_map<Vertex_index, IVertex>;
  using V_iedge_map    = typename Mesh::template Property_map<Vertex_index, IEdge>;
  using V_bool_map     = typename Mesh::template Property_map<Vertex_index, bool>;
  using E_iedge_map    = typename Mesh::template Property_map<Edge_index, IEdge>;
  using F_index_map    = typename Mesh::template Property_map<Face_index, std::vector<std::size_t> >;
  using F_uint_map     = typename Mesh::template Property_map<Face_index, unsigned int>;
  using V_original_map = typename Mesh::template Property_map<Vertex_index, bool>;
  using V_time_map     = typename Mesh::template Property_map<Vertex_index, std::vector<FT> >;

private:
  struct Data {
    bool is_bbox;
    Point_3 centroid;
    Plane_3 plane;
    Mesh mesh;
    V_vector_map direction;
    V_ivertex_map v_ivertex_map;
    V_iedge_map v_iedge_map;
    V_bool_map v_active_map;
    E_iedge_map e_iedge_map;
    F_index_map input_map;
    F_uint_map k_map;
    V_original_map v_original_map;
    V_time_map v_time_map;
    std::set<IEdge> unique_iedges;
    std::vector<IEdge> iedges;
    std::vector<Segment_2> isegments;
    std::vector<Bbox_2> ibboxes;
    unsigned int k;
  };

  std::shared_ptr<Data> m_data;

public:
  Support_plane() :
  m_data(std::make_shared<Data>()) {
    add_property_maps();
  }

  template<typename PointRange>
  Support_plane(const PointRange& polygon, const bool is_bbox) :
  m_data(std::make_shared<Data>()) {

    std::vector<Point_3> points;
    points.reserve(polygon.size());
    for (const auto& point : polygon) {
      points.push_back(Point_3(
        static_cast<FT>(point.x()),
        static_cast<FT>(point.y()),
        static_cast<FT>(point.z())));
    }
    const std::size_t n = points.size();
    CGAL_assertion(n == polygon.size());

    // Newell's method.
    FT cx = FT(0), cy = FT(0), cz = FT(0);
    Vector_3 normal = CGAL::NULL_VECTOR;
    for (std::size_t i = 0; i < n; ++i) {
      const std::size_t ip = (i + 1) % n;
      const auto& pa = points[i];
      const auto& pb = points[ip];
      const FT x = normal.x() + (pa.y() - pb.y()) * (pa.z() + pb.z());
      const FT y = normal.y() + (pa.z() - pb.z()) * (pa.x() + pb.x());
      const FT z = normal.z() + (pa.x() - pb.x()) * (pa.y() + pb.y());
      normal = Vector_3(x, y, z);
      cx += pa.x();
      cy += pa.y();
      cz += pa.z();
    }
    CGAL_assertion_msg(normal != CGAL::NULL_VECTOR, "ERROR: BBOX IS FLAT!");
    CGAL_assertion(n != 0);
    cx /= static_cast<FT>(n);
    cy /= static_cast<FT>(n);
    cz /= static_cast<FT>(n);

    m_data->k = 0;
    m_data->plane = Plane_3(points[0], KSR::normalize(normal));
    m_data->centroid = Point_3(cx, cy, cz);
    m_data->is_bbox = is_bbox;
    add_property_maps();
  }

  void add_property_maps() {

    m_data->direction      = m_data->mesh.template add_property_map<Vertex_index, Vector_2>(
      "v:direction", CGAL::NULL_VECTOR).first;

    m_data->v_ivertex_map  = m_data->mesh.template add_property_map<Vertex_index, IVertex>(
      "v:ivertex", Intersection_graph::null_ivertex()).first;

    m_data->v_iedge_map    = m_data->mesh.template add_property_map<Vertex_index, IEdge>(
      "v:iedge", Intersection_graph::null_iedge()).first;

    m_data->v_active_map   = m_data->mesh.template add_property_map<Vertex_index, bool>(
      "v:active", true).first;

    m_data->e_iedge_map    = m_data->mesh.template add_property_map<Edge_index, IEdge>(
      "e:iedge", Intersection_graph::null_iedge()).first;

    m_data->input_map      = m_data->mesh.template add_property_map<Face_index, std::vector<std::size_t> >(
      "f:input", std::vector<std::size_t>()).first;

    m_data->k_map          = m_data->mesh.template add_property_map<Face_index, unsigned int>(
      "f:k", 0).first;

    m_data->v_original_map = m_data->mesh.template add_property_map<Vertex_index, bool>(
      "v:original", false).first;

    // TODO: I can have a similar vector to push all ivertices/events of the polygon vertex
    // to keep track of the path it traversed. Later, we can return this path.
    std::vector<FT> time_vector(1, FT(0));
    m_data->v_time_map     = m_data->mesh.template add_property_map<Vertex_index, std::vector<FT> >(
      "v:time", time_vector).first;
  }

  template<typename IG, typename SP>
  void convert(const IG& ig, SP& sp) {

    using CFT       = typename SP::Kernel::FT;
    using CPoint_2  = typename SP::Kernel::Point_2;
    using CPoint_3  = typename SP::Kernel::Point_3;
    using CPlane_3  = typename SP::Kernel::Plane_3;
    using CVector_2 = typename SP::Kernel::Vector_2;

    // using Converter = CGAL::Cartesian_converter<Kernel, typename SP::Kernel>;
    // Converter converter;

    const auto& vmap = ig.vmap();
    const auto& emap = ig.emap();

    std::set<CPoint_2> pts;
    std::map<Vertex_index, Vertex_index> map_vi;
    sp.data().k = m_data->k;
    sp.data().is_bbox = m_data->is_bbox;
    sp.data().centroid = CPoint_3(
      static_cast<CFT>(CGAL::to_double(m_data->centroid.x())),
      static_cast<CFT>(CGAL::to_double(m_data->centroid.y())),
      static_cast<CFT>(CGAL::to_double(m_data->centroid.z())));
    // sp.data().plane = converter(m_data->plane);
    sp.data().plane = CPlane_3(
      static_cast<CFT>(CGAL::to_double(m_data->plane.a())),
      static_cast<CFT>(CGAL::to_double(m_data->plane.b())),
      static_cast<CFT>(CGAL::to_double(m_data->plane.c())),
      static_cast<CFT>(CGAL::to_double(m_data->plane.d())));
    for (const auto& vertex : m_data->mesh.vertices()) {
      // const auto converted = converter(m_data->mesh.point(vertex));
      const CPoint_2 converted = CPoint_2(
        static_cast<CFT>(CGAL::to_double(m_data->mesh.point(vertex).x())),
        static_cast<CFT>(CGAL::to_double(m_data->mesh.point(vertex).y())));
      const bool is_inserted = pts.insert(converted).second;
      const auto vi = sp.data().mesh.add_vertex();
      map_vi[vertex] = vi;

      if (is_inserted) {
        sp.data().mesh.point(vi) = converted;
      } else {
        sp.data().mesh.point(vi) = converted;

        // using CFT = typename SP::Kernel::FT;
        // const CFT b1 = CFT(9) / CFT(10);
        // const CFT b2 = CFT(1) / CFT(10);

        // const auto pi = this->prev(vertex);
        // const auto pc = converter(m_data->mesh.point(pi));
        // const auto ni = this->next(vertex);
        // const auto nc = converter(m_data->mesh.point(ni));

        // if (nc != converted) {
        //   const auto x = b1 * converted.x() + b2 * nc.x();
        //   const auto y = b1 * converted.y() + b2 * nc.y();
        //   const CPoint_2 new_point(x, y);
        //   sp.data().mesh.point(vi) = new_point;
        //   // std::cout << "or: " << to_3d(Point_2(converted.x(), converted.y())) << std::endl;
        //   // std::cout << "nc: " << to_3d(Point_2(new_point.x(), new_point.y())) << std::endl;
        // } else if (pc != converted) {
        //   const auto x = b1 * converted.x() + b2 * pc.x();
        //   const auto y = b1 * converted.y() + b2 * pc.y();
        //   const CPoint_2 new_point(x, y);
        //   sp.data().mesh.point(vi) = new_point;
        //   // std::cout << "or: " << to_3d(Point_2(converted.x(), converted.y())) << std::endl;
        //   // std::cout << "pc: " << to_3d(Point_2(new_point.x(), new_point.y())) << std::endl;
        // } else {
        //   CGAL_assertion_msg(false, "ERROR: WE HAVE THREE EQUAL POINTS!");
        // }
      }
    }
    CGAL_assertion(sp.data().mesh.number_of_vertices() == m_data->mesh.number_of_vertices());

    std::map<Face_index, Face_index> map_fi;
    std::vector<Vertex_index> mapped_vertices;
    for (const auto& face : m_data->mesh.faces()) {
      const auto vertices = CGAL::vertices_around_face(
        m_data->mesh.halfedge(face), m_data->mesh);

      mapped_vertices.clear();
      mapped_vertices.reserve(vertices.size());
      for (const auto vertex : vertices) {
        mapped_vertices.push_back(map_vi.at(vertex));
      }
      CGAL_assertion(mapped_vertices.size() == vertices.size());
      const auto fi = sp.data().mesh.add_face(mapped_vertices);
      map_fi[face] = fi;
    }
    CGAL_assertion(sp.data().mesh.number_of_faces() == m_data->mesh.number_of_faces());

    for (const auto& vertex : m_data->mesh.vertices()) {
      const auto vi = map_vi.at(vertex);
      // sp.data().direction[vi] = converter(m_data->direction[vertex]);
      sp.data().direction[vi] = CVector_2(
        static_cast<CFT>(CGAL::to_double(m_data->direction[vertex].x())),
        static_cast<CFT>(CGAL::to_double(m_data->direction[vertex].y())));

      const auto ivertex = m_data->v_ivertex_map[vertex];
      if (ivertex != IG::null_ivertex()) {
        sp.data().v_ivertex_map[vi] = vmap.at(ivertex);
      } else {
        sp.data().v_ivertex_map[vi] = ivertex;
      }

      const auto iedge = m_data->v_iedge_map[vertex];
      if (iedge != IG::null_iedge()) {
        sp.data().v_iedge_map[vi] = emap.at(iedge);
      } else {
        sp.data().v_iedge_map[vi] = iedge;
      }

      sp.data().v_active_map[vi]   = m_data->v_active_map[vertex];
      sp.data().v_original_map[vi] = m_data->v_original_map[vertex];

      // sp.data().v_time_map[vi] = converter(m_data->v_time_map[vertex]);
      // sp.data().v_time_map[vi] = static_cast<CFT>(CGAL::to_double(m_data->v_time_map[vertex]));

      sp.data().v_time_map[vi].clear();
      sp.data().v_time_map[vi].reserve(m_data->v_time_map[vertex].size());
      for (const auto vtime : m_data->v_time_map[vertex]) {
        sp.data().v_time_map[vi].push_back(static_cast<CFT>(CGAL::to_double(vtime)));
      }
      CGAL_assertion(
        sp.data().v_time_map[vi].size() == m_data->v_time_map[vertex].size());
    }

    for (const auto& edge : m_data->mesh.edges()) {
      const auto source = m_data->mesh.source(m_data->mesh.halfedge(edge));
      const auto target = m_data->mesh.target(m_data->mesh.halfedge(edge));

      const auto s = map_vi[source];
      const auto t = map_vi[target];
      const auto he = sp.data().mesh.halfedge(s, t);
      const auto ei = sp.data().mesh.edge(he);

      const auto iedge = m_data->e_iedge_map[edge];
      if (iedge != IG::null_iedge()) {
        sp.data().e_iedge_map[ei] = emap.at(iedge);
      } else {
        sp.data().e_iedge_map[ei] = iedge;
      }
    }

    for (const auto& face : m_data->mesh.faces()) {
      const auto fi = map_fi.at(face);
      sp.data().input_map[fi] = m_data->input_map[face];
      sp.data().k_map[fi]     = m_data->k_map[face];
    }

    sp.data().unique_iedges.clear();
    for (const auto& iedge : m_data->unique_iedges) {
      CGAL_assertion(iedge != IG::null_iedge());
      sp.data().unique_iedges.insert(emap.at(iedge));
    }
  }

  Data& data() { return *m_data; }

  void clear_pfaces() {
    m_data->mesh.clear();
    add_property_maps();
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
    input_vec.push_back(KSR::no_element());
    return vertices;
  }

  template<typename Pair>
  std::size_t add_input_polygon(
    const std::vector<Pair>& points,
    const Point_2& centroid,
    const std::vector<std::size_t>& input_indices) {

    CGAL_assertion(is_simple_polygon(points));
    CGAL_assertion(is_convex_polygon(points));
    CGAL_assertion(is_valid_polygon(points));

    std::vector<Vertex_index> vertices;
    const std::size_t n = points.size();
    CGAL_assertion(n >= 3);
    vertices.reserve(n);

    FT sum_length = FT(0);
    std::vector<Vector_2> directions;
    directions.reserve(n);

    for (const auto& pair : points) {
      const auto& point = pair.first;
      directions.push_back(Vector_2(centroid, point));
      const FT length = static_cast<FT>(
        CGAL::sqrt(CGAL::to_double(CGAL::abs(directions.back() * directions.back()))));
      sum_length += length;
    }
    CGAL_assertion(directions.size() == n);
    sum_length /= static_cast<FT>(n);

    for (std::size_t i = 0; i < n; ++i) {
      const auto& point = points[i].first;
      const auto vi = m_data->mesh.add_vertex(point);
      m_data->direction[vi] = directions[i] / sum_length;
      m_data->v_original_map[vi] = true;
      vertices.push_back(vi);
    }

    const auto fi = m_data->mesh.add_face(vertices);
    CGAL_assertion(fi != Mesh::null_face());
    auto& input_vec = m_data->input_map[fi];
    CGAL_assertion(input_vec.empty());
    for (const std::size_t input_index : input_indices) {
      input_vec.push_back(input_index);
    }
    return static_cast<std::size_t>(fi);
  }

  template<typename Pair>
  bool is_valid_polygon(const std::vector<Pair>& polygon) const {

    const FT ptol = KSR::point_tolerance<FT>();
    for (std::size_t i = 0; i < polygon.size(); ++i) {
      const std::size_t ip = (i + 1) % polygon.size();
      const auto& p = polygon[i].first;
      const auto& q = polygon[ip].first;
      const FT distance = KSR::distance(p, q);
      const bool is_equal_zero = (distance < ptol);
      CGAL_assertion_msg(!is_equal_zero,
      "ERROR: WE HAVE EQUAL POINTS IN THE INPUT POLYGON!");
      if (is_equal_zero) return false;
    }
    return true;
  }

  template<typename Pair>
  bool is_simple_polygon(const std::vector<Pair>& points) const {
    std::vector<Point_2> polygon;
    polygon.reserve(points.size());
    for (const auto& pair : points)
      polygon.push_back(pair.first);
    CGAL_assertion(polygon.size() == points.size());
    return CGAL::is_simple_2(polygon.begin(), polygon.end());
  }

  template<typename Pair>
  bool is_convex_polygon(const std::vector<Pair>& points) const {
    std::vector<Point_2> polygon;
    polygon.reserve(points.size());
    for (const auto& pair : points)
      polygon.push_back(pair.first);
    CGAL_assertion(polygon.size() == points.size());
    return CGAL::is_convex_2(polygon.begin(), polygon.end());
  }

  const Plane_3& plane() const { return m_data->plane; }
  const Point_3& centroid() const { return m_data->centroid; }
  bool is_bbox() const { return m_data->is_bbox; }

  const Mesh& mesh() const { return m_data->mesh; }
  Mesh& mesh() { return m_data->mesh; }

  const Point_2& get_point(const Vertex_index& vi) const {
    return m_data->mesh.point(vi);
  }

  void set_point(const Vertex_index& vi, const Point_2& point) {
    m_data->mesh.point(vi) = point;
  }

  void set_last_event_time(const Vertex_index& vi, const FT time) {
    // TODO: If we do not need the full vector, remove it.
    m_data->v_time_map[vi].push_back(time);
  }

  const FT last_event_time(const Vertex_index& vi, const FT /* curr_time */) const {

    // FT last_time = FT(-1);
    // const FT tol = KSR::tolerance<FT>();
    // CGAL_assertion(m_data->v_time_map[vi].size() > 0);

    // // std::cout << "----" << std::endl;
    // for (const FT vtime : m_data->v_time_map[vi]) {
    //   // std::cout << "vtime: " << vtime << std::endl;
    //   const FT time_diff = CGAL::abs(curr_time - vtime);
    //   if (time_diff < tol) continue;
    //   last_time = vtime;
    // }
    // CGAL_assertion(last_time >= FT(0));
    // return last_time;

    return m_data->v_time_map[vi].back();
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

  const Point_2 point_2(const Vertex_index& vi, const FT time) const {
    CGAL_assertion(time >= FT(0));
    return m_data->mesh.point(vi) + time * m_data->direction[vi];
  }

  const Point_3 point_3(const Vertex_index& vi, const FT time) const {
    return to_3d(point_2(vi, time));
  }

  const Segment_2 segment_2(const Edge_index& ei, const FT time) const {

    return Segment_2(
      point_2(m_data->mesh.source(m_data->mesh.halfedge(ei)), time),
      point_2(m_data->mesh.target(m_data->mesh.halfedge(ei)), time));
  }

  const Segment_3 segment_3(const Edge_index& ei, const FT time) const {
    return Segment_3(
      point_3(m_data->mesh.source(m_data->mesh.halfedge(ei)), time),
      point_3(m_data->mesh.target(m_data->mesh.halfedge(ei)), time));
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

  const Vector_2& direction(const Vertex_index& vi) const { return m_data->direction[vi]; }
  Vector_2& direction(const Vertex_index& vi) { return m_data->direction[vi]; }

  const FT speed(const Vertex_index& vi) const {
    return static_cast<FT>(CGAL::sqrt(
      CGAL::to_double(CGAL::abs(m_data->direction[vi].squared_length()))));
  }

  const std::vector<std::size_t>& input(const Face_index& fi) const { return m_data->input_map[fi]; }
  std::vector<std::size_t>& input(const Face_index& fi) { return m_data->input_map[fi]; }

  bool is_original(const Vertex_index& vi) const { return m_data->v_original_map[vi]; }

  const unsigned int& k() const { return m_data->k; }
  unsigned int& k() { return m_data->k; }

  const unsigned int& k(const Face_index& /* fi */) const {
    return m_data->k;
    // return m_data->k_map[fi];
  }
  unsigned int& k(const Face_index& /* fi */) {
    return m_data->k;
    // return m_data->k_map[fi];
  }

  bool is_active(const Vertex_index& vi) const { return m_data->v_active_map[vi]; }
  void set_active(const Vertex_index& vi, const bool value) { m_data->v_active_map[vi] = value; }

  bool is_frozen(const Vertex_index& vi) const {
    return (m_data->direction[vi] == CGAL::NULL_VECTOR);
  }

  const std::set<IEdge>& unique_iedges() const { return m_data->unique_iedges; }
  std::set<IEdge>& unique_iedges() { return m_data->unique_iedges; }

  const std::vector<IEdge>& iedges() const { return m_data->iedges; }
  std::vector<IEdge>& iedges() { return m_data->iedges; }

  const std::vector<Segment_2>& isegments() const { return m_data->isegments; }
  std::vector<Segment_2>& isegments() { return m_data->isegments; }

  const std::vector<Bbox_2>& ibboxes() const { return m_data->ibboxes; }
  std::vector<Bbox_2>& ibboxes() { return m_data->ibboxes; }

  const Point_2 to_2d(const Point_3& point) const {
    return m_data->plane.to_2d(point);
  }

  const Line_2 to_2d(const Line_3& line) const {
    return Line_2(
      m_data->plane.to_2d(line.point()),
      m_data->plane.to_2d(line.point() + line.to_vector()));
  }

  const Segment_2 to_2d(const Segment_3& segment) const {
    return Segment_2(
      m_data->plane.to_2d(segment.source()),
      m_data->plane.to_2d(segment.target()));
  }

  const Vector_3 to_3d(const Vector_2& vec) const {
    return Vector_3(
      m_data->plane.to_3d(Point_2(FT(0), FT(0))),
      m_data->plane.to_3d(Point_2(FT(0), FT(0)) + vec));
  }

  const Point_3 to_3d(const Point_2& point) const {
    return m_data->plane.to_3d(point);
  }

  const Edge_index edge(const Vertex_index& v0, const Vertex_index& v1) {

    // std::cout << int(v0) << " : " << int(v1) << std::endl;
    // std::cout << int(m_data->mesh.halfedge(v0, v1)) << std::endl;
    // std::cout << int(m_data->mesh.halfedge(v1, v0)) << std::endl;
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

  const Vertex_index duplicate_vertex(const Vertex_index& v) {
    const auto vi = m_data->mesh.add_vertex(m_data->mesh.point(v));
    m_data->direction[vi]     = m_data->direction[v];
    m_data->v_ivertex_map[vi] = m_data->v_ivertex_map[v];
    m_data->v_iedge_map[vi]   = m_data->v_iedge_map[v];
    return vi;
  }

  void remove_vertex(const Vertex_index& vi) {
    m_data->mesh.remove_vertex(vi);
  }

  const Edge_index split_vertex(const Vertex_index& vi) {
    return m_data->mesh.edge(
      CGAL::Euler::split_vertex(
        m_data->mesh.halfedge(vi),
        m_data->mesh.opposite(m_data->mesh.next(m_data->mesh.halfedge(vi))),
        m_data->mesh));
  }

  const Vertex_index split_edge(
    const Vertex_index& v0, const Vertex_index& v1) {

    return m_data->mesh.target(
      CGAL::Euler::split_edge(m_data->mesh.halfedge(v0, v1), m_data->mesh));
  }
};

template<typename Kernel>
bool operator==(const Support_plane<Kernel>& a, const Support_plane<Kernel>& b) {

  if (a.is_bbox() || b.is_bbox()) {
    return false;
  }

  using FT = typename Kernel::FT;
  const auto& planea = a.plane();
  const auto& planeb = b.plane();

  const auto va = planea.orthogonal_vector();
  const auto vb = planeb.orthogonal_vector();

  // Are the planes parallel?
  // const FT vtol = KSR::vector_tolerance<FT>();
  // const FT aval = CGAL::abs(va * vb);
  // if (aval < vtol) {
  //   return false;
  // }

  const FT vtol = FT(5); // degrees // TODO: We should put it as a parameter.
  const FT aval = KSR::angle_3d(va, vb);
  if (aval >= vtol) {
    return false;
  }

  // std::cout << "aval: " << aval << " : " << vtol << std::endl;

  // Are the planes coplanar?
  // const FT ptol = KSR::point_tolerance<FT>();
  // const auto pa = planea.point();
  // const auto pb = planeb.projection(pa);
  // const FT bval = KSR::distance(pa, pb);
  // TODO: Should we rotate the planes here before computing the distance?

  const FT ptol = FT(5) / FT(10); // TODO: We should put it as a parameter.
  const auto pa1 = a.centroid();
  const auto pb1 = planeb.projection(pa1);
  const auto pb2 = b.centroid();
  const auto pa2 = planea.projection(pb2);

  const FT bval1 = KSR::distance(pa1, pb1);
  const FT bval2 = KSR::distance(pa2, pb2);
  const FT bval = (CGAL::min)(bval1, bval2);

  // if (bval < ptol) {
  //   std::cout << "2 " << pa << " " << pb << std::endl;
  //   std::cout << "bval: " << bval << " : " << ptol << std::endl;
  // }

  // std::cout << "bval: " << bval << " : " << ptol << std::endl;
  if (bval >= ptol) return false;
  // std::cout << "- found coplanar planes" << std::endl;
  return true;
}

} // namespace KSR_3
} // namespace CGAL

#endif // CGAL_KSR_3_SUPPORT_LINE_H
