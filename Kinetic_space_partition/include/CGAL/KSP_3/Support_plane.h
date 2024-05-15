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

  using V_vector_map = typename Mesh::template Property_map<Vertex_index, Vector_2>;
  using V_ivertex_map = typename Mesh::template Property_map<Vertex_index, IVertex>;
  using V_iedge_map = typename Mesh::template Property_map<Vertex_index, IEdge>;
  using V_bool_map = typename Mesh::template Property_map<Vertex_index, bool>;
  using E_iedge_map = typename Mesh::template Property_map<Edge_index, IEdge>;
  using F_index_map = typename Mesh::template Property_map<Face_index, std::vector<std::size_t> >;
  using F_uint_map = typename Mesh::template Property_map<Face_index, unsigned int>;
  using F_bool_map = typename Mesh::template Property_map<Face_index, bool>;
  using V_original_map = typename Mesh::template Property_map<Vertex_index, bool>;
  using V_time_map = typename Mesh::template Property_map<Vertex_index, std::vector<FT> >;

  struct Face_event {
    Face_event() {}
    Face_event(std::size_t sp_idx, FT time, IEdge edge, IFace face) : support_plane(sp_idx), time(time), crossed_edge(edge), face(face) {}
    std::size_t support_plane;
    FT time;
    FT intersection_bary;
    IEdge crossed_edge;
    IFace face; // The face that does not yet belong to the region.
  };

  struct Data {
    Data() : mesh(),
      v_ivertex_map(mesh.template add_property_map<Vertex_index, IVertex>("v:ivertex", Intersection_graph::null_ivertex()).first),
      v_iedge_map(mesh.template add_property_map<Vertex_index, IEdge>("v:iedge", Intersection_graph::null_iedge()).first),
      e_iedge_map(mesh.template add_property_map<Edge_index, IEdge>("e:iedge", Intersection_graph::null_iedge()).first),
      input_map(mesh.template add_property_map<Face_index, std::vector<std::size_t> >("f:input", std::vector<std::size_t>()).first),
      f_initial_map(mesh.template add_property_map<Face_index, bool >("f:initial", false).first),
      v_original_map(mesh.template add_property_map<Vertex_index, bool>("v:original", false).first) {}

    bool is_bbox;
    Point_2 centroid;
    Plane_3 plane;
    typename Intersection_kernel::Plane_3 exact_plane;
    Mesh mesh;

    V_ivertex_map v_ivertex_map;
    V_iedge_map v_iedge_map;
    E_iedge_map e_iedge_map;
    F_index_map input_map;
    F_bool_map f_initial_map;
    V_original_map v_original_map;
    std::map<IEdge, std::pair<IFace, IFace> > iedge2ifaces;
    std::set<IFace> ifaces; // All ifaces in the support plane
    std::vector<IFace> initial_ifaces; // IFaces which intersect with the input polygon and are thus part of the mesh before the propagation starts.
    std::vector<Face_index> initial_pfaces;
    std::map<IVertex, Vertex_index> ivertex2pvertex;
    IEdge_set unique_iedges;
    std::set<std::size_t> crossed_lines;

    std::vector<IEdge> iedges;
    std::vector<Point_2> original_vertices;
    std::vector<Vector_2> original_vectors;
    std::vector<Direction_2> original_directions;
    std::vector<typename Intersection_kernel::Ray_2> original_rays;

    FT distance_tolerance;
    FT angle_tolerance;

    std::size_t actual_input_polygon;

    int k;
  };

private:

  std::shared_ptr<Data> m_data;

public:
  Support_plane() : m_data(std::make_shared<Data>()) {}

  template<typename PointRange>
  Support_plane(const PointRange& polygon, const bool is_bbox, typename Intersection_kernel::Plane_3 plane) :
    m_data(std::make_shared<Data>()) {

    std::vector<Point_3> points;
    points.reserve(polygon.size());
    for (const auto& point : polygon) {
      points.push_back(Point_3(
        static_cast<FT>(point.x()),
        static_cast<FT>(point.y()),
        static_cast<FT>(point.z())));
    }
    CGAL_assertion(points.size() == polygon.size());

    From_exact from_exact;

    m_data->k = 0;
    m_data->plane = from_exact(plane);
    m_data->exact_plane = plane;
    m_data->is_bbox = is_bbox;
    m_data->distance_tolerance = 0;
    m_data->angle_tolerance = 0;
    m_data->actual_input_polygon = static_cast<std::size_t>(- 1);

    std::vector<Triangle_2> tris(points.size() - 2);
    for (std::size_t i = 2; i < points.size(); i++) {
      tris[i - 2] = Triangle_2(to_2d(points[0]), to_2d(points[i - 1]), to_2d(points[i]));
    }

    m_data->centroid = CGAL::centroid(tris.begin(), tris.end(), CGAL::Dimension_tag<2>());

    add_property_maps();
  }

  template<typename PointRange>
  Support_plane(const PointRange& polygon, const bool is_bbox) :
    m_data(std::make_shared<Data>()) {
    To_exact to_exact;

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

    Vector_3 normal = CGAL::NULL_VECTOR;
    for (std::size_t i = 0; i < n; ++i) {
      const std::size_t ip = (i + 1) % n;
      const auto& pa = points[i];
      const auto& pb = points[ip];
      const FT x = normal.x() + (pa.y() - pb.y()) * (pa.z() + pb.z());
      const FT y = normal.y() + (pa.z() - pb.z()) * (pa.x() + pb.x());
      const FT z = normal.z() + (pa.x() - pb.x()) * (pa.y() + pb.y());
      normal = Vector_3(x, y, z);
    }
    CGAL_assertion_msg(normal != CGAL::NULL_VECTOR, "ERROR: BBOX IS FLAT!");
    CGAL_assertion(n != 0);

    m_data->k = 0;
    m_data->plane = Plane_3(points[0], KSP::internal::normalize(normal));
    m_data->exact_plane = to_exact(m_data->plane);
    m_data->is_bbox = is_bbox;
    m_data->distance_tolerance = 0;
    m_data->angle_tolerance = 0;
    m_data->actual_input_polygon = -1;

    std::vector<Triangle_2> tris(points.size() - 2);
    for (std::size_t i = 2; i < points.size(); i++) {
      tris[i - 2] = Triangle_2(to_2d(points[0]), to_2d(points[i - 1]), to_2d(points[i]));
    }

    m_data->centroid = CGAL::centroid(tris.begin(), tris.end(), CGAL::Dimension_tag<2>());

    add_property_maps();
  }

  Support_plane(const std::vector<typename Intersection_kernel::Point_3>& polygon, const bool is_bbox) :
    m_data(std::make_shared<Data>()) {
    From_exact from_exact;

    std::vector<Point_3> points;
    points.reserve(polygon.size());
    for (const auto& point : polygon) {
      points.push_back(Point_3(
        from_exact(point.x()),
        from_exact(point.y()),
        from_exact(point.z())));
    }
    CGAL_assertion_code(const std::size_t n = points.size();)
    CGAL_assertion(n == polygon.size());
    CGAL_assertion(n != 0);

    m_data->k = 0;
    m_data->exact_plane = typename Intersection_kernel::Plane_3(polygon[0], polygon[1], polygon[2]);
    m_data->plane = from_exact(m_data->exact_plane);
    m_data->is_bbox = is_bbox;
    m_data->distance_tolerance = 0;
    m_data->angle_tolerance = 0;
    m_data->actual_input_polygon = static_cast<std::size_t>(- 1);

    std::vector<Triangle_2> tris(points.size() - 2);
    for (std::size_t i = 2; i < points.size(); i++) {
      tris[i - 2] = Triangle_2(to_2d(points[0]), to_2d(points[i - 1]), to_2d(points[i]));
    }

    m_data->centroid = CGAL::centroid(tris.begin(), tris.end(), CGAL::Dimension_tag<2>());

    add_property_maps();
  }

  void add_property_maps() {

    m_data->v_ivertex_map = m_data->mesh.template add_property_map<Vertex_index, IVertex>("v:ivertex", Intersection_graph::null_ivertex()).first;
    m_data->v_iedge_map = m_data->mesh.template add_property_map<Vertex_index, IEdge>("v:iedge", Intersection_graph::null_iedge()).first;
    m_data->e_iedge_map = m_data->mesh.template add_property_map<Edge_index, IEdge>("e:iedge", Intersection_graph::null_iedge()).first;
    m_data->input_map = m_data->mesh.template add_property_map<Face_index, std::vector<std::size_t> >("f:input", std::vector<std::size_t>()).first;
    m_data->v_original_map = m_data->mesh.template add_property_map<Vertex_index, bool>("v:original", false).first;
    m_data->f_initial_map = m_data->mesh.template add_property_map<Face_index, bool >("f:initial", false).first;
  }

  void link_property_maps() {
    m_data->v_ivertex_map = m_data->mesh.template property_map<Vertex_index, IVertex>("v:ivertex").first;
    m_data->v_iedge_map = m_data->mesh.template property_map<Vertex_index, IEdge>("v:iedge").first;
    m_data->e_iedge_map = m_data->mesh.template property_map<Edge_index, IEdge>("e:iedge").first;
    m_data->input_map = m_data->mesh.template property_map<Face_index, std::vector<std::size_t> >("f:input").first;
    m_data->v_original_map = m_data->mesh.template property_map<Vertex_index, bool>("v:original").first;
    m_data->f_initial_map = m_data->mesh.template property_map<Face_index, bool >("f:initial").first;
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

  void get_border(Intersection_graph& igraph, std::vector<IEdge>& border) {
    border.clear();
    auto m = mesh();

    Vertex_index s = Mesh::null_vertex();

    for (auto v : m_data->mesh.vertices()) {
      if (m_data->mesh.is_border(v)) {
        s = v;
        break;
      }
    }

    if (s == Mesh::null_vertex()) {
      std::cout << "Support plane does not have border vertices" << std::endl;
      return;
    }

    auto h = m.halfedge(s);
    if (!m.is_border(h))
      h = m.opposite(h);

    auto n = h;
    IVertex last = ivertex(s);
    do {
      n = m.next(n);
      IVertex current = ivertex(m.target(n));
      border.push_back(igraph.edge(last, current));
      last = current;
    } while (n != h && n != Mesh::null_halfedge());

    if (n == Mesh::null_halfedge()) {
      std::cout << " NULL_HALFEDGE!";
    }
    //std::cout << "edges: " << border.size() << std::endl;
  }

  void get_border(Intersection_graph& igraph, const Face_index& fi, std::vector<IEdge>& border) {
    border.clear();
    auto m = mesh();

    auto first = m.halfedge(fi);
    auto h = first;
    do {
      auto o = m.opposite(h);

      if (m.is_border(o))
        border.push_back(igraph.edge(ivertex(m.target(h)), ivertex(m.target(o))));

      h = m.next(h);
    } while (h != first && h != Mesh::null_halfedge());

    if (h == Mesh::null_halfedge()) {
      std::cout << " NULL_HALFEDGE!";
    }
    //std::cout << "edges: " << border.size() << std::endl;
  }

  Data& data() { return *m_data; }

  FT distance_tolerance() const {
    return m_data->distance_tolerance;
  }

  FT angle_tolerance() const {
    return m_data->angle_tolerance;
  }

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
    input_vec.push_back(std::size_t(-1));
    return vertices;
  }

  template<typename Pair>
  std::size_t add_input_polygon(
    const std::vector<Pair>& points,
    const std::vector<std::size_t>& input_indices) {

    CGAL_assertion(is_simple_polygon(points));
    CGAL_assertion(is_convex_polygon(points));
    CGAL_assertion(is_valid_polygon(points));

    To_exact to_exact;

    CGAL_assertion(points.size() >= 3);
    std::vector<Triangle_2> tris(points.size() - 2);
    for (std::size_t i = 2; i < points.size(); i++)
      tris[i - 2] = Triangle_2(points[0].first, points[i - 1].first, points[i].first);

    m_data->centroid = CGAL::centroid(tris.begin(), tris.end(), CGAL::Dimension_tag<2>());

    std::vector<Vertex_index> vertices;
    const std::size_t n = points.size();
    CGAL_assertion(n >= 3);
    vertices.reserve(n);
    m_data->original_vertices.resize(n);
    m_data->original_vectors.resize(n);
    m_data->original_directions.resize(n);
    m_data->original_rays.resize(n);

    FT sum_length = FT(0);
    std::vector<Vector_2> directions;
    directions.reserve(n);

    std::vector<std::pair<std::size_t, Direction_2> > dir_vec;

    for (const auto& pair : points) {
      const auto& point = pair.first;
      directions.push_back(Vector_2(m_data->centroid, point));
      const FT length = static_cast<FT>(
        CGAL::sqrt(CGAL::to_double(CGAL::abs(directions.back() * directions.back()))));
      sum_length += length;
    }
    CGAL_assertion(directions.size() == n);
    sum_length /= static_cast<FT>(n);

    dir_vec.reserve(n);
    for (std::size_t i = 0; i < n; i++)
      dir_vec.push_back(std::pair<std::size_t, Direction_2>(i, directions[i]));

    std::sort(dir_vec.begin(), dir_vec.end(),
      [&](const std::pair<std::size_t, Direction_2>& a,
        const std::pair<std::size_t, Direction_2>& b) -> bool {
          return a.second < b.second;
      });

    for (std::size_t i = 0; i < n; ++i) {
      const auto& point = points[dir_vec[i].first].first;
      const auto vi = m_data->mesh.add_vertex(point);
      m_data->original_vertices[i] = point;
      m_data->original_vectors[i] = directions[dir_vec[i].first] / sum_length;
      m_data->original_directions[i] = Direction_2(directions[dir_vec[i].first]);
      m_data->original_rays[i] = typename Intersection_kernel::Ray_2(to_exact(point), to_exact(m_data->original_directions[i]));
      m_data->v_original_map[vi] = true;
      vertices.push_back(vi);
    }

    for (std::size_t i = 0; i < m_data->original_directions.size(); i++) {
      for (std::size_t j = 0; j < m_data->original_directions.size(); j++) {
        if (j < i)
          assert(m_data->original_directions[j] < m_data->original_directions[i]);
        if (j > i)
          assert(m_data->original_directions[i] < m_data->original_directions[j]);
      }
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

  bool has_crossed_line(std::size_t line) const {
    return m_data->crossed_lines.find(line) != m_data->crossed_lines.end();
  }

  void set_crossed_line(std::size_t line) {
    m_data->crossed_lines.insert(line);
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
      const bool is_equal_zero = (KSP::internal::distance(p, q) == 0);
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
  const typename Intersection_kernel::Plane_3& exact_plane() const { return m_data->exact_plane; }
  const Point_2& centroid() const { return m_data->centroid; }
  bool is_bbox() const { return m_data->is_bbox; }
  std::map<IVertex, Vertex_index>& ivertex2pvertex() { return m_data->ivertex2pvertex; }

  const Mesh& mesh() const { return m_data->mesh; }
  Mesh& mesh() { return m_data->mesh; }

  const Point_2& get_point(const Vertex_index& vi) const {
    return m_data->mesh.point(vi);
  }

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

  std::size_t has_ifaces(IEdge edge) const {
    auto it = m_data->iedge2ifaces.find(edge);
    if (it == m_data->iedge2ifaces.end())
      return 0;
    if (it->second.second != Intersection_graph::null_iface())
      return 2;
    else
      return 1;
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

  const Point_2 point_2(const Vertex_index& vi) const {
    return m_data->mesh.point(vi);
  }

  const Point_3 point_3(const Vertex_index& vi) const {
    return to_3d(m_data->mesh.point(vi));
  }

  const Segment_2 segment_2(const Edge_index& ei) const {
    return Segment_2(m_data->mesh.point(m_data->mesh.source(m_data->mesh.halfedge(ei))), m_data->mesh.point(m_data->mesh.target(m_data->mesh.halfedge(ei))));
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

  const Vector_2 original_edge_direction(std::size_t v1, std::size_t v2) const {
    const Vector_2 edge = m_data->original_vertices[v1] - m_data->original_vertices[v2];
    Vector_2 orth = Vector_2(-edge.y(), edge.x());
    orth = (1.0 / (CGAL::sqrt(orth * orth))) * orth;
    FT s1 = orth * m_data->original_vectors[v1];
    FT s2 = orth * m_data->original_vectors[v2];

    if (abs(s1 - s2) > 0.0001)
      std::cout << "edge speed seems inconsistent" << std::endl;

    return s1 * orth;
  }

  const FT speed(const Vertex_index& vi) const {
    return static_cast<FT>(CGAL::sqrt(
      CGAL::to_double(CGAL::abs(m_data->direction[vi].squared_length()))));
  }

  const std::vector<std::size_t>& input(const Face_index& fi) const { return m_data->input_map[fi]; }
  std::vector<std::size_t>& input(const Face_index& fi) { return m_data->input_map[fi]; }

  bool is_original(const Vertex_index& vi) const { return m_data->v_original_map[vi]; }

  bool is_initial(const Face_index& fi) const { return m_data->f_initial_map[fi]; }

  void set_initial(const Face_index& fi) { m_data->f_initial_map[fi] = true; }

  const int& k() const { return m_data->k; }
  int& k() { return m_data->k; }

  const std::set<IFace>& ifaces() const { return m_data->ifaces; }

  const IEdge_set& unique_iedges() const { return m_data->unique_iedges; }
  IEdge_set& unique_iedges() { return m_data->unique_iedges; }

  const std::vector<IEdge>& iedges() const { return m_data->iedges; }
  std::vector<IEdge>& iedges() { return m_data->iedges; }

  const Point_2 to_2d(const Point_3& point) const {
    return m_data->plane.to_2d(point);
  }

  const Vector_2 to_2d(const Vector_3& vec) const {
    return Vector_2(
      m_data->plane.to_2d(Point_3(0, 0, 0)),
      m_data->plane.to_2d(Point_3(0, 0, 0) + vec));
  }

  const typename Intersection_kernel::Point_2 to_2d(const typename Intersection_kernel::Point_3& point) const {
    return m_data->exact_plane.to_2d(point);
  }

  const Line_2 to_2d(const Line_3& line) const {
    return Line_2(
      m_data->plane.to_2d(line.point()),
      m_data->plane.to_2d(line.point() + line.to_vector()));
  }

  const typename Intersection_kernel::Line_2 to_2d(const typename Intersection_kernel::Line_3& line) const {
    return typename Intersection_kernel::Line_2(
      m_data->exact_plane.to_2d(line.point()),
      m_data->exact_plane.to_2d(line.point() + line.to_vector()));
  }

  const Segment_2 to_2d(const Segment_3& segment) const {
    return Segment_2(
      m_data->plane.to_2d(segment.source()),
      m_data->plane.to_2d(segment.target()));
  }

  const typename Intersection_kernel::Segment_2 to_2d(const typename Intersection_kernel::Segment_3& segment) const {
    return typename Intersection_kernel::Segment_2(
      m_data->exact_plane.to_2d(segment.source()),
      m_data->exact_plane.to_2d(segment.target()));
  }

  const Vector_3 to_3d(const Vector_2& vec) const {
    return Vector_3(
      m_data->plane.to_3d(Point_2(FT(0), FT(0))),
      m_data->plane.to_3d(Point_2(FT(0), FT(0)) + vec));
  }

  const Point_3 to_3d(const Point_2& point) const {
    return m_data->plane.to_3d(point);
  }

  const typename Intersection_kernel::Point_3 to_3d(const typename Intersection_kernel::Point_2& point) const {
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
/*

  const Vertex_index duplicate_vertex(const Vertex_index& v) {
    // TODO: We cannot take it by reference because it fails for EPECK
    // when called from  front_and_back_34() in Data_structure.
    const auto pp = m_data->mesh.point(v);
    const auto vi = m_data->mesh.add_vertex(pp);
    m_data->direction[vi] = m_data->direction[v];
    m_data->v_ivertex_map[vi] = m_data->v_ivertex_map[v];
    m_data->v_iedge_map[vi] = m_data->v_iedge_map[v];
    return vi;
  }*/

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

template<typename GeomTraits, typename IntersectionKernel>
bool operator==(const Support_plane<GeomTraits, IntersectionKernel>& a, const Support_plane<GeomTraits, IntersectionKernel>& b) {

  if (a.is_bbox() || b.is_bbox()) {
    return false;
  }

  using FT = typename GeomTraits::FT;
  const auto& planea = a.plane();
  const auto& planeb = b.plane();

  const auto va = planea.orthogonal_vector();
  const auto vb = planeb.orthogonal_vector();

  // Are the planes parallel?

  FT aval = approximate_angle(va, vb);
  CGAL_assertion(aval >= FT(0) && aval <= FT(180));
  if (aval >= FT(90))
    aval = FT(180) - aval;

  if (aval >= a.angle_tolerance()) {
    return false;
  }

  const auto pa1 = a.to_3d(a.centroid());
  const auto pb1 = planeb.projection(pa1);
  const auto pb2 = b.to_3d(b.centroid());
  const auto pa2 = planea.projection(pb2);

  const FT bval1 = KSP::internal::distance(pa1, pb1);
  const FT bval2 = KSP::internal::distance(pa2, pb2);
  const FT bval = (CGAL::max)(bval1, bval2);
  CGAL_assertion(bval >= FT(0));

  if (bval >= a.distance_tolerance())
    return false;

  return true;
}

#endif //DOXYGEN_RUNNING

} // namespace internal
} // namespace KSP_3
} // namespace CGAL

#endif // CGAL_KSP_3_SUPPORT_LINE_H
