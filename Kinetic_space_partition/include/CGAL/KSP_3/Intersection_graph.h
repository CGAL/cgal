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

#ifndef CGAL_KSP_3_INTERSECTION_GRAPH_H
#define CGAL_KSP_3_INTERSECTION_GRAPH_H

#include <CGAL/license/Kinetic_space_partition.h>

// Boost includes.
#include <boost/graph/adjacency_list.hpp>

// CGAL includes.
#include <CGAL/Cartesian_converter.h>
#include <CGAL/Polygon_2.h>

// Internal includes.
#include <CGAL/KSP/utils.h>

namespace CGAL {
namespace KSP_3 {
namespace internal {

#ifdef DOXYGEN_RUNNING
#else

template<typename GeomTraits, typename IntersectionKernel>
class Intersection_graph {

public:
  using Kernel = GeomTraits;
  using Intersection_kernel = IntersectionKernel;

  using IkFT = typename Intersection_kernel::FT;
  using Point_2 = typename Intersection_kernel::Point_2;
  using Point_3 = typename Intersection_kernel::Point_3;
  using Segment_3 = typename Intersection_kernel::Segment_3;
  using Line_3 = typename Intersection_kernel::Line_3;
  using Polygon_2 = typename CGAL::Polygon_2<Intersection_kernel>;

  struct Vertex_property {
    Point_3 point;
    Vertex_property() {}
    Vertex_property(const Point_3& point) : point(point) {}
  };

  using Kinetic_interval = std::vector<std::pair<IkFT, IkFT> >;

  struct Edge_property {
    std::size_t line;
    std::size_t order;
    std::map<std::size_t, std::pair<std::size_t, std::size_t> > faces; // For each intersecting support plane there is one pair of adjacent faces (or less if the edge is on the bbox)
    std::set<std::size_t> planes;
    std::map<std::size_t, Kinetic_interval> intervals; // Maps support plane index to the kinetic interval. std::pair<FT, FT> is the barycentric coordinate and intersection time.
    Edge_property() : line(std::size_t(-1)), order(edge_counter++) { }

    Edge_property(const Edge_property& e) = default;

    const Edge_property& operator=(const Edge_property& other) {
      line = other.line;
      faces = other.faces;
      planes = other.planes;
      intervals = other.intervals;

      return *this;
    }
  private:
    static std::size_t edge_counter;
  };

  using Kinetic_interval_iterator = typename std::map<std::size_t, Kinetic_interval>::const_iterator;

  using Graph = boost::adjacency_list<
    boost::setS, boost::vecS, boost::undirectedS,
    Vertex_property, Edge_property>;

  using Vertex_descriptor = typename boost::graph_traits<Graph>::vertex_descriptor;
  using Edge_descriptor = typename boost::graph_traits<Graph>::edge_descriptor;
  using Face_descriptor = std::size_t;

  struct lex {
    bool operator()(const Edge_descriptor& a, const Edge_descriptor& b) const {
      Edge_property* pa = (Edge_property*)a.get_property();
      Edge_property* pb = (Edge_property*)b.get_property();
      return pa->order < pb->order;
    }
  };

  using IEdge_set = std::set<Edge_descriptor, lex>;

  struct Face_property {
    Face_property() : support_plane(static_cast<std::size_t>(-1)), part_of_partition(false) {}
    Face_property(std::size_t support_plane_idx) : support_plane(support_plane_idx), part_of_partition(false) {}
    std::size_t support_plane;
    bool part_of_partition;
    CGAL::Polygon_2<Intersection_kernel> poly;
    std::vector<Point_2> pts;
    std::vector<Edge_descriptor> edges;
    std::vector<Vertex_descriptor> vertices;
    bool is_part(Edge_descriptor a, Edge_descriptor b) {
      std::size_t aidx = std::size_t(-1);
      for (std::size_t i = 0; i < edges.size(); i++) {
        if (edges[i] == a) {
          aidx = i;
          break;
        }
      }

      if (aidx == std::size_t(-1))
        return false;

      if (edges[(aidx + 1) % edges.size()] == b || edges[(aidx + edges.size() - 1) % edges.size()] == b)
        return true;

      return false;
    }
  };

private:
  Graph m_graph;
  std::vector<Line_3> m_lines;
  std::size_t m_nb_lines_on_bbox;
  std::map<Point_3, Vertex_descriptor> m_map_points;
  std::map<std::vector<std::size_t>, Vertex_descriptor> m_map_vertices;
  std::vector<Face_property> m_ifaces;

  std::vector<bool> m_initial_part_of_partition;
  std::vector<std::map<std::size_t, Kinetic_interval> > m_initial_intervals;

public:
  Intersection_graph() :
    m_nb_lines_on_bbox(0)
  { }

  void clear() {
    m_graph.clear();
    m_map_points.clear();
    m_map_vertices.clear();
  }

  std::size_t number_of_vertices() const {
    return static_cast<std::size_t>(boost::num_vertices(m_graph));
  }

  std::size_t number_of_faces() const {
    return m_ifaces.size();
  }

  static Vertex_descriptor null_ivertex() {
    return boost::graph_traits<Graph>::null_vertex();
  }

  static Edge_descriptor null_iedge() {
    return Edge_descriptor(null_ivertex(), null_ivertex(), nullptr);
  }

  static Face_descriptor null_iface() {
    return std::size_t(-1);
  }

  std::size_t add_line(const Line_3& line) {
    m_lines.push_back(line);
    return m_lines.size() - 1;
  }

  std::size_t nb_lines() const { return m_lines.size(); }

  const std::pair<Vertex_descriptor, bool> add_vertex(const Point_3& point) {
    const auto pair = m_map_points.insert(std::make_pair(point, Vertex_descriptor()));
    const auto is_inserted = pair.second;
    if (is_inserted) {
      pair.first->second = boost::add_vertex(m_graph);
      m_graph[pair.first->second].point = point;
    }
    return std::make_pair(pair.first->second, is_inserted);
  }

  const std::pair<Vertex_descriptor, bool> add_vertex(
    const Point_3& point, const std::vector<std::size_t>& intersected_planes) {
    const auto pair = m_map_vertices.insert(std::make_pair(intersected_planes, Vertex_descriptor()));
    const auto is_inserted = pair.second;
    if (is_inserted) {
      pair.first->second = boost::add_vertex(m_graph);
      m_graph[pair.first->second].point = point;
    }
    return std::make_pair(pair.first->second, is_inserted);
  }

  const std::pair<Edge_descriptor, bool> add_edge(
    const Vertex_descriptor& source, const Vertex_descriptor& target,
    const std::size_t support_plane_idx) {
    const auto out = boost::add_edge(source, target, m_graph);
    m_graph[out.first].planes.insert(support_plane_idx);

    return out;
  }

  template<typename IndexContainer>
  const std::pair<Edge_descriptor, bool> add_edge(
    const Vertex_descriptor& source, const Vertex_descriptor& target,
    const IndexContainer& support_planes_idx) {
    const auto out = boost::add_edge(source, target, m_graph);
    for (const auto support_plane_idx : support_planes_idx) {
      m_graph[out.first].planes.insert(support_plane_idx);
    }
    return out;
  }

  const std::pair<Edge_descriptor, bool> add_edge(const Point_3& source, const Point_3& target) {
    return add_edge(add_vertex(source).first, add_vertex(target).first);
  }

  std::size_t add_face(std::size_t support_plane_idx) {
    m_ifaces.push_back(Face_property(support_plane_idx));
    return std::size_t(m_ifaces.size() - 1);
  }

  bool add_face(std::size_t sp_idx, const Edge_descriptor& edge, const Face_descriptor& idx) {
    auto pair = m_graph[edge].faces.insert(std::make_pair(sp_idx, std::pair<Face_descriptor, Face_descriptor>(-1, -1)));
    if (pair.first->second.first == static_cast<std::size_t>(-1)) {
      pair.first->second.first = idx;
      return true;
    }
    else if (pair.first->second.second == static_cast<std::size_t>(-1)) {
      pair.first->second.second = idx;
      return true;
    }
    return false;
  }

  void get_faces(std::size_t sp_idx, const Edge_descriptor& edge, std::pair<Face_descriptor, Face_descriptor>& pair) const {
    auto it = m_graph[edge].faces.find(sp_idx);
    if (it != m_graph[edge].faces.end())
      pair = it->second;
  }

  const Face_property& face(Face_descriptor idx) const {
    CGAL_assertion(idx < m_ifaces.size());
    return m_ifaces[idx];
  }

  Face_property& face(Face_descriptor idx) {
    CGAL_assertion(idx < m_ifaces.size());
    return m_ifaces[idx];
  }

  const Edge_property& edge(Edge_descriptor idx) const {
    return m_graph[idx];
  }

  void set_line(const Edge_descriptor& edge, const std::size_t line_idx) {
    m_graph[edge].line = line_idx;
  }

  std::size_t line(const Edge_descriptor& edge) const {
    return m_graph[edge].line;
  }

  const Line_3& line(std::size_t line_idx) const {
    return m_lines[line_idx];
  }

  bool line_is_on_bbox(std::size_t line_idx) const {
    return line_idx < m_nb_lines_on_bbox;
  }

  bool line_is_bbox_edge(std::size_t line_idx) const {
    return line_idx < 12;
  }

  bool iedge_is_on_bbox(Edge_descriptor e) {
    return line(e) < m_nb_lines_on_bbox;
  }

  void finished_bbox() {
    m_nb_lines_on_bbox = m_lines.size();
  }

  void initialization_done() {
    auto e = edges();
    m_initial_intervals.resize(e.size());

    std::size_t idx = 0;
    for (const auto& edge : e) {
      m_initial_intervals[idx++] = m_graph[edge].intervals;
    }

    m_initial_part_of_partition.resize(m_ifaces.size());
    for (idx = 0; idx < m_ifaces.size(); idx++)
      m_initial_part_of_partition[idx] = m_ifaces[idx].part_of_partition;
  }

  void reset_to_initialization() {
    auto e = edges();
    CGAL_assertion(e.size() == m_initial_intervals.size());
    std::size_t idx = 0;

    for (auto edge : e) {
      m_graph[edge].intervals = m_initial_intervals[idx];
    }

    CGAL_assertion(m_ifaces.size() == m_initial_part_of_partition.size());
    for (idx = 0; idx < m_ifaces.size(); idx++)
      m_ifaces[idx].part_of_partition = m_initial_part_of_partition[idx];
  }

  const std::pair<Edge_descriptor, Edge_descriptor>
    split_edge(const Edge_descriptor& edge, const Vertex_descriptor& vertex) {

    const auto source = boost::source(edge, m_graph);
    const auto target = boost::target(edge, m_graph);
    const auto prop = m_graph[edge];
    boost::remove_edge(edge, m_graph);

    bool is_inserted;
    Edge_descriptor sedge;
    std::tie(sedge, is_inserted) = boost::add_edge(source, vertex, m_graph);
    if (!is_inserted) {
      std::cerr << "WARNING: " << segment_3(edge) << " " << point_3(vertex) << std::endl;
    }
    CGAL_assertion(is_inserted);
    m_graph[sedge] = prop;

    Edge_descriptor tedge;
    std::tie(tedge, is_inserted) = boost::add_edge(vertex, target, m_graph);
    if (!is_inserted) {
      std::cerr << "WARNING: " << segment_3(edge) << " " << point_3(vertex) << std::endl;
    }
    CGAL_assertion(is_inserted);
    m_graph[tedge] = prop;

    return std::make_pair(sedge, tedge);
  }

  decltype(auto) vertices() const {
    return CGAL::make_range(boost::vertices(m_graph));
  }

  decltype(auto) edges() const {
    return CGAL::make_range(boost::edges(m_graph));
  }

  std::vector<Face_descriptor>& faces() {
    return m_ifaces;
  }

  const std::vector<Face_descriptor>& faces() const {
    return m_ifaces;
  }

  const Vertex_descriptor source(const Edge_descriptor& edge) const {
    return boost::source(edge, m_graph);
  }

  const Vertex_descriptor target(const Edge_descriptor& edge) const {
    return boost::target(edge, m_graph);
  }

  bool is_edge(const Vertex_descriptor& source, const Vertex_descriptor& target) const {
    return boost::edge(source, target, m_graph).second;
  }

  const Edge_descriptor edge(const Vertex_descriptor& source, const Vertex_descriptor& target) const {
    return boost::edge(source, target, m_graph).first;
  }

  decltype(auto) incident_edges(const Vertex_descriptor& vertex) const {
    return CGAL::make_range(boost::out_edges(vertex, m_graph));
  }

  const std::set<std::size_t>& intersected_planes(const Edge_descriptor& edge) const {
    return m_graph[edge].planes;
  }

  std::set<std::size_t>& intersected_planes(const Edge_descriptor& edge) {
    return m_graph[edge].planes;
  }

  const std::pair<Kinetic_interval_iterator, Kinetic_interval_iterator> kinetic_intervals(const Edge_descriptor& edge) {
    return std::pair<Kinetic_interval_iterator, Kinetic_interval_iterator>(m_graph[edge].intervals.begin(), m_graph[edge].intervals.end());
  }
  Kinetic_interval& kinetic_interval(const Edge_descriptor& edge, std::size_t sp_idx) {
    return m_graph[edge].intervals[sp_idx];
  }

  const Point_3& point_3(const Vertex_descriptor& vertex) const {
    return m_graph[vertex].point;
  }

  const Segment_3 segment_3(const Edge_descriptor& edge) const {
    return Segment_3(
      m_graph[boost::source(edge, m_graph)].point,
      m_graph[boost::target(edge, m_graph)].point);
  }

  const Line_3 line_3(const Edge_descriptor& edge) const {
    return Line_3(
      m_graph[boost::source(edge, m_graph)].point,
      m_graph[boost::target(edge, m_graph)].point);
  }
};

template<typename GeomTraits, typename IntersectionKernel> std::size_t Intersection_graph<GeomTraits, IntersectionKernel>::Edge_property::edge_counter = 0;

#endif //DOXYGEN_RUNNING

} // namespace internal
} // namespace KSP_3
} // namespace CGAL

#endif // CGAL_KSP_3_INTERSECTION_GRAPH_H
