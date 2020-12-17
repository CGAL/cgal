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

#ifndef CGAL_KSR_3_INTERSECTION_GRAPH_H
#define CGAL_KSR_3_INTERSECTION_GRAPH_H

// #include <CGAL/license/Kinetic_shape_reconstruction.h>

// Boost includes.
#include <boost/graph/adjacency_list.hpp>

// CGAL includes.
#include <CGAL/Cartesian_converter.h>

// Internal includes.
#include <CGAL/KSR/utils.h>

namespace CGAL {
namespace KSR_3 {

template<typename GeomTraits>
class Intersection_graph {

public:
  using Kernel = GeomTraits;

  using FT        = typename Kernel::FT;
  using Point_3   = typename Kernel::Point_3;
  using Segment_3 = typename Kernel::Segment_3;
  using Line_3    = typename Kernel::Line_3;

  struct Vertex_property {
    Point_3 point;
    bool active;
    Vertex_property() :
    active(true)
    { }
    Vertex_property(const Point_3& point) :
    point(point),
    active(true)
    { }
  };

  struct Edge_property {
    KSR::size_t line;
    KSR::Idx_set planes;
    bool active;
    Edge_property() :
    line(KSR::no_element()),
    active(true)
    { }
  };

  using Graph = boost::adjacency_list<
    boost::setS, boost::vecS, boost::undirectedS,
    Vertex_property, Edge_property>;

  using Vertex_descriptor = typename boost::graph_traits<Graph>::vertex_descriptor;
  using Edge_descriptor   = typename boost::graph_traits<Graph>::edge_descriptor;

private:
  Graph m_graph;
  KSR::size_t m_nb_lines;
  std::map<Point_3, Vertex_descriptor> m_map_points;
  std::map<KSR::Idx_vector, Vertex_descriptor> m_map_vertices;
  std::map<Vertex_descriptor, Vertex_descriptor> m_vmap;
  std::map<Edge_descriptor, Edge_descriptor> m_emap;

public:
  Intersection_graph() :
  m_nb_lines(0)
  { }

  void clear() {
    m_graph.clear();
    m_nb_lines = 0;
    m_map_points.clear();
    m_map_vertices.clear();
  }

  const std::size_t number_of_vertices() const {
    return static_cast<std::size_t>(boost::num_vertices(m_graph));
  }

  const std::size_t number_of_edges() const {
    return static_cast<std::size_t>(boost::num_edges(m_graph));
  }

  template<typename IG>
  void convert(IG& ig) {

    using Converter = CGAL::Cartesian_converter<Kernel, typename IG::Kernel>;

    Converter converter;
    ig.set_nb_lines(m_nb_lines);

    const auto vpair = boost::vertices(m_graph);
    const auto vertex_range = CGAL::make_range(vpair);
    for (const auto vertex : vertex_range) {
      const auto vd = boost::add_vertex(ig.graph());
      ig.graph()[vd].point  = converter(m_graph[vertex].point);
      ig.graph()[vd].active = m_graph[vertex].active;
      CGAL_assertion(m_graph[vertex].active);
      m_vmap[vertex] = vd;
    }
    CGAL_assertion(boost::num_vertices(ig.graph()) == boost::num_vertices(m_graph));

    const auto epair = boost::edges(m_graph);
    const auto edge_range = CGAL::make_range(epair);
    for (const auto edge : edge_range) {
      const auto ed = boost::add_edge(
        boost::source(edge, m_graph), boost::target(edge, m_graph), ig.graph()).first;

      CGAL_assertion(m_graph[edge].line >= 0);
      ig.graph()[ed].line   = m_graph[edge].line;

      CGAL_assertion(m_graph[edge].planes.size() >= 1);
      ig.graph()[ed].planes = m_graph[edge].planes;

      CGAL_assertion(m_graph[edge].active);
      ig.graph()[ed].active = m_graph[edge].active;

      m_emap[edge] = ed;
    }
    CGAL_assertion(boost::num_edges(ig.graph()) == boost::num_edges(m_graph));

    // for (const auto& mp : m_map_points) {
    //   ig.mapped_points()[converter(mp.first)] = m_vmap.at(mp.second);
    // }
    // for (const auto& mv : m_map_vertices) {
    //   ig.mapped_vertices()[mv.first] = m_vmap.at(mv.second);
    // }
  }

  const std::map<Vertex_descriptor, Vertex_descriptor>& vmap() const {
    return m_vmap;
  }

  const std::map<Edge_descriptor, Edge_descriptor>& emap() const {
    return m_emap;
  }

  static Vertex_descriptor null_ivertex() {
    return boost::graph_traits<Graph>::null_vertex();
  }

  static Edge_descriptor null_iedge() {
    return Edge_descriptor(null_ivertex(), null_ivertex(), nullptr);
  }

  const KSR::size_t add_line() { return ( m_nb_lines++ ); }
  const KSR::size_t nb_lines() const { return m_nb_lines; }
  void set_nb_lines(const KSR::size_t value) { m_nb_lines = value; }
  Graph& graph() { return m_graph; }

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
    const Point_3& point, const KSR::Idx_vector& intersected_planes) {

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
    const KSR::size_t support_plane_idx) {

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

  const std::pair<Edge_descriptor, bool> add_edge(
    const Point_3& source, const Point_3& target) {
    return add_edge(add_vertex(source).first, add_vertex(target).first);
  }

  void set_line(const Edge_descriptor& edge, const KSR::size_t line_idx) {
    m_graph[edge].line = line_idx;
  }

  const KSR::size_t line(const Edge_descriptor& edge) const { return m_graph[edge].line; }

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

  decltype(auto) vertices() const { return CGAL::make_range(boost::vertices(m_graph)); }
  decltype(auto) edges() const { return CGAL::make_range(boost::edges(m_graph)); }

  const Vertex_descriptor source(const Edge_descriptor& edge) const { return boost::source(edge, m_graph); }
  const Vertex_descriptor target(const Edge_descriptor& edge) const { return boost::target(edge, m_graph); }

  const bool is_edge(const Vertex_descriptor& source, const Vertex_descriptor& target) const {
    return boost::edge(source, target, m_graph).second;
  }

  decltype(auto) incident_edges(const Vertex_descriptor& vertex) const {
    return CGAL::make_range(boost::out_edges(vertex, m_graph));
  }

  const KSR::Idx_set& intersected_planes(const Edge_descriptor& edge) const { return m_graph[edge].planes; }
  KSR::Idx_set& intersected_planes(const Edge_descriptor& edge) { return m_graph[edge].planes; }

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
      m_graph[source(edge, m_graph)].point,
      m_graph[target(edge, m_graph)].point);
  }

  const bool& is_active(const Vertex_descriptor& vertex) const { return m_graph[vertex].active; }
  bool& is_active(const Vertex_descriptor& vertex) { return m_graph[vertex].active; }
  const bool& is_active(const Edge_descriptor& edge) const { return m_graph[edge].active; }
  bool& is_active(const Edge_descriptor& edge) { return m_graph[edge].active; }
};

} // namespace KSR_3
} // namespace CGAL

#endif // CGAL_KSR_3_INTERSECTION_GRAPH_H
