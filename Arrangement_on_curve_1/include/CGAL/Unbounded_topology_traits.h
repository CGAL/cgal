// Copyright (c) 2026 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Efi Fogel         <efifogel@gmail.com>

#ifndef CGAL_UNBOUNDED_TOPOLOGY_TRAITS_H
#define CGAL_UNBOUNDED_TOPOLOGY_TRAITS_H

#include <list>
#include <type_traits>
#include <boost/property_map/property_map.hpp>

namespace CGAL {
namespace Arrangement_on_curve_1 {

// ============================================================================
// DATA CONTAINERS (Leveraging Empty Class Optimization)
// ============================================================================

// Primary template for handling user extension data attributes
template <typename Data>
struct Data_container {
  Data m_data{};
  Data_container() = default;
  Data_container(const Data& data) : m_data(data) {}
  Data& data() { return m_data; }
  const Data& data() const { return m_data; }
};

// Explicit specialization for 'void' to leverage Empty Class Optimization (ECO)
template <>
struct Data_container<void> {
  // Completely empty: takes up 0 bytes as a base class
};

// ============================================================================
// UNBOUNDED TOPOLOGY TRAITS
// ============================================================================

template <typename Point_1, typename VertexData = void, typename EdgeData = void>
class Unbounded_topology_traits {
public:
  struct Vertex;
  struct Edge;

  using Vertex_list = std::list<Vertex>;
  using Edge_list   = std::list<Edge>;

  using Vertex_descriptor = typename Vertex_list::iterator;
  using Edge_descriptor   = typename Edge_list::iterator;
  using Vertex_const_descriptor = typename Vertex_list::const_iterator;
  using Edge_const_descriptor   = typename Edge_list::const_iterator;

  // Topological entities inherit from data containers to achieve zero-overhead when void
  struct Vertex : public Data_container<VertexData> {
    Point_1 m_point;
    Edge_descriptor m_left;
    Edge_descriptor m_right;
    Vertex(const Point_1& p) : m_point(p) {}
  };

  struct Edge : public Data_container<EdgeData> {
    Vertex_descriptor m_left_v;
    Vertex_descriptor m_right_v;
    bool m_has_left = false;
    bool m_has_right = false;
  };

  // ============================================================================
  // PROPERTY MAPS
  // ============================================================================

  // 1. Geometric Point Map
  class Vertex_point_map {
  public:
    using key_type = Vertex_descriptor;
    using value_type = Point_1;
    using reference = const Point_1&;
    using category = boost::lvalue_property_map_tag;

    friend reference get(const Vertex_point_map&, key_type k) { return k->m_point; }
    friend void put(const Vertex_point_map&, key_type k, const Point_1& p) { k->m_point = p; }
  };

  // 2. Vertex User Data Map
  class Vertex_data_map {
  public:
    using key_type = Vertex_descriptor;
    using value_type = VertexData;
    using reference = std::conditional_t<std::is_void_v<VertexData>, void, VertexData&>;
    using category = boost::lvalue_property_map_tag;

    friend reference get(const Vertex_data_map&, key_type k) {
      if constexpr (std::is_void_v<VertexData>) return;
      else return k->data();
    }
    friend void put(const Vertex_data_map&, key_type k, const VertexData& val) {
      if constexpr (! std::is_void_v<VertexData>) k->data() = val;
    }
  };

  // 3. Edge User Data Map
  class Edge_data_map {
  public:
    using key_type = Edge_descriptor;
    using value_type = EdgeData;
    using reference = std::conditional_t<std::is_void_v<EdgeData>, void, EdgeData&>;
    using category = boost::lvalue_property_map_tag;

    friend reference get(const Edge_data_map&, key_type k) {
      if constexpr (std::is_void_v<EdgeData>) return;
      else return k->data();
    }
    friend void put(const Edge_data_map&, key_type k, const EdgeData& val) {
      if constexpr (! std::is_void_v<EdgeData>) k->data() = val;
    }
  };

private:
  Vertex_list m_vertices;
  Edge_list m_edges;

public:
  Unbounded_topology_traits() = default;

  bool is_empty() const { return m_vertices.empty(); }
  size_t number_of_vertices() const { return m_vertices.size(); }
  size_t number_of_edges() const { return m_edges.size(); }

  Vertex_point_map vertex_point_map() const { return Vertex_point_map(); }
  Vertex_data_map  vertex_data_map() const { return Vertex_data_map(); }
  Edge_data_map    edge_data_map() const { return Edge_data_map(); }

  const Vertex_list& vertices() const { return m_vertices; }
  const Edge_list& edges() const { return m_edges; }

  Edge_descriptor left_edge(Vertex_descriptor v) const { return v->m_left; }
  Edge_descriptor right_edge(Vertex_descriptor v) const { return v->m_right; }
  Vertex_descriptor left_vertex(Edge_descriptor e) const { return e->m_left_v; }
  Vertex_descriptor right_vertex(Edge_descriptor e) const { return e->m_right_v; }

  bool has_left_vertex(Edge_descriptor e) const { return e->m_has_left; }
  bool has_right_vertex(Edge_descriptor e) const { return e->m_has_right; }

  // ============================================================================
  // LOW-LEVEL STORAGE PRIMITIVES
  // ============================================================================
  Vertex_descriptor create_vertex(const Point_1& p) {
    m_vertices.emplace_back(p);
    return std::prev(m_vertices.end());
  }

  Vertex_descriptor create_vertex(Vertex_descriptor position, const Point_1& p) {
    return m_vertices.emplace(position, p);
  }

  Vertex_descriptor create_vertex_front(const Point_1& p) {
    m_vertices.emplace_front(p);
    return m_vertices.begin();
  }

  Edge_descriptor create_edge() {
    m_edges.emplace_back();
    return std::prev(m_edges.end());
  }

  void set_left_edge(Vertex_descriptor v, Edge_descriptor e) { v->m_left = e; }
  void set_right_edge(Vertex_descriptor v, Edge_descriptor e) { v->m_right = e; }

  void set_left_vertex(Edge_descriptor e, Vertex_descriptor v) { e->m_left_v = v; e->m_has_left = true; }
  void set_right_vertex(Edge_descriptor e, Vertex_descriptor v) { e->m_right_v = v; e->m_has_right = true; }

  void clear_left_vertex(Edge_descriptor e) { e->m_has_left = false; }
  void clear_right_vertex(Edge_descriptor e) { e->m_has_right = false; }
};

} // namespace Arrangement_on_curve_1
} // namespace CGAL

#endif
