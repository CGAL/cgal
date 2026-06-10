// Copyright (c) 2026 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Efi Fogel         <efifogel@gmail.com>

#ifndef CGAL_ARRANGEMENT_ON_CURVE_1_UNBOUNDED_TOPOLOGY_TRAITS_H
#define CGAL_ARRANGEMENT_ON_CURVE_1_UNBOUNDED_TOPOLOGY_TRAITS_H

#include <list>
#include <type_traits>

#include <boost/property_map/property_map.hpp>
#include <boost/iterator/iterator_adaptor.hpp>
#include <boost/range/iterator_range.hpp>

namespace CGAL {
namespace Arrangement_on_curve_1 {

// ============================================================================
// DATA CONTAINERS (Leveraging Empty Class Optimization)
// ============================================================================

template <typename Data>
struct Data_container {
  Data m_data{};
  Data_container() = default;
  Data_container(const Data& data) : m_data(data) {}
  Data& data() { return m_data; }
  const Data& data() const { return m_data; }
  void set_data(const Data& data) { m_data = data; }
};

// Empty base class optimization (0 bytes)
template <>
struct Data_container<void> {};

// Helper traits to prevent invalid reference-to-void compilation errors
template <typename T> struct Map_reference_traits { using type = T&; };
template <> struct Map_reference_traits<void> { using type = void; };

template <typename T> struct Map_const_reference_traits { using type = const T&; };
template <> struct Map_const_reference_traits<void> { using type = void; };

template <typename T> struct Map_param_traits { using type = const T&; };
template <> struct Map_param_traits<void> { using type = int; }; // Safe dummy fallback type

// ============================================================================
// UNBOUNDED TOPOLOGY TRAITS
// ============================================================================

template <typename Point_1, typename VertexData = void, typename EdgeData = void>
class Unbounded_topology_traits {
public:
  struct Vertex;
  struct Edge;

  using Vertex_list = std::list<Vertex>;
  using Edge_list = std::list<Edge>;

  using Vertex_descriptor = typename Vertex_list::iterator;
  using Edge_descriptor = typename Edge_list::iterator;
  using Vertex_const_descriptor = typename Vertex_list::const_iterator;
  using Edge_const_descriptor = typename Edge_list::const_iterator;

#if 0
  // ============================================================================
  // DESCRIPTOR ITERATOR
  //
  // An iterator adaptor whose operator* returns the underlying list
  // const_iterator (i.e. the descriptor) rather than the list node itself.
  // This allows vertices() / edges() to yield descriptor values that callers
  // can pass directly to property maps, left_vertex(), right_edge(), etc.
  // ============================================================================
  template <typename ListConstIterator>
  class Descriptor_iterator {
  public:
    using iterator_category = std::bidirectional_iterator_tag;
    using value_type = ListConstIterator;
    using difference_type = std::ptrdiff_t;
    using pointer = const ListConstIterator*;
    using reference = ListConstIterator; // returns the descriptor by value

    Descriptor_iterator() = default;
    explicit Descriptor_iterator(ListConstIterator it) : m_it(it) {}

    // Dereferencing yields the const_iterator itself — the descriptor.
    reference operator*() const { return m_it; }
    pointer operator->() const { m_tmp = m_it; return &m_tmp; }

    Descriptor_iterator& operator++() { ++m_it; return *this; }
    Descriptor_iterator operator++(int) { auto tmp = *this; ++m_it; return tmp; }
    Descriptor_iterator& operator--() { --m_it; return *this; }
    Descriptor_iterator operator--(int) { auto tmp = *this; --m_it; return tmp; }

    bool operator==(const Descriptor_iterator& o) const { return m_it == o.m_it; }
    bool operator!=(const Descriptor_iterator& o) const { return m_it != o.m_it; }

  private:
    ListConstIterator m_it;
    mutable ListConstIterator m_tmp; // backing store for operator->
  };
#else
  template <typename ListConstIterator>
  class Descriptor_iterator :
    public boost::iterator_adaptor<Descriptor_iterator<ListConstIterator>, // Derived class (CRTP)
                                   ListConstIterator,          // Base iterator to wrap
                                   ListConstIterator,          // value_type (the descriptor itself)
                                   boost::use_default,         // iterator_category (deduced as bidirectional)
                                   ListConstIterator> {        // reference type (forces operator* to return by value)
    friend class boost::iterator_core_access;

  public:
    Descriptor_iterator() = default;
    explicit Descriptor_iterator(ListConstIterator it) : Descriptor_iterator::iterator_adaptor_(it) {}

  private:
    // Direct core access hook: intercepts dereferencing.
    // boost::iterator_adaptor automatically implements operator* and operator-> using this.
    ListConstIterator dereference() const { return this->base(); }
  };
#endif

  using Vertex_descriptor_iterator = Descriptor_iterator<Vertex_const_descriptor>;
  using Edge_descriptor_iterator = Descriptor_iterator<Edge_const_descriptor>;

  // Range types returned by vertices() and edges()
  using Vertex_descriptor_range = boost::iterator_range<Vertex_descriptor_iterator>;
  using Edge_descriptor_range = boost::iterator_range<Edge_descriptor_iterator>;

  struct Vertex : public Data_container<VertexData> {
    Vertex(const Point_1& p) : m_point(p) {}
    const Point_1& point() const { return m_point; }

    Edge_descriptor left() const { return m_left; }
    Edge_descriptor right() const { return m_right; }

    Point_1 m_point;
    Edge_descriptor m_left;   // edge immediately to the left of this vertex
    Edge_descriptor m_right;  // edge immediately to the right of this vertex
  };

  struct Edge : public Data_container<EdgeData> {
    Vertex_descriptor left() const { return m_left; }
    Vertex_descriptor right() const { return m_right; }
    bool has_left() const { return m_has_left; }
    bool has_right() const { return m_has_right; }

    Vertex_descriptor m_left;
    Vertex_descriptor m_right;
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
    using key_const_type = Vertex_const_descriptor;
    using value_type = Point_1;
    using reference = const Point_1&;
    using category = boost::lvalue_property_map_tag;

    friend reference get(const Vertex_point_map&, key_const_type k) { return k->m_point; }
    friend void put(const Vertex_point_map&, key_type k, const Point_1& p) { k->m_point = p; }
  };

  // 2. Vertex User Data Map
  class Vertex_data_map {
  public:
    using key_const_type = Vertex_const_descriptor;
    using key_type = Vertex_descriptor;
    using value_type = VertexData;
    using reference = typename Map_reference_traits<VertexData>::type;
    using const_reference = typename Map_const_reference_traits<VertexData>::type;
    using category = boost::lvalue_property_map_tag;

    friend reference get(const Vertex_data_map&, key_type k) {
      if constexpr (std::is_void_v<VertexData>) return;
      else return k->data();
    }

    friend const_reference get(const Vertex_data_map&, key_const_type k) {
      if constexpr (std::is_void_v<VertexData>) return;
      else return k->data();
    }

    friend void put(const Vertex_data_map&, key_type k, typename Map_param_traits<VertexData>::type val)
    { if constexpr (! std::is_void_v<VertexData>) k->set_data(val); }
  };

  // 3. Edge User Data Map
  class Edge_data_map {
  public:
    using key_const_type = Edge_const_descriptor;
    using key_type = Edge_descriptor;
    using value_type = EdgeData;
    using reference = typename Map_reference_traits<EdgeData>::type;
    using const_reference = typename Map_const_reference_traits<EdgeData>::type;
    using category = boost::lvalue_property_map_tag;

    friend reference get(const Edge_data_map&, key_type k) {
      if constexpr (std::is_void_v<EdgeData>) return;
      else return k->data();
    }

    friend const_reference get(const Edge_data_map&, key_const_type k) {
      if constexpr (std::is_void_v<EdgeData>) return;
      else return k->data();
    }

    friend void put(const Edge_data_map&, key_type k, typename Map_param_traits<EdgeData>::type val)
    { if constexpr (! std::is_void_v<EdgeData>) k->set_data(val); }
  };

private:
  Vertex_list m_vertices;
  Edge_list m_edges;
  Edge_descriptor m_unbounded_left_edge;
  Edge_descriptor m_unbounded_right_edge;
public:
  // Default constructor
  Unbounded_topology_traits() {
    m_edges.emplace_back();
    m_unbounded_left_edge = m_edges.begin();
    m_unbounded_right_edge = m_edges.begin();
  }

  // ============================================================================
  // QUERIES
  // ============================================================================

  bool is_empty() const { return m_vertices.empty(); }
  size_t number_of_vertices() const { return m_vertices.size(); }
  size_t number_of_edges() const { return m_edges.size(); }

  Vertex_point_map vertex_point_map() const { return Vertex_point_map(); }
  Vertex_data_map vertex_data_map() const { return Vertex_data_map(); }
  Edge_data_map edge_data_map() const { return Edge_data_map(); }

  // Returns const references to the sorted vertex/edge lists.
  const Vertex_list& raw_vertices() const { return m_vertices; }
  const Edge_list& raw_edges() const { return m_edges; }

  // Returns a range of vertex descriptors (Vertex_const_descriptor values).
  // Iterating over this range yields descriptors usable as keys into
  // vertex_point_map(), vertex_data_map(), left_edge(), right_edge(), etc.
  Vertex_descriptor_range vertices() const {
    return boost::make_iterator_range(Vertex_descriptor_iterator(m_vertices.cbegin()),
                                      Vertex_descriptor_iterator(m_vertices.cend()));
  }

  // Returns a range of edge descriptors (Edge_const_descriptor values).
  // Iterating over this range yields descriptors usable as keys into
  // edge_data_map(), has_left_vertex(), left_vertex(), right_vertex(), etc.
  Edge_descriptor_range edges() const {
    return boost::make_iterator_range(Edge_descriptor_iterator(m_edges.cbegin()),
                                      Edge_descriptor_iterator(m_edges.cend()));
  }

  // Mutable iterators over the vertex list
  Vertex_descriptor vertices_begin() { return m_vertices.begin(); }
  Vertex_descriptor vertices_end() { return m_vertices.end(); }

  // UNBOUNDED EDGE ACCESSORS
  Edge_descriptor unbounded_left_edge() { return m_unbounded_left_edge; }
  Edge_const_descriptor unbounded_left_edge() const { return m_unbounded_left_edge; }
  Edge_descriptor unbounded_right_edge() { return m_unbounded_right_edge; }
  Edge_const_descriptor unbounded_right_edge() const { return m_unbounded_right_edge; }

  Edge_descriptor left_edge(Vertex_descriptor v) { return v->m_left; }
  Edge_descriptor right_edge(Vertex_descriptor v) { return v->m_right; }
  Vertex_descriptor left_vertex(Edge_descriptor e) { return e->m_left; }
  Vertex_descriptor right_vertex(Edge_descriptor e) { return e->m_right; }

  Edge_const_descriptor left_edge(Vertex_const_descriptor v) const { return v->m_left; }
  Edge_const_descriptor right_edge(Vertex_const_descriptor v) const { return v->m_right; }
  Vertex_const_descriptor left_vertex(Edge_const_descriptor e) const { return e->m_left; }
  Vertex_const_descriptor right_vertex(Edge_const_descriptor e) const { return e->m_right; }

  bool has_left_vertex(Edge_const_descriptor e) const { return e->m_has_left; }
  bool has_right_vertex(Edge_const_descriptor e) const { return e->m_has_right; }

  // ============================================================================
  // LOW-LEVEL STORAGE PRIMITIVES
  // ============================================================================

  // Append a new vertex to the back of the list (used for the first/rightmost insert).
  Vertex_descriptor create_vertex(const Point_1& p) {
    m_vertices.emplace_back(p);
    return std::prev(m_vertices.end());
  }

  // Append a new edge to the back of the list
  Edge_descriptor create_edge() {
    m_edges.emplace_back();
    return std::prev(m_edges.end());
  }

  void set_left_edge(Vertex_descriptor v, Edge_descriptor e) { v->m_left = e; }
  void set_right_edge(Vertex_descriptor v, Edge_descriptor e) { v->m_right = e; }

  void set_left_vertex(Edge_descriptor e, Vertex_descriptor v) { e->m_left = v; e->m_has_left = true; }
  void set_right_vertex(Edge_descriptor e, Vertex_descriptor v) { e->m_right = v; e->m_has_right = true; }

  // Package-level setters for the cached boundaries
  void set_unbounded_left_edge(Edge_descriptor e)  { m_unbounded_left_edge = e; }
  void set_unbounded_right_edge(Edge_descriptor e) { m_unbounded_right_edge = e; }

  void clear_left_vertex(Edge_descriptor e) { e->m_has_left  = false; }
  void clear_right_vertex(Edge_descriptor e) { e->m_has_right = false; }

  void erase_vertex(Vertex_descriptor v) { m_vertices.erase(v); }
  void erase_edge(Edge_descriptor e) { m_edges.erase(e); }
};

} // namespace Arrangement_on_curve_1
} // namespace CGAL

#endif
