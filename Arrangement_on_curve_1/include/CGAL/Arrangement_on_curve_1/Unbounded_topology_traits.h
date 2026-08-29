// Copyright (c) 2026 Tel-Aviv University (Israel).
// All rights reserved.
// This file is part of CGAL (www.cgal.org).
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
// Author(s): Efi Fogel         <efifogel@gmail.com>

#ifndef CGAL_ARRANGEMENT_ON_CURVE_1_UNBOUNDED_TOPOLOGY_TRAITS_H
#define CGAL_ARRANGEMENT_ON_CURVE_1_UNBOUNDED_TOPOLOGY_TRAITS_H

#include <list>
#include <vector>
#include <memory>
#include <type_traits>
#include <unordered_map>

#include <boost/property_map/property_map.hpp>
#include <boost/iterator/iterator_facade.hpp>
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
//
// Template parameters:
//   Point_1      - the geometric point type stored at each vertex.
//   VertexData   - optional user data attached to each vertex (void = none).
//   EdgeData     - optional user data attached to each edge (void = none).
//   UseVector    - when true, store vertices and edges in std::vector instead
//                  of std::list.  std::vector gives cache-friendly traversal
//                  and enables binary search in locate(), but descriptors
//                  change from iterators (list mode) to integer indices
//                  (vector mode) because vector iterators are invalidated on
//                  every insertion that causes reallocation.
//   Allocator    - allocator for the internal containers.
// ============================================================================
template <typename Point_1, typename VertexData = void, typename EdgeData = void,
          bool UseVector = false, typename Allocator = std::allocator<char>>
class Unbounded_topology_traits {
public:
  struct Vertex;
  struct Edge;

  // --------------------------------------------------------------------------
  // Allocator rebinds
  // --------------------------------------------------------------------------
  using Allocator_type = Allocator;
  using Vertex_allocator = typename std::allocator_traits<Allocator>::template rebind_alloc<Vertex>;
  using Edge_allocator = typename std::allocator_traits<Allocator>::template rebind_alloc<Edge>;

  // --------------------------------------------------------------------------
  // Container selection.
  //
  // UseVector = false (default): std::list
  //   - Iterators are stable across all insertions and most erasures.
  //   - Descriptors are list iterators; cross-references inside nodes are
  //     stored as iterators.
  //   - locate() is O(n) linear scan.
  //
  // UseVector = true: std::vector
  //   - Cache-friendly storage; random-access iterators enable binary search.
  //   - Iterators are invalidated by push_back (reallocation), so descriptors
  //     are plain std::size_t indices instead.  Indices are unaffected by
  //     reallocation.
  //   - erase() is O(n) and shifts subsequent elements; all stored indices
  //     above the erased slot are patched automatically (see erase_vertex /
  //     erase_edge).
  //   - locate() is O(log n) when binary search is enabled.
  // --------------------------------------------------------------------------
  static constexpr bool use_vector = UseVector;

  using Vertex_container =
    std::conditional_t<UseVector, std::vector<Vertex, Vertex_allocator>, std::list<Vertex, Vertex_allocator>>;

  using Edge_container =
    std::conditional_t<UseVector, std::vector<Edge, Edge_allocator>, std::list<Edge, Edge_allocator>>;

  // Descriptor types:
  //   list mode   → stable bidirectional iterators
  //   vector mode → std::size_t indices (immune to reallocation)
  using Vertex_descriptor = std::conditional_t<UseVector, std::size_t, typename Vertex_container::iterator>;
  using Edge_descriptor = std::conditional_t<UseVector, std::size_t, typename Edge_container::iterator>;
  using Vertex_const_descriptor = std::conditional_t<UseVector, std::size_t, typename Vertex_container::const_iterator>;
  using Edge_const_descriptor = std::conditional_t<UseVector, std::size_t, typename Edge_container::const_iterator>;

  using Size = std::size_t;

  // --------------------------------------------------------------------------
  // Vertex and Edge node types.
  // Cross-reference fields (m_left, m_right) use the descriptor type matching
  // the chosen container.
  // --------------------------------------------------------------------------
  struct Vertex : public Data_container<VertexData> {
    Vertex() = default;
    Vertex(const Point_1& p) : m_point(p) {}
    const Point_1& point() const { return m_point; }

    Edge_descriptor left() const { return m_left; }
    Edge_descriptor right() const { return m_right; }

    Point_1 m_point;
    Edge_descriptor m_left {};
    Edge_descriptor m_right {};
  };

  struct Edge : public Data_container<EdgeData> {
    Vertex_descriptor left() const { return m_left; }
    Vertex_descriptor right() const { return m_right; }
    bool has_left() const { return m_has_left; }
    bool has_right() const { return m_has_right; }

    Vertex_descriptor m_left {};
    Vertex_descriptor m_right {};
    bool m_has_left = false;
    bool m_has_right = false;
  };

  // --------------------------------------------------------------------------
  // DESCRIPTOR ITERATORS
  //
  // vertices() and edges() return ranges whose value_type IS the descriptor
  // (not a reference to the node), so callers can pass the yielded values
  // directly to property maps, left_edge(), right_vertex(), etc.
  //
  // List mode  : wrap list::const_iterator in a boost::iterator_adaptor whose
  //              operator* returns the iterator itself (the descriptor).
  //
  // Vector mode: expose a simple random-access counting iterator from 0 to
  //              container.size(); dereferencing it yields the index.
  // --------------------------------------------------------------------------

  // --- List-mode: iterator-as-descriptor adaptor ---
  template <typename ListConstIterator>
  class List_descriptor_iterator :
    public boost::iterator_adaptor<List_descriptor_iterator<ListConstIterator>, ListConstIterator, ListConstIterator,
                                   boost::bidirectional_traversal_tag, ListConstIterator> {
    friend class boost::iterator_core_access;
  public:
    List_descriptor_iterator() = default;
    explicit List_descriptor_iterator(ListConstIterator it) : List_descriptor_iterator::iterator_adaptor_(it) {}

  private:
    ListConstIterator dereference() const { return this->base(); }
  };

  // --- Vector-mode: index counting iterator ---
  class Index_descriptor_iterator :
    public boost::iterator_facade<Index_descriptor_iterator,
                                  std::size_t,                          // value type
                                  boost::random_access_traversal_tag,   // traversal category
                                  std::size_t,                          // reference type (returned by value)
                                  std::ptrdiff_t                        // difference type
                                  > {
  public:
    Index_descriptor_iterator() = default;
    explicit Index_descriptor_iterator(std::size_t idx) : m_idx(idx) {}

  private:
    friend class boost::iterator_core_access;

    std::size_t dereference() const { return m_idx; }
    bool equal(const Index_descriptor_iterator& other) const { return m_idx == other.m_idx; }
    void increment() { ++m_idx; }
    void decrement() { --m_idx; }
    void advance(std::ptrdiff_t n) { m_idx += static_cast<std::size_t>(n); }

    std::ptrdiff_t distance_to(const Index_descriptor_iterator& other) const
    { return static_cast<std::ptrdiff_t>(other.m_idx) - static_cast<std::ptrdiff_t>(m_idx); }

  private:
    std::size_t m_idx = 0;
  };

  using Vertex_descriptor_iterator =
    std::conditional_t<UseVector, Index_descriptor_iterator, List_descriptor_iterator<Vertex_const_descriptor>>;

  using Edge_descriptor_iterator =
    std::conditional_t<UseVector, Index_descriptor_iterator, List_descriptor_iterator<Edge_const_descriptor>>;

  using Vertex_descriptor_range = boost::iterator_range<Vertex_descriptor_iterator>;
  using Edge_descriptor_range = boost::iterator_range<Edge_descriptor_iterator>;

  // --------------------------------------------------------------------------
  // PROPERTY MAPS
  //
  // The maps now carry a non-owning pointer to the container so that in vector
  // mode they can index into it.  In list mode the pointer is unused but
  // harmless (the iterator already encodes the node address).
  // --------------------------------------------------------------------------

  // 1. Geometric Point Map
  class Vertex_point_map {
  public:
    using key_type = Vertex_descriptor;
    using key_const_type = Vertex_const_descriptor;
    using value_type = Point_1;
    using reference = const Point_1&;
    using category = boost::lvalue_property_map_tag;

    explicit Vertex_point_map(Vertex_container* verts) : m_verts(verts) {}

    friend reference get(const Vertex_point_map& m, key_const_type k) {
      if constexpr (UseVector) return (*m.m_verts)[k].m_point;
      else return k->m_point;
    }
    friend void put(const Vertex_point_map& m, key_type k, const Point_1& p) {
      if constexpr (UseVector) (*m.m_verts)[k].m_point = p;
      else k->m_point = p;
    }
  private:
    Vertex_container* m_verts;
  };

  // 2. Vertex User Data Map
  class Vertex_data_map {
  public:
    using key_const_type  = Vertex_const_descriptor;
    using key_type = Vertex_descriptor;
    using value_type = VertexData;
    using reference = typename Map_reference_traits<VertexData>::type;
    using const_reference = typename Map_const_reference_traits<VertexData>::type;
    using category = boost::lvalue_property_map_tag;

    explicit Vertex_data_map(Vertex_container* verts) : m_verts(verts) {}

    /*! Single templated get() handles BOTH key_type and key_const_type
     * (iterator, const_iterator, or std::size_t) without overload conflicts!
     */
    template <typename Key>
    friend decltype(auto) get(const Vertex_data_map& m, Key k) {
      if constexpr (std::is_void_v<VertexData>) return;
      else if constexpr (UseVector) return (*m.m_verts)[k].data();
      else return k->data();
    }

    /*!
     */
    friend void put(const Vertex_data_map& m, key_type k, typename Map_param_traits<VertexData>::type val) {
      if constexpr (std::is_void_v<VertexData>) return;
      else if constexpr (UseVector) (*m.m_verts)[k].set_data(val);
      else k->set_data(val);
    }
  private:
    Vertex_container* m_verts;
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

    explicit Edge_data_map(Edge_container* edges) : m_edges(edges) {}

    /*! Single templated get() handles BOTH key_type and key_const_type
     * (iterator, const_iterator, or std::size_t) without overload conflicts!
     */
    template <typename Key>
    friend decltype(auto) get(const Edge_data_map& m, Key k) {
      if constexpr (std::is_void_v<EdgeData>) return;
      else if constexpr (UseVector) return (*m.m_edges)[k].data();
      else return k->data();
    }

    /*!
     */
    friend void put(const Edge_data_map& m, key_type k, typename Map_param_traits<EdgeData>::type val) {
      if constexpr (! std::is_void_v<EdgeData>) {
        if constexpr (UseVector) (*m.m_edges)[k].set_data(val);
        else k->set_data(val);
      }
    }
  private:
    Edge_container* m_edges;
  };

  // --------------------------------------------------------------------------
  // get_allocator()
  // --------------------------------------------------------------------------
  // m_vertices.get_allocator() returns Vertex_allocator (= rebind_alloc<Vertex>),
  // not Allocator_type.  Convert back via the cross-rebind converting
  // constructor, which is always available for the Allocator concept.
  Allocator_type get_allocator() const noexcept { return Allocator_type(m_vertices.get_allocator()); }

private:
  Vertex_container m_vertices;
  Edge_container m_edges;
  Edge_descriptor m_unbounded_left_edge  {};
  Edge_descriptor m_unbounded_right_edge {};

  // --------------------------------------------------------------------------
  // Helper: dereference a descriptor into a mutable or const node reference.
  // Abstracts over the iterator-vs-index difference so all navigation methods
  // can be written once.
  // --------------------------------------------------------------------------
  Vertex& v_ref(Vertex_descriptor v) { if constexpr (UseVector) return m_vertices[v]; else return *v; }
  const Vertex& v_ref(Vertex_const_descriptor v) const { if constexpr (UseVector) return m_vertices[v]; else return *v; }
  Edge& e_ref(Edge_descriptor e) { if constexpr (UseVector) return m_edges[e]; else return *e; }
  const Edge& e_ref(Edge_const_descriptor e) const { if constexpr (UseVector) return m_edges[e]; else return *e; }

public:
  // --------------------------------------------------------------------------
  // CONSTRUCTORS / SPECIAL MEMBERS
  // --------------------------------------------------------------------------

  // Default constructor
  Unbounded_topology_traits() {
    m_edges.emplace_back();
    if constexpr (UseVector) {
      m_unbounded_left_edge = 0;
      m_unbounded_right_edge = 0;
    }
    else {
      m_unbounded_left_edge = m_edges.begin();
      m_unbounded_right_edge = m_edges.begin();
    }
  }

  // Allocator-aware constructor
  explicit Unbounded_topology_traits(const Allocator_type& alloc) :
    m_vertices(Vertex_allocator(alloc)),
    m_edges(Edge_allocator(alloc)) {
    m_edges.emplace_back();
    if constexpr (UseVector) {
      m_unbounded_left_edge = 0;
      m_unbounded_right_edge = 0;
    }
    else {
      m_unbounded_left_edge = m_edges.begin();
      m_unbounded_right_edge = m_edges.begin();
    }
  }

  // --------------------------------------------------------------------------
  // Move constructor and move assignment — compiler-generated is correct.
  //
  // List mode : std::list move transfers heap nodes without reallocation;
  //             stored iterators remain valid in the moved-to object.
  // Vector mode: std::vector move transfers the heap buffer; stored indices
  //              are plain integers and remain valid.
  // --------------------------------------------------------------------------
  Unbounded_topology_traits(Unbounded_topology_traits&&) noexcept = default;

  // --------------------------------------------------------------------------
  // Move assignment operator
  // --------------------------------------------------------------------------
  Unbounded_topology_traits& operator=(Unbounded_topology_traits&&) noexcept = default;

  ~Unbounded_topology_traits() = default;

  // --------------------------------------------------------------------------
  // Copy constructor.
  //
  // List mode : std::list copy produces new nodes in the same order, but all
  //             stored iterators still point into the source's lists.
  //             patch_cross_references() rebuilds them.
  //
  // Vector mode: descriptors are plain std::size_t values; a value copy is
  //              trivially correct.  No patching needed.
  // --------------------------------------------------------------------------
  Unbounded_topology_traits(const Unbounded_topology_traits& other) :
    m_vertices(other.m_vertices),
    m_edges(other.m_edges),
    m_unbounded_left_edge(other.m_unbounded_left_edge),
    m_unbounded_right_edge(other.m_unbounded_right_edge)
  { if constexpr (! UseVector) patch_cross_references(other); }

  // --------------------------------------------------------------------------
  // Copy assignment operator
  // --------------------------------------------------------------------------
  Unbounded_topology_traits& operator=(const Unbounded_topology_traits& other) {
    if (this != &other) {
      Unbounded_topology_traits tmp(other);
      swap(tmp);
    }
    return *this;
  }

  /*! swaps
   */
  void swap(Unbounded_topology_traits& other) noexcept {
    std::swap(m_vertices, other.m_vertices);
    std::swap(m_edges, other.m_edges);
    std::swap(m_unbounded_left_edge, other.m_unbounded_left_edge);
    std::swap(m_unbounded_right_edge, other.m_unbounded_right_edge);
  }

  /*! clears
   */
  void clear() {
    m_vertices.clear();
    m_edges.clear();
    m_edges.emplace_back();
    if constexpr (UseVector) {
      m_unbounded_left_edge = 0;
      m_unbounded_right_edge = 0;
    }
    else {
      m_unbounded_left_edge = m_edges.begin();
      m_unbounded_right_edge = m_edges.begin();
    }
  }

  // --------------------------------------------------------------------------
  // QUERIES
  // --------------------------------------------------------------------------

  bool empty() const { return m_vertices.empty(); }
  Size number_of_vertices() const { return m_vertices.size(); }
  Size number_of_edges() const { return m_edges.size(); }

  Vertex_point_map vertex_point_map() const { return Vertex_point_map(const_cast<Vertex_container*>(&m_vertices)); }
  Vertex_data_map vertex_data_map() const { return Vertex_data_map(const_cast<Vertex_container*>(&m_vertices)); }
  Edge_data_map edge_data_map() const { return Edge_data_map(const_cast<Edge_container*>(&m_edges)); }

  const Vertex_container& raw_vertices() const { return m_vertices; }
  const Edge_container& raw_edges() const { return m_edges; }

  Vertex_descriptor_range vertices() const {
    if constexpr (UseVector) {
      return boost::make_iterator_range(Index_descriptor_iterator(0), Index_descriptor_iterator(m_vertices.size()));
    }
    else {
      return boost::make_iterator_range(List_descriptor_iterator<Vertex_const_descriptor>(m_vertices.cbegin()),
                                        List_descriptor_iterator<Vertex_const_descriptor>(m_vertices.cend()));
    }
  }

  Edge_descriptor_range edges() const {
    if constexpr (UseVector) {
      return boost::make_iterator_range(Index_descriptor_iterator(0), Index_descriptor_iterator(m_edges.size()));
    }
    else {
      return boost::make_iterator_range(List_descriptor_iterator<Edge_const_descriptor>(m_edges.cbegin()),
                                        List_descriptor_iterator<Edge_const_descriptor>(m_edges.cend()));
    }
  }

  // Mutable begin/end for locate_impl (linear scan path)
  Vertex_descriptor vertices_begin() {
    if constexpr (UseVector) return 0;
    else return m_vertices.begin();
  }
  Vertex_descriptor vertices_end() {
    if constexpr (UseVector) return m_vertices.size();
    else return m_vertices.end();
  }

  // --------------------------------------------------------------------------
  // UNBOUNDED EDGE ACCESSORS
  // --------------------------------------------------------------------------
  Edge_descriptor unbounded_left_edge() { return m_unbounded_left_edge;  }
  Edge_const_descriptor unbounded_left_edge() const { return m_unbounded_left_edge;  }
  Edge_descriptor unbounded_right_edge() { return m_unbounded_right_edge; }
  Edge_const_descriptor unbounded_right_edge() const { return m_unbounded_right_edge; }

  // --------------------------------------------------------------------------
  // NAVIGATION
  // --------------------------------------------------------------------------
  Edge_descriptor left_edge(Vertex_descriptor v) { return v_ref(v).m_left; }
  Edge_descriptor right_edge(Vertex_descriptor v) { return v_ref(v).m_right; }
  Vertex_descriptor left_vertex(Edge_descriptor e) { return e_ref(e).m_left; }
  Vertex_descriptor right_vertex(Edge_descriptor e) { return e_ref(e).m_right; }

  Edge_const_descriptor left_edge(Vertex_const_descriptor v) const { return v_ref(v).m_left; }
  Edge_const_descriptor right_edge(Vertex_const_descriptor v) const { return v_ref(v).m_right; }
  Vertex_const_descriptor left_vertex (Edge_const_descriptor e) const { return e_ref(e).m_left; }
  Vertex_const_descriptor right_vertex(Edge_const_descriptor e) const { return e_ref(e).m_right; }

  bool has_left_vertex(Edge_const_descriptor e) const { return e_ref(e).m_has_left;  }
  bool has_right_vertex(Edge_const_descriptor e) const { return e_ref(e).m_has_right; }

  // Direct point access used by binary_search_vertex (avoids property map overhead).
  const Point_1& vertex_point(Vertex_const_descriptor v) const { return v_ref(v).m_point; }

  // --------------------------------------------------------------------------
  // LOW-LEVEL STORAGE PRIMITIVES
  // --------------------------------------------------------------------------

  Vertex_descriptor create_vertex(const Point_1& p) {
    m_vertices.emplace_back(p);
    if constexpr (UseVector) return m_vertices.size() - 1;
    else return std::prev(m_vertices.end());
  }

  Edge_descriptor create_edge() {
    m_edges.emplace_back();
    if constexpr (UseVector) return m_edges.size() - 1;
    else return std::prev(m_edges.end());
  }

  void set_left_edge(Vertex_descriptor v, Edge_descriptor e) { v_ref(v).m_left  = e; }
  void set_right_edge(Vertex_descriptor v, Edge_descriptor e) { v_ref(v).m_right = e; }
  void set_left_vertex(Edge_descriptor e, Vertex_descriptor v) { e_ref(e).m_left  = v; e_ref(e).m_has_left  = true; }
  void set_right_vertex(Edge_descriptor e, Vertex_descriptor v) { e_ref(e).m_right = v; e_ref(e).m_has_right = true; }

  void set_unbounded_left_edge(Edge_descriptor e) { m_unbounded_left_edge  = e; }
  void set_unbounded_right_edge(Edge_descriptor e) { m_unbounded_right_edge = e; }

  void clear_left_vertex(Edge_descriptor e) { e_ref(e).m_has_left  = false; }
  void clear_right_vertex(Edge_descriptor e) { e_ref(e).m_has_right = false; }

  // erase_vertex
  //
  // List mode: O(1). Iterators to other nodes remain valid.
  // Vector mode: O(1) when popping from the back; O(n) when erasing interior slots.
  void erase_vertex(Vertex_descriptor v) {
    if constexpr (! UseVector) {
      m_vertices.erase(v);
      return;
    }

    // Fast path: O(1) removal if erasing the last element
    if (v == m_vertices.size() - 1) {
      m_vertices.pop_back();
      return;
    }

    // Slow path: O(n) interior erase requiring index patching
    m_vertices.erase(m_vertices.begin() + static_cast<std::ptrdiff_t>(v));
    for (auto& edge : m_edges) {
      if (edge.m_has_left  && edge.m_left > v) --edge.m_left;
      if (edge.m_has_right && edge.m_right > v) --edge.m_right;
    }
  }

  // erase_edge
  //
  // List mode: O(1).
  // Vector mode: O(1) when popping from the back; O(n) when erasing interior slots.
  void erase_edge(Edge_descriptor e) {
    if constexpr (!UseVector) {
      m_edges.erase(e);
      return;
    }

    // Fast path: O(1) removal if erasing the last element
    if (e == m_edges.size() - 1) {
      m_edges.pop_back();
      if (m_unbounded_left_edge == e)  --m_unbounded_left_edge;
      if (m_unbounded_right_edge == e) --m_unbounded_right_edge;
      return;
    }

    // Slow path: O(n) interior erase requiring index patching
    m_edges.erase(m_edges.begin() + static_cast<std::ptrdiff_t>(e));
    for (auto& vert : m_vertices) {
      if (vert.m_left  > e) --vert.m_left;
      if (vert.m_right > e) --vert.m_right;
    }
    if (m_unbounded_left_edge > e) --m_unbounded_left_edge;
    if (m_unbounded_right_edge > e) --m_unbounded_right_edge;
  }

private:
  // --------------------------------------------------------------------------
  // patch_cross_references() — list mode only
  // --------------------------------------------------------------------------
  void patch_cross_references(const Unbounded_topology_traits& other) {
    using Edge_map_value = std::pair<const Edge_const_descriptor,   Edge_descriptor>;
    using Vertex_map_value = std::pair<const Vertex_const_descriptor, Vertex_descriptor>;
    using Edge_map_alloc = typename std::allocator_traits<Allocator>::template rebind_alloc<Edge_map_value>;
    using Vertex_map_alloc = typename std::allocator_traits<Allocator>::template rebind_alloc<Vertex_map_value>;

    std::unordered_map<Edge_const_descriptor, Edge_descriptor, Iterator_hash, std::equal_to<Edge_const_descriptor>,
                       Edge_map_alloc>
      edge_map(0, Iterator_hash{}, std::equal_to<Edge_const_descriptor>{}, Edge_map_alloc(get_allocator()));
    edge_map.reserve(m_edges.size());
    {
      auto src = other.m_edges.cbegin();
      auto dst = m_edges.begin();
      for (; src != other.m_edges.cend(); ++src, ++dst)
        edge_map.emplace(src, dst);
    }

    std::unordered_map<Vertex_const_descriptor, Vertex_descriptor, Iterator_hash, std::equal_to<Vertex_const_descriptor>,
                       Vertex_map_alloc>
      vertex_map(0, Iterator_hash{}, std::equal_to<Vertex_const_descriptor>{}, Vertex_map_alloc(get_allocator()));
    vertex_map.reserve(m_vertices.size());
    {
      auto src = other.m_vertices.cbegin();
      auto dst = m_vertices.begin();
      for (; src != other.m_vertices.cend(); ++src, ++dst)
        vertex_map.emplace(src, dst);
    }

    {
      auto src_v = other.m_vertices.cbegin();
      auto dst_v = m_vertices.begin();
      for (; src_v != other.m_vertices.cend(); ++src_v, ++dst_v) {
        dst_v->m_left  = edge_map.at(src_v->m_left);
        dst_v->m_right = edge_map.at(src_v->m_right);
      }
    }
    {
      auto src_e = other.m_edges.cbegin();
      auto dst_e = m_edges.begin();
      for (; src_e != other.m_edges.cend(); ++src_e, ++dst_e) {
        if (src_e->m_has_left)  dst_e->m_left  = vertex_map.at(src_e->m_left);
        if (src_e->m_has_right) dst_e->m_right = vertex_map.at(src_e->m_right);
      }
    }

    // Step 4: patch the cached unbounded-edge descriptors
    m_unbounded_left_edge = edge_map.at(other.m_unbounded_left_edge);
    m_unbounded_right_edge = edge_map.at(other.m_unbounded_right_edge);
  }

  // Minimal hash for iterators: use the address of the node they point to.
  struct Iterator_hash {
    template <typename It>
    std::size_t operator()(It it) const noexcept { return std::hash<const void*>{}(static_cast<const void*>(&*it)); }
  };

  void reserve_vertices(std::size_t n) { if constexpr (UseVector) m_vertices.reserve(n); }
  void reserve_edges(std::size_t n) { if constexpr (UseVector) m_edges.reserve(n); }
};

} // namespace Arrangement_on_curve_1
} // namespace CGAL

#endif
