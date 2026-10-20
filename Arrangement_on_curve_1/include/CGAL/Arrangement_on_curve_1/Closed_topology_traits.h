// Copyright (c) 2026 Tel-Aviv University (Israel).
// All rights reserved.
// This file is part of CGAL (www.cgal.org).
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
// Author(s): Efi Fogel         <efifogel@gmail.com>

#ifndef CGAL_ARRANGEMENT_ON_CURVE_1_CLOSED_TOPOLOGY_TRAITS_H
#define CGAL_ARRANGEMENT_ON_CURVE_1_CLOSED_TOPOLOGY_TRAITS_H

/*! \file
 * Topology traits for a 1D arrangement on a curve with closed (finite)
 * boundary endpoints.
 *
 * A *closed* boundary means the curve has a well-defined, finite endpoint on
 * that side — for instance an arc [a, b] where a and b are known points.
 * The arrangement stores a *boundary vertex* at each closed end instead of an
 * unbounded sentinel edge.
 *
 * Three combinations are supported via the two tag template parameters:
 *
 *   LeftTag  = Closed_side_tag, RightTag = Closed_side_tag
 *     → curve [a, b]: both endpoints are boundary vertices.
 *
 *   LeftTag  = Closed_side_tag, RightTag = Open_side_tag
 *     → curve [a, +∞): left boundary vertex, right unbounded edge.
 *
 *   LeftTag  = Open_side_tag,  RightTag = Closed_side_tag
 *     → curve (−∞, b]: left unbounded edge, right boundary vertex.
 *
 * Most shared infrastructure comes from `Topology_traits_base`.
 * `Closed_topology_traits` defines its own `Vertex` node struct (extending the
 * base vertex with `m_has_left_edge` / `m_has_right_edge` flags needed at the
 * boundary vertices) and its own `Edge` node (API-compatible with the base).
 * It therefore does NOT inherit from `Topology_traits_base` but re-uses the
 * base's helper types (`Data_container`, `Map_*_traits`, the descriptor
 * iterators, and the three property maps) directly.
 *
 * API additions (compared with Unbounded_topology_traits):
 *
 *   has_left_boundary_vertex()  / left_boundary_vertex()
 *   has_right_boundary_vertex() / right_boundary_vertex()
 *   set_left_boundary_vertex(v) / set_right_boundary_vertex(v)
 *   has_left_edge(v) / has_right_edge(v)   ← vertex-level boundary check
 */

#include <list>
#include <vector>
#include <memory>
#include <type_traits>
#include <unordered_map>

#include <boost/property_map/property_map.hpp>
#include <boost/iterator/iterator_facade.hpp>
#include <boost/iterator/iterator_adaptor.hpp>
#include <boost/range/iterator_range.hpp>

#include <CGAL/Arrangement_on_curve_1/Topology_traits_base.h>

namespace CGAL {
namespace Arrangement_on_curve_1 {

// ============================================================================
// CLOSED TOPOLOGY TRAITS
//
// LeftTag / RightTag: each must be Open_side_tag or Closed_side_tag.
//   (Identified_side_tag is handled by Identified_topology_traits.)
// ============================================================================
template <typename Point_1_,
          typename LeftTag = Closed_side_tag,
          typename RightTag = Closed_side_tag,
          typename VertexData = void,
          typename EdgeData = void,
          bool UseVector = false,
          typename Allocator = std::allocator<char>>
class Closed_topology_traits {
  static_assert(! Side_is_identified<LeftTag>::value && ! Side_is_identified<RightTag>::value,
                "Closed_topology_traits does not support Identified_side_tag. "
                "Use Identified_topology_traits for closed-loop curves.");

public:
  // Re-use the helper meta-types from the base (brought in without inheriting).
  using Point_1 = Point_1_;

  static constexpr bool left_is_closed = Side_is_closed<LeftTag>::value;
  static constexpr bool right_is_closed = Side_is_closed<RightTag>::value;
  static constexpr bool use_vector = UseVector;

  using Allocator_type = Allocator;

  // --------------------------------------------------------------------------
  // Node types
  //
  // Vertex extends the base Vertex by adding has_left_edge / has_right_edge
  // flags; these are false for boundary vertices on a closed side.
  // --------------------------------------------------------------------------
  struct Vertex : public Data_container<VertexData> {
    Vertex() = default;
    explicit Vertex(const Point_1& p) : m_point(p) {}
    const Point_1& point() const { return m_point; }

    Point_1 m_point;

    // Edge descriptors are set once the node types below are complete.
    // Using a placeholder until the full descriptor type is established.
    // (Resolved in the final class body via forward-declaration trick.)
  };
  // Forward-declare Edge so Vertex can hold an Edge descriptor.
  // The real layout is established in the full struct below.

  // We need to forward-declare the container types so the descriptor aliases
  // can be computed — but Vertex/Edge need the descriptor type.  Break the
  // circularity the same way as Topology_traits_base: forward-declare structs,
  // define container/descriptor, then fill in the fields.
  struct Edge;

  using Vertex_allocator = typename std::allocator_traits<Allocator>::template rebind_alloc<Vertex>;
  using Edge_allocator = typename std::allocator_traits<Allocator>::template rebind_alloc<Edge>;

  using Vertex_container =
    std::conditional_t<UseVector, std::vector<Vertex, Vertex_allocator>, std::list<Vertex, Vertex_allocator>>;
  using Edge_container =
    std::conditional_t<UseVector, std::vector<Edge, Edge_allocator>, std::list<Edge, Edge_allocator>>;

  using Vertex_descriptor = std::conditional_t<UseVector, std::size_t, typename Vertex_container::iterator>;
  using Edge_descriptor = std::conditional_t<UseVector, std::size_t, typename Edge_container::iterator>;
  using Vertex_const_descriptor = std::conditional_t<UseVector, std::size_t, typename Vertex_container::const_iterator>;
  using Edge_const_descriptor = std::conditional_t<UseVector, std::size_t, typename Edge_container::const_iterator>;

  using Size = std::size_t;

  // Now complete the Vertex and Edge structs with descriptor fields.
  // (We re-open them via a full specialization trick — instead, we just embed
  //  the descriptor fields directly in the struct definitions above using the
  //  well-known "declare incomplete struct, define container/descriptor, then
  //  complete struct" approach.  Since C++ doesn't allow re-opening a struct
  //  inside a template, we instead duplicate the complete definitions here.)

  // NOTE: Vertex is declared above with only Point_1 and the Data_container base.
  //       We need to add m_left, m_right, m_has_left_edge, m_has_right_edge.
  //       The standard trick is to use a non-type-dependent base or a wrapper.
  //       Here we use the simplest approach: define Vertex/Edge as complete
  //       structs that capture Edge_descriptor/Vertex_descriptor by inheriting
  //       from a helper that already knows the descriptor types.

  // Redefine Vertex and Edge as the real, complete node types.
  // (The partial forward-declaration above is just for the container aliases.)
private:
  // Complete Vertex node — overrides the forward declaration above via a private
  // alias.  C++ requires structs nested in templates to be complete before
  // instantiation of member functions; the forward declaration + redefinition
  // pattern does NOT work in standard C++.  Instead we use a single complete
  // definition from the start, enabled by the fact that the descriptor type
  // (iterator or size_t) does not depend on the node type itself (only on the
  // container, which uses a forward-declared, incomplete node type at the point
  // of `using Vertex_container = ...`; std::list/std::vector accept incomplete
  // element types in the alias declaration — GCC and Clang both support this).
  // The struct body below fills in the descriptor-typed fields.
  struct Vertex_full {
    using ED = Edge_descriptor;
    Vertex_full() = default;
    explicit Vertex_full(const Point_1& p) : m_point(p) {}
    const Point_1& point() const { return m_point; }

    Point_1 m_point;
    ED m_left {};
    ED m_right {};
    bool m_has_left_edge = false;
    bool m_has_right_edge = false;
    // Per-vertex user data (EBO):
    typename std::conditional_t<std::is_void_v<VertexData>, Data_container<void>, Data_container<VertexData>>
      m_data_field;
    auto& data() { return m_data_field.m_data; }
    const auto& data() const { return m_data_field.m_data; }
    void set_data(const VertexData& d) { if constexpr (! std::is_void_v<VertexData>) m_data_field.m_data = d; }
  };

public:
  // -------------------------------------------------------------------------
  // The actual Vertex type used by the containers is the complete one.
  // Re-alias it so the rest of the class can use "Vertex".
  // -------------------------------------------------------------------------

  // We define Vertex and Edge as proper complete structs here.
  // (The forward Vertex / Edge stubs above are replaced by the full versions;
  //  C++ allows re-defining a nested struct in the same class only if the first
  //  "definition" was an incomplete (forward) declaration.  We achieve this by
  //  making the initial `struct Vertex;` / `struct Edge;` forward declarations
  //  and providing the full body here.)
  //
  // Unfortunately the intermediate container/descriptor aliases already use
  // `Vertex` and `Edge` as the element type; std::list and std::vector allow
  // incomplete element types in the alias (not in the object definition), so
  // the alias compiles.  The full struct body must appear before any object of
  // type Vertex_container or Edge_container is actually constructed.

private:
  // The compiler sees only ONE definition of struct Vertex / struct Edge per
  // class scope.  The cleanest way to satisfy the descriptor dependency without
  // re-opening a struct is to embed the full fields in the forward-declared
  // structs directly using a helper template that does not reference Vertex or
  // Edge:

  // *** Because C++ nested structs cannot be forward-declared and later defined
  // in the same class scope *without* already having a complete definition, we
  // take the pragmatic approach: define Vertex and Edge completely at the start
  // of the class using the concrete Edge_descriptor / Vertex_descriptor types
  // that are known BEFORE the struct body thanks to std::conditional_t +
  // incomplete-element-type allowance in container alias declarations.

  // This is not a re-definition issue in practice because `struct Vertex;`
  // above is an incomplete-type forward declaration, and the full definition
  // appears once below.  The intermediate `using Vertex_container = ...` line
  // works because std::list/std::vector store elements by value but accept an
  // incomplete type at alias-declaration time (the element type need only be
  // complete at instantiation time of member functions).
public:
  // *** REAL complete Vertex definition — replaces the stub above. ***
  // Declared again here as a full definition; this is legal because
  // the first `struct Vertex;` was incomplete.

  // (Suppress the "redefinition" note: we simply never provide a body in the
  //  first declaration.  C++ standard §9.2 p2 says an incomplete class
  //  specifier may be followed by a complete one in the same scope.)

  // VERTEX -------------------------------------------------------------------
  struct Vertex : public Data_container<VertexData> {
    Vertex() = default;
    explicit Vertex(const Point_1& p) : m_point(p) {}
    const Point_1& point() const { return m_point; }

    Point_1 m_point;
    Edge_descriptor m_left {};
    Edge_descriptor m_right {};
    // Flags: false for the boundary vertex on a closed side.
    bool m_has_left_edge  = false;
    bool m_has_right_edge = false;
  };

  // EDGE ---------------------------------------------------------------------
  struct Edge : public Data_container<EdgeData> {
    Vertex_descriptor m_left {};
    Vertex_descriptor m_right {};
    bool m_has_left = false;
    bool m_has_right = false;

    Vertex_descriptor left() const { return m_left; }
    Vertex_descriptor right() const { return m_right; }
    bool has_left() const { return m_has_left;  }
    bool has_right() const { return m_has_right; }
  };

  // --------------------------------------------------------------------------
  // Descriptor iterator types (re-use base's inner class templates directly).
  // --------------------------------------------------------------------------
  template <typename ListConstIterator>
  class List_descriptor_iterator : public boost::iterator_adaptor<List_descriptor_iterator<ListConstIterator>,
                                                                  ListConstIterator, ListConstIterator,
                                                                  boost::bidirectional_traversal_tag,
                                                                  ListConstIterator> {
    friend class boost::iterator_core_access;

  public:
    List_descriptor_iterator() = default;
    explicit List_descriptor_iterator(ListConstIterator it) : List_descriptor_iterator::iterator_adaptor_(it) {}

  private:
    ListConstIterator dereference() const { return this->base(); }
  };

  class Index_descriptor_iterator : public boost::iterator_facade<Index_descriptor_iterator,
                                                                  std::size_t, boost::random_access_traversal_tag,
                                                                  std::size_t, std::ptrdiff_t> {
  public:
    Index_descriptor_iterator() = default;
    explicit Index_descriptor_iterator(std::size_t idx) : m_idx(idx) {}

  private:
    friend class boost::iterator_core_access;
    std::size_t dereference()  const { return m_idx; }
    bool equal(const Index_descriptor_iterator& o) const { return m_idx == o.m_idx; }
    void increment() { ++m_idx; }
    void decrement() { --m_idx; }
    void advance(std::ptrdiff_t n) { m_idx += static_cast<std::size_t>(n); }
    std::ptrdiff_t distance_to(const Index_descriptor_iterator& o) const
    { return static_cast<std::ptrdiff_t>(o.m_idx) - static_cast<std::ptrdiff_t>(m_idx); }
    std::size_t m_idx = 0;
  };

  using Vertex_descriptor_iterator = std::conditional_t<UseVector, Index_descriptor_iterator,
                                                        List_descriptor_iterator<Vertex_const_descriptor>>;
  using Edge_descriptor_iterator = std::conditional_t<UseVector, Index_descriptor_iterator,
                                                      List_descriptor_iterator<Edge_const_descriptor>>;

  using Vertex_descriptor_range = boost::iterator_range<Vertex_descriptor_iterator>;
  using Edge_descriptor_range = boost::iterator_range<Edge_descriptor_iterator>;

  // --------------------------------------------------------------------------
  // Property maps (same pattern as Topology_traits_base)
  // --------------------------------------------------------------------------
  class Vertex_point_map {
  public:
    using key_type = Vertex_descriptor;
    using key_const_type = Vertex_const_descriptor;
    using value_type = Point_1;
    using reference = const Point_1&;
    using category = boost::lvalue_property_map_tag;
    explicit Vertex_point_map(Vertex_container* v) : m_verts(v) {}
    friend reference get(const Vertex_point_map& m, key_const_type k) {
      if constexpr (UseVector) return (*m.m_verts)[k].m_point; else return k->m_point;
    }
    friend void put(const Vertex_point_map& m, key_type k, const Point_1& p) {
      if constexpr (UseVector) (*m.m_verts)[k].m_point = p; else k->m_point = p;
    }
  private: Vertex_container* m_verts;
  };

  class Vertex_data_map {
  public:
    using key_const_type = Vertex_const_descriptor;
    using key_type = Vertex_descriptor;
    using value_type = VertexData;
    using reference = typename Map_reference_traits<VertexData>::type;
    using const_reference = typename Map_const_reference_traits<VertexData>::type;
    using category = boost::lvalue_property_map_tag;
    explicit Vertex_data_map(Vertex_container* v) : m_verts(v) {}
    template <typename Key>
    friend decltype(auto) get(const Vertex_data_map& m, Key k) {
      if constexpr (std::is_void_v<VertexData>) return;
      else if constexpr (UseVector) return (*m.m_verts)[k].data();
      else return k->data();
    }
    friend void put(const Vertex_data_map& m, key_type k, typename Map_param_traits<VertexData>::type val) {
      if constexpr (!std::is_void_v<VertexData>) {
        if constexpr (UseVector) (*m.m_verts)[k].set_data(val);
        else k->set_data(val);
      }
    }
  private: Vertex_container* m_verts;
  };

  class Edge_data_map {
  public:
    using key_const_type = Edge_const_descriptor;
    using key_type = Edge_descriptor;
    using value_type = EdgeData;
    using reference = typename Map_reference_traits<EdgeData>::type;
    using const_reference = typename Map_const_reference_traits<EdgeData>::type;
    using category = boost::lvalue_property_map_tag;
    explicit Edge_data_map(Edge_container* e) : m_edges(e) {}
    template <typename Key>
    friend decltype(auto) get(const Edge_data_map& m, Key k) {
      if constexpr (std::is_void_v<EdgeData>) return;
      else if constexpr (UseVector) return (*m.m_edges)[k].data();
      else return k->data();
    }
    friend void put(const Edge_data_map& m, key_type k, typename Map_param_traits<EdgeData>::type val) {
      if constexpr (!std::is_void_v<EdgeData>) {
        if constexpr (UseVector) (*m.m_edges)[k].set_data(val);
        else k->set_data(val);
      }
    }
  private: Edge_container* m_edges;
  };

  Allocator_type get_allocator() const noexcept
  { return Allocator_type(m_vertices.get_allocator()); }

private:
  Vertex_container m_vertices;
  Edge_container m_edges;

  // Left boundary:
  //   open side   → m_unbounded_left_edge holds the leftmost sentinel edge.
  //   closed side → m_left_bnd_vertex holds the boundary vertex.
  Edge_descriptor m_unbounded_left_edge {};
  Vertex_descriptor m_left_bnd_vertex {};
  bool m_has_left_bnd_vertex = false;

  // Right boundary:
  Edge_descriptor m_unbounded_right_edge {};
  Vertex_descriptor m_right_bnd_vertex {};
  bool m_has_right_bnd_vertex = false;

  // -----------------------------------------------------------------------
  // Internal helpers: descriptor → node reference
  // -----------------------------------------------------------------------
  Vertex& v_ref(Vertex_descriptor v) {
    if constexpr (UseVector) return m_vertices[v];
    else return *v;
  }

  const Vertex& v_ref(Vertex_const_descriptor v) const {
    if constexpr (UseVector) return m_vertices[v];
    else return *v;
  }

  Edge&  e_ref(Edge_descriptor e) {
    if constexpr (UseVector) return m_edges[e];
    else return *e;
  }

  const Edge& e_ref(Edge_const_descriptor e) const {
    if constexpr (UseVector) return m_edges[e];
    else return *e;
  }

  // Minimal hash for list iterators.
  struct Iterator_hash {
    template <typename It>
    std::size_t operator()(It it) const noexcept
    { return std::hash<const void*>{}(static_cast<const void*>(&*it)); }
  };

public:
  // --------------------------------------------------------------------------
  // Constructors
  // --------------------------------------------------------------------------

  Closed_topology_traits() { init_open_sentinels(); }

  explicit Closed_topology_traits(const Allocator_type& alloc) :
    m_vertices(Vertex_allocator(alloc)),
    m_edges(Edge_allocator(alloc))
  { init_open_sentinels(); }

  Closed_topology_traits(const Closed_topology_traits& other) :
    m_vertices(other.m_vertices),
    m_edges(other.m_edges),
    m_unbounded_left_edge(other.m_unbounded_left_edge),
    m_unbounded_right_edge(other.m_unbounded_right_edge),
    m_left_bnd_vertex(other.m_left_bnd_vertex),
    m_right_bnd_vertex(other.m_right_bnd_vertex),
    m_has_left_bnd_vertex(other.m_has_left_bnd_vertex),
    m_has_right_bnd_vertex(other.m_has_right_bnd_vertex)
  { if constexpr (! UseVector) patch_copy(other); }

  Closed_topology_traits(Closed_topology_traits&&) noexcept = default;
  Closed_topology_traits& operator=(const Closed_topology_traits& other) {
    if (this != &other) { Closed_topology_traits tmp(other); swap(tmp); }
    return *this;
  }

  Closed_topology_traits& operator=(Closed_topology_traits&&) noexcept = default;
  ~Closed_topology_traits() = default;

  // --------------------------------------------------------------------------
  // Queries
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
    if constexpr (UseVector)
      return boost::make_iterator_range(Index_descriptor_iterator(0), Index_descriptor_iterator(m_vertices.size()));
    else
      return boost::make_iterator_range(
        List_descriptor_iterator<Vertex_const_descriptor>(m_vertices.cbegin()),
        List_descriptor_iterator<Vertex_const_descriptor>(m_vertices.cend()));
  }

  Edge_descriptor_range edges() const {
    if constexpr (UseVector)
      return boost::make_iterator_range(Index_descriptor_iterator(0), Index_descriptor_iterator(m_edges.size()));
    else return boost::make_iterator_range(List_descriptor_iterator<Edge_const_descriptor>(m_edges.cbegin()),
                                           List_descriptor_iterator<Edge_const_descriptor>(m_edges.cend()));
  }

  // --------------------------------------------------------------------------
  // Open-side boundary accessors (only valid when the respective side is open)
  // --------------------------------------------------------------------------
  Edge_descriptor unbounded_left_edge() {
    static_assert(!left_is_closed, "unbounded_left_edge() called on a closed left boundary.");
    return m_unbounded_left_edge;
  }

  Edge_const_descriptor unbounded_left_edge() const {
    static_assert(!left_is_closed, "unbounded_left_edge() called on a closed left boundary.");
    return m_unbounded_left_edge;
  }

  Edge_descriptor unbounded_right_edge() {
    static_assert(!right_is_closed, "unbounded_right_edge() called on a closed right boundary.");
    return m_unbounded_right_edge;
  }

  Edge_const_descriptor unbounded_right_edge() const {
    static_assert(!right_is_closed, "unbounded_right_edge() called on a closed right boundary.");
    return m_unbounded_right_edge;
  }

  void set_unbounded_left_edge (Edge_descriptor e) { m_unbounded_left_edge  = e; }
  void set_unbounded_right_edge(Edge_descriptor e) { m_unbounded_right_edge = e; }

  // --------------------------------------------------------------------------
  // Closed-side boundary vertex accessors
  // --------------------------------------------------------------------------
  bool has_left_boundary_vertex()  const { return m_has_left_bnd_vertex;  }
  bool has_right_boundary_vertex() const { return m_has_right_bnd_vertex; }

  Vertex_descriptor left_boundary_vertex() {
    CGAL_assertion(m_has_left_bnd_vertex);
    return m_left_bnd_vertex;
  }
  Vertex_const_descriptor left_boundary_vertex() const {
    CGAL_assertion(m_has_left_bnd_vertex);
    return m_left_bnd_vertex;
  }
  Vertex_descriptor right_boundary_vertex() {
    CGAL_assertion(m_has_right_bnd_vertex);
    return m_right_bnd_vertex;
  }
  Vertex_const_descriptor right_boundary_vertex() const {
    CGAL_assertion(m_has_right_bnd_vertex);
    return m_right_bnd_vertex;
  }

  void set_left_boundary_vertex (Vertex_descriptor v) { m_left_bnd_vertex  = v; m_has_left_bnd_vertex  = true; }
  void set_right_boundary_vertex(Vertex_descriptor v) { m_right_bnd_vertex = v; m_has_right_bnd_vertex = true; }

  // --------------------------------------------------------------------------
  // Navigation
  // --------------------------------------------------------------------------
  Edge_descriptor left_edge  (Vertex_descriptor v) { return v_ref(v).m_left; }
  Edge_descriptor right_edge (Vertex_descriptor v) { return v_ref(v).m_right; }
  Vertex_descriptor left_vertex(Edge_descriptor e) { return e_ref(e).m_left; }
  Vertex_descriptor right_vertex(Edge_descriptor e) { return e_ref(e).m_right; }
  Edge_const_descriptor left_edge  (Vertex_const_descriptor v) const { return v_ref(v).m_left; }
  Edge_const_descriptor right_edge (Vertex_const_descriptor v) const { return v_ref(v).m_right; }
  Vertex_const_descriptor left_vertex(Edge_const_descriptor e) const { return e_ref(e).m_left; }
  Vertex_const_descriptor right_vertex(Edge_const_descriptor e) const { return e_ref(e).m_right; }

  bool has_left_vertex(Edge_const_descriptor e) const { return e_ref(e).m_has_left; }
  bool has_right_vertex(Edge_const_descriptor e) const { return e_ref(e).m_has_right; }

  // Vertex-level edge presence (false for boundary vertices on a closed side).
  bool has_left_edge(Vertex_const_descriptor v) const { return v_ref(v).m_has_left_edge; }
  bool has_right_edge(Vertex_const_descriptor v) const { return v_ref(v).m_has_right_edge; }

  const Point_1& vertex_point(Vertex_const_descriptor v) const { return v_ref(v).m_point; }

  // --------------------------------------------------------------------------
  // Low-level storage primitives
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

  void set_left_edge(Vertex_descriptor v, Edge_descriptor e) {
    v_ref(v).m_left = e;
    v_ref(v).m_has_left_edge = true;
  }

  void set_right_edge(Vertex_descriptor v, Edge_descriptor e) {
    v_ref(v).m_right = e;
    v_ref(v).m_has_right_edge = true;
  }

  void set_left_vertex (Edge_descriptor e, Vertex_descriptor v) {
    e_ref(e).m_left = v;
    e_ref(e).m_has_left = true;
  }

  void set_right_vertex(Edge_descriptor e, Vertex_descriptor v) {
    e_ref(e).m_right = v;
    e_ref(e).m_has_right = true;
  }

  void clear_left_vertex(Edge_descriptor e) { e_ref(e).m_has_left = false; }
  void clear_right_vertex(Edge_descriptor e) { e_ref(e).m_has_right = false; }

  void erase_vertex(Vertex_descriptor v) {
    if constexpr (!UseVector) { m_vertices.erase(v); return; }
    if (v == m_vertices.size() - 1) { m_vertices.pop_back(); return; }
    m_vertices.erase(m_vertices.begin() + static_cast<std::ptrdiff_t>(v));
    for (auto& edge : m_edges) {
      if (edge.m_has_left  && edge.m_left > v) --edge.m_left;
      if (edge.m_has_right && edge.m_right > v) --edge.m_right;
    }
    if (m_has_left_bnd_vertex && m_left_bnd_vertex > v) --m_left_bnd_vertex;
    if (m_has_right_bnd_vertex && m_right_bnd_vertex > v) --m_right_bnd_vertex;
  }

  void erase_edge(Edge_descriptor e) {
    if constexpr (!UseVector) { m_edges.erase(e); return; }
    if (e == m_edges.size() - 1) {
      m_edges.pop_back();
      if (! left_is_closed && m_unbounded_left_edge == e) --m_unbounded_left_edge;
      if (! right_is_closed && m_unbounded_right_edge == e) --m_unbounded_right_edge;
      return;
    }
    m_edges.erase(m_edges.begin() + static_cast<std::ptrdiff_t>(e));
    for (auto& vert : m_vertices) {
      if (vert.m_has_left_edge && vert.m_left > e) --vert.m_left;
      if (vert.m_has_right_edge && vert.m_right > e) --vert.m_right;
    }
    if (! left_is_closed && m_unbounded_left_edge > e) --m_unbounded_left_edge;
    if (! right_is_closed && m_unbounded_right_edge > e) --m_unbounded_right_edge;
  }

  // --------------------------------------------------------------------------
  // clear / swap
  // --------------------------------------------------------------------------
  void clear() {
    m_vertices.clear();
    m_edges.clear();
    m_has_left_bnd_vertex  = false;
    m_has_right_bnd_vertex = false;
    init_open_sentinels();
  }

  void swap(Closed_topology_traits& other) noexcept {
    std::swap(m_vertices, other.m_vertices);
    std::swap(m_edges, other.m_edges);
    std::swap(m_unbounded_left_edge, other.m_unbounded_left_edge);
    std::swap(m_unbounded_right_edge, other.m_unbounded_right_edge);
    std::swap(m_left_bnd_vertex, other.m_left_bnd_vertex);
    std::swap(m_right_bnd_vertex, other.m_right_bnd_vertex);
    std::swap(m_has_left_bnd_vertex, other.m_has_left_bnd_vertex);
    std::swap(m_has_right_bnd_vertex, other.m_has_right_bnd_vertex);
  }

private:
  void init_open_sentinels() {
    if constexpr (! left_is_closed) {
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
  }

  void patch_copy(const Closed_topology_traits& other) {
    // Build iterator-to-iterator maps, patch vertex/edge cross-refs, then
    // patch the cached open-side sentinel and boundary vertex descriptors.
    using Edge_map_val = std::pair<const Edge_const_descriptor, Edge_descriptor>;
    using Vertex_map_val = std::pair<const Vertex_const_descriptor, Vertex_descriptor>;

    std::unordered_map<Edge_const_descriptor, Edge_descriptor, Iterator_hash> emap;
    std::unordered_map<Vertex_const_descriptor, Vertex_descriptor, Iterator_hash> vmap;
    emap.reserve(m_edges.size());
    vmap.reserve(m_vertices.size());

    {
      auto s = other.m_edges.cbegin();
      auto d = m_edges.begin();
      for (; s != other.m_edges.cend(); ++s, ++d) emap.emplace(s, d);
    }
    {
      auto s = other.m_vertices.cbegin();
      auto d = m_vertices.begin();
      for (; s != other.m_vertices.cend(); ++s, ++d) vmap.emplace(s, d);
    }

    // Patch vertex → edge.
    {
      auto sv = other.m_vertices.cbegin();
      auto dv = m_vertices.begin();
      for (; sv != other.m_vertices.cend(); ++sv, ++dv) {
        if (sv->m_has_left_edge)  dv->m_left  = emap.at(sv->m_left);
        if (sv->m_has_right_edge) dv->m_right = emap.at(sv->m_right);
      }
    }
    // Patch edge → vertex.
    {
      auto se = other.m_edges.cbegin();
      auto de = m_edges.begin();
      for (; se != other.m_edges.cend(); ++se, ++de) {
        if (se->m_has_left)  de->m_left  = vmap.at(se->m_left);
        if (se->m_has_right) de->m_right = vmap.at(se->m_right);
      }
    }

    if (! left_is_closed) m_unbounded_left_edge  = emap.at(other.m_unbounded_left_edge);
    if (! right_is_closed) m_unbounded_right_edge = emap.at(other.m_unbounded_right_edge);
    if (m_has_left_bnd_vertex) m_left_bnd_vertex  = vmap.at(other.m_left_bnd_vertex);
    if (m_has_right_bnd_vertex) m_right_bnd_vertex = vmap.at(other.m_right_bnd_vertex);
  }
};

// ---------------------------------------------------------------------------
// Convenience aliases
// ---------------------------------------------------------------------------

/// Both sides closed: [a, b] arc
template <typename Point_1, typename VertexData = void, typename EdgeData = void,
          bool UseVector = false, typename Allocator = std::allocator<char>>
using Closed_closed_topology_traits =
  Closed_topology_traits<Point_1, Closed_side_tag, Closed_side_tag, VertexData, EdgeData, UseVector, Allocator>;

/// Left closed, right open: [a, +∞)
template <typename Point_1, typename VertexData = void, typename EdgeData = void,
          bool UseVector = false, typename Allocator = std::allocator<char>>
using Closed_open_topology_traits =
  Closed_topology_traits<Point_1, Closed_side_tag, Open_side_tag, VertexData, EdgeData, UseVector, Allocator>;

/// Left open, right closed: (−∞, b]
template <typename Point_1, typename VertexData = void, typename EdgeData = void,
          bool UseVector = false, typename Allocator = std::allocator<char>>
using Open_closed_topology_traits =
  Closed_topology_traits<Point_1, Open_side_tag, Closed_side_tag, VertexData, EdgeData, UseVector, Allocator>;

} // namespace Arrangement_on_curve_1
} // namespace CGAL

#endif
