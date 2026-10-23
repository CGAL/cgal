// Copyright (c) 2026 Tel-Aviv University (Israel).
// All rights reserved.
// This file is part of CGAL (www.cgal.org).
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
// Author(s): Efi Fogel         <efifogel@gmail.com>

#ifndef CGAL_ARRANGEMENT_ON_CURVE_1_TOPOLOGY_TRAITS_BASE_H
#define CGAL_ARRANGEMENT_ON_CURVE_1_TOPOLOGY_TRAITS_BASE_H

/*! \file
 * Base class shared by all topology traits in Arrangement_on_curve_1.
 *
 * Topology_traits_base factors out everything that the five concrete topology
 * traits, `Open_topology_traits`, `Closed_topology_traits`,
 * `Left_closed_topology_traits`, `Right_closed_topology_traits`, and
 * `Identified_topology_traits`, have in common:
 *
 * - Data_container and its void-specialization (empty-base optimization)
 * - Map helper traits (reference, const_reference, parameter)
 * - The Vertex and Edge node structs (with optional user data)
 * - Container and descriptor type definitions
 * - List-mode and vector-mode descriptor iterators
 * - All three property maps (Vertex_point_map, Vertex_data_map, Edge_data_map)
 * - Core navigation (left/right_edge, left/right_vertex, has_left/right_vertex, vertex_point)
 * - The left and right boundary edges (left_boundary_edge, set_left_boundary_edge,
 *   right_boundary_edge, set_right_boundary_edge)
 * - Low-level storage primitives (create_vertex, create_edge, set_*, erase_vertex, erase_edge)
 * - Construction (the spanning edge), copying (list-mode iterators are
 *   redirected into the copy), copy assignment, clear, and swap
 * - Allocator access
 *
 * Each derived class adds only its topology-specific state and logic:
 * - `Open_topology_traits`---nothing
 * - `Closed_topology_traits`---the left and right boundary vertex descriptors
 * - `Left_closed_topology_traits`---the left boundary vertex descriptor
 * - `Right_closed_topology_traits`---the right boundary vertex descriptor
 * - `Identified_topology_traits`---the identification vertex descriptor
 */

#include <list>
#include <vector>
#include <memory>
#include <type_traits>
#include <unordered_map>

#include <CGAL/assertions.h>

#include <boost/property_map/property_map.hpp>
#include <boost/iterator/iterator_facade.hpp>
#include <boost/iterator/iterator_adaptor.hpp>
#include <boost/range/iterator_range.hpp>

namespace CGAL {
namespace Arrangement_on_curve_1 {

// ============================================================================
// DATA CONTAINERS (empty-base optimization)
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

// Specialization for void; zero bytes, no members.
template <>
struct Data_container<void> {};

// Helper traits to avoid reference-to-void compilation errors.
template <typename T> struct Map_reference_traits { using type = T&; };
template <> struct Map_reference_traits<void> { using type = void; };

template <typename T> struct Map_const_reference_traits { using type = const T&; };
template <> struct Map_const_reference_traits<void> { using type = void; };

template <typename T> struct Map_param_traits { using type = const T&; };
template <>  struct Map_param_traits<void> { using type = int; }; // safe dummy

// ============================================================================
// Topology_traits_base
//
// Template parameters (same as the concrete topology traits):
//   Point_1    - geometric point type stored at each vertex.
//   VertexData - optional user data per vertex (void = none).
//   EdgeData   - optional user data per edge   (void = none).
//   UseVector  - true -> std::vector + index descriptors;
//                false -> std::list + iterator descriptors.
//   Allocator  - allocator for internal containers.
// ============================================================================
template <typename Point_1_,
          typename VertexData = void,
          typename EdgeData = void,
          bool UseVector = false,
          typename Allocator_ = std::allocator<char>>
class Topology_traits_base {
public:
  using Point_1 = Point_1_;
  using Vertex_data = VertexData;
  using Edge_data = EdgeData;
  using Allocator = Allocator_;

  //! The standard stl type.
  using allocator_type = Allocator;

  static constexpr bool use_vector = UseVector;

  // Forward-declare cell types so container/descriptor aliases can reference them.
  struct Vertex;
  struct Edge;

protected:
  // --------------------------------------------------------------------------
  // Allocator rebinds
  // --------------------------------------------------------------------------
  using Vertex_allocator = typename std::allocator_traits<Allocator>::template rebind_alloc<Vertex>;
  using Edge_allocator = typename std::allocator_traits<Allocator>::template rebind_alloc<Edge>;

public:
  // --------------------------------------------------------------------------
  // Container and descriptor types
  //
  // UseVector = false (default): std::list
  //   Descriptors are stable bidirectional iterators.
  //
  // UseVector = true: std::vector
  //   Descriptors are plain std::size_t indices, immune to reallocation.
  //   erase() patches all stored indices automatically.
  // --------------------------------------------------------------------------
  using Vertex_container = std::conditional_t<UseVector,
                                              std::vector<Vertex, Vertex_allocator>,
                                              std::list  <Vertex, Vertex_allocator>>;

  using Edge_container = std::conditional_t<UseVector,
                                            std::vector<Edge, Edge_allocator>,
                                            std::list  <Edge, Edge_allocator>>;

  using Vertex_descriptor = std::conditional_t<UseVector, std::size_t, typename Vertex_container::iterator>;
  using Edge_descriptor = std::conditional_t<UseVector, std::size_t, typename Edge_container::iterator>;
  using Vertex_const_descriptor = std::conditional_t<UseVector, std::size_t, typename Vertex_container::const_iterator>;
  using Edge_const_descriptor = std::conditional_t<UseVector, std::size_t, typename Edge_container::const_iterator>;

  using Size = std::size_t;

  // --------------------------------------------------------------------------
  // Vertex
  // --------------------------------------------------------------------------
  struct Vertex : public Data_container<Vertex_data> {
    Vertex() = default;
    Vertex(const Point_1& p) : m_point(p) {}
    const Point_1& point() const { return m_point; }

    Edge_descriptor left() const { return m_left; }
    Edge_descriptor right() const { return m_right; }

    Point_1 m_point;
    Edge_descriptor m_left {};
    Edge_descriptor m_right {};
  };

  // --------------------------------------------------------------------------
  // Edge
  // --------------------------------------------------------------------------
  struct Edge : public Data_container<Edge_data> {
    Vertex_descriptor left() const { return m_left; }
    Vertex_descriptor right() const { return m_right; }
    bool has_left() const { return m_has_left;  }
    bool has_right() const { return m_has_right; }

    Vertex_descriptor m_left {};
    Vertex_descriptor m_right {};
    bool m_has_left  = false;
    bool m_has_right = false;
  };

  // --------------------------------------------------------------------------
  // Descriptor iterators
  //
  // vertices() / edges() return ranges whose value_type IS the descriptor,
  // so callers can pass yielded values directly to property maps and navigation.
  // --------------------------------------------------------------------------

  // List-mode: wraps list::const_iterator; operator* returns the iterator itself.
  template <typename ListConstIterator> class List_descriptor_iterator :
    public boost::iterator_adaptor<List_descriptor_iterator<ListConstIterator>,
                                   ListConstIterator,
                                   ListConstIterator,
                                   boost::bidirectional_traversal_tag,
                                   ListConstIterator> {
    friend class boost::iterator_core_access;

  public:
    List_descriptor_iterator() = default;
    explicit List_descriptor_iterator(ListConstIterator it) : List_descriptor_iterator::iterator_adaptor_(it) {}

  private:
    ListConstIterator dereference() const { return this->base(); }
  };

  // Vector-mode: counting iterator from 0 to container.size().
  class Index_descriptor_iterator :
    public boost::iterator_facade<Index_descriptor_iterator,
                                  std::size_t,
                                  boost::random_access_traversal_tag,
                                  std::size_t,
                                  std::ptrdiff_t> {
    friend class boost::iterator_core_access;

  public:
    Index_descriptor_iterator() = default;
    explicit Index_descriptor_iterator(std::size_t idx) : m_idx(idx) {}

  private:
    std::size_t dereference()  const { return m_idx; }
    bool equal(const Index_descriptor_iterator& o) const { return m_idx == o.m_idx; }
    void increment() { ++m_idx; }
    void decrement() { --m_idx; }
    void advance(std::ptrdiff_t n) { m_idx += static_cast<std::size_t>(n); }

    std::ptrdiff_t distance_to(const Index_descriptor_iterator& o) const
    { return static_cast<std::ptrdiff_t>(o.m_idx) - static_cast<std::ptrdiff_t>(m_idx); }

    std::size_t m_idx = 0;
  };

  using Vertex_descriptor_iterator = std::conditional_t<UseVector,
                                                        Index_descriptor_iterator,
                                                        List_descriptor_iterator<Vertex_const_descriptor>>;

  using Edge_descriptor_iterator = std::conditional_t<UseVector,
                                                      Index_descriptor_iterator,
                                                      List_descriptor_iterator<Edge_const_descriptor>>;

  using Vertex_descriptor_range = boost::iterator_range<Vertex_descriptor_iterator>;
  using Edge_descriptor_range = boost::iterator_range<Edge_descriptor_iterator>;

  // --------------------------------------------------------------------------
  // Property maps
  // --------------------------------------------------------------------------

  // 1. Geometric point map
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

  // 2. Vertex user-data map
  class Vertex_data_map {
  public:
    using key_type = Vertex_descriptor;
    using key_const_type  = Vertex_const_descriptor;
    using value_type = Vertex_data;
    using reference = typename Map_reference_traits<Vertex_data>::type;
    using const_reference = typename Map_const_reference_traits<Vertex_data>::type;
    using category = boost::lvalue_property_map_tag;

    explicit Vertex_data_map(Vertex_container* verts) : m_verts(verts) {}

    template <typename Key>
    friend decltype(auto) get(const Vertex_data_map& m, Key k) {
      if constexpr (std::is_void_v<Vertex_data>) return;
      else if constexpr (UseVector) return (*m.m_verts)[k].data();
      else return k->data();
    }
    friend void put(const Vertex_data_map& m, key_type k,
                    typename Map_param_traits<Vertex_data>::type val) {
      if constexpr (!std::is_void_v<Vertex_data>) {
        if constexpr (UseVector) (*m.m_verts)[k].set_data(val);
        else k->set_data(val);
      }
    }

  private:
    Vertex_container* m_verts;
  };

  // 3. Edge user-data map
  class Edge_data_map {
  public:
    using key_type = Edge_descriptor;
    using key_const_type  = Edge_const_descriptor;
    using value_type = Edge_data;
    using reference = typename Map_reference_traits<Edge_data>::type;
    using const_reference = typename Map_const_reference_traits<Edge_data>::type;
    using category = boost::lvalue_property_map_tag;

    explicit Edge_data_map(Edge_container* edges) : m_edges(edges) {}

    template <typename Key>
    friend decltype(auto) get(const Edge_data_map& m, Key k) {
      if constexpr (std::is_void_v<Edge_data>) return;
      else if constexpr (UseVector) return (*m.m_edges)[k].data();
      else return k->data();
    }
    friend void put(const Edge_data_map& m, key_type k,
                    typename Map_param_traits<Edge_data>::type val) {
      if constexpr (!std::is_void_v<Edge_data>) {
        if constexpr (UseVector) (*m.m_edges)[k].set_data(val);
        else k->set_data(val);
      }
    }

  private:
    Edge_container* m_edges;
  };

  // --------------------------------------------------------------------------
  // Allocator access
  // --------------------------------------------------------------------------
  Allocator get_allocator() const noexcept { return Allocator(m_vertices.get_allocator()); }

protected:
  // --------------------------------------------------------------------------
  // Protected containers owned by derived classes through this base.
  // --------------------------------------------------------------------------
  Vertex_container m_vertices;
  Edge_container   m_edges;

  // The leftmost edge; in an identified arrangement, the edge whose left end
  // is at the identification point.
  Edge_descriptor m_left_boundary_edge {};

  // The rightmost edge; in an identified arrangement, the edge whose right end
  // is at the identification point.  Without the identification vertex, it is
  // the left boundary edge, which contains the identification point.
  Edge_descriptor m_right_boundary_edge {};

  // --------------------------------------------------------------------------
  // Protected helpers: descriptor -> node reference
  // --------------------------------------------------------------------------
  Vertex& v_ref(Vertex_descriptor v) {
    if constexpr (UseVector) return m_vertices[v];
    else return *v;
  }

  const Vertex& v_ref(Vertex_const_descriptor v) const {
    if constexpr (UseVector) return m_vertices[v];
    else return *v;
  }

  Edge& e_ref(Edge_descriptor e) {
    if constexpr (UseVector) return m_edges[e];
    else return *e;
  }

  const Edge& e_ref(Edge_const_descriptor e) const {
    if constexpr (UseVector) return m_edges[e];
    else return *e;
  }

public:
  // --------------------------------------------------------------------------
  // Constructors
  //
  // An empty arrangement has one edge that spans the entire parameter space
  // and no vertices.  Boundary vertices, if any, are created later by the
  // arrangement.
  // --------------------------------------------------------------------------
  Topology_traits_base() { init_spanning_edge(); }

  explicit Topology_traits_base(const Allocator& alloc) :
    m_vertices(Vertex_allocator(alloc)),
    m_edges(Edge_allocator(alloc))
  { init_spanning_edge(); }

  Topology_traits_base(Topology_traits_base&&) noexcept = default;
  Topology_traits_base& operator=(Topology_traits_base&&) noexcept = default;

  // Copy constructor: in list mode, the iterators stored inside the copied
  // nodes, and the boundary edges, are redirected into the copy.  Derived
  // classes that cache vertex descriptors redirect them with copied_vertex().
  Topology_traits_base(const Topology_traits_base& other) :
    m_vertices(other.m_vertices),
    m_edges(other.m_edges),
    m_left_boundary_edge(other.m_left_boundary_edge),
    m_right_boundary_edge(other.m_right_boundary_edge) {
    if constexpr (! UseVector) patch_cross_references(other);
  }

  Topology_traits_base& operator=(const Topology_traits_base& other) {
    if (this != &other) {
      Topology_traits_base tmp(other);
      swap(tmp);
    }
    return *this;
  }

  ~Topology_traits_base() = default;

  // --------------------------------------------------------------------------
  // removes all vertices and edges, and recreates the single spanning edge as
  // the left and right boundary edge.  Derived classes that cache descriptors hide this
  // function with one that calls it and resets their cache.
  // --------------------------------------------------------------------------
  void clear() {
    m_vertices.clear();
    m_edges.clear();
    init_spanning_edge();
  }

  // --------------------------------------------------------------------------
  // swaps the containers and the boundary edges.  Derived classes that
  // cache descriptors hide this function with one that calls it and swaps
  // their cache.
  // --------------------------------------------------------------------------
  void swap(Topology_traits_base& other) noexcept {
    std::swap(m_vertices, other.m_vertices);
    std::swap(m_edges, other.m_edges);
    std::swap(m_left_boundary_edge, other.m_left_boundary_edge);
    std::swap(m_right_boundary_edge, other.m_right_boundary_edge);
  }

  // --------------------------------------------------------------------------
  // Queries
  // --------------------------------------------------------------------------
  bool empty() const { return m_vertices.empty(); }
  Size number_of_vertices() const { return m_vertices.size();  }
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
      return boost::make_iterator_range(List_descriptor_iterator<Vertex_const_descriptor>(m_vertices.cbegin()),
                                        List_descriptor_iterator<Vertex_const_descriptor>(m_vertices.cend()));
  }

  Edge_descriptor_range edges() const {
    if constexpr (UseVector)
      return boost::make_iterator_range(Index_descriptor_iterator(0), Index_descriptor_iterator(m_edges.size()));
    else
      return boost::make_iterator_range(List_descriptor_iterator<Edge_const_descriptor>(m_edges.cbegin()),
                                        List_descriptor_iterator<Edge_const_descriptor>(m_edges.cend()));
  }

  // Mutable begin/end used by locate_impl linear scan.
  Vertex_descriptor vertices_begin() {
    if constexpr (UseVector) return 0;
    else return m_vertices.begin();
  }

  Vertex_descriptor vertices_end() {
    if constexpr (UseVector) return m_vertices.size();
    else return m_vertices.end();
  }

  // --------------------------------------------------------------------------
  // Navigation
  // --------------------------------------------------------------------------
  Edge_descriptor left_edge(Vertex_descriptor v) { return v_ref(v).m_left;  }
  Edge_descriptor right_edge(Vertex_descriptor v) { return v_ref(v).m_right; }
  Vertex_descriptor left_vertex(Edge_descriptor e) { return e_ref(e).m_left;  }
  Vertex_descriptor right_vertex(Edge_descriptor e) { return e_ref(e).m_right; }

  Edge_const_descriptor left_edge(Vertex_const_descriptor v) const { return v_ref(v).m_left;  }
  Edge_const_descriptor right_edge(Vertex_const_descriptor v) const { return v_ref(v).m_right; }
  Vertex_const_descriptor left_vertex(Edge_const_descriptor e) const { return e_ref(e).m_left;  }
  Vertex_const_descriptor right_vertex(Edge_const_descriptor e) const { return e_ref(e).m_right; }

  bool has_left_vertex (Edge_const_descriptor e) const { return e_ref(e).m_has_left;  }
  bool has_right_vertex(Edge_const_descriptor e) const { return e_ref(e).m_has_right; }

  const Point_1& vertex_point(Vertex_const_descriptor v) const { return v_ref(v).m_point; }

  // --------------------------------------------------------------------------
  // Boundary edge accessors
  //
  // left_boundary_edge() and right_boundary_edge() are always valid; both are
  // set to the spanning edge at construction.
  // --------------------------------------------------------------------------
  Edge_descriptor left_boundary_edge() { return m_left_boundary_edge; }
  Edge_const_descriptor left_boundary_edge() const { return m_left_boundary_edge; }

  void set_left_boundary_edge(Edge_descriptor e) { m_left_boundary_edge = e; }

  Edge_descriptor right_boundary_edge() { return m_right_boundary_edge; }
  Edge_const_descriptor right_boundary_edge() const { return m_right_boundary_edge; }

  void set_right_boundary_edge(Edge_descriptor e) { m_right_boundary_edge = e; }

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

  void set_left_edge(Vertex_descriptor v, Edge_descriptor e) { v_ref(v).m_left = e; }

  void set_right_edge(Vertex_descriptor v, Edge_descriptor e) { v_ref(v).m_right = e; }

  void set_left_vertex(Edge_descriptor e, Vertex_descriptor v) {
    e_ref(e).m_left = v;
    e_ref(e).m_has_left  = true;
  }

  void set_right_vertex(Edge_descriptor e, Vertex_descriptor v) {
    e_ref(e).m_right = v;
    e_ref(e).m_has_right = true;
  }

  void clear_left_vertex (Edge_descriptor e) { e_ref(e).m_has_left  = false; }
  void clear_right_vertex(Edge_descriptor e) { e_ref(e).m_has_right = false; }

  // erases the vertex `v`.  O(1) in list mode; in vector mode, O(1) for the
  // last vertex and O(n) otherwise, as every stored vertex index greater than
  // `v` is decremented.  Derived classes that cache vertex descriptors (e.g.
  // Closed_topology_traits) hide this function with one that patches their
  // cache.
  void erase_vertex(Vertex_descriptor v) {
    if constexpr (! UseVector) m_vertices.erase(v);
    else {
      if (v == m_vertices.size() - 1) { m_vertices.pop_back(); return; }
      m_vertices.erase(m_vertices.begin() + static_cast<std::ptrdiff_t>(v));
      for (auto& edge : m_edges) {
        if (edge.m_has_left  && edge.m_left  > v) --edge.m_left;
        if (edge.m_has_right && edge.m_right > v) --edge.m_right;
      }
    }
  }

  // erases the edge `e`.  In vector mode, decrements every stored edge index
  // greater than `e`, including the boundary edges.
  // Pre: `e` is neither the left nor the right boundary edge.
  void erase_edge(Edge_descriptor e) {
    CGAL_precondition(e != m_left_boundary_edge);
    CGAL_precondition(e != m_right_boundary_edge);
    if constexpr (! UseVector) m_edges.erase(e);
    else {
      if (e == m_edges.size() - 1) { m_edges.pop_back(); return; }
      m_edges.erase(m_edges.begin() + static_cast<std::ptrdiff_t>(e));
      for (auto& vert : m_vertices) {
        if (vert.m_left  > e) --vert.m_left;
        if (vert.m_right > e) --vert.m_right;
      }
      if (m_left_boundary_edge > e) --m_left_boundary_edge;
      if (m_right_boundary_edge > e) --m_right_boundary_edge;
    }
  }

private:
  // Minimal hash for list iterators: hash the address of the pointed-to node.
  struct Iterator_hash {
    template <typename It>
    std::size_t operator()(It it) const noexcept
    { return std::hash<const void*>{}(static_cast<const void*>(&*it)); }
  };

  // --------------------------------------------------------------------------
  // patch_cross_references() list mode only
  //
  // Called by the copy constructor after the member containers have been
  // value-copied (so all stored iterators still point into the source).
  // Rebuilds every stored iterator, including the boundary edges, to point
  // into the newly copied containers.
  // --------------------------------------------------------------------------
  struct Patch_maps {
    std::unordered_map<Edge_const_descriptor, Edge_descriptor, Iterator_hash> edge_map;
    std::unordered_map<Vertex_const_descriptor, Vertex_descriptor, Iterator_hash> vertex_map;
  };

  void patch_cross_references(const Topology_traits_base& other) {
    Patch_maps maps;
    maps.edge_map.reserve(m_edges.size());
    {
      auto src = other.m_edges.cbegin();
      auto dst = m_edges.begin();
      for (; src != other.m_edges.cend(); ++src, ++dst)
        maps.edge_map.emplace(src, dst);
    }
    maps.vertex_map.reserve(m_vertices.size());
    {
      auto src = other.m_vertices.cbegin();
      auto dst = m_vertices.begin();
      for (; src != other.m_vertices.cend(); ++src, ++dst)
        maps.vertex_map.emplace(src, dst);
    }
    // Patch vertex -> edge and edge -> vertex cross-references.
    {
      auto sv = other.m_vertices.cbegin();
      auto dv = m_vertices.begin();
      for (; sv != other.m_vertices.cend(); ++sv, ++dv) {
        dv->m_left  = maps.edge_map.at(sv->m_left);
        dv->m_right = maps.edge_map.at(sv->m_right);
      }
    }
    {
      auto se = other.m_edges.cbegin();
      auto de = m_edges.begin();
      for (; se != other.m_edges.cend(); ++se, ++de) {
        if (se->m_has_left)  de->m_left  = maps.vertex_map.at(se->m_left);
        if (se->m_has_right) de->m_right = maps.vertex_map.at(se->m_right);
      }
    }
    m_left_boundary_edge = maps.edge_map.at(other.m_left_boundary_edge);
    m_right_boundary_edge = maps.edge_map.at(other.m_right_boundary_edge);
  }

protected:

  // --------------------------------------------------------------------------
  // obtains the vertex of this copy that corresponds to the vertex `v` of
  // `other`, the source of the copy, by its position; list mode only.  Used by
  // the copy constructors of derived classes that cache vertex descriptors.
  // --------------------------------------------------------------------------
  Vertex_descriptor copied_vertex(const Topology_traits_base& other, Vertex_const_descriptor v) {
    static_assert(! UseVector, "copied_vertex() is needed only in list mode.");
    return std::next(m_vertices.begin(), std::distance(other.m_vertices.cbegin(), v));
  }

  // --------------------------------------------------------------------------
  // creates the single edge that spans the entire parameter space and sets it
  // as the left and right boundary edge.  No vertices are created here regardless of the
  // boundary tags; boundary vertices are created later by the arrangement when
  // the user inserts a point at a closed boundary.  Called from the
  // constructors and clear().
  // --------------------------------------------------------------------------
  void init_spanning_edge() {
    m_edges.emplace_back();
    if constexpr (UseVector) m_left_boundary_edge = 0;
    else m_left_boundary_edge = m_edges.begin();
    m_right_boundary_edge = m_left_boundary_edge;
  }

  void reserve_vertices(std::size_t n) { if constexpr (UseVector) m_vertices.reserve(n); }
  void reserve_edges(std::size_t n) { if constexpr (UseVector) m_edges.reserve(n); }
};

} // namespace Arrangement_on_curve_1
} // namespace CGAL

#endif
