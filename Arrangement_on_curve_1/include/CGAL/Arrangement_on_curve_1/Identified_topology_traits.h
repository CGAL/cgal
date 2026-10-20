// Copyright (c) 2026 Tel-Aviv University (Israel).
// All rights reserved.
// This file is part of CGAL (www.cgal.org).
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
// Author(s): Efi Fogel         <efifogel@gmail.com>

#ifndef CGAL_ARRANGEMENT_ON_CURVE_1_IDENTIFIED_TOPOLOGY_TRAITS_H
#define CGAL_ARRANGEMENT_ON_CURVE_1_IDENTIFIED_TOPOLOGY_TRAITS_H

/*! \file
 * Topology traits for a 1D arrangement on a closed loop (identified boundary).
 *
 * When both sides of the parameter space are identified (tag:
 * `Identified_side_tag`), the curve is a closed loop — for example a circle,
 * a direction-space S¹, or any periodic 1D domain.
 *
 * Structure:
 *
 *   • An empty arrangement contains one edge that spans the entire parameter
 *     space and no vertices.  This is the same invariant as the unbounded and
 *     bounded topology types; the difference is that on a closed loop the two
 *     ends of that edge are identified (they represent the same point).
 *
 *   • Vertices and edges together form a single doubly-linked cyclic list:
 *       …  ←e₀→ v₀ ←e₁→ v₁ ←e₂→ v₂ ←e₀→ …
 *     where e₀.right = v₀, v₀.left = e₀, v₀.right = e₁, e₁.left = v₀, …
 *     and the last edge's right wraps back to v₀.
 *
 *   • `left_boundary_edge()`: the edge whose left end is at the identification
 *     point (the canonical entry point for cyclic traversal).  Set to the
 *     spanning edge at construction and updated after insertions and deletions.
 *
 *   • `split_edge(e, p)`: splits edge e into (e, e_new), inserts vertex v
 *     between them, and rewires both endpoints cyclically.
 *
 *   • `remove(v)`: merges the left and right edges of v back into a single
 *     edge.  Removing the last vertex restores the single spanning edge.
 *
 * All shared storage infrastructure is inherited from `Topology_traits_base`.
 * `Identified_topology_traits` adds `m_left_boundary_edge`, the
 * identification-vertex cache, and `init_spanning_edge()`.
 *
 * Template parameters:
 *   GeometryTraits_1 - a model of `AocIdentifiedBoundaryTraits_1`.
 *   VertexData       - optional user data per vertex (void = none).
 *   EdgeData         - optional user data per edge   (void = none).
 *   UseVector        - true → std::vector + index descriptors;
 *                      false → std::list  + iterator descriptors (default).
 *   Allocator        - allocator for internal containers.
 */

#include <CGAL/Arrangement_on_curve_1/Topology_traits_base.h>

namespace CGAL {
namespace Arrangement_on_curve_1 {

template <typename GeometryTraits_1,
          typename VertexData = void,
          typename EdgeData = void,
          bool UseVector = false,
          typename Allocator = std::allocator<char>>
class Identified_topology_traits :
    public Topology_traits_base<typename GeometryTraits_1::Point_1, VertexData, EdgeData, UseVector, Allocator> {
  using Base = Topology_traits_base<typename GeometryTraits_1::Point_1, VertexData, EdgeData, UseVector, Allocator>;

public:
  using Geometry_traits_1 = GeometryTraits_1;
  using Point_1 = typename Geometry_traits_1::Point_1;

  using typename Base::Allocator_type;
  using typename Base::Vertex_container;
  using typename Base::Edge_container;
  using typename Base::Vertex_descriptor;
  using typename Base::Edge_descriptor;
  using typename Base::Vertex_const_descriptor;
  using typename Base::Edge_const_descriptor;
  using typename Base::Vertex_descriptor_iterator;
  using typename Base::Edge_descriptor_iterator;
  using typename Base::Vertex_descriptor_range;
  using typename Base::Edge_descriptor_range;
  using typename Base::Size;
  using typename Base::Vertex;
  using typename Base::Edge;
  using typename Base::Vertex_point_map;
  using typename Base::Vertex_data_map;
  using typename Base::Edge_data_map;

  static constexpr bool use_vector = Base::use_vector;

  // --------------------------------------------------------------------------
  // Constructors
  //
  // An empty identified arrangement has one spanning edge and no vertices.
  // --------------------------------------------------------------------------
  Identified_topology_traits() { init_spanning_edge(); }

  explicit Identified_topology_traits(const Allocator_type& alloc) :
    Base(alloc) { init_spanning_edge(); }

  Identified_topology_traits(Identified_topology_traits&&) noexcept = default;
  Identified_topology_traits& operator=(Identified_topology_traits&&) noexcept = default;

  Identified_topology_traits(const Identified_topology_traits& other) :
    Base(other),
    m_left_boundary_edge(other.m_left_boundary_edge),
    m_identification_vertex(other.m_identification_vertex),
    m_has_identification_vertex(other.m_has_identification_vertex) {
    if constexpr (! UseVector) {
      auto maps = this->patch_cross_references(other);
      m_left_boundary_edge = maps.edge_map.at(other.m_left_boundary_edge);
      if (m_has_identification_vertex) m_identification_vertex = maps.vertex_map.at(other.m_identification_vertex);
    }
  }

  Identified_topology_traits& operator=(const Identified_topology_traits& other) {
    if (this != &other) {
      Identified_topology_traits tmp(other);
      swap(tmp);
    }
    return *this;
  }

  ~Identified_topology_traits() = default;

  // --------------------------------------------------------------------------
  // Identification vertex accessors.
  //
  // The identification vertex is the unique vertex that lies on the identified
  // boundary (both sides of the parameter space map to the same point).
  // has_identification_vertex() returns false until the arrangement creates
  // the vertex by inserting the identification point, and may return false
  // again after clear() or after the vertex is erased.
  // --------------------------------------------------------------------------
  bool has_identification_vertex() const { return m_has_identification_vertex; }

  Vertex_descriptor identification_vertex() {
    CGAL_assertion(m_has_identification_vertex);
    return m_identification_vertex;
  }

  Vertex_const_descriptor identification_vertex() const {
    CGAL_assertion(m_has_identification_vertex);
    return m_identification_vertex;
  }

  void set_identification_vertex(Vertex_descriptor v) {
    m_identification_vertex = v;
    m_has_identification_vertex = true;
  }

  // --------------------------------------------------------------------------
  // create_vertex — creates a vertex and, when the point lies on the
  // identification, registers it as the identification vertex.
  //
  // GeometryTraits_1 must provide is_on_x_identification_1_object().
  // --------------------------------------------------------------------------
  Vertex_descriptor create_vertex(const Point_1& p, const GeometryTraits_1& gt) {
    Vertex_descriptor v = Base::create_vertex(p);
    if (gt.is_on_x_identification_1_object()(p)) set_identification_vertex(v);
    return v;
  }

  // --------------------------------------------------------------------------
  // left_boundary_edge — the edge whose left end is at the identification
  // point (canonical entry point for cyclic traversal).
  // Always valid (the spanning edge is created at construction).
  // --------------------------------------------------------------------------
  Edge_descriptor left_boundary_edge() { return m_left_boundary_edge; }
  Edge_const_descriptor left_boundary_edge() const { return m_left_boundary_edge; }
  void set_left_boundary_edge(Edge_descriptor e) { m_left_boundary_edge = e; }

  // --------------------------------------------------------------------------
  // erase_edge — patches m_left_boundary_edge in vector mode.
  // --------------------------------------------------------------------------
  void erase_edge(Edge_descriptor e) {
    if constexpr (! UseVector) {
      this->m_edges.erase(e);
      return;
    }
    if (e == this->m_edges.size() - 1) {
      this->m_edges.pop_back();
      if (m_left_boundary_edge == e) --m_left_boundary_edge;
      return;
    }
    this->erase_edge_base(e);
    if (m_left_boundary_edge > e) --m_left_boundary_edge;
  }

  // --------------------------------------------------------------------------
  // erase_vertex — clears the identification vertex cache if it is removed,
  // then patches the stored index in vector mode.
  // --------------------------------------------------------------------------
  void erase_vertex(Vertex_descriptor v) {
    if (m_has_identification_vertex && m_identification_vertex == v) m_has_identification_vertex = false;
    this->Base::erase_vertex(v);
    if constexpr (UseVector) {
      if (m_has_identification_vertex && m_identification_vertex > v) --m_identification_vertex;
    }
  }

  // --------------------------------------------------------------------------
  // clear / swap
  // --------------------------------------------------------------------------
  void clear() {
    this->m_vertices.clear();
    this->m_edges.clear();
    m_has_identification_vertex = false;
    init_spanning_edge();
  }

  void swap(Identified_topology_traits& other) noexcept {
    std::swap(this->m_vertices, other.m_vertices);
    std::swap(this->m_edges, other.m_edges);
    std::swap(m_left_boundary_edge, other.m_left_boundary_edge);
    std::swap(m_identification_vertex, other.m_identification_vertex);
    std::swap(m_has_identification_vertex, other.m_has_identification_vertex);
  }

private:
  Edge_descriptor m_left_boundary_edge {};
  Vertex_descriptor m_identification_vertex {};
  bool m_has_identification_vertex = false;

  // --------------------------------------------------------------------------
  // init_spanning_edge — creates the single edge that spans the entire
  // parameter space and sets it as the left boundary edge.
  // Called from constructors and clear().
  // --------------------------------------------------------------------------
  void init_spanning_edge() {
    this->m_edges.emplace_back();
    if constexpr (UseVector) m_left_boundary_edge = 0;
    else m_left_boundary_edge = this->m_edges.begin();
  }
};

} // namespace Arrangement_on_curve_1
} // namespace CGAL

#endif
