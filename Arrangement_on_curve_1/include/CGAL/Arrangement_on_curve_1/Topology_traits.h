// Copyright (c) 2026 Tel-Aviv University (Israel).
// All rights reserved.
// This file is part of CGAL (www.cgal.org).
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
// Author(s): Efi Fogel         <efifogel@gmail.com>

#ifndef CGAL_ARRANGEMENT_ON_CURVE_1_TOPOLOGY_TRAITS_H
#define CGAL_ARRANGEMENT_ON_CURVE_1_TOPOLOGY_TRAITS_H

/*! \file
 * Generic topology traits for a 1D arrangement on a master curve with open or
 * closed (but not identified) boundary conditions.
 *
 * `Topology_traits` derives from `Topology_traits_base` and handles all
 * combinations of `Open_side_tag` and `Closed_side_tag` on the two sides of
 * the parameter space.  `Identified_side_tag` is not supported; use
 * `Identified_topology_traits` for closed-loop curves.
 */

#include <CGAL/assertions.h>
#include <CGAL/Arrangement_on_curve_1/Topology_traits_base.h>
#include <CGAL/Arrangement_on_curve_1/Side_tags.h>

namespace CGAL {
namespace Arrangement_on_curve_1 {

template <typename GeometryTraits_1,
          typename LeftTag = typename GeometryTraits_1::Left_side_category,
          typename RightTag = typename GeometryTraits_1::Right_side_category,
          typename VertexData = void,
          typename EdgeData = void,
          bool UseVector = false,
          typename Allocator_ = std::allocator<char>>
class Topology_traits :
  public Topology_traits_base<typename GeometryTraits_1::Point_1, VertexData, EdgeData, UseVector, Allocator_> {
  static_assert(! Side_is_identified<LeftTag>::value && ! Side_is_identified<RightTag>::value,
                "Topology_traits does not support Identified_side_tag. "
                "Use Identified_topology_traits for closed-loop (identified) curves.");
public:
  using Geometry_traits_1 = GeometryTraits_1;
  using Vertex_data = VertexData;
  using Edge_data = EdgeData;
  using Allocator = Allocator_;

  //! The standard stl type.
  using allocator_type = Allocator;

private:
  using Base = Topology_traits_base<typename Geometry_traits_1::Point_1, Vertex_data, Edge_data, UseVector, Allocator>;

public:
  using Base::Base;

  using typename Base::Point_1;
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
  static constexpr bool left_is_closed = Side_is_closed<LeftTag>::value;
  static constexpr bool right_is_closed = Side_is_closed<RightTag>::value;

  // --------------------------------------------------------------------------
  // Constructors
  // --------------------------------------------------------------------------

  // Default constructor — valid for all tag combinations.
  // Creates a single edge spanning the parameter space and no vertices.
  // Boundary vertices, if any, are created later by the arrangement.
  Topology_traits() { init_spanning_edge(); }

  // Allocator-aware default constructor.
  explicit Topology_traits(const Allocator& alloc) : Base(alloc) { init_spanning_edge(); }

  // --------------------------------------------------------------------------
  // Move constructor / assignment
  // --------------------------------------------------------------------------
  Topology_traits(Topology_traits&&) noexcept = default;
  Topology_traits& operator=(Topology_traits&&) noexcept = default;

  // --------------------------------------------------------------------------
  // Copy constructor
  // --------------------------------------------------------------------------
  Topology_traits(const Topology_traits& other) :
    Base(other),
    m_left_boundary_edge(other.m_left_boundary_edge),
    m_left_boundary_vertex(other.m_left_boundary_vertex),
    m_right_boundary_vertex(other.m_right_boundary_vertex),
    m_has_left_boundary_vertex(other.m_has_left_boundary_vertex),
    m_has_right_boundary_vertex(other.m_has_right_boundary_vertex) {
    if constexpr (! UseVector) {
      auto maps = this->patch_cross_references(other);
      m_left_boundary_edge = maps.edge_map.at(other.m_left_boundary_edge);
      if (m_has_left_boundary_vertex) m_left_boundary_vertex  = maps.vertex_map.at(other.m_left_boundary_vertex);
      if (m_has_right_boundary_vertex) m_right_boundary_vertex = maps.vertex_map.at(other.m_right_boundary_vertex);
    }
  }

  Topology_traits& operator=(const Topology_traits& other) {
    if (this != &other) { Topology_traits tmp(other); swap(tmp); }
    return *this;
  }

  ~Topology_traits() = default;

  // --------------------------------------------------------------------------
  // create_vertex — creates a vertex and stores its point.
  //
  // The base implementation (open/open topology) ignores the geometry traits
  // entirely: no boundary registration is possible when both sides are open.
  // Derived classes (Closed_topology_traits, Left_closed_topology_traits,
  // Right_closed_topology_traits) override this to register the vertex as the
  // left or right boundary vertex when the point lies on a closed boundary.
  // --------------------------------------------------------------------------
  Vertex_descriptor create_vertex(const Point_1& p, const GeometryTraits_1& /*gt*/) { return Base::create_vertex(p); }

  // --------------------------------------------------------------------------
  // Left edge accessor
  // --------------------------------------------------------------------------
  Edge_descriptor left_boundary_edge() { return m_left_boundary_edge; }
  Edge_const_descriptor left_boundary_edge() const { return m_left_boundary_edge; }

  void set_left_boundary_edge(Edge_descriptor e) { m_left_boundary_edge = e; }

  // --------------------------------------------------------------------------
  // Boundary vertex accessors — meaningful only when the corresponding side is
  // closed.  has_left/right_boundary_vertex() returns false until the
  // arrangement creates the vertex; it may return false again after clear() or
  // after the vertex is erased.
  // --------------------------------------------------------------------------
  bool has_left_boundary_vertex() const { return m_has_left_boundary_vertex; }
  bool has_right_boundary_vertex() const { return m_has_right_boundary_vertex; }

  Vertex_descriptor left_boundary_vertex() {
    CGAL_assertion(m_has_left_boundary_vertex);
    return m_left_boundary_vertex;
  }

  Vertex_const_descriptor left_boundary_vertex() const {
    CGAL_assertion(m_has_left_boundary_vertex);
    return m_left_boundary_vertex;
  }

  Vertex_descriptor right_boundary_vertex() {
    CGAL_assertion(m_has_right_boundary_vertex);
    return m_right_boundary_vertex;
  }

  Vertex_const_descriptor right_boundary_vertex() const {
    CGAL_assertion(m_has_right_boundary_vertex);
    return m_right_boundary_vertex;
  }

  void set_left_boundary_vertex(Vertex_descriptor v) {
    m_left_boundary_vertex  = v;
    m_has_left_boundary_vertex  = true;
  }

  void set_right_boundary_vertex(Vertex_descriptor v) {
    m_right_boundary_vertex = v;
    m_has_right_boundary_vertex = true;
  }

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

  // erase_vertex — patches boundary vertex descriptors in vector mode.
  void erase_vertex(Vertex_descriptor v) {
    if (m_has_left_boundary_vertex && m_left_boundary_vertex == v) m_has_left_boundary_vertex = false;
    if (m_has_right_boundary_vertex && m_right_boundary_vertex == v) m_has_right_boundary_vertex = false;
    this->Base::erase_vertex(v);
    if constexpr (UseVector) {
      if (m_has_left_boundary_vertex && m_left_boundary_vertex > v) --m_left_boundary_vertex;
      if (m_has_right_boundary_vertex && m_right_boundary_vertex > v) --m_right_boundary_vertex;
    }
  }

  // --------------------------------------------------------------------------
  // clear / swap
  // --------------------------------------------------------------------------
  void clear() {
    this->m_vertices.clear();
    this->m_edges.clear();
    m_has_left_boundary_vertex  = false;
    m_has_right_boundary_vertex = false;
    init_spanning_edge();
  }

  void swap(Topology_traits& other) noexcept {
    std::swap(this->m_vertices, other.m_vertices);
    std::swap(this->m_edges, other.m_edges);
    std::swap(m_left_boundary_edge, other.m_left_boundary_edge);
    std::swap(m_left_boundary_vertex,  other.m_left_boundary_vertex);
    std::swap(m_right_boundary_vertex, other.m_right_boundary_vertex);
    std::swap(m_has_left_boundary_vertex,  other.m_has_left_boundary_vertex);
    std::swap(m_has_right_boundary_vertex, other.m_has_right_boundary_vertex);
  }

private:
  Edge_descriptor m_left_boundary_edge {};
  Vertex_descriptor m_left_boundary_vertex {};
  Vertex_descriptor m_right_boundary_vertex {};
  bool m_has_left_boundary_vertex = false;
  bool m_has_right_boundary_vertex = false;

  // --------------------------------------------------------------------------
  // init_spanning_edge
  //
  // Creates the single edge that spans the entire parameter space.
  // No vertices are created here regardless of the boundary tags; boundary
  // vertices are created later by the arrangement when the user inserts a
  // point at a closed boundary.
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
