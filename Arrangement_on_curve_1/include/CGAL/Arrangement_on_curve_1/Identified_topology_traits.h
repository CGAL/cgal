// Copyright (c) 2026 Tel-Aviv University (Israel).
// All rights reserved.
// This file is part of CGAL (www.cgal.org).
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
// Author(s): Efi Fogel         <efifogel@gmail.com>

#ifndef CGAL_ARRANGEMENT_ON_CURVE_1_IDENTIFIED_TOPOLOGY_TRAITS_H
#define CGAL_ARRANGEMENT_ON_CURVE_1_IDENTIFIED_TOPOLOGY_TRAITS_H

/*! Topology traits for a 1D arrangement on a cyclic curve (identified boundary).
 *
 * When both sides of the parameter space are identified (tag:
 * `Identified_side_tag`), the curve is cyclic; for example a circle,
 * a direction-space \f$\mathbb{S}^1\f$, or any periodic 1D domain.
 *
 * Structure:
 *
 * - An empty arrangement contains one edge that spans the entire parameter space and no vertices.
 *
 * - Vertices and edges together form a single doubly-linked cyclic list.
 *
 * - `left_boundary_edge()`: the edge whose left end is at the identification point
 *                           (the canonical entry point for cyclic traversal). Set to the
 *                           spanning edge at construction and updated after insertions and deletions.
 *
 * - `right_boundary_edge()`: the edge whose right end is at the identification point, or, if
 *                            the identification vertex does not exist, `left_boundary_edge()`.
 *
 * - `split_edge(e, p)`: splits edge `e` into (e, e_new), inserts vertex `v` between them,
 *                       and rewires both endpoints cyclically.
 *
 * - `remove(v)`: merges the left and right edges of `v` back into a single edge.
 *                Removing the last vertex restores the single spanning edge.
 *
 * All shared storage infrastructure, including the boundary edges, is
 * inherited from `Topology_traits_base`.  `Identified_topology_traits` adds the
 * identification-vertex cache.
 *
 * Template parameters:
 *   GeometryTraits_1 - a model of `AocIdentifiedBoundaryTraits_1`.
 *   VertexData       - optional user data per vertex (void = none).
 *   EdgeData         - optional user data per edge   (void = none).
 *   UseVector        - true -> std::vector + index descriptors;
 *                      false -> std::list + iterator descriptors (default).
 *   Allocator        - allocator for internal containers.
 */

#include <CGAL/Arrangement_on_curve_1/Topology_traits_base.h>

namespace CGAL {
namespace Arrangement_on_curve_1 {

template <typename GeometryTraits_1,
          typename VertexData = void,
          typename EdgeData = void,
          bool UseVector = false,
          typename Allocator_ = std::allocator<char>>
class Identified_topology_traits :
  public Topology_traits_base<typename GeometryTraits_1::Point_1, VertexData, EdgeData, UseVector, Allocator_> {
public:
  using Geometry_traits_1 = GeometryTraits_1;
  using Vertex_data = VertexData;
  using Edge_data = EdgeData;
  using Allocator = Allocator_;

  //! The standard stl type.
  using allocator_type = Allocator;

private:
  using Base = Topology_traits_base<typename GeometryTraits_1::Point_1, VertexData, EdgeData, UseVector, Allocator>;

public:
  using Point_1 = typename Geometry_traits_1::Point_1;

  using Base::Base;

  using typename Base::Vertex_descriptor;
  using typename Base::Vertex_const_descriptor;

  // --------------------------------------------------------------------------
  // Constructors
  // --------------------------------------------------------------------------
  Identified_topology_traits() = default;

  Identified_topology_traits(Identified_topology_traits&&) noexcept = default;
  Identified_topology_traits& operator=(Identified_topology_traits&&) noexcept = default;

  // In list mode, the cached identification vertex descriptor still points
  // into the source after the base copy; it is redirected to the copied vertex.
  Identified_topology_traits(const Identified_topology_traits& other) :
    Base(other),
    m_identification_vertex(other.m_identification_vertex),
    m_has_identification_vertex(other.m_has_identification_vertex) {
    if constexpr (! UseVector) {
      if (m_has_identification_vertex)
        m_identification_vertex = this->copied_vertex(other, other.m_identification_vertex);
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
  // creates a vertex and, when the point lies on the identification, registers
  // it as the identification vertex.  GeometryTraits_1 must provide
  // is_on_x_identification_1_object().
  // --------------------------------------------------------------------------
  Vertex_descriptor create_vertex(const Point_1& p, const GeometryTraits_1& gt) {
    Vertex_descriptor v = Base::create_vertex(p);
    if (gt.is_on_x_identification_1_object()(p)) set_identification_vertex(v);
    return v;
  }

  // --------------------------------------------------------------------------
  // clears the identification vertex cache if it is removed, then patches the
  // stored index in vector mode.
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
    Base::clear();
    m_has_identification_vertex = false;
  }

  void swap(Identified_topology_traits& other) noexcept {
    Base::swap(other);
    std::swap(m_identification_vertex, other.m_identification_vertex);
    std::swap(m_has_identification_vertex, other.m_has_identification_vertex);
  }

private:
  Vertex_descriptor m_identification_vertex {};
  bool m_has_identification_vertex = false;
};

} // namespace Arrangement_on_curve_1
} // namespace CGAL

#endif
