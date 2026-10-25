// Copyright (c) 2026 Tel-Aviv University (Israel).
// All rights reserved.
// This file is part of CGAL (www.cgal.org).
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
// Author(s): Efi Fogel         <efifogel@gmail.com>

#ifndef CGAL_ARRANGEMENT_ON_CURVE_1_CLOSED_TOPOLOGY_TRAITS_H
#define CGAL_ARRANGEMENT_ON_CURVE_1_CLOSED_TOPOLOGY_TRAITS_H

#include <CGAL/assertions.h>
#include <CGAL/Arrangement_on_curve_1/Topology_traits_base.h>
#include <CGAL/Arrangement_on_curve_1/Parameter_space.h>

namespace CGAL {
namespace Arrangement_on_curve_1 {

/*! \ingroup PkgArrangementOnCurve1TopologyTraitsClasses
 *
 * \brief Topology traits for an arrangement on a closed master curve.
 *
 * Derives from `Topology_traits_base`, maintains the left
 * and right boundary vertices, and overrides `create_vertex` to register the
 * vertex as the left or right boundary vertex when the point lies on the
 * closed left or right boundary, using `GeometryTraits_1::Parameter_space_in_x_1`.
 */
template <typename GeometryTraits_1,
          typename VertexData = void,
          typename EdgeData = void,
          bool UseVector = false,
          typename Allocator_ = std::allocator<char>>
class Closed_topology_traits :
  public Topology_traits_base<typename GeometryTraits_1::Point_1, VertexData, EdgeData, UseVector, Allocator_> {
public:
  using Geometry_traits_1 = GeometryTraits_1;
  using Vertex_data = VertexData;
  using Edge_data = EdgeData;
  using Allocator = Allocator_;

  // The standard stl type
  using allocator_type = Allocator_;

private:
  using Base = Topology_traits_base<typename GeometryTraits_1::Point_1, VertexData, EdgeData, UseVector, Allocator>;

public:
  using Base::Base;
  using typename Base::Point_1;
  using typename Base::Vertex_descriptor;
  using typename Base::Vertex_const_descriptor;

  // --------------------------------------------------------------------------
  // Constructors
  // --------------------------------------------------------------------------
  Closed_topology_traits() = default;

  Closed_topology_traits(Closed_topology_traits&&) noexcept = default;
  Closed_topology_traits& operator=(Closed_topology_traits&&) noexcept = default;

  // In list mode, the cached boundary vertex descriptors still point into the
  // source after the base copy; they are redirected to the copied vertices.
  Closed_topology_traits(const Closed_topology_traits& other) :
    Base(other),
    m_left_boundary_vertex(other.m_left_boundary_vertex),
    m_has_left_boundary_vertex(other.m_has_left_boundary_vertex),
    m_right_boundary_vertex(other.m_right_boundary_vertex),
    m_has_right_boundary_vertex(other.m_has_right_boundary_vertex) {
    if constexpr (! UseVector) {
      if (m_has_left_boundary_vertex) m_left_boundary_vertex = this->copied_vertex(other, other.m_left_boundary_vertex);
      if (m_has_right_boundary_vertex)
        m_right_boundary_vertex = this->copied_vertex(other, other.m_right_boundary_vertex);
    }
  }

  Closed_topology_traits& operator=(const Closed_topology_traits& other) {
    if (this != &other) {
      Closed_topology_traits tmp(other);
      swap(tmp);
    }
    return *this;
  }

  ~Closed_topology_traits() = default;

  // creates a new vertex associated with `p` without registering it as a
  // boundary vertex; the caller registers it if needed.
  using Base::create_vertex;

  //! creates a new vertex associated with point `p`.
  Vertex_descriptor create_vertex(const Point_1& p, const GeometryTraits_1& gt) {
    Vertex_descriptor v = Base::create_vertex(p);
    auto ps = gt.parameter_space_in_x_1_object()(p);
    if (ps == Parameter_space::LEFT_BOUNDARY) set_left_boundary_vertex(v);
    else if (ps == Parameter_space::RIGHT_BOUNDARY) set_right_boundary_vertex(v);
    return v;
  }

  // --------------------------------------------------------------------------
  // Boundary vertex accessors.
  // has_*_boundary_vertex() returns false until the arrangement creates the
  // vertex; it may return false again after clear() or after the vertex is
  // erased.
  // --------------------------------------------------------------------------
  bool has_left_boundary_vertex() const { return m_has_left_boundary_vertex; }

  Vertex_descriptor left_boundary_vertex() {
    CGAL_assertion(m_has_left_boundary_vertex);
    return m_left_boundary_vertex;
  }

  Vertex_const_descriptor left_boundary_vertex() const {
    CGAL_assertion(m_has_left_boundary_vertex);
    return m_left_boundary_vertex;
  }

  void set_left_boundary_vertex(Vertex_descriptor v) {
    m_left_boundary_vertex = v;
    m_has_left_boundary_vertex = true;
  }

  bool has_right_boundary_vertex() const { return m_has_right_boundary_vertex; }

  Vertex_descriptor right_boundary_vertex() {
    CGAL_assertion(m_has_right_boundary_vertex);
    return m_right_boundary_vertex;
  }

  Vertex_const_descriptor right_boundary_vertex() const {
    CGAL_assertion(m_has_right_boundary_vertex);
    return m_right_boundary_vertex;
  }

  void set_right_boundary_vertex(Vertex_descriptor v) {
    m_right_boundary_vertex = v;
    m_has_right_boundary_vertex = true;
  }

  // patches the boundary vertex descriptors in vector mode.
  void erase_vertex(Vertex_descriptor v) {
    if (m_has_left_boundary_vertex && (m_left_boundary_vertex == v)) m_has_left_boundary_vertex = false;
    if (m_has_right_boundary_vertex && (m_right_boundary_vertex == v)) m_has_right_boundary_vertex = false;
    Base::erase_vertex(v);
    if constexpr (UseVector) {
      if (m_has_left_boundary_vertex && (m_left_boundary_vertex > v)) --m_left_boundary_vertex;
      if (m_has_right_boundary_vertex && (m_right_boundary_vertex > v)) --m_right_boundary_vertex;
    }
  }

  // --------------------------------------------------------------------------
  // clear / swap
  // --------------------------------------------------------------------------
  void clear() {
    Base::clear();
    m_has_left_boundary_vertex = false;
    m_has_right_boundary_vertex = false;
  }

  void swap(Closed_topology_traits& other) noexcept {
    Base::swap(other);
    std::swap(m_left_boundary_vertex, other.m_left_boundary_vertex);
    std::swap(m_has_left_boundary_vertex, other.m_has_left_boundary_vertex);
    std::swap(m_right_boundary_vertex, other.m_right_boundary_vertex);
    std::swap(m_has_right_boundary_vertex, other.m_has_right_boundary_vertex);
  }

private:
  Vertex_descriptor m_left_boundary_vertex {};
  bool m_has_left_boundary_vertex = false;
  Vertex_descriptor m_right_boundary_vertex {};
  bool m_has_right_boundary_vertex = false;
};

} // namespace Arrangement_on_curve_1
} // namespace CGAL

#endif
