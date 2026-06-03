// Copyright (c) 2026 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Efi Fogel         <efifogel@gmail.com>

#ifndef CGAL_ARRANGEMENT_ON_CURVE_1_H
#define CGAL_ARRANGEMENT_ON_CURVE_1_H

#include <variant>
#include <CGAL/basic.h>

namespace CGAL {
namespace Arrangement_on_curve_1 {

template <typename GeometryTraits_1, typename TopologyTraits>
class Arrangement_on_curve_1 {
public:
  using Geometry_traits_1 = GeometryTraits_1;
  using Topology_traits = TopologyTraits;
  using Point_1 = typename Geometry_traits_1::Point_1;

  using Vertex_descriptor = typename Topology_traits::Vertex_descriptor;
  using Edge_descriptor = typename Topology_traits::Edge_descriptor;

  using Location_result = std::variant<Vertex_descriptor, Edge_descriptor, void*>;
  using Const_location_result = std::variant<Vertex_descriptor, Edge_descriptor, void*>;

private:
  Geometry_traits_1 m_geometry_traits;
  Topology_traits m_topology_traits;

public:
  Arrangement_on_curve_1(const Geometry_traits_1& geom_tr = Geometry_traits_1(),
                         const Topology_traits& topo_tr = Topology_traits()) :
    m_geometry_traits(geom_tr), m_topology_traits(topo_tr) {}

  const Geometry_traits_1& geometry_traits_1() const { return m_geometry_traits; }
  const Topology_traits& topology_traits() const { return m_topology_traits; }
  Topology_traits& topology_traits() { return m_topology_traits; }

  bool is_empty() const { return m_topology_traits.is_empty(); }
  size_t number_of_vertices() const { return m_topology_traits.number_of_vertices(); }
  size_t number_of_edges() const { return m_topology_traits.number_of_edges(); }

  auto vertices() const { return m_topology_traits.vertices(); }
  auto edges() const { return m_topology_traits.edges(); }

  auto vertex_point_map() const { return m_topology_traits.vertex_point_map(); }
  auto vertex_data_map() const { return m_topology_traits.vertex_data_map(); }
  auto edge_data_map() const { return m_topology_traits.edge_data_map(); }

  Vertex_descriptor left_vertex(Edge_descriptor e) const { return m_topology_traits.left_vertex(e); }
  Vertex_descriptor right_vertex(Edge_descriptor e) const { return m_topology_traits.right_vertex(e); }
  bool has_left_vertex(Edge_descriptor e) const { return m_topology_traits.has_left_vertex(e); }
  bool has_right_vertex(Edge_descriptor e) const { return m_topology_traits.has_right_vertex(e); }

  // ============================================================================
  // HIGH-LEVEL TOPOLOGICAL OPERATIONS
  // ============================================================================
  Vertex_descriptor insert_empty(const Point_1& p) {
    Vertex_descriptor v = m_topology_traits.create_vertex(p);
    Edge_descriptor e_left = m_topology_traits.create_edge();
    Edge_descriptor e_right = m_topology_traits.create_edge();

    m_topology_traits.set_right_vertex(e_left, v);
    m_topology_traits.set_left_vertex(e_right, v);

    m_topology_traits.set_left_edge(v, e_left);
    m_topology_traits.set_right_edge(v, e_right);
    return v;
  }

  Vertex_descriptor insert_before(Vertex_descriptor v_first, const Point_1& p) {
    Vertex_descriptor v_new = m_topology_traits.create_vertex_front(p);
    Edge_descriptor e_left_unbounded = m_topology_traits.left_edge(v_first);
    m_topology_traits.set_right_vertex(e_left_unbounded, v_new);

    Edge_descriptor e_new = m_topology_traits.create_edge();
    m_topology_traits.set_left_vertex(e_new, v_new);
    m_topology_traits.set_right_vertex(e_new, v_first);

    m_topology_traits.set_left_edge(v_new, e_left_unbounded);
    m_topology_traits.set_right_edge(v_new, e_new);
    m_topology_traits.set_left_edge(v_first, e_new);
    return v_new;
  }

  Vertex_descriptor insert_after(Vertex_descriptor v_last, const Point_1& p) {
    Vertex_descriptor v_new = m_topology_traits.create_vertex(p);
    Edge_descriptor e_right_unbounded = m_topology_traits.right_edge(v_last);
    m_topology_traits.set_left_vertex(e_right_unbounded, v_new);

    Edge_descriptor e_new = m_topology_traits.create_edge();
    m_topology_traits.set_left_vertex(e_new, v_last);
    m_topology_traits.set_right_vertex(e_new, v_new);

    m_topology_traits.set_right_edge(v_last, e_new);
    m_topology_traits.set_left_edge(v_new, e_new);
    m_topology_traits.set_right_edge(v_new, e_right_unbounded);
    return v_new;
  }

  Vertex_descriptor split_edge(Edge_descriptor e, const Point_1& p) {
    if (! m_topology_traits.has_left_vertex(e)) return insert_before(m_topology_traits.right_vertex(e), p);
    if (! m_topology_traits.has_right_vertex(e)) return insert_after(m_topology_traits.left_vertex(e), p);

    Vertex_descriptor v_tgt = m_topology_traits.right_vertex(e);
    Vertex_descriptor v_new = m_topology_traits.create_vertex(v_tgt, p);

    Edge_descriptor e_new = m_topology_traits.create_edge();
    m_topology_traits.set_left_vertex(e_new, v_new);
    m_topology_traits.set_right_vertex(e_new, v_tgt);
    m_topology_traits.set_right_vertex(e, v_new);

    m_topology_traits.set_left_edge(v_new, e);
    m_topology_traits.set_right_edge(v_new, e_new);
    m_topology_traits.set_left_edge(v_tgt, e_new);
    return v_new;
  }

  void remove(Vertex_descriptor v) {
    // Extract the immediate topological neighbors of the vertex
    Edge_descriptor e_left  = m_topology_traits.left_edge(v);
    Edge_descriptor e_right = m_topology_traits.right_edge(v);

    bool has_left_neighbor  = m_topology_traits.has_left_vertex(e_left);
    bool has_right_neighbor = m_topology_traits.has_right_vertex(e_right);

    // Case 1: The arrangement collapses to empty space if this was the final vertex
    if (! has_left_neighbor && ! has_right_neighbor) {
      m_topology_traits.erase_vertex(v);
      m_topology_traits.erase_edge(e_left);
      m_topology_traits.erase_edge(e_right);
      return;
    }

    // Case 2: Removing the leftmost vertex shifts the unbounded track boundary rightward
    if (! has_left_neighbor) {
      Vertex_descriptor v_right_neighbor = m_topology_traits.right_vertex(e_right);

      // e_right becomes the new left-unbounded tracking edge
      m_topology_traits.clear_left_vertex(e_right);
      m_topology_traits.set_left_edge(v_right_neighbor, e_right);

      m_topology_traits.erase_vertex(v);
      m_topology_traits.erase_edge(e_left);
      return;
    }

    // Case 3: Removing the rightmost vertex shifts the unbounded track boundary leftward
    if (! has_right_neighbor) {
      Vertex_descriptor v_left_neighbor = m_topology_traits.left_vertex(e_left);

      // e_left becomes the new right-unbounded tracking edge
      m_topology_traits.clear_right_vertex(e_left);
      m_topology_traits.set_right_edge(v_left_neighbor, e_left);

      m_topology_traits.erase_vertex(v);
      m_topology_traits.erase_edge(e_right);
      return;
    }

    // Case 4: Internal Node removal. Stitch the gap by building a unified replacement edge.
    Vertex_descriptor v_left_neighbor = m_topology_traits.left_vertex(e_left);
    Vertex_descriptor v_right_neighbor = m_topology_traits.right_vertex(e_right);

    Edge_descriptor e_new = m_topology_traits.create_edge();
    m_topology_traits.set_left_vertex(e_new, v_left_neighbor);
    m_topology_traits.set_right_vertex(e_new, v_right_neighbor);

    m_topology_traits.set_right_edge(v_left_neighbor, e_new);
    m_topology_traits.set_left_edge(v_right_neighbor, e_new);

    // Clean up old topology structures from memory
    m_topology_traits.erase_vertex(v);
    m_topology_traits.erase_edge(e_left);
    m_topology_traits.erase_edge(e_right);
  }
};

} // namespace Arrangement_on_curve_1
} // namespace CGAL

#endif
