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

  using Vertex_const_descriptor = typename Topology_traits::Vertex_const_descriptor;
  using Edge_const_descriptor = typename Topology_traits::Edge_const_descriptor;

  // Location result: a point is either ON a vertex or INSIDE an edge
  using Location_result = std::variant<Vertex_descriptor, Edge_descriptor>;
  using Const_location_result = std::variant<Vertex_const_descriptor, Edge_const_descriptor>;

private:
  Geometry_traits_1 m_geometry_traits;
  Topology_traits m_topology_traits;

public:
  Arrangement_on_curve_1(const Geometry_traits_1& geom_tr) :
    m_geometry_traits(geom_tr),
    m_topology_traits()
  {}

  Arrangement_on_curve_1(const Geometry_traits_1& geom_tr, const Topology_traits& topo_tr) :
    m_geometry_traits(geom_tr),
    m_topology_traits(topo_tr)
  {}

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

  Edge_descriptor unbounded_edge() { return m_topology_traits.unbounded_edge(); }
  Edge_const_descriptor unbounded_edge() const { return m_topology_traits.unbounded_edge(); }

  Vertex_descriptor left_vertex(Edge_descriptor e) { return m_topology_traits.left_vertex(e); }
  Vertex_descriptor right_vertex(Edge_descriptor e) { return m_topology_traits.right_vertex(e); }
  Edge_descriptor left_edge(Vertex_descriptor v) { return m_topology_traits.left_edge(v); }
  Edge_descriptor right_edge(Vertex_descriptor v) { return m_topology_traits.right_edge(v); }

  Vertex_const_descriptor left_vertex(Edge_const_descriptor e) const { return m_topology_traits.left_vertex(e); }
  Vertex_const_descriptor right_vertex(Edge_const_descriptor e) const { return m_topology_traits.right_vertex(e); }
  Edge_const_descriptor left_edge(Vertex_const_descriptor v) const{ return m_topology_traits.left_edge(v); }
  Edge_const_descriptor right_edge(Vertex_const_descriptor v) const { return m_topology_traits.right_edge(v); }

  bool has_left_vertex(Edge_const_descriptor e) const { return m_topology_traits.has_left_vertex(e); }
  bool has_right_vertex(Edge_const_descriptor e) const { return m_topology_traits.has_right_vertex(e); }

  // ============================================================================
  // HIGH-LEVEL TOPOLOGICAL OPERATIONS
  // ============================================================================

  // Insert the very first vertex into an empty arrangement.
  // The single unbounded edge is split into: (-inf, v) and (v, +inf).
  Vertex_descriptor insert_empty(const Point_1& p) {
    auto& topo = m_topology_traits;

    // The arrangement must be empty: one unbounded edge exists.
    auto e_unbounded = topo.unbounded_edge();

    // Create the new vertex.
    Vertex_descriptor v = topo.create_vertex(p);

    // Create the right new edge: (v, +inf)
    Edge_descriptor e_right = topo.create_edge();

    // Reuse e_unbounded as the left edge: (-inf, v)
    // It already has no left vertex; set its right vertex to v.
    topo.set_right_vertex(e_unbounded, v);

    // The right edge has v as its left vertex, and no right vertex.
    topo.set_left_vertex(e_right, v);

    // Wire the vertex to its two adjacent edges.
    topo.set_left_edge(v, e_unbounded);
    topo.set_right_edge(v, e_right);

    return v;
  }

  // Insert a new vertex p strictly to the left of an existing first vertex v_first.
  // Splits the unbounded left edge of v_first.
  Vertex_descriptor insert_before(Vertex_descriptor v_first, const Point_1& p) {
    auto& topo = m_topology_traits;

    // e_left is the unbounded edge to the left of v_first: (-inf, v_first)
    Edge_descriptor e_left = topo.left_edge(v_first);

    // Create the new vertex.
    Vertex_descriptor v_new = topo.create_vertex(p);

    // Create a new edge to sit between v_new and v_first: (v_new, v_first)
    Edge_descriptor e_between = topo.create_edge();
    topo.set_left_vertex(e_between, v_new);
    topo.set_right_vertex(e_between, v_first);

    // The old left unbounded edge now terminates at v_new on its right: (-inf, v_new)
    topo.set_right_vertex(e_left, v_new);

    // v_new sits between e_left and e_between.
    topo.set_left_edge(v_new, e_left);
    topo.set_right_edge(v_new, e_between);

    // v_first's left edge is now e_between.
    topo.set_left_edge(v_first, e_between);

    return v_new;
  }

  // Insert a new vertex p strictly to the right of an existing last vertex v_last.
  // Splits the unbounded right edge of v_last.
  Vertex_descriptor insert_after(Vertex_descriptor v_last, const Point_1& p) {
    auto& topo = m_topology_traits;

    // e_right is the unbounded edge to the right of v_last: (v_last, +inf)
    Edge_descriptor e_right = topo.right_edge(v_last);

    // Create the new vertex.
    Vertex_descriptor v_new = topo.create_vertex(p);

    // Create a new edge between v_last and v_new: (v_last, v_new)
    Edge_descriptor e_between = topo.create_edge();
    topo.set_left_vertex(e_between, v_last);
    topo.set_right_vertex(e_between, v_new);

    // The old right unbounded edge now starts from v_new: (v_new, +inf)
    topo.set_left_vertex(e_right, v_new);

    // v_new is wired between e_between and e_right.
    topo.set_left_edge(v_new, e_between);
    topo.set_right_edge(v_new, e_right);

    // v_last's right edge is now e_between.
    topo.set_right_edge(v_last, e_between);

    return v_new;
  }

  // Split the edge e by inserting a new vertex p inside it.
  // e becomes the left sub-edge; a new edge becomes the right sub-edge.
  Vertex_descriptor split_edge(Edge_descriptor e, const Point_1& p) {
    auto& topo = m_topology_traits;

    // Create the new vertex.
    Vertex_descriptor v_new = topo.create_vertex(p);

    // Create a new right sub-edge.
    Edge_descriptor e_right = topo.create_edge();

    // If e had a right vertex, transfer it to e_right.
    if (topo.has_right_vertex(e)) {
      Vertex_descriptor v_old_right = topo.right_vertex(e);
      topo.set_right_vertex(e_right, v_old_right);
      topo.set_left_edge(v_old_right, e_right);
    }
    // e_right's left vertex is v_new.
    topo.set_left_vertex(e_right, v_new);

    // e (the left sub-edge) now ends at v_new.
    topo.set_right_vertex(e, v_new);

    // Wire v_new between e and e_right.
    topo.set_left_edge(v_new, e);
    topo.set_right_edge(v_new, e_right);

    return v_new;
  }

  // Remove vertex v from the arrangement.
  // Merge its two adjacent edges into one, keeping the left edge and removing the right.
  void remove(Vertex_descriptor v) {
    auto& topo = m_topology_traits;

    Edge_descriptor e_left  = topo.left_edge(v);
    Edge_descriptor e_right = topo.right_edge(v);

    // Extend e_left to cover the span formerly covered by e_right.
    if (topo.has_right_vertex(e_right)) {
      Vertex_descriptor v_right = topo.right_vertex(e_right);
      topo.set_right_vertex(e_left, v_right);
      topo.set_left_edge(v_right, e_left);
    }
    else {
      topo.clear_right_vertex(e_left);
    }

    // Remove the now-redundant right edge and the vertex.
    topo.erase_edge(e_right);
    topo.erase_vertex(v);
  }
};

} // namespace Arrangement_on_curve_1
} // namespace CGAL

#endif
