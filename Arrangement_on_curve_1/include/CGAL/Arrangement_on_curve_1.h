// Copyright (c) 2026 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s): Efi Fogel         <efif@post.tau.ac.il>

#ifndef CGAL_ARRANGEMENT_ON_CURVE_1_H
#define CGAL_ARRANGEMENT_ON_CURVE_1_H

#include <list>
#include <type_traits>
#include <variant>

#include <CGAL/basic.h>

namespace CGAL {
namespace Arrangement_on_curve_1 {

// Primary template for handling user extension data attributes
template <typename Data>
struct Data_container {
  Data m_data;
  Data_container() = default;
  Data_container(const Data& data) : m_data(data) {}
  Data& data() { return m_data; }
  const Data& data() const { return m_data; }
};

// Explicit specialization for 'void' to leverage empty class optimization
template <>
struct Data_container<void> {
  // Intentionally left empty: no fields are constructed when Data is void
};

template <typename Traits_, typename VertexData = void, typename EdgeData = void>
class Arrangement_on_curve_1 {
public:
  using Traits = Traits_;
  using Point_1 = typename Traits::Point_1;

  class Vertex;
  class Edge;

  using Vertex_list = std::list<Vertex>;
  using Edge_list = std::list<Edge>;

  using Vertex_handle = typename Vertex_list::iterator;
  using Edge_handle = typename Edge_list::iterator;
  using Vertex_const_handle = typename Vertex_list::const_iterator;
  using Edge_const_handle = typename Edge_list::const_iterator;

  using Location_result = std::variant<Vertex_handle, Edge_handle, void*>;
  using Const_location_result = std::variant<Vertex_const_handle, Edge_const_handle, void*>;

  class Vertex : public Data_container<VertexData> {
  private:
    Point_1 m_point;
    Edge_handle m_left_edge;
    Edge_handle m_right_edge;

  public:
    Vertex(const Point_1& pt) : m_point(pt) {}

    const Point_1& point() const { return m_point; }
    void set_point(const Point_1& pt) { m_point = pt; }

    Edge_handle left() { return m_left_edge; }
    Edge_handle right() { return m_right_edge; }
    Edge_const_handle left() const { return m_left_edge; }
    Edge_const_handle right() const { return m_right_edge; }

    void set_left(Edge_handle e) { m_left_edge = e; }
    void set_right(Edge_handle e) { m_right_edge = e; }
  };

  class Edge : public Data_container<EdgeData> {
  private:
    Vertex_handle m_source_v;
    Vertex_handle m_target_v;
    bool m_has_source = false;
    bool m_has_target = false;

  public:
    Edge() = default;

    Vertex_handle source() { return m_source_v; }
    Vertex_handle target() { return m_target_v; }
    Vertex_const_handle source() const { return m_source_v; }
    Vertex_const_handle target() const { return m_target_v; }

    bool has_source() const { return m_has_source; }
    bool has_target() const { return m_has_target; }

    void set_source(Vertex_handle v) { m_source_v = v; m_has_source = true; }
    void set_target(Vertex_handle v) { m_target_v = v; m_has_target = true; }
  };

private:
  Traits m_traits;
  Vertex_list m_vertices;
  Edge_list m_edges;

public:
  Arrangement_on_curve_1(const Traits& tr = Traits()) : m_traits(tr) {}

  const Traits& traits() const { return m_traits; }

  typename Vertex_list::iterator vertices_begin() { return m_vertices.begin(); }
  typename Vertex_list::iterator vertices_end() { return m_vertices.end(); }
  typename Vertex_list::const_iterator vertices_begin() const { return m_vertices.begin(); }
  typename Vertex_list::const_iterator vertices_end() const { return m_vertices.end(); }

  typename Edge_list::iterator edges_begin() { return m_edges.begin(); }
  typename Edge_list::iterator edges_end() { return m_edges.end(); }
  typename Edge_list::const_iterator edges_begin() const { return m_edges.begin(); }
  typename Edge_list::const_iterator edges_end() const { return m_edges.end(); }

  bool is_empty() const { return m_vertices.empty(); }
  size_t number_of_vertices() const { return m_vertices.size(); }
  size_t number_of_edges() const { return m_edges.size(); }

  // Handles inserting into empty space and sets up the initial infinite left/right edges
  Vertex_handle insert_empty(const Point_1& p) {
    m_vertices.emplace_back(p);
    Vertex_handle v = std::prev(m_vertices.end());

    m_edges.emplace_back(); // Infinite left edge
    Edge_handle e_left = std::prev(m_edges.end());
    e_left->set_target(v);

    m_edges.emplace_back(); // Infinite right edge
    Edge_handle e_right = std::prev(m_edges.end());
    e_right->set_source(v);

    v->set_left(e_left);
    v->set_right(e_right);
    return v;
  }

  // Inserts a new vertex to the left of the entire structure
  Vertex_handle insert_before(Vertex_handle v_first, const Point_1& p) {
    m_vertices.emplace_front(p);
    Vertex_handle v_new = m_vertices.begin();

    Edge_handle e_left_unbounded = v_first->left();
    e_left_unbounded->set_target(v_new);

    m_edges.emplace_back();
    Edge_handle e_new = std::prev(m_edges.end());
    e_new->set_source(v_new);
    e_new->set_target(v_first);

    v_new->set_left(e_left_unbounded);
    v_new->set_right(e_new);
    v_first->set_left(e_new);

    return v_new;
  }

  // Inserts a new vertex to the right of the entire structure
  Vertex_handle insert_after(Vertex_handle v_last, const Point_1& p) {
    m_vertices.emplace_back(p);
    Vertex_handle v_new = std::prev(m_vertices.end());

    Edge_handle e_right_unbounded = v_last->right();
    e_right_unbounded->set_source(v_new);

    m_edges.emplace_back();
    Edge_handle e_new = std::prev(m_edges.end());
    e_new->set_source(v_last);
    e_new->set_target(v_new);

    v_last->set_right(e_new);
    v_new->set_left(e_new);
    v_new->set_right(e_right_unbounded);

    return v_new;
  }

  Vertex_handle split_edge(Edge_handle e, const Point_1& p) {
    // If the edge lacks a source or a target, it's an unbounded outer edge
    if (!e->has_source()) return insert_before(e->target(), p);
    if (!e->has_target()) return insert_after(e->source(), p);

    // Standard internal edge split logic
    Vertex_handle v_source = e->source();
    Vertex_handle v_target = e->target();

    // Find correct linear position in std::list to keep vertices sorted
    typename Vertex_list::iterator place_it = v_target;
    m_vertices.emplace(place_it, p);
    Vertex_handle v_new = std::prev(place_it);

    m_edges.emplace_back();
    Edge_handle e_new = std::prev(m_edges.end());

    e_new->set_source(v_new);
    e_new->set_target(v_target);
    e->set_target(v_new);

    v_new->set_left(e);
    v_new->set_right(e_new);
    v_target->set_left(e_new);

    return v_new;
  }
};

} // namespace Arrangement_on_curve_1
} // namespace CGAL

#endif
