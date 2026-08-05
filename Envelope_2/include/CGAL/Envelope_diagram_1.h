// Copyright (c) 2011  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Ron Wein   <wein@post.tau.ac.il>
//            Efi Fogel  <efif@post.tau.ac.il>

#ifndef CGAL_ENVELOPE_DIAGRAM_1_H
#define CGAL_ENVELOPE_DIAGRAM_1_H

#include <CGAL/license/Envelope_2.h>

#include <list>
#include <memory>
#include <CGAL/basic.h>
#include <CGAL/memory.h>

#include <CGAL/Arrangement_on_curve_1/Arrangement_on_curve_1.h>
#include <CGAL/Arrangement_on_curve_1/Geom_traits_2_adaptor_1.h>
#include <CGAL/Arrangement_on_curve_1/Unbounded_topology_traits.h>

namespace CGAL {

/*! \class Envelope_diagram_1
 * A minimization (or maximization) diagram representing the lower (or upper)
 * envelope of a set of curves in the plane, inheriting directly from Arrangement_on_curve_1.
 */
template <typename Traits_, typename Allocator = CGAL_ALLOCATOR(int)>
class Envelope_diagram_1 :
  public Arrangement_on_curve_1::Arrangement_on_curve_1<
    Arrangement_on_curve_1::Geom_traits_2_adaptor_1<Traits_>,
    Arrangement_on_curve_1::Unbounded_topology_traits<
      typename Traits_::Point_2,
      std::list<typename Traits_::X_monotone_curve_2>,
      std::list<typename Traits_::X_monotone_curve_2>,
      Allocator, false, false>> {
public:
  using Traits_2 = Traits_;
  using Point_2 = typename Traits_2::Point_2;
  using X_monotone_curve_2 = typename Traits_2::X_monotone_curve_2;
  using Curve_container = std::list<X_monotone_curve_2>;
  using Curve_iterator = typename Curve_container::iterator;
  using Curve_const_iterator = typename Curve_container::const_iterator;
  using Size = std::size_t;

  using Geom_traits_1 = Arrangement_on_curve_1::Geom_traits_2_adaptor_1<Traits_2>;
  using Topol_traits =
    Arrangement_on_curve_1::Unbounded_topology_traits<Point_2, Curve_container, Curve_container, Allocator,
                                                      false, false>;

  using Base = Arrangement_on_curve_1::Arrangement_on_curve_1<Geom_traits_1, Topol_traits>;

  // Descriptors are the native handles used everywhere
  using Vertex_handle = typename Base::Vertex_descriptor;
  using Vertex_const_handle = typename Base::Vertex_const_descriptor;
  using Edge_handle = typename Base::Edge_descriptor;
  using Edge_const_handle = typename Base::Edge_const_descriptor;

public:
  /*! Default constructor.
   * Constructs an empty diagram with a default-constructed geometry traits.
   * Used when a diagram is needed before the traits is known (rare).
   */
  Envelope_diagram_1() : Base() {}

  /*! Constructor passing a 2D geometry traits shared pointer.
   * This is the preferred constructor; it avoids a heap allocation by
   * sharing the traits object rather than copying it.
   */
  Envelope_diagram_1(std::shared_ptr<const Traits_2> traits_2_ptr) :
    Base(std::make_shared<const Geom_traits_1>(traits_2_ptr))
  {}

  /*! Constructor passing a 2D geometry traits adaptor shared pointer.
   * This is the preferred constructor; it avoids a heap allocation by
   * sharing the traits object rather than copying it.
   */
  Envelope_diagram_1(std::shared_ptr<const Geom_traits_1> traits_2_ptr) : Base(traits_2_ptr) {}

  Envelope_diagram_1(const Envelope_diagram_1&) = default;
  Envelope_diagram_1& operator=(const Envelope_diagram_1&) = default;

  Envelope_diagram_1(Envelope_diagram_1&&) noexcept = default;
  Envelope_diagram_1& operator=(Envelope_diagram_1&&) noexcept = default;

  ~Envelope_diagram_1() = default;

  // --------------------------------------------------------------------------
  // DIAGRAM BOUNDARY ACCESSORS
  // --------------------------------------------------------------------------
  Edge_const_handle leftmost() const { return this->unbounded_left_edge(); }
  Edge_handle leftmost() { return this->unbounded_left_edge(); }

  Edge_const_handle rightmost() const { return this->unbounded_right_edge(); }
  Edge_handle rightmost() { return this->unbounded_right_edge(); }

  void set_leftmost(Edge_handle e) { this->topology_traits().set_unbounded_left_edge(e); }
  void set_rightmost(Edge_handle e) { this->topology_traits().set_unbounded_right_edge(e); }

  void clear() {
    // Reset the topology traits in-place (clears all vertex/edge lists and
    // reinitialises the single unbounded edge) without touching the shared
    // geometry traits pointer.
    // This is strictly faster than Base::operator=(Base(...)) which would
    // copy the shared_ptr (atomic increment + later decrement) and also
    // call make_shared<Geom_traits_1> in the default-constructor path.
    this->topology_traits() = Topol_traits{};
  }

  // --------------------------------------------------------------------------
  // CURVE DATA ACCESSORS & MUTATORS
  // --------------------------------------------------------------------------
  const Point_2& point(Vertex_const_handle v) const { return get(this->vertex_point_map(), v); }

  void set_point(Vertex_handle v, const Point_2& p) { put(this->vertex_point_map(), v, p); }

  Curve_container& vertex_curves(Vertex_handle v) { return get(this->vertex_data_map(), v); }
  const Curve_container& vertex_curves(Vertex_const_handle v) const { return get(this->vertex_data_map(), v); }

  Curve_container& edge_curves(Edge_handle e) { return get(this->edge_data_map(), e); }
  const Curve_container& edge_curves(Edge_const_handle e) const { return get(this->edge_data_map(), e); }

  Size number_of_vertex_curves(Vertex_const_handle v) const { return vertex_curves(v).size(); }
  Size number_of_edge_curves(Edge_const_handle e) const { return edge_curves(e).size(); }

  bool empty_vertex_curves(Vertex_const_handle v) const { return vertex_curves(v).empty(); }
  bool empty_edge_curves(Edge_const_handle e) const { return edge_curves(e).empty(); }

  const X_monotone_curve_2& vertex_curve(Vertex_const_handle v) const {
    CGAL_precondition(! empty_vertex_curves(v));
    return vertex_curves(v).front();
  }

  const X_monotone_curve_2& edge_curve(Edge_const_handle e) const {
    CGAL_precondition(! empty_edge_curves(e));
    return edge_curves(e).front();
  }

  void add_edge_curve(Edge_handle e, const X_monotone_curve_2& cv) { edge_curves(e).push_back(cv); }

  void add_vertex_curve(Vertex_handle v, const X_monotone_curve_2& cv) { vertex_curves(v).push_back(cv); }

  void add_edge_curves(Edge_handle e, Curve_const_iterator begin, Curve_const_iterator end) {
    auto& c = edge_curves(e);
    for (auto it = begin; it != end; ++it) c.push_back(*it);
  }

  void add_vertex_curves(Vertex_handle v, Curve_const_iterator begin, Curve_const_iterator end) {
    auto& c = vertex_curves(v);
    for (auto it = begin; it != end; ++it) c.push_back(*it);
  }

  void clear_vertex_curves(Vertex_handle v) { vertex_curves(v).clear(); }
  void clear_edge_curves(Edge_handle e) { edge_curves(e).clear(); }

  // --------------------------------------------------------------------------
  // TOPOLOGY MUTATION HELPER DELEGATES
  // --------------------------------------------------------------------------
  Vertex_handle new_vertex(const Point_2& p) { return this->topology_traits().create_vertex(p); }

  Edge_handle new_edge() { return this->topology_traits().create_edge(); }

  void delete_vertex(Vertex_handle v) { this->topology_traits().erase_vertex(v); }

  void delete_edge(Edge_handle e) { this->topology_traits().erase_edge(e); }
};

} // namespace CGAL

#endif
