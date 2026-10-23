// Copyright (c) 2026 Tel-Aviv University (Israel).
// All rights reserved.
// This file is part of CGAL (www.cgal.org).
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
// Author(s): Efi Fogel         <efifogel@gmail.com>

#ifndef CGAL_ARRANGEMENT_ON_CURVE_1_OVERLAY_H
#define CGAL_ARRANGEMENT_ON_CURVE_1_OVERLAY_H

#include <variant>
#include <iostream>
#include <iterator>
#include <type_traits>

#include <CGAL/Arrangement_on_curve_1/Arrangement_on_curve_1.h>
#include <CGAL/Arrangement_on_curve_1/Default_overlay_visitor.h>
#include <CGAL/property_map.h>

namespace CGAL {
namespace Arrangement_on_curve_1 {

// ============================================================================
// Overlay with overlay visitor
// ============================================================================
template <typename GeometryTraitsA, typename GeometryTraitsB, typename GeometryTraitsRes,
          typename TopologyTraitsA, typename TopologyTraitsB, typename TopologyTraitsRes,
          bool BinarySearchA, bool BinarySearchB, bool BinarySearchRes,
          typename OverlayVisitor>
void overlay(const Arrangement_on_curve_1<GeometryTraitsA, TopologyTraitsA, BinarySearchA>& arr_a,
             const Arrangement_on_curve_1<GeometryTraitsB, TopologyTraitsB, BinarySearchB>& arr_b,
             Arrangement_on_curve_1<GeometryTraitsRes, TopologyTraitsRes, BinarySearchRes>& arr_res,
             OverlayVisitor& visitor) {
  CGAL_precondition(arr_res.empty());

  // If the Result arrangement uses the exact same traits type as Input A,
  // we can share the pointer instance instead of keeping separate identical copies.
  if constexpr (std::is_same_v<GeometryTraitsA, GeometryTraitsRes>) {
    if (arr_res.empty() && arr_res.shared_geometry_traits_1() != arr_a.shared_geometry_traits_1()) {
      arr_res.reset_shared_geometry_traits(arr_a.shared_geometry_traits_1());
    }
  }
  else if constexpr (std::is_same_v<GeometryTraitsB, GeometryTraitsRes>) {
    if (arr_res.empty() && arr_res.shared_geometry_traits_1() != arr_b.shared_geometry_traits_1()) {
      arr_res.reset_shared_geometry_traits(arr_b.shared_geometry_traits_1());
    }
  }

  using ArrA = std::remove_cv_t<std::remove_reference_t<decltype(arr_a)>>;
  using ArrB = std::remove_cv_t<std::remove_reference_t<decltype(arr_b)>>;
  using ArrRes = std::remove_cv_t<std::remove_reference_t<decltype(arr_res)>>;
  static_assert(std::is_same_v<typename ArrA::Left_side_category, typename ArrRes::Left_side_category> &&
                std::is_same_v<typename ArrA::Right_side_category, typename ArrRes::Right_side_category> &&
                std::is_same_v<typename ArrB::Left_side_category, typename ArrRes::Left_side_category> &&
                std::is_same_v<typename ArrB::Right_side_category, typename ArrRes::Right_side_category>,
                "The input and result arrangements must have the same parameter space.");
  constexpr bool is_identified = ArrRes::is_identified;
  constexpr bool left_is_closed = ArrRes::left_is_closed;
  constexpr bool right_is_closed = ArrRes::right_is_closed;
  using Vertex_a = typename ArrA::Vertex_const_descriptor;
  using Vertex_b = typename ArrB::Vertex_const_descriptor;

  auto cmp_x = arr_res.geometry_traits_1().compare_x_1_object();

  auto lb_a = arr_a.left_boundary_edge();
  auto lb_b = arr_b.left_boundary_edge();
  auto lb_res = arr_res.left_boundary_edge();

  // The boundary vertices and the identification vertex are not visited by
  // the walk below; they are inserted separately.
  auto has_left_boundary_vertex = [](const auto& arr) -> bool {
    if constexpr (left_is_closed) return arr.topology_traits().has_left_boundary_vertex();
    else return false;
  };
  auto has_right_boundary_vertex = [](const auto& arr) -> bool {
    if constexpr (right_is_closed) return arr.topology_traits().has_right_boundary_vertex();
    else return false;
  };
  auto has_identification_vertex = [](const auto& arr) -> bool {
    if constexpr (is_identified) return arr.topology_traits().has_identification_vertex();
    else return false;
  };
  const bool lbv_in_a = has_left_boundary_vertex(arr_a);
  const bool lbv_in_b = has_left_boundary_vertex(arr_b);
  const bool rbv_in_a = has_right_boundary_vertex(arr_a);
  const bool rbv_in_b = has_right_boundary_vertex(arr_b);
  const bool id_in_a = has_identification_vertex(arr_a);
  const bool id_in_b = has_identification_vertex(arr_b);

  // Walk each input rightward starting from left_boundary_edge().  The walk
  // visits exactly the interior vertices, in increasing order, so the number
  // of steps is known in advance.  (A count is needed, because in an
  // identified arrangement every edge has a right vertex, and in a right-closed
  // arrangement the right boundary vertex is its own right edge's endpoint.)
  std::size_t n_a = arr_a.number_of_vertices() - (lbv_in_a ? 1 : 0) - (rbv_in_a ? 1 : 0) - (id_in_a ? 1 : 0);
  std::size_t n_b = arr_b.number_of_vertices() - (lbv_in_b ? 1 : 0) - (rbv_in_b ? 1 : 0) - (id_in_b ? 1 : 0);
  auto e_a = lb_a;
  auto e_b = lb_b;
  auto e_res = lb_res;

  // Inserts into arr_res the point of v_a, or of v_b if v_a does not exist,
  // using insert_fn, and reports the new vertex.  in_a and in_b indicate
  // whether v_a and v_b exist; when only one of them exists, the current edge
  // (e_a or e_b) of the other input contains the point.
  auto insert_and_report = [&](bool in_a, Vertex_a v_a, bool in_b, Vertex_b v_b, auto insert_fn) {
    auto v_res = in_a ? insert_fn(get(arr_a.vertex_point_map(), v_a)) : insert_fn(get(arr_b.vertex_point_map(), v_b));
    if (in_a && in_b) visitor.create_vertex_vv(v_a, v_b, v_res);
    else if (in_a) visitor.create_vertex_ve(v_a, e_b, v_res);
    else visitor.create_vertex_ev(e_a, v_b, v_res);
    return v_res;
  };

  // The left boundary vertex (if in either input) is inserted first.
  if constexpr (left_is_closed) {
    if (lbv_in_a || lbv_in_b) {
      Vertex_a v_a = lbv_in_a ? arr_a.topology_traits().left_boundary_vertex() : Vertex_a();
      Vertex_b v_b = lbv_in_b ? arr_b.topology_traits().left_boundary_vertex() : Vertex_b();
      insert_and_report(lbv_in_a, v_a, lbv_in_b, v_b,
                        [&](const auto& p) { return arr_res.insert_left_boundary_vertex(p); });
    }
  }

  // In an identified arrangement left_boundary_edge() of arr_res is the
  // wrap-around edge, which changes with every insertion; it is reported once
  // at the end.  Otherwise, it keeps its left part when split, so its final
  // pair of input edges is already known.
  if constexpr (! is_identified) visitor.create_edge_ee(e_a, e_b, e_res);

  // Inserts p (larger than all vertices inserted so far) into arr_res.
  // Pre: (e_a, e_b) are the input edges that contain the edge left of p.
  auto insert_res = [&](const auto& p) {
    if constexpr (is_identified) {
      // left_boundary_edge() of arr_res is the wrap-around edge (v_n, v_1); it
      // keeps the part (p, v_1), and the new edge (v_n, p) is to the left of p.
      if (arr_res.empty()) return arr_res.insert_empty(p);
      auto v_res = arr_res.insert_after(arr_res.left_vertex(lb_res), p);
      visitor.create_edge_ee(e_a, e_b, arr_res.left_edge(v_res));
      return v_res;
    }
    else {
      auto v_res = arr_res.split_edge(e_res, p);
      e_res = arr_res.right_edge(v_res);
      return v_res;
    }
  };

  // Reports the edge to the right of v_res.  In an identified arrangement this
  // is left_boundary_edge(), which is reported once at the end.
  // Pre: (e_a, e_b) are the input edges that contain it.
  auto report_right_edge = [&](auto v_res) {
    if constexpr (! is_identified) visitor.create_edge_ee(e_a, e_b, arr_res.right_edge(v_res));
  };

  while ((n_a > 0) || (n_b > 0)) {
    Comparison_result res;
    if (n_a == 0) res = LARGER;
    else if (n_b == 0) res = SMALLER;
    else res = cmp_x(get(arr_a.vertex_point_map(), arr_a.right_vertex(e_a)),
                     get(arr_b.vertex_point_map(), arr_b.right_vertex(e_b)));
    if (res == EQUAL) {
      auto v_a = arr_a.right_vertex(e_a);
      auto v_b = arr_b.right_vertex(e_b);
      auto v_res = insert_res(get(arr_a.vertex_point_map(), v_a));
      visitor.create_vertex_vv(v_a, v_b, v_res);
      e_a = arr_a.right_edge(v_a); --n_a;
      e_b = arr_b.right_edge(v_b); --n_b;
      report_right_edge(v_res);
    }
    else if (res == SMALLER) {
      auto v_a = arr_a.right_vertex(e_a);
      auto v_res = insert_res(get(arr_a.vertex_point_map(), v_a));
      visitor.create_vertex_ve(v_a, e_b, v_res);
      e_a = arr_a.right_edge(v_a); --n_a;
      report_right_edge(v_res);
    }
    else {
      auto v_b = arr_b.right_vertex(e_b);
      auto v_res = insert_res(get(arr_b.vertex_point_map(), v_b));
      visitor.create_vertex_ev(e_a, v_b, v_res);
      e_b = arr_b.right_edge(v_b); --n_b;
      report_right_edge(v_res);
    }
  }

  // The right boundary vertex (if in either input) is inserted last.  Now
  // (e_a, e_b) are the input edges that contain it, and e_res is the rightmost
  // edge of arr_res.  No new edge is created.
  if constexpr (right_is_closed) {
    if (rbv_in_a || rbv_in_b) {
      Vertex_a v_a = rbv_in_a ? arr_a.topology_traits().right_boundary_vertex() : Vertex_a();
      Vertex_b v_b = rbv_in_b ? arr_b.topology_traits().right_boundary_vertex() : Vertex_b();
      insert_and_report(rbv_in_a, v_a, rbv_in_b, v_b,
                        [&](const auto& p) { return arr_res.insert_right_boundary_vertex(e_res, p); });
    }
  }

  // The identification vertex (if in either input) is inserted last.  Now
  // (e_a, e_b) are the input edges that contain the edge left of it.
  if constexpr (is_identified) {
    if (id_in_a || id_in_b) {
      Vertex_a v_a = id_in_a ? arr_a.topology_traits().identification_vertex() : Vertex_a();
      Vertex_b v_b = id_in_b ? arr_b.topology_traits().identification_vertex() : Vertex_b();
      const bool was_empty = arr_res.empty();
      auto v_res = insert_and_report(id_in_a, v_a, id_in_b, v_b, [&](const auto& p) {
        return was_empty ? arr_res.insert_empty(p) : arr_res.insert_identification_vertex(p);
      });
      if (! was_empty) visitor.create_edge_ee(e_a, e_b, arr_res.left_edge(v_res));
    }
    // left_boundary_edge() of arr_res is final only now.  It is contained in
    // the left boundary edges of the inputs, whether or not they have an
    // identification vertex.
    visitor.create_edge_ee(lb_a, lb_b, lb_res);
  }
}

// ============================================================================
// Overlay without visitor; uses Default_overlay_visitor.
// ============================================================================
template <typename GeometryTraitsA, typename GeometryTraitsB, typename GeometryTraitsRes,
          typename TopologyTraitsA, typename TopologyTraitsB, typename TopologyTraitsRes,
          bool BinarySearchA, bool BinarySearchB, bool BinarySearchRes>
void overlay(const Arrangement_on_curve_1<GeometryTraitsA, TopologyTraitsA, BinarySearchA>& arr_a,
             const Arrangement_on_curve_1<GeometryTraitsB, TopologyTraitsB, BinarySearchB>& arr_b,
             Arrangement_on_curve_1<GeometryTraitsRes, TopologyTraitsRes, BinarySearchRes>& arr_res) {
  using Arr_a = Arrangement_on_curve_1<GeometryTraitsA, TopologyTraitsA, BinarySearchA>;
  using Arr_b = Arrangement_on_curve_1<GeometryTraitsB, TopologyTraitsB, BinarySearchB>;
  using Arr_res = Arrangement_on_curve_1<GeometryTraitsRes, TopologyTraitsRes, BinarySearchRes>;

  Default_overlay_visitor<Arr_a, Arr_b, Arr_res> obs(arr_a, arr_b, arr_res);
  overlay(arr_a, arr_b, arr_res, obs);
}

} // namespace Arrangement_on_curve_1
} // namespace CGAL

#endif
