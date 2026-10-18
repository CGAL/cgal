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
  constexpr bool identified = ArrRes::is_identified_curve;

  auto cmp_x = arr_res.geometry_traits_1().compare_x_1_object();
  auto is_on_id = [&](const auto& pt) -> bool {
    if constexpr (identified) return arr_res.geometry_traits_1().is_on_identification_1_object()(pt);
    else return false;
  };

  auto lb_a = arr_a.left_boundary_edge();
  auto lb_b = arr_b.left_boundary_edge();
  auto lb_res = arr_res.left_boundary_edge();

  // In an identified arrangement the identification vertex, if present, is
  // the left vertex of left_boundary_edge().
  auto has_id_vertex = [&](const auto& arr, auto lb) -> bool {
    if constexpr (identified)
      return arr.has_left_vertex(lb) && is_on_id(get(arr.vertex_point_map(), arr.left_vertex(lb)));
    else return false;
  };
  const bool id_in_a = has_id_vertex(arr_a, lb_a);
  const bool id_in_b = has_id_vertex(arr_b, lb_b);

  // Walk each input rightward starting from left_boundary_edge().  The walk
  // visits exactly the non-identification vertices, in increasing order, so
  // the number of steps is known in advance.  (A count is needed because in
  // an identified arrangement every edge has a right vertex.)
  std::size_t n_a = arr_a.number_of_vertices() - (id_in_a ? 1 : 0);
  std::size_t n_b = arr_b.number_of_vertices() - (id_in_b ? 1 : 0);
  auto e_a = lb_a;
  auto e_b = lb_b;
  auto e_res = lb_res;
  // In open / bounded arrangements the walk ends at an edge without a right vertex.
  auto more_a = [&]() -> bool {
    if constexpr (identified) return n_a > 0;
    else return arr_a.has_right_vertex(e_a);
  };
  auto more_b = [&]() -> bool {
    if constexpr (identified) return n_b > 0;
    else return arr_b.has_right_vertex(e_b);
  };

  visitor.create_edge_ee(e_a, e_b, e_res);

  // Inserts p (larger than all vertices inserted so far) into arr_res.
  // Pre: (e_a, e_b) are the input edges that contain the edge left of p.
  auto insert_res = [&](const auto& p) {
    if constexpr (identified) {
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

  while (more_a() || more_b()) {
    Comparison_result res;
    if (! more_a()) res = LARGER;
    else if (! more_b()) res = SMALLER;
    else res = cmp_x(get(arr_a.vertex_point_map(), arr_a.right_vertex(e_a)),
                     get(arr_b.vertex_point_map(), arr_b.right_vertex(e_b)));
    if (res == EQUAL) {
      auto v_a = arr_a.right_vertex(e_a);
      auto v_b = arr_b.right_vertex(e_b);
      auto v_res = insert_res(get(arr_a.vertex_point_map(), v_a));
      visitor.create_vertex_vv(v_a, v_b, v_res);
      e_a = arr_a.right_edge(v_a); --n_a;
      e_b = arr_b.right_edge(v_b); --n_b;
      visitor.create_edge_ee(e_a, e_b, arr_res.right_edge(v_res));
    }
    else if (res == SMALLER) {
      auto v_a = arr_a.right_vertex(e_a);
      auto v_res = insert_res(get(arr_a.vertex_point_map(), v_a));
      visitor.create_vertex_ve(v_a, e_b, v_res);
      e_a = arr_a.right_edge(v_a); --n_a;
      visitor.create_edge_ee(e_a, e_b, arr_res.right_edge(v_res));
    }
    else {
      auto v_b = arr_b.right_vertex(e_b);
      auto v_res = insert_res(get(arr_b.vertex_point_map(), v_b));
      visitor.create_vertex_ev(e_a, v_b, v_res);
      e_b = arr_b.right_edge(v_b); --n_b;
      visitor.create_edge_ee(e_a, e_b, arr_res.right_edge(v_res));
    }
  }

  // The identification vertex (if in either input) is inserted last.  Now
  // (e_a, e_b) are the input edges that contain the edge left of it.
  if constexpr (identified) {
    if (! id_in_a && ! id_in_b) return;
    const auto& p_id = id_in_a ?
      get(arr_a.vertex_point_map(), arr_a.left_vertex(lb_a)) : get(arr_b.vertex_point_map(), arr_b.left_vertex(lb_b));
    auto v_res = arr_res.empty() ? arr_res.insert_empty(p_id) : arr_res.insert_identification_vertex(p_id);
    if (! arr_res.empty() && arr_res.number_of_vertices() > 1)
      visitor.create_edge_ee(e_a, e_b, arr_res.left_edge(v_res));
    if (id_in_a && id_in_b) visitor.create_vertex_vv(arr_a.left_vertex(lb_a), arr_b.left_vertex(lb_b), v_res);
    else if (id_in_a) visitor.create_vertex_ve(arr_a.left_vertex(lb_a), e_b, v_res);
    else visitor.create_vertex_ev(e_a, arr_b.left_vertex(lb_b), v_res);
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
