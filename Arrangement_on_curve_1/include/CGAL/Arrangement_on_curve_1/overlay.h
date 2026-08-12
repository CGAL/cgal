// Copyright (c) 2026 Tel-Aviv University (Israel).
// All rights reserved.
// This file is part of CGAL (www.cgal.org).
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
// Author(s): Efi Fogel         <efifogel@gmail.com>

#ifndef CGAL_ARRANGEMENT_ON_CURVE_1_OVERLAY_H
#define CGAL_ARRANGEMENT_ON_CURVE_1_OVERLAY_H

#include <variant>
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

  auto comp = arr_res.geometry_traits_1().compare_x_1_object();

  auto e_a = arr_a.unbounded_left_edge();
  auto e_b = arr_b.unbounded_left_edge();
  auto e_res = arr_res.unbounded_left_edge();
  visitor.create_edge_ee(e_a, e_b, e_res);
  while (arr_a.has_right_vertex(e_a) && arr_b.has_right_vertex(e_b)) {
    auto v_a = arr_a.right_vertex(e_a);
    auto p_a = get(arr_a.vertex_point_map(), v_a);
    auto v_b = arr_b.right_vertex(e_b);
    auto p_b = get(arr_b.vertex_point_map(), v_b);
    auto cmp = comp(p_a, p_b);

    if (cmp == EQUAL) {
      auto v_res = arr_res.split_edge(e_res, p_a);
      visitor.create_vertex_vv(v_a, v_b, v_res);
      e_a = arr_a.right_edge(v_a);
      e_b = arr_b.right_edge(v_b);
      e_res = arr_res.right_edge(v_res);
      visitor.create_edge_ee(e_a, e_b, e_res);
      continue;
    }
    else if (cmp == SMALLER) {
      auto v_res = arr_res.split_edge(e_res, p_a);
      visitor.create_vertex_ve(v_a, e_b, v_res);
      e_a = arr_a.right_edge(v_a);
      e_res = arr_res.right_edge(v_res);
      visitor.create_edge_ee(e_a, e_b, e_res);
      continue;
    }
    CGAL_assertion(cmp == LARGER);
    auto v_res = arr_res.split_edge(e_res, p_b);
    visitor.create_vertex_ev(e_a, v_b, v_res);
    e_b = arr_b.right_edge(v_b);
    e_res = arr_res.right_edge(v_res);
    visitor.create_edge_ee(e_a, e_b, e_res);
  }

  while (arr_a.has_right_vertex(e_a)) {
    auto v_a = arr_a.right_vertex(e_a);
    auto p_a = get(arr_a.vertex_point_map(), v_a);
    auto v_res = arr_res.split_edge(e_res, p_a);
    visitor.create_vertex_ve(v_a, e_b, v_res);
    e_a = arr_a.right_edge(v_a);
    e_res = arr_res.right_edge(v_res);
    visitor.create_edge_ee(e_a, e_b, e_res);
  }

  while (arr_b.has_right_vertex(e_b)) {
    auto v_b = arr_b.right_vertex(e_b);
    auto p_b = get(arr_b.vertex_point_map(), v_b);
    auto v_res = arr_res.split_edge(e_res, p_b);
    visitor.create_vertex_ev(e_a, v_b, v_res);
    e_b = arr_b.right_edge(v_b);
    e_res = arr_res.right_edge(v_res);
    visitor.create_edge_ee(e_a, e_b, e_res);
  }
}

// ============================================================================
// Overlay without visitor; uses Default_overlay_visitor)
// ============================================================================
template <typename GeometryTraitsA, typename GeometryTraitsB, typename GeometryTraitsRes,
          typename TopologyTraitsA, typename TopologyTraitsB, typename TopologyTraitsRes,
          bool BinarySearchA, bool BinarySearchB, bool BinarySearchRes>
void overlay(const Arrangement_on_curve_1<GeometryTraitsA, TopologyTraitsA, BinarySearchA>& arr_a,
             const Arrangement_on_curve_1<GeometryTraitsB, TopologyTraitsB, BinarySearchB>& arr_b,
             Arrangement_on_curve_1<GeometryTraitsRes, TopologyTraitsRes, BinarySearchRes>& arr_res) {
  using Arr_a = Arrangement_on_curve_1<GeometryTraitsA, TopologyTraitsA>;
  using Arr_b = Arrangement_on_curve_1<GeometryTraitsB, TopologyTraitsB>;
  using Arr_res = Arrangement_on_curve_1<GeometryTraitsRes, TopologyTraitsRes>;

  Default_overlay_visitor<Arr_a, Arr_b, Arr_res> obs(arr_a, arr_b, arr_res);
  overlay(arr_a, arr_b, arr_res, obs);
}

} // namespace Arrangement_on_curve_1
} // namespace CGAL

#endif
