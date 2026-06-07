// Copyright (c) 2026 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Efi Fogel         <efifogel@gmail.com>

#ifndef CGAL_ARRANGEMENT_ON_CURVE_1_OVERLAY_H
#define CGAL_ARRANGEMENT_ON_CURVE_1_OVERLAY_H

#include <variant>
#include <iterator>

#include <CGAL/Arrangement_on_curve_1.h>
#include <CGAL/Arrangement_on_curve_1_functions.h>
#include <CGAL/property_map.h>

namespace CGAL {
namespace Arrangement_on_curve_1 {

// [locate() and insert() implementations remain unchanged and are omitted for brevity]

// ============================================================================
// 1. OVERLAY (Fully Decoupled Signatures)
// ============================================================================
template <typename GeometryTraitsA, typename GeometryTraitsB, typename GeometryTraitsRes,
          typename TopologyTraitsA, typename TopologyTraitsB, typename TopologyTraitsRes>
void overlay(const Arrangement_on_curve_1<GeometryTraitsA, TopologyTraitsA>& arr1,
             const Arrangement_on_curve_1<GeometryTraitsB, TopologyTraitsB>& arr2,
             Arrangement_on_curve_1<GeometryTraitsRes, TopologyTraitsRes>& result) {
  auto comp = result.geometry_traits_1().compare_x_1_object();

  auto pmap1 = arr1.vertex_point_map();
  auto pmap2 = arr2.vertex_point_map();

  auto r1 = arr1.vertices();
  auto r2 = arr2.vertices();

  auto vit1 = r1.begin();
  auto vit2 = r2.begin();

  while (vit1 != r1.end() && vit2 != r2.end()) {
    typename GeometryTraitsRes::Point_1 p1 = get(pmap1, vit1);
    typename GeometryTraitsRes::Point_1 p2 = get(pmap2, vit2);

    Comparison_result res = comp(p1, p2);
    if (res == EQUAL) {
      insert(result, p1);
      ++vit1;
      ++vit2;
    }
    else
    if (res == SMALLER) {
      insert(result, p1);
      ++vit1;
    }
    else {
      insert(result, p2);
      ++vit2;
    }
  }

  while (vit1 != r1.end()) {
    insert(result, typename GeometryTraitsRes::Point_1(get(pmap1, vit1)));
    ++vit1;
  }
  while (vit2 != r2.end()) {
    insert(result, typename GeometryTraitsRes::Point_1(get(pmap2, vit2)));
    ++vit2;
  }
}

// ============================================================================
// 2. OVERLAY WITH OVERLAY OBSERVER (User-defined Data Consolidation)
// ============================================================================
template <typename GeometryTraitsA, typename GeometryTraitsB, typename GeometryTraitsRes,
          typename TopologyTraitsA, typename TopologyTraitsB, typename TopologyTraitsRes,
          typename OverlayObserver>
void overlay(const Arrangement_on_curve_1<GeometryTraitsA, TopologyTraitsA>& arr1,
             const Arrangement_on_curve_1<GeometryTraitsB, TopologyTraitsB>& arr2,
             Arrangement_on_curve_1<GeometryTraitsRes, TopologyTraitsRes>& result,
             OverlayObserver& observer) {
  auto comp = result.geometry_traits_1().compare_x_1_object();

  auto pmap1 = arr1.vertex_point_map();
  auto pmap2 = arr2.vertex_point_map();

  auto r1 = arr1.vertices();
  auto r2 = arr2.vertices();

  auto vit1 = r1.begin();
  auto vit2 = r2.begin();

  auto topo1 = arr1.topology_traits();
  auto topo2 = arr2.topology_traits();

  using Point_1 = typename GeometryTraitsRes::Point_1;

  using Edge_descriptor_a = typename TopologyTraitsA::Edge_descriptor;
  using Edge_descriptor_b = typename TopologyTraitsB::Edge_descriptor;
  using Vertex_descriptor_r = typename TopologyTraitsRes::Vertex_descriptor;
  using Edge_descriptor_r = typename TopologyTraitsRes::Edge_descriptor;

  // Track the current active edge of each arrangement.
  // Initially, before the first vertex, we are in the leftmost unbounded edge.
  Edge_descriptor_a active_e1 = (vit1 != r1.end()) ? topo1.left_edge(vit1) : Edge_descriptor_a();
  Edge_descriptor_b active_e2 = (vit2 != r2.end()) ? topo2.left_edge(vit2) : Edge_descriptor_b();

  Vertex_descriptor_r last_v_res = Vertex_descriptor_r();
  bool first_vertex = true;

  auto trigger_edge_overlay = [&](Edge_descriptor_a e1, Edge_descriptor_b e2, Vertex_descriptor_r curr_v_res) {
    if (! first_vertex) {
      auto topo_res = result.topology_traits();
      Edge_descriptor_r e_res = topo_res.left_edge(curr_v_res);
      observer.create_edge(e1, e2, e_res);
    }
  };

  while (vit1 != r1.end() && vit2 != r2.end()) {
    Point_1 p1 = get(pmap1, vit1);
    Point_1 p2 = get(pmap2, vit2);

    Comparison_result res = comp(p1, p2);
    if (res == EQUAL) {
      auto v_res = insert(result, p1);
      observer.create_vertex(vit1, vit2, v_res);
      trigger_edge_overlay(active_e1, active_e2, v_res);

      active_e1 = topo1.right_edge(vit1);
      active_e2 = topo2.right_edge(vit2);
      last_v_res = v_res; first_vertex = false;
      ++vit1; ++vit2;
    }
    else if (res == SMALLER) {
      auto v_res = insert(result, p1);
      // Vertex from A splits the active edge of B
      observer.create_vertex(vit1, active_e2, v_res);
      trigger_edge_overlay(active_e1, active_e2, v_res);

      active_e1 = topo1.right_edge(vit1);
      last_v_res = v_res; first_vertex = false;
      ++vit1;
    }
    else {
      auto v_res = insert(result, p2);
      // Vertex from B splits the active edge of A
      observer.create_vertex(active_e1, vit2, v_res);
      trigger_edge_overlay(active_e1, active_e2, v_res);

      active_e2 = topo2.right_edge(vit2);
      last_v_res = v_res; first_vertex = false;
      ++vit2;
    }
  }

  while (vit1 != r1.end()) {
    Point_1 p = get(pmap1, vit1);
    auto v_res = insert(result, p);
    observer.create_vertex(vit1, active_e2, v_res);
    trigger_edge_overlay(active_e1, active_e2, v_res);

    active_e1 = topo1.right_edge(vit1);
    last_v_res = v_res; first_vertex = false;
    ++vit1;
  }

  while (vit2 != r2.end()) {
    Point_1 p = get(pmap2, vit2);
    auto v_res = insert(result, p);
    observer.create_vertex(active_e1, vit2, v_res);
    trigger_edge_overlay(active_e1, active_e2, v_res);

    active_e2 = topo2.right_edge(vit2);
    last_v_res = v_res; first_vertex = false;
    ++vit2;
  }

  // Handle the final right-unbounded overlapping interval edge
  if (! result.is_empty()) {
    auto topo_res = result.topology_traits();
    auto last_vit = std::prev(result.vertices().end());
    Edge_descriptor_r final_e_res = topo_res.right_edge(last_vit);
    observer.create_edge(active_e1, active_e2, final_e_res);
  }
}

} // namespace Arrangement_on_curve_1
} // namespace CGAL

#endif
