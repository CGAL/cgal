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

  while (vit1 != r1.end() && vit2 != r2.end()) {
    typename GeometryTraitsRes::Point_1 p1 = get(pmap1, vit1);
    typename GeometryTraitsRes::Point_1 p2 = get(pmap2, vit2);

    Comparison_result res = comp(p1, p2);
    if (res == EQUAL) {
      auto v_res = insert(result, p1);
      observer.create_vertex(vit1, vit2, v_res);
      ++vit1;
      ++vit2;
    }
    else
    if (res == SMALLER) {
      auto v_res = insert(result, p1);
      observer.create_vertex_from_a(vit1, v_res);
      ++vit1;
    }
    else {
      auto v_res = insert(result, p2);
      observer.create_vertex_from_b(vit2, v_res);
      ++vit2;
    }
  }

  while (vit1 != r1.end()) {
    typename GeometryTraitsRes::Point_1 p = get(pmap1, vit1);
    auto v_res = insert(result, p);
    observer.create_vertex_from_a(vit1, v_res);
    ++vit1;
  }
  while (vit2 != r2.end()) {
    typename GeometryTraitsRes::Point_1 p = get(pmap2, vit2);
    auto v_res = insert(result, p);
    observer.create_vertex_from_b(vit2, v_res);
    ++vit2;
  }

  // Note: For a 1D arrangement, edges are completely implied by the sequential
  // vertex chain ordering. Edge modification attributes can be parsed by the
  // observer if necessary by checking the resulting arrangement boundaries.
}

} // namespace Arrangement_on_curve_1
} // namespace CGAL

#endif
