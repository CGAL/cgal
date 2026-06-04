// Copyright (c) 2026 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Efi Fogel         <efifogel@gmail.com>

#ifndef CGAL_ARRANGEMENT_ON_CURVE_1_FUNCTIONS_H
#define CGAL_ARRANGEMENT_ON_CURVE_1_FUNCTIONS_H

#include <variant>
#include <iterator>

#include <CGAL/Arrangement_on_curve_1.h>
#include <CGAL/property_map.h>

namespace CGAL {
namespace Arrangement_on_curve_1 {

// ==========================================
// LOCATE (Const)
// ==========================================
template <typename GeometryTraits, typename TopologyTraits>
typename Arrangement_on_curve_1<GeometryTraits, TopologyTraits>::Const_location_result
locate(const Arrangement_on_curve_1<GeometryTraits, TopologyTraits>& arr,
       const typename GeometryTraits::Point_1& p) {
  if (arr.is_empty()) return nullptr;

  auto comp = arr.geometry_traits_1().compare_x_1_object();
  auto pmap = arr.vertex_point_map();
  auto topo = arr.topology_traits();
  auto v_range = arr.vertices();

  for (auto vit = v_range.begin(); vit != v_range.end(); ++vit) {
    Comparison_result res = comp(p, get(pmap, vit));
    if (res == EQUAL) return vit;
    if (res == SMALLER) return topo.left_edge(vit);
  }

  auto last_vit = std::prev(v_range.end());
  return topo.right_edge(last_vit);
}

// ==========================================
// LOCATE (Mutable)
// ==========================================
template <typename GeometryTraits, typename TopologyTraits>
typename Arrangement_on_curve_1<GeometryTraits, TopologyTraits>::Location_result
locate(Arrangement_on_curve_1<GeometryTraits, TopologyTraits>& arr,
       const typename GeometryTraits::Point_1& p) {
  if (arr.is_empty()) return nullptr;

  auto comp = arr.geometry_traits_1().compare_x_1_object();
  auto pmap = arr.vertex_point_map();
  auto topo = arr.topology_traits();
  auto v_range = arr.vertices();

  for (auto vit = v_range.begin(); vit != v_range.end(); ++vit) {
    Comparison_result res = comp(p, get(pmap, vit));
    if (res == EQUAL) return vit;
    if (res == SMALLER) return topo.left_edge(vit);
  }

  auto last_vit = std::prev(v_range.end());
  return topo.right_edge(last_vit);
}

  // ==========================================
// INSERT
// ==========================================
template <typename GeometryTraits, typename TopologyTraits>
typename Arrangement_on_curve_1<GeometryTraits, TopologyTraits>::Vertex_descriptor
insert(Arrangement_on_curve_1<GeometryTraits, TopologyTraits>& arr, const typename GeometryTraits::Point_1& p) {
  using Arr = Arrangement_on_curve_1<GeometryTraits, TopologyTraits>;
  auto location = locate(arr, p);

  if (auto v_ptr = std::get_if<typename Arr::Vertex_descriptor>(&location)) return *v_ptr;
  if (auto e_ptr = std::get_if<typename Arr::Edge_descriptor>(&location)) return arr.split_edge(*e_ptr, p);

  return arr.insert_empty(p);
}

} // namespace Arrangement_on_curve_1
} // namespace CGAL

#endif
