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

#ifndef CGAL_ARRANGEMENT_ON_CURVE_1_FUNCTIONS_H
#define CGAL_ARRANGEMENT_ON_CURVE_1_FUNCTIONS_H

#include <variant>

#include <CGAL/Arrangement_on_curve_1.h>

namespace CGAL {
namespace Arrangement_on_curve_1 {

// ==========================================
// LOCATE (Const)
// ==========================================
template <typename Traits, typename VertexData, typename EdgeData>
typename Arrangement_on_curve_1<Traits, VertexData, EdgeData>::Const_location_result
locate(const Arrangement_on_curve_1<Traits, VertexData, EdgeData>& arr, const typename Traits::Point_1& p) {
  auto traits = arr.traits();
  auto comp = traits.compare_x_1_object();

  if (arr.is_empty()) return nullptr;

  for (auto vit = arr.vertices_begin(); vit != arr.vertices_end(); ++vit) {
    Comparison_result res = comp(p, vit->point());
    if (res == EQUAL) return vit;
    if (res == SMALLER) return vit->left();
  }

  auto last_vit = std::prev(arr.vertices_end());
  return last_vit->right();
}

// ==========================================
// LOCATE (Mutable)
// ==========================================
template <typename Traits, typename VertexData, typename EdgeData>
typename Arrangement_on_curve_1<Traits, VertexData, EdgeData>::Location_result
locate(Arrangement_on_curve_1<Traits, VertexData, EdgeData>& arr, const typename Traits::Point_1& p) {
  auto traits = arr.traits();
  auto comp = traits.compare_x_1_object();

  if (arr.is_empty()) return nullptr;

  for (auto vit = arr.vertices_begin(); vit != arr.vertices_end(); ++vit) {
    Comparison_result res = comp(p, vit->point());
    if (res == EQUAL) return vit;
    if (res == SMALLER) return vit->left();
  }

  auto last_vit = std::prev(arr.vertices_end());
  return last_vit->right();
}

// ==========================================
// INSERT
// ==========================================
template <typename Traits, typename VertexData, typename EdgeData>
typename Arrangement_on_curve_1<Traits, VertexData, EdgeData>::Vertex_handle
insert(Arrangement_on_curve_1<Traits, VertexData, EdgeData>& arr, const typename Traits::Point_1& p) {
  using Arr = Arrangement_on_curve_1<Traits, VertexData, EdgeData>;

  auto location = locate(arr, p);

  if (auto v_ptr = std::get_if<typename Arr::Vertex_handle>(&location)) return *v_ptr;
  if (auto e_ptr = std::get_if<typename Arr::Edge_handle>(&location)) return arr.split_edge(*e_ptr, p);

  return arr.insert_empty(p);
}

// ==========================================
// OVERLAY
// ==========================================
template <typename Traits, typename VertexData, typename EdgeData>
void overlay(const Arrangement_on_curve_1<Traits, VertexData, EdgeData>& arr1,
             const Arrangement_on_curve_1<Traits, VertexData, EdgeData>& arr2,
             Arrangement_on_curve_1<Traits, VertexData, EdgeData>& result) {
  auto comp = arr1.traits().compare_x_1_object();
  auto vit1 = arr1.vertices_begin();
  auto vit2 = arr2.vertices_begin();

  while (vit1 != arr1.vertices_end() && vit2 != arr2.vertices_end()) {
    Comparison_result res = comp(vit1->point(), vit2->point());
    if (res == EQUAL) {
      insert(result, vit1->point());
      ++vit1;
      ++vit2;
    }
    else
    if (res == SMALLER) {
      insert(result, vit1->point());
      ++vit1;
    }
    else {
      insert(result, vit2->point());
      ++vit2;
    }
  }

  while (vit1 != arr1.vertices_end()) {
    insert(result, vit1->point());
    ++vit1;
  }
  while (vit2 != arr2.vertices_end()) {
    insert(result, vit2->point());
    ++vit2;
  }
}

} // namespace Arrangement_on_curve_1
} // namespace CGAL

#endif
