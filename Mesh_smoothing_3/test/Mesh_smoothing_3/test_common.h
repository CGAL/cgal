// Copyright (c) 2026  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial

#ifndef CGAL_MESH_SMOOTHING_3_TEST_COMMON_H
#define CGAL_MESH_SMOOTHING_3_TEST_COMMON_H

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_smoothing_3/boundary_aware_mesh_smoothing.h>
#include <CGAL/Mesh_smoothing_3/projectors.h>
#include <CGAL/Mesh_smoothing_3/Smoothing_status.h>
#include <CGAL/Simplicial_mesh_cell_base_3.h>
#include <CGAL/Simplicial_mesh_vertex_base_3.h>
#include <CGAL/Triangulation_3.h>
#include <CGAL/Triangulation_data_structure_3.h>
#include <CGAL/enum.h>

#include <boost/property_map/property_map.hpp>

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <map>
#include <tuple>
#include <utility>
#include <vector>

namespace ms3_test {

using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_3 = K::Point_3;
using Subdomain_index = int;
using Surface_patch_index = unsigned;
using Curve_index = unsigned;
using Corner_index = unsigned;

using Cb = CGAL::Simplicial_mesh_cell_base_3<K, Subdomain_index, Surface_patch_index>;
using Vb = CGAL::Simplicial_mesh_vertex_base_3<K, Subdomain_index, Surface_patch_index,
                                              Curve_index, Corner_index>;
using Tds = CGAL::Triangulation_data_structure_3<Vb, Cb, CGAL::Sequential_tag>;
using Tr = CGAL::Triangulation_3<K, Tds>;
using C3t3 = CGAL::Mesh_complex_3_in_triangulation_3<Tr>;
using Vertex_handle = C3t3::Vertex_handle;
using Cell_handle = C3t3::Cell_handle;
using Facet = C3t3::Facet;
using Edge = C3t3::Edge;

inline bool same_point(Point_3 const& a, Point_3 const& b, bool exact = true)
{
  if (exact) return a == b;
  auto d = a - b;
  return d.squared_length() <= 1e-20;
}

inline bool finite_point(Point_3 const& p)
{
  return std::isfinite(CGAL::to_double(p.x())) &&
         std::isfinite(CGAL::to_double(p.y())) &&
         std::isfinite(CGAL::to_double(p.z()));
}

inline int local_index(Cell_handle c, Vertex_handle v)
{
  for (int i = 0; i < 4; ++i) {
    if (c->vertex(i) == v) return i;
  }
  assert(false);
  return -1;
}

inline Facet facet_opposite(Cell_handle c, Vertex_handle v)
{
  return Facet(c, local_index(c, v));
}

inline Edge find_edge(Tr const& tr, Vertex_handle a, Vertex_handle b)
{
  auto fallback = *tr.finite_edges_begin();
  for (const auto& e : tr.finite_edges()) {
    auto v0 = e.first->vertex(e.second);
    auto v1 = e.first->vertex(e.third);
    if ((v0 == a && v1 == b) || (v0 == b && v1 == a)) return e;
  }
  assert(false);
  return fallback;
}

inline Cell_handle find_cell_with_vertex(Tr const& tr, Vertex_handle v)
{
  auto fallback = *tr.finite_cell_handles().begin();
  for (const auto& c : tr.finite_cell_handles()) {
    for (int i = 0; i < 4; ++i) {
      if (c->vertex(i) == v) return c;
    }
  }
  assert(false);
  return fallback;
}


struct Fixture
{
  C3t3 c3t3;
  std::array<Vertex_handle, 5> v{};
  Cell_handle complex_cell;
  Cell_handle other_cell;
  Facet shared_facet;
  Facet complex_facet_2;
  Facet unmarked_other_facet;
  Edge curve_edge_1;
  Edge curve_edge_2;
  Edge unmarked_edge;
};

inline Fixture make_fixture()
{
  Fixture f;
  auto& tr = f.c3t3.triangulation();

  f.v[0] = tr.insert(Point_3(0.0, 0.0, 0.0));
  f.v[1] = tr.insert(Point_3(1.3, 0.1, 0.0));
  f.v[2] = tr.insert(Point_3(0.2, 1.1, 0.0));
  f.v[3] = tr.insert(Point_3(0.2, 0.1, 1.4));
  f.v[4] = tr.insert(Point_3(0.15, 0.2, -0.9));

  f.complex_cell = find_cell_with_vertex(tr, f.v[3]);
  f.other_cell = find_cell_with_vertex(tr, f.v[4]);

  f.c3t3.add_to_complex(f.complex_cell, 3);

  f.shared_facet = facet_opposite(f.complex_cell, f.v[3]);        // ABC
  f.complex_facet_2 = facet_opposite(f.complex_cell, f.v[0]);     // BCD
  f.unmarked_other_facet = facet_opposite(f.other_cell, f.v[4]);  // ABC on other cell

  f.curve_edge_1 = find_edge(tr, f.v[1], f.v[3]);                 // BD
  f.curve_edge_2 = find_edge(tr, f.v[2], f.v[3]);                 // CD
  f.unmarked_edge = find_edge(tr, f.v[0], f.v[1]);               // AB

  f.c3t3.add_to_complex(f.shared_facet, 11);
  f.c3t3.add_to_complex(f.complex_facet_2, 12);
  f.c3t3.add_to_complex(f.curve_edge_1, 21);
  f.c3t3.add_to_complex(f.curve_edge_2, 22);
  f.c3t3.add_to_complex(f.v[0], 31);
  return f;
}

template <typename C3t3T>
struct Recording_cts
{
  using Geom_traits = typename C3t3T::Triangulation::Geom_traits;
  using Tangent_space = CGAL::Mesh_smoothing_3::Tangent_space<Geom_traits>;
  using Patch_face = std::pair<typename C3t3T::Surface_patch_index, typename C3t3T::Facet>;
  using Curve_edge = std::pair<typename C3t3T::Curve_index, typename C3t3T::Edge>;

  mutable std::vector<Patch_face> patch_queries;
  mutable std::vector<Curve_edge> curve_queries;
  std::map<typename C3t3T::Surface_patch_index, Tangent_space> patch_results;
  std::map<typename C3t3T::Curve_index, Tangent_space> curve_results;

  Tangent_space patch_face_projection_plane(Patch_face patch_face, std::vector<typename Geom_traits::Point_3> const&) const
  {
    patch_queries.push_back(patch_face);
    return patch_results.at(patch_face.first);
  }

  Tangent_space curve_edge_projection_line(Curve_edge curve_edge, std::array<typename Geom_traits::Point_3, 2> const&) const
  {
    curve_queries.push_back(curve_edge);
    return curve_results.at(curve_edge.first);
  }
};

template <typename C3t3T>
inline std::vector<Point_3> finite_vertices(C3t3T const& c3t3)
{
  std::vector<Point_3> pts;
  for (auto v : c3t3.triangulation().finite_vertex_handles()) {
    pts.push_back(v->point());
  }
  return pts;
}

template <typename C3t3T>
inline void assert_all_finite(C3t3T const& c3t3)
{
  for (auto v : c3t3.triangulation().finite_vertex_handles()) {
    assert(finite_point(v->point()));
  }
}

} // namespace ms3_test

#endif
