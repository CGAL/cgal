// Copyright (c) 2026 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Efi Fogel         <efifogel@gmail.com>

#include <iostream>
#include <string>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Arrangement_on_curve_1/Arrangement_on_curve_1.h>
#include <CGAL/Arrangement_on_curve_1/insert.h>
#include <CGAL/Arrangement_on_curve_1/overlay.h>
#include <CGAL/Arrangement_on_curve_1/Unbounded_topology_traits.h>
#include <CGAL/Arrangement_on_curve_1/Line_3_traits_1.h>
#include <CGAL/Arrangement_on_curve_1/Copy_overlay_observer.h>

using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;
using Traits = CGAL::Arrangement_on_curve_1::Line_3_traits_1<Kernel>;
using Topo = CGAL::Arrangement_on_curve_1::Unbounded_topology_traits<typename Traits::Point_1, std::string, void>;
using Arrangement = CGAL::Arrangement_on_curve_1::Arrangement_on_curve_1<Traits, Topo>;
using Observer = CGAL::Arrangement_on_curve_1::Copy_overlay_observer<Arrangement, Arrangement, Arrangement>;

int main() {
  // Define an infinite master 3D line passing through (0,0,0) with direction vector (1,1,1)
  Kernel::Point_3 origin(0, 0, 0);
  Kernel::Point_3 direction_pt(1, 1, 1);
  Kernel::Line_3 master_line(origin, direction_pt);

  Traits traits(master_line);

  // Initialize the two source structures and the result container
  Arrangement arr1(traits);
  Arrangement arr2(traits);
  Arrangement res_arr(traits);

  // Define some 3D points resting strictly along the line trajectory
  Kernel::Point_3 pA(1, 1, 1);
  Kernel::Point_3 pB(3, 3, 3);
  Kernel::Point_3 pC(5, 5, 5);

  std::cout << "Populating the first 1D arrangement (arr1)...\n";
  auto v1_A = CGAL::Arrangement_on_curve_1::insert(arr1, pA);
  auto v1_B = CGAL::Arrangement_on_curve_1::insert(arr1, pB);

  auto v1_data = arr1.vertex_data_map();
  put(v1_data, v1_A, std::string("Station Alpha"));
  put(v1_data, v1_B, std::string("Station Beta"));

  std::cout << "Populating the second 1D arrangement (arr2)...\n";
  auto v2_B = CGAL::Arrangement_on_curve_1::insert(arr2, pB);
  auto v2_C = CGAL::Arrangement_on_curve_1::insert(arr2, pC);

  auto v2_data = arr2.vertex_data_map();
  put(v2_data, v2_B, std::string("Station Beta (Updated)"));
  put(v2_data, v2_C, std::string("Station Gamma"));

  // Instantiate the concrete data copying observer
  std::cout << "\nInstantiating the copying observer interface...\n";
  Observer observer(arr1, arr2, res_arr);

  // Compute the topological overlay
  std::cout << "Computing the arrangement overlay along the 3D line track...\n";
  CGAL::Arrangement_on_curve_1::overlay(arr1, arr2, res_arr, observer);

  // Traverse and inspect the merged arrangement sequence
  std::cout << "\nResulting Overlay Arrangement Sequence:\n";
  auto res_pnt_map = res_arr.vertex_point_map();
  auto res_data_map = res_arr.vertex_data_map();

  for (auto v : res_arr.vertices()) {
    std::cout << "  Vertex 3D Point: (" << get(res_pnt_map, v)
              << ") | Copied Attribute: \"" << get(res_data_map, v) << "\"\n";
  }

  std::cout << "\nStatistics:\n";
  std::cout << "  Number of Vertices: " << res_arr.number_of_vertices() << "\n";
  std::cout << "  Number of Edges:    " << res_arr.number_of_edges() << "\n";

  return 0;
}
