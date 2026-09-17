// Copyright (c) 2026 Shubham Padkonde.
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Alpha_shape_3.h>
#include <CGAL/Alpha_shape_cell_base_3.h>
#include <CGAL/Alpha_shape_vertex_base_3.h>
#include <CGAL/Regular_triangulation_3.h>
#include <cassert>
#include <vector>

template <typename K>
void test_clear()
{
  typedef CGAL::Alpha_shape_vertex_base_3<K,
    CGAL::Regular_triangulation_vertex_base_3<K> > Vb;
  typedef CGAL::Alpha_shape_cell_base_3<K,
    CGAL::Regular_triangulation_cell_base_3<K> > Cb;
  typedef CGAL::Triangulation_data_structure_3<Vb, Cb> Tds;
  typedef CGAL::Alpha_shape_3<CGAL::Regular_triangulation_3<K, Tds> > Shape;
  typedef typename K::Point_3 Point;
  typedef typename K::Weighted_point_3 Weighted_point;

  const std::vector<Weighted_point> points = {
    Weighted_point(Point(0, 0, 0), 1),
    Weighted_point(Point(1.5, 0, 0), 1),
    Weighted_point(Point(0, 0, 10), 1),
    Weighted_point(Point(0, 1.5, 10), 1)
  };

  Shape shape;
  shape.set_mode(Shape::GENERAL);
  shape.make_alpha_shape(points.begin(), points.end());
  const auto edge_count = shape.get_edge_alpha_map()->size();
  assert(edge_count != 0);
  const auto alpha_count = shape.number_of_alphas();

  shape.clear();
  assert(shape.get_edge_alpha_map()->empty());
  assert(shape.number_of_vertices() == 0);
  assert(shape.number_of_alphas() == 0);

  // Rebuilding also calls clear(), so exercise consecutive builds as well.
  for(int build = 0; build < 3; ++build)
  {
    shape.make_alpha_shape(points.begin(), points.end());
    assert(shape.get_edge_alpha_map()->size() == edge_count);
    assert(shape.number_of_alphas() == alpha_count);
    assert(shape.number_of_vertices() == points.size());
  }
}

int main()
{
  test_clear<CGAL::Exact_predicates_exact_constructions_kernel>();
  test_clear<CGAL::Exact_predicates_inexact_constructions_kernel>();
}
