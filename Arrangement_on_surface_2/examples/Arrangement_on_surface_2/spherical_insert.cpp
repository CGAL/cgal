//! \file examples/Arrangement_on_surface_2/spherical_insert.cpp
// Constructing an arrangement of arcs of great circles.

#include <list>
#include <cmath>
#include <cstdio>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Arrangement_on_surface_2.h>
#include <CGAL/Arr_geodesic_arc_on_sphere_traits_2.h>
#include <CGAL/Arr_spherical_topology_traits_2.h>

#include "arr_geodesic.h"
#include "arr_print.h"

int main() {
  Geom_traits traits;
  auto ctr_p = traits.construct_point_2_object();
  auto ctr_cv = traits.construct_curve_2_object();
  Arrangement arr(&traits);

  Point p1 = ctr_p(0, 0, -1), p3 = ctr_p(0, -1, 0), p5 = ctr_p(-1, 0, 0);
  Point p2 = ctr_p(0, 0,  1), p4 = ctr_p(0,  1, 0), p6 = ctr_p( 1, 0, 0);
  Curve arcs[] = {
    ctr_cv(p6, p1), ctr_cv(p6, p2), ctr_cv(p4, p1), ctr_cv(p4, p2),
    ctr_cv(p5, p1), ctr_cv(p5, p2), ctr_cv(p3, p1), ctr_cv(p3, p2),
    ctr_cv(p6, p4), ctr_cv(p6, p3), ctr_cv(p5, p4), ctr_cv(p5, p3) };
  CGAL::insert(arr, arcs, arcs + sizeof(arcs)/sizeof(Curve));
  print_arrangement_size(arr);
  return 0;
}
