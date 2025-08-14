#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Arrangement_on_surface_2.h>
#include <CGAL/Arr_geodesic_arc_on_sphere_traits_2.h>
#include <CGAL/Arr_spherical_topology_traits_2.h>
#include <CGAL/draw_arrangement_2.h>

#include "arr_geodesic.h"
#include "arr_print.h"

int main() {
  Geom_traits traits;
  auto ctr_p = traits.construct_point_2_object();
  auto ctr_cv = traits.construct_curve_2_object();
  Arrangement arr(&traits);
  // identification curve in two parts
  CGAL::insert(arr, ctr_cv(ctr_p(0, 0, -1), ctr_p(-1, 0, 0)));
  CGAL::insert(arr, ctr_cv(ctr_p(-1, 0, 0), ctr_p(0, 0, 1)));

  {
    // outer face
    Point p1 = ctr_p(0, 1, 0), p2 = ctr_p(0, 0, 1), p3 = ctr_p(-1, 0, 0);
    auto arcs = {ctr_cv(p1, p2), ctr_cv(p1, p3), ctr_cv(p2, p3)};
    CGAL::insert(arr, arcs.begin(), arcs.end());
  }
  {
    // lake
    Point p1 = ctr_p(0.45, -0.45, -0.75), p2 = ctr_p(0.45, -0.85, -0.30), p3 = ctr_p(0.85, -0.35, -0.35);
    auto arcs = {ctr_cv(p1, p2), ctr_cv(p2, p3), ctr_cv(p3, p1)};
    CGAL::insert(arr, arcs.begin(), arcs.end());
  }

  print_arrangement_size(arr);
  CGAL::draw(arr, "spherical lakes");
  return 0;
}
