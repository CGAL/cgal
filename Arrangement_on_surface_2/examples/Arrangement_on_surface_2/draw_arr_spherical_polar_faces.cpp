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

  CGAL::insert(arr, ctr_cv(ctr_p(0, 0, 1), ctr_p(-1, 0, 0)));
  Point p1 = ctr_p(1, 0, CGAL_PI / 4.0), p2 = ctr_p(0, 1, CGAL_PI / 4.0), p3 = ctr_p(-1, 0, CGAL_PI / 4.0),
        p4 = ctr_p(0, -1, CGAL_PI / 4.0);
  Curve arcs[] = {ctr_cv(p1, p2), ctr_cv(p2, p3), ctr_cv(p3, p4), ctr_cv(p4, p1)};
  CGAL::insert(arr, arcs, arcs + sizeof(arcs) / sizeof(Curve));

  print_arrangement_size(arr);
  CGAL::draw(arr, "spherical faces containing one of the pole");
  return 0;
}
