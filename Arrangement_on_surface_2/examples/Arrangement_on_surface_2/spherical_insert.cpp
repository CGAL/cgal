//! \file examples/Arrangement_on_surface_2/spherical_insert.cpp
// Constructing an arrangement of arcs of great circles.

#include <list>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Arrangement_on_surface_2.h>
#include <CGAL/Arr_geodesic_arc_on_sphere_traits_2.h>
#include <CGAL/Arr_spherical_topology_traits_2.h>

#include "arr_print.h"

typedef CGAL::Exact_predicates_exact_constructions_kernel         Kernel;
typedef CGAL::Arr_geodesic_arc_on_sphere_traits_2<Kernel>         Geom_traits;
typedef Geom_traits::Point_2                                      Point;
typedef Geom_traits::Curve_2                                      Curve;
typedef CGAL::Arr_spherical_topology_traits_2<Geom_traits>        Topol_traits;
typedef CGAL::Arrangement_on_surface_2<Geom_traits, Topol_traits> Arrangement;

int main() {
  // Construct the arrangement from 12 geodesic arcs.
  Geom_traits traits;
  Arrangement arr(&traits);

  auto ctr_p = traits.construct_point_2_object();
  auto ctr_cv = traits.construct_curve_2_object();

  // Observe that the identification curve is a meridian that contains the
  // point (-11, 7, 0). The curve (-1,0,0),(0,1,0) intersects the identification
  // curve.

  std::list<Curve> arcs;

  arcs.push_back(ctr_cv(ctr_p(1, 0, 0), ctr_p(0, 0, -1)));
  arcs.push_back(ctr_cv(ctr_p(1, 0, 0), ctr_p(0, 0, 1)));
  arcs.push_back(ctr_cv(ctr_p(0, 1, 0), ctr_p(0, 0, -1)));
  arcs.push_back(ctr_cv(ctr_p(0, 1, 0), ctr_p(0, 0, 1)));
  arcs.push_back(ctr_cv(ctr_p(-1, 0, 0), ctr_p(0, 0, -1)));
  arcs.push_back(ctr_cv(ctr_p(-1, 0, 0), ctr_p(0, 0, 1)));
  arcs.push_back(ctr_cv(ctr_p(0, -1, 0), ctr_p(0, 0, -1)));
  arcs.push_back(ctr_cv(ctr_p(0, -1, 0), ctr_p(0, 0, 1)));
  arcs.push_back(ctr_cv(ctr_p(1, 0, 0), ctr_p(0, 1, 0)));
  arcs.push_back(ctr_cv(ctr_p(1, 0, 0), ctr_p(0, -1, 0)));
  arcs.push_back(ctr_cv(ctr_p(-1, 0, 0), ctr_p(0, 1, 0)));
  arcs.push_back(ctr_cv(ctr_p(-1, 0, 0), ctr_p(0, -1, 0)));

  CGAL::insert(arr, arcs.begin(), arcs.end());
  print_arrangement_size(arr);          // print the arrangement size
  // print_arrangement(arr);

  return 0;
}
