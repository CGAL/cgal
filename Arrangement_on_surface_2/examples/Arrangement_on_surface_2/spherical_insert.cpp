//! \file examples/Arrangement_on_surface_2/spherical_insert.cpp
// Constructing an arrangement of arcs of great circles.

#include <list>
#include <cmath>
#include <cstdio>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Arrangement_on_surface_2.h>
#include <CGAL/Arr_geodesic_arc_on_sphere_traits_2.h>
#include <CGAL/Arr_spherical_topology_traits_2.h>
#include <CGAL/draw_arrangement_2.h>

#include "arr_print.h"

using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;
using Geom_traits = CGAL::Arr_geodesic_arc_on_sphere_traits_2<Kernel>;
using Point = Geom_traits::Point_2;
using X_monotone_curve = Geom_traits::X_monotone_curve_2;
using Topol_traits = CGAL::Arr_spherical_topology_traits_2<Geom_traits>;
using Arrangement = CGAL::Arrangement_on_surface_2<Geom_traits, Topol_traits>;

int main(int argc, char* argv[]) {
  // Construct the arrangement from 12 geodesic arcs.
  Geom_traits traits;
  Arrangement arr(&traits);

  auto ctr_p = traits.construct_point_2_object();
  auto ctr_xcv = traits.construct_x_monotone_curve_2_object();

  // Observe that the identification curve is a meridian that contains the
  // point (-11, 7, 0). The curve (-1,0,0),(0,1,0) intersects the identification
  // curve.

  std::list<X_monotone_curve> arcs;

  arcs.push_back(ctr_xcv(ctr_p(1, 0, 0), ctr_p(0, 0, -1)));
  arcs.push_back(ctr_xcv(ctr_p(1, 0, 0), ctr_p(0, 0, 1)));
  arcs.push_back(ctr_xcv(ctr_p(0, 1, 0), ctr_p(0, 0, -1)));
  arcs.push_back(ctr_xcv(ctr_p(0, 1, 0), ctr_p(0, 0, 1)));
  arcs.push_back(ctr_xcv(ctr_p(-1, 0, 0), ctr_p(0, 0, -1)));
  arcs.push_back(ctr_xcv(ctr_p(-1, 0, 0), ctr_p(0, 0, 1)));
  arcs.push_back(ctr_xcv(ctr_p(0, -1, 0), ctr_p(0, 0, -1)));
  arcs.push_back(ctr_xcv(ctr_p(0, -1, 0), ctr_p(0, 0, 1)));
  arcs.push_back(ctr_xcv(ctr_p(1, 0, 0), ctr_p(0, 1, 0)));
  arcs.push_back(ctr_xcv(ctr_p(1, 0, 0), ctr_p(0, -1, 0)));
  arcs.push_back(ctr_xcv(ctr_p(-1, 0, 0), ctr_p(0, 1, 0)));
  arcs.push_back(ctr_xcv(ctr_p(-1, 0, 0), ctr_p(0, -1, 0)));

  // Construct meridians
  std::size_t num_meridians(16);
  if (argc > 1) sscanf(argv[1], "%zu", &num_meridians);

  if (num_meridians > 0) {
    auto np = ctr_p(0, 0, 1);
    auto sp = ctr_p(0, 0, -1);
    using FT = Kernel::FT;
    using Direction_3 = Kernel::Direction_3;
    auto delta = 2.0 * CGAL_PI / num_meridians;
    FT z(0);
    for (std::size_t i = 0; i < num_meridians; ++i) {
      auto alpha = delta * i;
      FT x(std::sin(alpha));
      FT y(std::cos(alpha));
      Direction_3 normal(x, y, z);
      arcs.push_back(ctr_xcv(sp, np, normal));
    }
  }

  CGAL::insert(arr, arcs.begin(), arcs.end());
  print_arrangement_size(arr);
  // print_arrangement(arr);
  CGAL::draw(arr, "Aos", true);

  return 0;
}
