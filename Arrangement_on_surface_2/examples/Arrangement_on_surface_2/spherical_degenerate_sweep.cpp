//! \file examples/Arrangement_on_surface_2/spherical_degenerate_sweep.cpp
// Using the global aggregated insertion function.

// #define CGAL_SL_VERBOSE 0
// #define CGAL_ARRANGEMENT_ON_SURFACE_INSERT_VERBOSE 1
// #define CGAL_ARR_CONSTRUCTION_SL_VISITOR_VERBOSE 1

#include <cassert>
#include <vector>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Arrangement_on_surface_2.h>
#include <CGAL/Arr_geodesic_arc_on_sphere_traits_2.h>
#include <CGAL/Arr_spherical_topology_traits_2.h>

#include "arr_print.h"

using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;

#if 0
using Geom_traits_2 = CGAL::Arr_geodesic_arc_on_sphere_traits_2<Kernel, -8, 6>;
#elif 0
using Geom_traits_2 = CGAL::Arr_geodesic_arc_on_sphere_traits_2<Kernel, -11, 7>;
#else
using Geom_traits_2 = CGAL::Arr_geodesic_arc_on_sphere_traits_2<Kernel, -1, 0>;
#endif

using Point_2 = Geom_traits_2::Point_2;
using X_monotone_curve_2 = Geom_traits_2::X_monotone_curve_2;
using Topol_traits_2 = CGAL::Arr_spherical_topology_traits_2<Geom_traits_2>;
using Arrangement_2 =
  CGAL::Arrangement_on_surface_2<Geom_traits_2, Topol_traits_2>;
using Vertex_handle = Arrangement_2::Vertex_handle;

int main() {
  Geom_traits_2 traits;
  auto ctr_p = traits.construct_point_2_object();
  auto ctr_xcv = traits.construct_x_monotone_curve_2_object();

  std::vector<Point_2> points;
  std::vector<X_monotone_curve_2> xcvs;

  CGAL::IO::set_pretty_mode(std::cout);

  Point_2 sp = ctr_p(0, 0, -1);
  Point_2 np = ctr_p(0, 0, 1);
  points.push_back(sp);
  points.push_back(np);

  Point_2 p1 = ctr_p(-1,  0, 0);
  Point_2 p2 = ctr_p( 0, -1, 0);
  Point_2 p3 = ctr_p( 1,  0, 0);
  points.push_back(p1);
  points.push_back(p2);
  points.push_back(p3);

  X_monotone_curve_2 xcv_sp1 = ctr_xcv(sp, p1);
  X_monotone_curve_2 xcv_sp2 = ctr_xcv(sp, p2);
  X_monotone_curve_2 xcv_sp3 = ctr_xcv(sp, p3);
  assert(xcv_sp1.is_vertical());
  assert(xcv_sp2.is_vertical());
  assert(xcv_sp3.is_vertical());
  xcvs.push_back(xcv_sp1);      // 0
  xcvs.push_back(xcv_sp2);      // 1
  xcvs.push_back(xcv_sp3);      // 2

  X_monotone_curve_2 xcv_12 = ctr_xcv(p1, p2);
  X_monotone_curve_2 xcv_23 = ctr_xcv(p2, p3);
  assert(!xcv_12.is_vertical());
  assert(!xcv_23.is_vertical());
  xcvs.push_back(xcv_12);       // 3
  xcvs.push_back(xcv_23);       // 4

  X_monotone_curve_2 xcv_np1 = ctr_xcv(np, p1);
  X_monotone_curve_2 xcv_np2 = ctr_xcv(np, p2);
  X_monotone_curve_2 xcv_np3 = ctr_xcv(np, p3);
  assert(xcv_np1.is_vertical());
  assert(xcv_np2.is_vertical());
  assert(xcv_np3.is_vertical());
  xcvs.push_back(xcv_np1);      // 5
  xcvs.push_back(xcv_np2);      // 6
  xcvs.push_back(xcv_np3);      // 7

  unsigned subsetsp = (1 << points.size());
  std::cout << "#subsets points: " << subsetsp << std::endl;

  unsigned subsets = (1 << xcvs.size());
  std::cout << "#subsets curves: " << subsets << std::endl;

  std::cout << "total combinations: " << (subsetsp)*(subsets)
            << std::endl<< std::endl;

  for (unsigned up = 0; up < subsetsp; up++) {
    std::vector< Point_2 > points_sub;
    for (unsigned ep = 0; ep <= points.size(); ep++) {
      if (up & (1 << ep)) {
        std::cout << "take pt: "  << ep << std::endl;
        points_sub.push_back(points[ep]);
      }
    }

    for (unsigned u = 0; u < subsets; u++) {
      std::vector< X_monotone_curve_2 > xcvs_sub;
      for (unsigned e = 0; e <= xcvs.size(); e++) {
        if (u & (1 << e)) {
          std::cout << "take xcv: "  << e << std::endl;
          xcvs_sub.push_back(xcvs[e]);
        }
      }

      std::cout << "subsetpoints #" << up << " has size: "
                << points_sub.size() << std::endl;
      std::cout << "subsetcurves #" << u << " has size: "
                << xcvs_sub.size() << std::endl;

#if 1
      Arrangement_2 arr;
      std::cout << "inserting "
                << xcvs_sub.size() << " x-monotone curves and "
                << points_sub.size() << " isolated points."
                << std::endl;

      // TODO why is this signature not available as "insert(...)"
      CGAL::insert_empty(arr, xcvs_sub.begin(), xcvs_sub.end(),
                         points_sub.begin(), points_sub.end());

      print_arrangement_size(arr);          // print the arrangement size

      std::cout << "======================================================="
                << std::endl << std::endl << std::endl;
#endif
      //std::cout << "arr: " << arr << std::endl;
      std::cout << std::endl;
    }
  }

  return 0;
}
