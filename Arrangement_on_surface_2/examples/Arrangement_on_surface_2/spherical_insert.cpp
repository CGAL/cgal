// #define CGAL_SS_VERBOSE 1

//! \file examples/Arrangement_on_surface_2/aggregated_insertion.cpp
// Using the global aggregated insertion functions.

#include <vector>
#include <list>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Arrangement_on_surface_2.h>
#include <CGAL/Arr_geodesic_arc_on_sphere_traits_2.h>
#include <CGAL/Arr_spherical_topology_traits_2.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel    Kernel;
typedef CGAL::Arr_geodesic_arc_on_sphere_traits_2<Kernel, -11, 7>    Geom_traits_2;
typedef Geom_traits_2::Point_2                               Point_2;
typedef Geom_traits_2::X_monotone_curve_2                    X_monotone_curve_2;
typedef CGAL::Arr_spherical_topology_traits_2<Geom_traits_2> Topol_traits_2;
typedef CGAL::Arrangement_on_surface_2<Geom_traits_2, Topol_traits_2>
                                                             Arrangement_2;

int main()
{
  // Construct the arrangement of five intersecting segments.
  Geom_traits_2 traits;
  Arrangement_2 arr(&traits);

  Geom_traits_2::Construct_point_2 ctr_p = traits.construct_point_2_object();
  Geom_traits_2::Construct_x_monotone_curve_2 ctr_xcv =
    traits.construct_x_monotone_curve_2_object();

  std::vector<Point_2> points;

  points.push_back(ctr_p(0, -2, 0));   // 0
  points.push_back(ctr_p(2, 0, -2));   // 1
  points.push_back(ctr_p(2, 0, 2));    // 2
  points.push_back(ctr_p(0, 2, 0));    // 3
  points.push_back(ctr_p(0, 0, 1));    // 4
  points.push_back(ctr_p(-2, 0, 2));   // 5
  points.push_back(ctr_p(-11, 7, 4));  // 6
  points.push_back(ctr_p(-2, 2, 0));   // 7
  points.push_back(ctr_p(-2, -2, 0));  // 8
  points.push_back(ctr_p(0, 0, -1));   // 9
  points.push_back(ctr_p(-2, 0, -2));  // 10
  points.push_back(ctr_p(-11, 7, -4)); // 11
  ////
  points.push_back(ctr_p(0, 0, -1));   // 12
  points.push_back(ctr_p(-1, -1, 0));  // 13
  points.push_back(ctr_p(-1, 1, 0));   // 14
  points.push_back(ctr_p(-11, 7, 0));  // 15
  points.push_back(ctr_p(2, 0, 2));    // 16

  // Observe that by default the identification curve is a meridian that
  // contains the point (-1, 0, 0). The following curves do not intersect
  // the identification curve; thus, they are all x-monotone.

  std::list<X_monotone_curve_2> arcs;

  // arcs.push_back(ctr_xcv(ctr_p(1, 0, 0), ctr_p(0, 1, 0)));
  // arcs.push_back(ctr_xcv(ctr_p(0, 1, 0), ctr_p(0, 0, 1)));
  // arcs.push_back(ctr_xcv(ctr_p(0, 0, 1), ctr_p(1, 0, 0)));
  // arcs.push_back(ctr_xcv(ctr_p(0, -1, 0), ctr_p(0, 0, -1)));
  // arcs.push_back(ctr_xcv(ctr_p(0, -1, 0), ctr_p(0, 0, 1)));
  // arcs.push_back(ctr_xcv(ctr_p(0, -1, 0), ctr_p(1, 0, 0)));
  // arcs.push_back(ctr_xcv(ctr_p(0, 1, 0), ctr_p(-1, 0, 0)));
  // arcs.push_back(ctr_xcv(ctr_p(-1, 0, 0), ctr_p(0, -1, 0)));
  // arcs.push_back(ctr_xcv(ctr_p(0, 0, -1), ctr_p(1, 0, 0)));
  // arcs.push_back(ctr_xcv(ctr_p(0, 1, 0), ctr_p(0, 0, -1)));
  // arcs.push_back(ctr_xcv(ctr_p(-1, 0, 0), ctr_p(0, 0, 1)));
  // arcs.push_back(ctr_xcv(ctr_p(-1, 0, 0), ctr_p(0, 0, -1)));

  // arcs.push_back(ctr_xcv(points[1], points[0]));
  // arcs.push_back(ctr_xcv(points[2], points[0]));
  arcs.push_back(ctr_xcv(points[2], points[1]));
  // arcs.push_back(ctr_xcv(points[3], points[2]));
  // arcs.push_back(ctr_xcv(points[3], points[1]));
  // arcs.push_back(ctr_xcv(points[4], points[2]));
  // arcs.push_back(ctr_xcv(points[5], points[4]));
  // arcs.push_back(ctr_xcv(points[6], points[5]));
  // arcs.push_back(ctr_xcv(points[7], points[6]));
  // arcs.push_back(ctr_xcv(points[3], points[7]));
  // arcs.push_back(ctr_xcv(points[8], points[0]));
  // arcs.push_back(ctr_xcv(points[5], points[8]));
  arcs.push_back(ctr_xcv(points[9], points[1]));
  // arcs.push_back(ctr_xcv(points[10], points[9]));
  // arcs.push_back(ctr_xcv(points[8], points[10]));
  // arcs.push_back(ctr_xcv(points[11], points[7]));
  // arcs.push_back(ctr_xcv(points[10], points[11]));
  // ////
  // arcs.push_back(ctr_xcv(points[13], points[12]));
  // arcs.push_back(ctr_xcv(points[14], points[12]));
  // arcs.push_back(ctr_xcv(points[15], points[14]));
  // arcs.push_back(ctr_xcv(points[15], points[13]));
  // arcs.push_back(ctr_xcv(points[16], points[14]));
  // arcs.push_back(ctr_xcv(points[13], points[16]));
  // arcs.push_back(ctr_xcv(points[16], points[12]));

#if 1
  // insert_non_intersecting_curves (arr, arcs.begin(), arcs.end());
  // insert_x_monotone_curves (arr, arcs.begin(), arcs.end());
  CGAL::insert(arr, arcs.begin(), arcs.end());
#else
  std::list<X_monotone_curve_2>::iterator it;
  for (it = arcs.begin(); it != arcs.end(); ++it) {
    // insert_non_intersecting_curve(arr, *it);
    CGAL::insert(arr, *it);
  }
#endif

  // Print the size of the arrangement.
  std::cout << "The arrangement size:" << std::endl
            << "   V = " << arr.number_of_vertices()
            << ",  E = " << arr.number_of_edges()
            << ",  F = " << arr.number_of_faces() << std::endl;

  return 0;
}
