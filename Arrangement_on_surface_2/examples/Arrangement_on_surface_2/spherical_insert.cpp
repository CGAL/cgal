//! \file examples/Arrangement_on_surface_2/aggregated_insertion.cpp
// Using the global aggregated insertion functions.

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Arrangement_on_surface_2.h>
#include <CGAL/Arr_geodesic_arc_on_sphere_traits_2.h>
#include <CGAL/Arr_spherical_topology_traits_2.h>
#include <list>

typedef CGAL::Exact_predicates_exact_constructions_kernel    Kernel;
// typedef CGAL::Arr_geodesic_arc_on_sphere_traits_2<Kernel, -1, 0>
typedef CGAL::Arr_geodesic_arc_on_sphere_traits_2<Kernel, -8, 6>
                                                             Geom_traits_2;
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

  std::list<X_monotone_curve_2> arcs;

  Geom_traits_2::Construct_point_2 ctr_p = traits.construct_point_2_object();
  Geom_traits_2::Construct_x_monotone_curve_2 ctr_xcv =
    traits.construct_x_monotone_curve_2_object();

  arcs.push_back(ctr_xcv(ctr_p(1, 0, 0), ctr_p(0, 1, 0)));
  arcs.push_back(ctr_xcv(ctr_p(0, 1, 0), ctr_p(0, 0, 1)));
  arcs.push_back(ctr_xcv(ctr_p(0, 0, 1), ctr_p(1, 0, 0)));
  arcs.push_back(ctr_xcv(ctr_p(0, -1, 0), ctr_p(0, 0, -1)));
  arcs.push_back(ctr_xcv(ctr_p(0, -1, 0), ctr_p(0, 0, 1)));
  arcs.push_back(ctr_xcv(ctr_p(0, -1, 0), ctr_p(1, 0, 0)));
  arcs.push_back(ctr_xcv(ctr_p(0, 1, 0), ctr_p(-1, 0, 0)));
  arcs.push_back(ctr_xcv(ctr_p(-1, 0, 0), ctr_p(0, -1, 0)));
  arcs.push_back(ctr_xcv(ctr_p(0, 0, -1), ctr_p(1, 0, 0)));
  arcs.push_back(ctr_xcv(ctr_p(0, 1, 0), ctr_p(0, 0, -1)));

#if 0
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
