//! \file examples/Arrangement_on_surface_2/spherical_insert.cpp
// Using the global aggregated insertion functions.

#include <vector>
#include <list>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Arrangement_on_surface_2.h>
#include <CGAL/Arr_geodesic_arc_on_sphere_traits_2.h>
#include <CGAL/Arr_spherical_topology_traits_2.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel    Kernel;
typedef CGAL::Arr_geodesic_arc_on_sphere_traits_2<Kernel, -11, 7>
                                                             Geom_traits_2;
typedef Geom_traits_2::Point_2                               Point_2;
typedef Geom_traits_2::Curve_2                               Curve_2;
typedef CGAL::Arr_spherical_topology_traits_2<Geom_traits_2> Topol_traits_2;
typedef CGAL::Arrangement_on_surface_2<Geom_traits_2, Topol_traits_2>
                                                             Arrangement_2;

int main()
{
  // Construct the arrangement from 12 geodesic arcs.
  Geom_traits_2 traits;
  Arrangement_2 arr(&traits);

  Geom_traits_2::Construct_point_2 ctr_p = traits.construct_point_2_object();
  Geom_traits_2::Construct_curve_2 ctr_cv = traits.construct_curve_2_object();

  std::vector<Point_2> points;

  // Observe that the identification curve is a meridian that contains the
  // point (-11, 7, 0). The curve (-1,0,0),(0,1,0) intersects the identification
  // curve.

  std::list<Curve_2> arcs;

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

  // Print the size of the arrangement.
  std::cout << "The arrangement size:" << std::endl
            << "   V = " << arr.number_of_vertices()
            << ",  E = " << arr.number_of_edges()
            << ",  F = " << arr.number_of_faces() << std::endl;

  return 0;
}
