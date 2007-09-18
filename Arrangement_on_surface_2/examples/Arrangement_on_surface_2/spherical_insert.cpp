//! \file examples/Arrangement_on_surface_2/aggregated_insertion.cpp
// Using the global aggregated insertion functions.

#include <CGAL/Gmpq.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Arrangement_on_surface_2.h>
#include <CGAL/Arr_great_circular_arc_on_sphere_traits_2.h>
#include <CGAL/Arr_tracing_traits_2.h>
#include <CGAL/Arr_spherical_topology_traits_2.h>
#include <list>

typedef CGAL::Gmpq                                            Number_type;
typedef CGAL::Cartesian<Number_type>                          Kernel;
typedef CGAL::Arr_great_circular_arc_on_sphere_traits_2<Kernel>              Geom_traits_base;
typedef CGAL::Arr_tracing_traits_2<Geom_traits_base>          Geom_traits_2;
typedef Geom_traits_2::Point_2                                Point_2;
typedef Geom_traits_2::X_monotone_curve_2                     Spherical_arc_3;
typedef CGAL::Arr_spherical_topology_traits_2<Geom_traits_2>  Topol_traits_2;
typedef CGAL::Arrangement_on_surface_2<Geom_traits_2, Topol_traits_2>
                                                              Arrangement_2;

int main ()
{
  // Construct the arrangement of five intersecting segments.
  Arrangement_2 arr;
  std::list<Spherical_arc_3> arcs;

  arcs.push_back(Spherical_arc_3(Point_2(1, 0, 0), Point_2(0, 1, 0)));
  arcs.push_back(Spherical_arc_3(Point_2(0, 1, 0), Point_2(0, 0, 1)));
  arcs.push_back(Spherical_arc_3(Point_2(0, 0, 1), Point_2(1, 0, 0)));
  arcs.push_back(Spherical_arc_3(Point_2(0, -1, 0), Point_2(0, 0, -1)));
  arcs.push_back(Spherical_arc_3(Point_2(0, -1, 0), Point_2(0, 0, 1)));
  arcs.push_back(Spherical_arc_3(Point_2(0, -1, 0), Point_2(1, 0, 0)));
  arcs.push_back(Spherical_arc_3(Point_2(0, 1, 0), Point_2(-1, 0, 0)));
  arcs.push_back(Spherical_arc_3(Point_2(-1, 0, 0), Point_2(0, -1, 0)));
  arcs.push_back(Spherical_arc_3(Point_2(0, 0, -1), Point_2(1, 0, 0)));
  arcs.push_back(Spherical_arc_3(Point_2(0, 1, 0), Point_2(0, 0, -1)));
  
#if 0
  // insert_non_intersecting_curves (arr, arcs.begin(), arcs.end());
  // insert_x_monotone_curves (arr, arcs.begin(), arcs.end());
  insert_curves (arr, arcs.begin(), arcs.end());
#else
  std::list<Spherical_arc_3>::iterator it;
  for (it = arcs.begin(); it != arcs.end(); ++it) {
    // insert_non_intersecting_curve(arr, *it);
    // insert_curve(arr, *it);
    insert_x_monotone_curve(arr, *it);
  }
#endif
  
  // Print the size of the arrangement.
  std::cout << "The arrangement size:" << std::endl
            << "   V = " << arr.number_of_vertices()
            << ",  E = " << arr.number_of_edges() 
            << ",  F = " << arr.number_of_faces() << std::endl;

  return 0;
}
