//! \file examples/Arrangement_on_surface_2/spherical_degenerate_sweep.cpp
// Using the global aggregated insertion function.

//#define CGAL_SL_VERBOSE 1

//#define CGAL_ARRANGEMENT_ON_SURFACE_INSERT_VERBOSE 1

#include <vector>

#include <CGAL/config.h>
#include <CGAL/Exact_rational.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Arrangement_on_surface_2.h>
#include <CGAL/Arr_geodesic_arc_on_sphere_traits_2.h>
#include <CGAL/Arr_spherical_topology_traits_2.h>

typedef CGAL::Exact_rational                               Number_type;
typedef CGAL::Cartesian<Number_type>                       Kernel;
typedef CGAL::Arr_geodesic_arc_on_sphere_traits_2<Kernel>  Geom_traits_2;
typedef Geom_traits_2::Point_2                             Point_2;
typedef Geom_traits_2::X_monotone_curve_2                  X_monotone_curve_2;
typedef CGAL::Arr_spherical_topology_traits_2<Geom_traits_2> Topol_traits_2;
typedef CGAL::Arrangement_on_surface_2<Geom_traits_2, Topol_traits_2>
                                                           Arrangement_2;
typedef Arrangement_2::Vertex_handle                       Vertex_handle;

int main()
{
  std::vector< Point_2 > points;
  std::vector< X_monotone_curve_2 > xcvs;

  CGAL::set_pretty_mode(std::cout);

  Arrangement_2 arr;
  Point_2 sp(0, 0, -1);
  Point_2 np(0, 0, 1);
  points.push_back(sp);
  points.push_back(np);

  Point_2 p1(-1,  0, 0);
  Point_2 p2( 0, -1, 0);
  Point_2 p3( 1,  0, 0);
  points.push_back(p1);
  points.push_back(p2);
  points.push_back(p3);

  X_monotone_curve_2 xcv_sp1(sp, p1);
  X_monotone_curve_2 xcv_sp2(sp, p2);
  X_monotone_curve_2 xcv_sp3(sp, p3);
  CGAL_assertion(xcv_sp1.is_vertical());
  CGAL_assertion(xcv_sp2.is_vertical());
  CGAL_assertion(xcv_sp3.is_vertical());
  xcvs.push_back(xcv_sp1);
  xcvs.push_back(xcv_sp2);
  xcvs.push_back(xcv_sp3);

  X_monotone_curve_2 xcv_12(p1, p2);
  X_monotone_curve_2 xcv_23(p2, p3);
  CGAL_assertion(!xcv_12.is_vertical());
  CGAL_assertion(!xcv_23.is_vertical());
  xcvs.push_back(xcv_12);
  xcvs.push_back(xcv_23);

  X_monotone_curve_2 xcv_np1(np, p1);
  X_monotone_curve_2 xcv_np2(np, p2);
  X_monotone_curve_2 xcv_np3(np, p3);
  CGAL_assertion(xcv_np1.is_vertical());
  CGAL_assertion(xcv_np2.is_vertical());
  CGAL_assertion(xcv_np3.is_vertical());
  xcvs.push_back(xcv_np1);
  xcvs.push_back(xcv_np2);
  xcvs.push_back(xcv_np3);

  std::cout << "inserting "
            << std::distance(xcvs.begin(), xcvs.end()) << " x-monotone curves and "
            << std::distance(points.begin(), points.end()) << " isolated points."
            << std::endl;

  // TODO why is this signature not available as "insert(...)"
  CGAL::insert_empty(arr, xcvs.begin(), xcvs.end(), points.begin(), points.end());

  // Print the size of the arrangement.
  std::cout << "The arrangement size:" << std::endl
            << "   V = " << arr.number_of_vertices()
            << ",  E = " << arr.number_of_edges()
            << ",  F = " << arr.number_of_faces() << std::endl;

  std::cout << "arr: " << arr << std::endl;

  return 0;
}
