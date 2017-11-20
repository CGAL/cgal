//! \file examples/Arrangement_on_surface_2/aggregated_insertion.cpp
// Using the global aggregated insertion functions.

// #define CGAL_IDENTIFICATION_XY CGAL_X_MINUS_11_Y_7
// #define CGAL_ARRANGEMENT_ON_SURFACE_INSERT_VERBOSE 1

#include <list>

#include <CGAL/basic.h>
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
  Geom_traits_2 traits;
  Geom_traits_2::Construct_point_2 ctr_p = traits.construct_point_2_object();
  Geom_traits_2::Construct_x_monotone_curve_2 ctr_xcv =
    traits.construct_x_monotone_curve_2_object();
  Arrangement_2 arr(&traits);
  Point_2 sp = ctr_p(0, 0, -1);
  Point_2 np = ctr_p(0, 0, 1);

  Point_2 p1 = ctr_p(-1,  0, -1);
  Point_2 p2 = ctr_p(-1,  0,  1);
  X_monotone_curve_2 xcv_sp_p2 = ctr_xcv(sp, p2);
  X_monotone_curve_2 xcv_np_p1 = ctr_xcv(np, p1);
  // std::cout << "Inserting " << xcv_sp_p2 << std::endl;
  insert(arr, xcv_sp_p2);
  // std::cout << "Inserting " << xcv_np_p1 << std::endl;
  insert(arr, xcv_np_p1);

  Point_2 q1 = ctr_p(-1,  -1, -1);
  Point_2 q2 = ctr_p(-1,  -1,  1);
  X_monotone_curve_2 xcv_sp_q2 = ctr_xcv(sp, q2);
  X_monotone_curve_2 xcv_np_q1 = ctr_xcv(np, q1);
  insert(arr, xcv_sp_q2);
  insert(arr, xcv_np_q1);

  X_monotone_curve_2 xcv_p1_q1 = ctr_xcv(p1, q1);
  X_monotone_curve_2 xcv_p2_q2 = ctr_xcv(p2, q2);
  insert(arr, xcv_p1_q1);
  insert(arr, xcv_p2_q2);

  // Print the size of the arrangement.
  std::cout << "The arrangement size:" << std::endl
            << "   V = " << arr.number_of_vertices()
            << ",  E = " << arr.number_of_edges()
            << ",  F = " << arr.number_of_faces() << std::endl;

  // std::cout << "arr: " << arr << std::endl;

  return 0;
}
