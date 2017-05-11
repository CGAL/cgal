//! \file examples/Arrangement_on_surface_2/aggregated_insertion.cpp
// Using the global aggregated insertion functions.

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
  Vertex_handle spv = arr.insert_in_face_interior(sp, arr.reference_face());
  Vertex_handle npv = arr.insert_in_face_interior(np, arr.reference_face());

#if 1
  Point_2 p1 = ctr_p(-1,  0, 0);
#else
  Point_2 p1 = ctr_p(-1,  -1, 0);
#endif
  Point_2 p2 = ctr_p( 0, -1, 0);
  Point_2 p3 = ctr_p( 1,  0, 0);

  Vertex_handle v1 = arr.insert_in_face_interior(p1, arr.reference_face());
  Vertex_handle v2 = arr.insert_in_face_interior(p2, arr.reference_face());
  Vertex_handle v3 = arr.insert_in_face_interior(p3, arr.reference_face());

  X_monotone_curve_2 xcv_sp1 = ctr_xcv(sp, p1);
  X_monotone_curve_2 xcv_sp2 = ctr_xcv(sp, p2);
  X_monotone_curve_2 xcv_sp3 = ctr_xcv(sp, p3);
  arr.insert_at_vertices(xcv_sp1, spv, v1);
  arr.insert_at_vertices(xcv_sp2, spv, v2);
  arr.insert_at_vertices(xcv_sp3, spv, v3);

  X_monotone_curve_2 xcv_np1 = ctr_xcv(np, p1);
  X_monotone_curve_2 xcv_np2 = ctr_xcv(np, p2);
  X_monotone_curve_2 xcv_np3 = ctr_xcv(np, p3);
  arr.insert_at_vertices(xcv_np1, npv, v1);
  arr.insert_at_vertices(xcv_np2, npv, v2);
  arr.insert_at_vertices(xcv_np3, npv, v3);

  X_monotone_curve_2 xcv_12 = ctr_xcv(p1, p2);
  X_monotone_curve_2 xcv_23 = ctr_xcv(p2, p3);
  arr.insert_at_vertices(xcv_12, v1, v2);
  arr.insert_at_vertices(xcv_23, v2, v3);

  // Print the size of the arrangement.
  std::cout << "The arrangement size:" << std::endl
            << "   V = " << arr.number_of_vertices()
            << ",  E = " << arr.number_of_edges()
            << ",  F = " << arr.number_of_faces() << std::endl;

  // std::cout << "arr: " << arr << std::endl;

  return 0;
}
