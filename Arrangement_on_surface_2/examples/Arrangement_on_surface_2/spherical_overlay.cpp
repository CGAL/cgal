//! \file examples/Arrangement_on_surface_2/spherical_overlay.cpp
// Overlay of two arrangements embedded on the sphere.

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Arrangement_on_surface_2.h>
#include <CGAL/Arr_geodesic_arc_on_sphere_traits_2.h>
#include <CGAL/Arr_spherical_topology_traits_2.h>
#include <CGAL/Arr_overlay_2.h>
#include <CGAL/Arr_default_overlay_traits.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel    Kernel;
typedef CGAL::Arr_geodesic_arc_on_sphere_traits_2<Kernel>    Geom_traits_2;
typedef Geom_traits_2::Point_2                               Point_2;
typedef Geom_traits_2::X_monotone_curve_2                    X_monotone_curve_2;
typedef CGAL::Arr_spherical_topology_traits_2<Geom_traits_2> Topol_traits_2;
typedef CGAL::Arrangement_on_surface_2<Geom_traits_2, Topol_traits_2>
                                                             Arrangement_2;
typedef CGAL::Arr_default_overlay_traits<Arrangement_2>      Overlay_traits;

int main()
{
  Geom_traits_2 traits;
  Geom_traits_2::Construct_point_2 ctr_p = traits.construct_point_2_object();
  Geom_traits_2::Construct_x_monotone_curve_2 ctr_xcv =
    traits.construct_x_monotone_curve_2_object();

  Kernel::Direction_3 dir1(0, 0, 1);
  X_monotone_curve_2 g11 = ctr_xcv(ctr_p(-1, 0, 0), ctr_p(0, -1, 0), dir1);
  X_monotone_curve_2 g12 = ctr_xcv(ctr_p(0, -1, 0), ctr_p(-1, 0, 0), dir1);

  Arrangement_2 arr1(&traits);
  CGAL::insert(arr1, g11);
  CGAL::insert(arr1, g12);
  std::cout << "No. of vertices: " << arr1.number_of_vertices() << std::endl;

  Kernel::Direction_3 dir2(0, 0, -1);
  X_monotone_curve_2 g21 = ctr_xcv(ctr_p(-1, 0, 0), ctr_p(0, 1, 0), dir2);
  X_monotone_curve_2 g22 = ctr_xcv(ctr_p(0, 1, 0), ctr_p(-1, 0, 0), dir2);

  Arrangement_2 arr2(&traits);
  CGAL::insert(arr2, g21);
  CGAL::insert(arr2, g22);
  std::cout << "No. of vertices: " << arr2.number_of_vertices() << std::endl;

  Arrangement_2 overlay_arr;
  Overlay_traits overlay_traits;
  overlay(arr1, arr2, overlay_arr, overlay_traits);
  std::cout << "No. of vertices: " << overlay_arr.number_of_vertices()
            << std::endl;

  return 0;
}
