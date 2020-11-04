//! \file examples/Arrangement_on_surface_2/spherical_overlay.cpp
// Overlay of two arrangements embedded on the sphere.

#include <CGAL/basic.h>
#include <CGAL/Arr_overlay_2.h>
#include <CGAL/Arr_default_overlay_traits.h>

#include "arr_geodesic_on_sphere.h"

typedef CGAL::Arr_default_overlay_traits<Arrangement>         Overlay_traits;

int main() {
  Geom_traits traits;
  auto ctr_p = traits.construct_point_2_object();
  auto ctr_xcv = traits.construct_x_monotone_curve_2_object();

  Kernel::Direction_3 dir1(0, 0, 1);
  X_monotone_curve g11 = ctr_xcv(ctr_p(-1, 0, 0), ctr_p(0, -1, 0), dir1);
  X_monotone_curve g12 = ctr_xcv(ctr_p(0, -1, 0), ctr_p(-1, 0, 0), dir1);

  Arrangement arr1(&traits);
  CGAL::insert(arr1, g11);
  CGAL::insert(arr1, g12);
  std::cout << "No. of vertices: " << arr1.number_of_vertices() << std::endl;

  Kernel::Direction_3 dir2(0, 0, -1);
  X_monotone_curve g21 = ctr_xcv(ctr_p(-1, 0, 0), ctr_p(0, 1, 0), dir2);
  X_monotone_curve g22 = ctr_xcv(ctr_p(0, 1, 0), ctr_p(-1, 0, 0), dir2);

  Arrangement arr2(&traits);
  CGAL::insert(arr2, g21);
  CGAL::insert(arr2, g22);
  std::cout << "No. of vertices: " << arr2.number_of_vertices() << std::endl;

  Arrangement overlay_arr;
  Overlay_traits overlay_traits;
  overlay(arr1, arr2, overlay_arr, overlay_traits);
  std::cout << "No. of vertices: " << overlay_arr.number_of_vertices()
            << std::endl;

  return 0;
}
