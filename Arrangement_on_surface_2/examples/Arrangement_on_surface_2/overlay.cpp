//! \file examples/Arrangement_on_surface_2/overlay.cpp
// A simple overlay of two arrangements.

#include <CGAL/Cartesian.h>
#include <CGAL/Exact_rational.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_overlay_2.h>
#include <CGAL/Arr_default_overlay_traits.h>

typedef CGAL::Cartesian<CGAL::Exact_rational>            Kernel;
typedef CGAL::Arr_segment_traits_2<Kernel>               Traits_2;
typedef Traits_2::Point_2                                Point_2;
typedef Traits_2::X_monotone_curve_2                     Segment_2;
typedef CGAL::Arrangement_2<Traits_2>                    Arrangement_2;
typedef CGAL::Arr_default_overlay_traits<Arrangement_2>  Overlay_traits;

int main ()
{
  // Construct the first arrangement, containing a square-shaped face.
  Arrangement_2          arr1;

  Segment_2      s1 (Point_2(2, 2), Point_2(6, 2));
  Segment_2      s2 (Point_2(6, 2), Point_2(6, 6));
  Segment_2      s3 (Point_2(6, 6), Point_2(2, 6));
  Segment_2      s4 (Point_2(2, 6), Point_2(2, 2));

  insert_non_intersecting_curve (arr1, s1);
  insert_non_intersecting_curve (arr1, s2);
  insert_non_intersecting_curve (arr1, s3);
  insert_non_intersecting_curve (arr1, s4);

  // Construct the second arrangement, containing a rhombus-shaped face.
  Arrangement_2          arr2;

  Segment_2      t1 (Point_2(4, 1), Point_2(7, 4));
  Segment_2      t2 (Point_2(7, 4), Point_2(4, 7));
  Segment_2      t3 (Point_2(4, 7), Point_2(1, 4));
  Segment_2      t4 (Point_2(1, 4), Point_2(4, 1));

  insert_non_intersecting_curve (arr2, t1);
  insert_non_intersecting_curve (arr2, t2);
  insert_non_intersecting_curve (arr2, t3);
  insert_non_intersecting_curve (arr2, t4);

  // Compute the overlay of the two arrangements.
  Arrangement_2          overlay_arr;
  Overlay_traits         overlay_traits;

  overlay (arr1, arr2, overlay_arr, overlay_traits);

  // Print the size of the overlaid arrangement.
  std::cout << "The overlaid arrangement size:" << std::endl
            << "   V = " << overlay_arr.number_of_vertices()
            << ",  E = " << overlay_arr.number_of_edges()
            << ",  F = " << overlay_arr.number_of_faces() << std::endl;

  return 0;
}
