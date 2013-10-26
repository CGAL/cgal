//! \file examples/Arrangement_on_surface_2/polylines.cpp
// Constructing an arrangement of polylines.

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arr_polyline_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <vector>
#include <list>

#include "arr_print.h"
/*
  Define the Arrangement traits class to be used. You can either use some user
  defined kernel and Segment_traits_2 or the defaults.
 */

// Instantiate the traits class using a user-defined kernel
// and Segment_traits_2.
typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef CGAL::Arr_segment_traits_2<Kernel>                Segment_traits_2;
typedef CGAL::Arr_polyline_traits_2<Segment_traits_2>     Geom_traits_2;

// Identical instantiation can be achieved using the default Kernel:
// typedef CGAL::Arr_polyline_traits_2<>                    Geom_traits_2;


typedef Geom_traits_2::Point_2                            Point_2;
typedef Geom_traits_2::Segment_2                          Segment_2;
typedef Geom_traits_2::Curve_2                            Polyline_2;
typedef CGAL::Arrangement_2<Geom_traits_2>                Arrangement_2;

int main()
{
  Geom_traits_2 traits;
  Arrangement_2 arr(&traits);

  Geom_traits_2::Construct_curve_2 polyline_construct =
    traits.construct_curve_2_object();

  Point_2 points1[5];
  points1[0] = Point_2(0, 0);
  points1[1] = Point_2(2, 4);
  points1[2] = Point_2(3, 0);
  points1[3] = Point_2(4, 4);
  points1[4] = Point_2(6, 0);
  Polyline_2 pi1 = polyline_construct(&points1[0], &points1[5]);

  std::list<Point_2> points2;
  points2.push_back(Point_2(1, 3));
  points2.push_back(Point_2(0, 2));
  points2.push_back(Point_2(1, 0));
  points2.push_back(Point_2(2, 1));
  points2.push_back(Point_2(3, 0));
  points2.push_back(Point_2(4, 1));
  points2.push_back(Point_2(5, 0));
  points2.push_back(Point_2(6, 2));
  points2.push_back(Point_2(5, 3));
  points2.push_back(Point_2(4, 2));
  Polyline_2 pi2 = polyline_construct(points2.begin(), points2.end());

  std::vector<Segment_2> segs;
  segs.push_back(Segment_2(Point_2(0, 2), Point_2(1, 2)));
  segs.push_back(Segment_2(Point_2(1, 2), Point_2(3, 6)));
  segs.push_back(Segment_2(Point_2(3, 6), Point_2(5, 2)));
  Polyline_2 pi3 = polyline_construct(segs.begin(), segs.end());

  insert(arr, pi1);
  insert(arr, pi2);
  insert(arr, pi3);

  print_arrangement(arr);
  return 0;
}
