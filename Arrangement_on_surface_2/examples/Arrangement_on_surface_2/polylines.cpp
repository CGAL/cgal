//! \file examples/Arrangement_on_surface_2/polylines.cpp
// Constructing an arrangement of polylines.

#include <vector>
#include <list>

#include <CGAL/draw_arrangement_2.h>

#include "arr_polylines.h"
#include "arr_print.h"

int main() {
  Traits traits;
  Arrangement arr(&traits);

  auto polyline_construct = traits.construct_curve_2_object();

  Point points1[5];
  points1[0] = Point(0, 0);
  points1[1] = Point(2, 4);
  points1[2] = Point(3, 0);
  points1[3] = Point(4, 4);
  points1[4] = Point(6, 0);
  auto pi1 = polyline_construct(&points1[0], &points1[5]);

  std::list<Point> points2;
  points2.push_back(Point(1, 3));
  points2.push_back(Point(0, 2));
  points2.push_back(Point(1, 0));
  points2.push_back(Point(2, 1));
  points2.push_back(Point(3, 0));
  points2.push_back(Point(4, 1));
  points2.push_back(Point(5, 0));
  points2.push_back(Point(6, 2));
  points2.push_back(Point(5, 3));
  points2.push_back(Point(4, 2));
  auto pi2 = polyline_construct(points2.begin(), points2.end());

  std::vector<Segment> segs;
  segs.push_back(Segment(Point(0, 2), Point(1, 2)));
  segs.push_back(Segment(Point(1, 2), Point(3, 6)));
  segs.push_back(Segment(Point(3, 6), Point(5, 2)));
  auto pi3 = polyline_construct(segs.begin(), segs.end());

  insert(arr, pi1);
  insert(arr, pi2);
  insert(arr, pi3);
  print_arrangement_size(arr);          // print the arrangement size
  CGAL::draw(arr);
  return 0;
}
