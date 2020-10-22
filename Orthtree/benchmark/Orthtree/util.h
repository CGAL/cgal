#ifndef OCTREE_UTIL_H
#define OCTREE_UTIL_H

#include <CGAL/Point_set_3.h>
#include <CGAL/Point_set_3/IO.h>
#include <CGAL/point_generators_3.h>

template<class Kernel>
CGAL::Point_set_3<typename Kernel::Point_3> generate(size_t num_points = 1) {

  typedef typename Kernel::Point_3 Point;
  typedef CGAL::Point_set_3<Point> Point_set;

  // Create an empty point set
  Point_set points;
  points.reserve(num_points);

  // Fill the point set with random points
  CGAL::Random_points_on_sphere_3<Point> generator;
  for (size_t i = 0; i < num_points; ++i)
    points.insert(*generator++);

  return points;
}

#endif //OCTREE_UTIL_H
