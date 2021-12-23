#ifndef OCTREE_UTIL_H
#define OCTREE_UTIL_H

#include <CGAL/Point_set_3.h>
#include <CGAL/Point_set_3/IO.h>
#include <CGAL/point_generators_3.h>

#include <chrono>

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

template<typename Time_unit>
Time_unit bench(const std::function<void(void)> &f, size_t repetitions = 1) {

  // Start the timer
  auto start = std::chrono::high_resolution_clock::now();

  // Run the function being benchmarked as many times as requested
  for (int i = 0; i < repetitions; ++i)
    f();

  // End the timer
  auto end = std::chrono::high_resolution_clock::now();

  // Return the elapsed time
  return std::chrono::duration_cast<Time_unit>(end - start) / repetitions;
}

#endif //OCTREE_UTIL_H
