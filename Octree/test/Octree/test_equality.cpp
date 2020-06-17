
#define CGAL_TRACE_STREAM std::cerr

#include <iostream>
#include <CGAL/Octree.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Point_set_3.h>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point;
typedef CGAL::Point_set_3<Point> Point_set;
typedef CGAL::Octree
        <Kernel, Point_set, typename Point_set::Point_map>
        Octree;

int test_identical_trees() {

  int failures = 0;

  // Create a simple point set
  Point_set points;
  points.insert({-1, -1, -1});
  points.insert({1, -1, -1});
  points.insert({-1, 1, -1});
  points.insert({1, 1, -1});
  points.insert({-1, -1, 1});
  points.insert({1, -1, 1});
  points.insert({-1, 1, 1});
  points.insert({1, 1, 1});

  // Create a pair of trees from the same point set
  auto point_map = points.point_map();
  Octree a(points, point_map);
  Octree b(points, point_map);

  // Refine both trees using the same criteria
  a.refine(10, 1);
  b.refine(10, 1);

  // Check if those trees are considered equal
  std::cout << "test equality: identical octrees" << std::endl;
  if (a == b) {
    std::cout << "[PASS]" << std::endl;
  } else {
    failures++;
    std::cout << "[FAIL]" << std::endl;
  }

  return failures;
}

int test_different_trees() {

  int failures = 0;

  // Create a simple point set
  Point_set points;
  points.insert({-1, -1, -1});
  points.insert({1, -1, -1});
  points.insert({-1, 1, -1});
  points.insert({1, 1, -1});
  points.insert({-1, -1, 1});
  points.insert({1, -1, 1});
  points.insert({-1, 1, 1});
  points.insert({1, 1, 1});

  // Create a pair of trees from the same point set
  auto point_map = points.point_map();
  Octree a(points, point_map);
  Octree b(points, point_map);

  // Refine both trees using different criteria
  a.refine(10, 1);
  b.refine(10, 4);

  // Check if those trees are considered equal
  std::cout << "test equality: different octrees" << std::endl;
  if (!(a == b)) {
    std::cout << "[PASS]" << std::endl;
  } else {
    failures++;
    std::cout << "[FAIL]" << std::endl;
  }

  return failures;
}

int main(void) {

  int failures = 0;

  failures += test_identical_trees();
  failures += test_different_trees();

  if (0 == failures)
    return EXIT_SUCCESS;
  return failures;
}