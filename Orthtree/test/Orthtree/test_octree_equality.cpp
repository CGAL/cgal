
#define CGAL_TRACE_STREAM std::cerr

#include <CGAL/Octree.h>
#include <CGAL/Point_set_3.h>
#include <CGAL/Simple_cartesian.h>
#include <iostream>

using Kernel = CGAL::Simple_cartesian<double>;
using Point = Kernel::Point_3;
using Point_set = CGAL::Point_set_3<Point>;
using Octree = CGAL::Octree<Kernel, Point_set, typename Point_set::Point_map>;

void test_identical_trees() {

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
  Octree a(points, points.point_map());
  Octree b(points, points.point_map());

  // Refine both trees using the same criteria
  a.refine(10, 1);
  b.refine(10, 1);

  // Check if those trees are considered equal
  assert(a == b);

}

void test_identical_contents_different_criteria() {

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
  Octree a(points, points.point_map());
  Octree b(points, points.point_map());

  // Refine both trees using different criteria
  a.refine(10, 1);
  b.refine(10, 9);

  // Check if those trees are considered equal
  assert(a != b);
}

void test_different_contents_identical_criteria() {

  // Create a couple of simple point sets
  Point_set points_a;
  points_a.insert({-1, -1, -1});
  points_a.insert({1, -1, -1});
  points_a.insert({-1, 1, -1});
  points_a.insert({1, 1, -1});
  points_a.insert({-1, -1, 1});
  points_a.insert({1, -1, 1});
  points_a.insert({-1, 1, 1});
  points_a.insert({1, 1, 1});

  Point_set points_b;
  points_b.insert({-1, -1, -1});
  points_b.insert({1, -1, -1});
  points_b.insert({-1, 1, -1});
  points_b.insert({1, 1, -1});
  points_b.insert({-1, -1, 1});
  points_b.insert({1, -1, 1});
  points_b.insert({-1, 1, 1});
  points_b.insert({1, 1, 2});

  // Create a pair of trees from the different point sets
  Octree a(points_a, points_a.point_map());
  Octree b(points_b, points_b.point_map());

  // Refine both trees using the same criteria
  a.refine(10, 1);
  b.refine(10, 1);

  // Check if those trees are considered equal
  assert(a != b);
}


int main(void) {


  test_identical_trees();

  test_identical_contents_different_criteria();

  test_different_contents_identical_criteria();

  return EXIT_SUCCESS;
}
