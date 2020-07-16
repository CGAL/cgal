
#define CGAL_TRACE_STREAM std::cerr

#include <iostream>
#include <CGAL/Octree.h>
#include <CGAL/Octree/IO.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Point_set_3.h>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point;
typedef CGAL::Point_set_3<Point> Point_set;
typedef CGAL::Octree::Octree
        <Point_set, typename Point_set::Point_map>
        Octree;

void test_identical_trees() {

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
  assert(a == b);

}

void test_identical_contents_different_criteria() {

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
  auto point_map_a = points_a.point_map();
  Octree a(points_a, point_map_a);
  auto point_map_b = points_b.point_map();
  Octree b(points_b, point_map_b);

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