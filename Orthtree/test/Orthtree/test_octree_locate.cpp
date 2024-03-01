
#define CGAL_TRACE_STREAM std::cerr

#include <CGAL/Octree.h>
#include <CGAL/Point_set_3.h>

#include <CGAL/Simple_cartesian.h>
#include <iostream>
#include <cassert>

using Kernel = CGAL::Simple_cartesian<double>;
using Point = Kernel::Point_3;
using FT = Kernel::FT;
using Point_set = CGAL::Point_set_3<Point>;
using Octree = CGAL::Octree<Kernel, Point_set, typename Point_set::Point_map>;

void test_1_point() {

  // Define the dataset
  Point_set points;
  points.insert({-1, -1, -1});

  // Create the octree
  Octree octree(points, points.point_map());
  octree.refine(10, 1);

  // Because there's only the root node, any point should be placed in it
  assert(octree.root() == octree.locate(Point(-1, -1, -1)));

  // These points would be placed outside the root node
//  assert(octree.root() == octree.locate({0, 0, 0}));
//  assert(octree.root() == octree.locate({1000, 0, -1000}));

}

void test_8_points() {

  // Define the dataset
  Point_set points;
  points.insert({-1, -1, -1});
  points.insert({1, -1, -1});
  points.insert({-1, 1, -1});
  points.insert({1, 1, -1});
  points.insert({-1, -1, 1});
  points.insert({1, -1, 1});
  points.insert({-1, 1, 1});
  points.insert({1, 1, 1});

  // Create the octree
  Octree octree(points, points.point_map());
  octree.refine(10, 1);

  // Existing points should end up in the same place
  assert(octree.node(0) == octree.locate({-1, -1, -1}));
  assert(octree.node(1) == octree.locate({1, -1, -1}));
  assert(octree.node(2) == octree.locate({-1, 1, -1}));
  assert(octree.node(3) == octree.locate({1, 1, -1}));
  assert(octree.node(4) == octree.locate({-1, -1, 1}));
  assert(octree.node(5) == octree.locate({1, -1, 1}));
  assert(octree.node(6) == octree.locate({-1, 1, 1}));
  assert(octree.node(7) == octree.locate({1, 1, 1}));

  // Points adjacent to the existing points should also end up in the same place
  assert(octree.node(0) == octree.locate({-0.99, -0.99, -0.99}));
  assert(octree.node(1) == octree.locate({0.99, -0.99, -0.99}));
  assert(octree.node(2) == octree.locate({-0.99, 0.99, -0.99}));
  assert(octree.node(3) == octree.locate({0.99, 0.99, -0.99}));
  assert(octree.node(4) == octree.locate({-0.99, -0.99, 0.99}));
  assert(octree.node(5) == octree.locate({0.99, -0.99, 0.99}));
  assert(octree.node(6) == octree.locate({-0.99, 0.99, 0.99}));
  assert(octree.node(7) == octree.locate({0.99, 0.99, 0.99}));

}

void test_10_points() {

  // Define the dataset
  Point_set points;
  points.insert({-1, -1, -1});
  points.insert({1, -1, -1});
  points.insert({-1, 1, -1});
  points.insert({1, 1, -1});
  points.insert({-1, -1, 1});
  points.insert({1, -1, 1});
  points.insert({-1, 1, 1});
  points.insert({1, 1, 1});
  points.insert({0.875, 1, -1});
  points.insert({-1, -0.75, 1});

  // Create the octree
  Octree octree(points, points.point_map());
  octree.refine(10, 1);

  // Existing points should end up in the same place
  assert(octree.node(0) == octree.locate({-1, -1, -1}));
  assert(octree.node(1) == octree.locate({1, -1, -1}));
  assert(octree.node(2) == octree.locate({-1, 1, -1}));
  assert(octree.node(3, 3, 3, 3, 3) == octree.locate({1, 1, -1}));
  assert(octree.node(4, 4, 4) == octree.locate({-1, -1, 1}));
  assert(octree.node(5) == octree.locate({1, -1, 1}));
  assert(octree.node(6) == octree.locate({-1, 1, 1}));
  assert(octree.node(7) == octree.locate({1, 1, 1}));

  // Points adjacent to the existing points might end up in different places
  assert(octree.node(0) == octree.locate({-0.99, -0.99, -0.99}));
  assert(octree.node(1) == octree.locate({0.99, -0.99, -0.99}));
  assert(octree.node(2) == octree.locate({-0.99, 0.99, -0.99}));
  assert(octree.node(3, 3, 3, 3, 3) == octree.locate({0.99, 0.99, -0.99}));
  assert(octree.node(4, 4, 4) == octree.locate({-0.99, -0.99, 0.99}));
  assert(octree.node(5) == octree.locate({0.99, -0.99, 0.99}));
  assert(octree.node(6) == octree.locate({-0.99, 0.99, 0.99}));
  assert(octree.node(7) == octree.locate({0.99, 0.99, 0.99}));

}

int main(void) {

  test_1_point();
  test_8_points();
  test_10_points();

  return EXIT_SUCCESS;
}
