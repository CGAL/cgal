
#define CGAL_TRACE_STREAM std::cerr

#include <iostream>
#include <CGAL/Octree.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Point_set_3.h>

#include <cassert>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point;
typedef Kernel::FT FT;
typedef CGAL::Point_set_3<Point> Point_set;
typedef CGAL::Octree<Kernel, Point_set, typename Point_set::Point_map>
Octree;

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
  assert(octree.root()[0] == octree.locate({-1, -1, -1}));
  assert(octree.root()[1] == octree.locate({1, -1, -1}));
  assert(octree.root()[2] == octree.locate({-1, 1, -1}));
  assert(octree.root()[3] == octree.locate({1, 1, -1}));
  assert(octree.root()[4] == octree.locate({-1, -1, 1}));
  assert(octree.root()[5] == octree.locate({1, -1, 1}));
  assert(octree.root()[6] == octree.locate({-1, 1, 1}));
  assert(octree.root()[7] == octree.locate({1, 1, 1}));

  // Points adjacent to the existing points should also end up in the same place
  assert(octree.root()[0] == octree.locate({-1.1, -1.1, -1.1}));
  assert(octree.root()[1] == octree.locate({1.1, -1.1, -1.1}));
  assert(octree.root()[2] == octree.locate({-1.1, 1.1, -1.1}));
  assert(octree.root()[3] == octree.locate({1.1, 1.1, -1.1}));
  assert(octree.root()[4] == octree.locate({-1.1, -1.1, 1.1}));
  assert(octree.root()[5] == octree.locate({1.1, -1.1, 1.1}));
  assert(octree.root()[6] == octree.locate({-1.1, 1.1, 1.1}));
  assert(octree.root()[7] == octree.locate({1.1, 1.1, 1.1}));

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
  assert(octree.root()[0] == octree.locate({-1, -1, -1}));
  assert(octree.root()[1] == octree.locate({1, -1, -1}));
  assert(octree.root()[2] == octree.locate({-1, 1, -1}));
  assert(octree.root()[3][3][3] == octree.locate({1, 1, -1}));
  assert(octree.root()[4][4][4] == octree.locate({-1, -1, 1}));
  assert(octree.root()[5] == octree.locate({1, -1, 1}));
  assert(octree.root()[6] == octree.locate({-1, 1, 1}));
  assert(octree.root()[7] == octree.locate({1, 1, 1}));

  // Points adjacent to the existing points might end up in different places
  assert(octree.root()[0] == octree.locate({-1.1, -1.1, -1.1}));
  assert(octree.root()[1] == octree.locate({1.1, -1.1, -1.1}));
  assert(octree.root()[2] == octree.locate({-1.1, 1.1, -1.1}));
  assert(octree.root()[3][3][3] == octree.locate({1.1, 1.1, -1.1}));
  assert(octree.root()[4][4][4] == octree.locate({-1.1, -1.1, 1.1}));
  assert(octree.root()[5] == octree.locate({1.1, -1.1, 1.1}));
  assert(octree.root()[6] == octree.locate({-1.1, 1.1, 1.1}));
  assert(octree.root()[7] == octree.locate({1.1, 1.1, 1.1}));

}

int main(void) {

  test_1_point();
  test_8_points();
  test_10_points();

  return EXIT_SUCCESS;
}
