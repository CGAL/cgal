
#define CGAL_TRACE_STREAM std::cerr

#include <iostream>
#include <CGAL/Octree.h>
#include <CGAL/Octree/IO.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Point_set_3.h>

#include <cassert>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point;
typedef CGAL::Point_set_3<Point> Point_set;
typedef CGAL::Octree::Octree
        <Point_set, typename Point_set::Point_map>
        Octree;

void test_1_point() {

  // Define the dataset
  Point_set points;
  points.insert({-1, -1, -1});

  // Create the octree
  Octree octree(points, points.point_map());
  octree.refine(10, 1);

  // Check that the topology matches
  Octree::Node single_node{};
  single_node.points() = octree.root().points();
  assert(single_node == octree.root());
  assert(0 == octree.max_depth_reached());

}

void test_2_points() {

  // Define the dataset
  Point_set points;
  points.insert({-1, -1, -1});
  points.insert({1, -1, -1});

  // Create the octree
  Octree octree(points, points.point_map());
  octree.refine(10, 1);

  // The octree should have been split once
  Octree::Node other{};
  other.split();
  assert(other == octree.root());
  std::cout << octree << std::endl;
  std::cout << octree.max_depth_reached() << std::endl;
  assert(1 == octree.max_depth_reached());

}

void test_4_points() {

  Point_set points;
  points.insert({1, 1, 1});
  points.insert({1, 1, 2});
  points.insert({1, 1, 3});
  points.insert({1, 1, 4});

  // Create the octree
  Octree octree(points, points.point_map());
  octree.refine(10, 1);

  // The octree should have been split once on the first level, and twice on the second
  Octree::Node other{};
  other.split();
  other[3].split();
  other[7].split();
  assert(other == octree.root());
  assert(2 == octree.max_depth_reached());
}

int main(void) {


  test_1_point();
  test_2_points();
  test_4_points();

  return EXIT_SUCCESS;
}
