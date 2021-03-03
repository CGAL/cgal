
#define CGAL_TRACE_STREAM std::cerr

#include <iostream>
#include <CGAL/Octree.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Point_set_3.h>

#include <cassert>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point;
typedef CGAL::Point_set_3<Point> Point_set;
typedef CGAL::Octree<Kernel, Point_set, typename Point_set::Point_map> Octree;
typedef Octree::Node Node;

void test_1_point() {

  // Define the dataset
  Point_set points;
  points.insert({-1, -1, -1});

  // Create the octree
  Octree octree(points, points.point_map());
  octree.refine(10, 1);

  // Check that the topology matches
  Node single_node = CGAL::Orthtrees::Node_access::create_node(Node(), 0);
  CGAL::Orthtrees::Node_access::points(single_node)
    = CGAL::Orthtrees::Node_access::points(octree.root());
  assert(Node::is_topology_equal(single_node, octree.root()));
  assert(0 == octree.depth());

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
  Node other = CGAL::Orthtrees::Node_access::create_node(Node(), 0);
  CGAL::Orthtrees::Node_access::split(other);
  assert(Node::is_topology_equal(other, octree.root()));
  assert(1 == octree.depth());

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
  Node other = CGAL::Orthtrees::Node_access::create_node(Node(), 0);
  CGAL::Orthtrees::Node_access::split(other);
  CGAL::Orthtrees::Node_access::split(other[3]);
  CGAL::Orthtrees::Node_access::split(other[7]);
  assert(Node::is_topology_equal(other, octree.root()));
  assert(2 == octree.depth());
}

int main(void) {


  test_1_point();
  test_2_points();
  test_4_points();

  return EXIT_SUCCESS;
}
