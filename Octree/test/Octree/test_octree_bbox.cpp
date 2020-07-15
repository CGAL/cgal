
#include <iostream>
#include <CGAL/Octree.h>
#include <CGAL/Octree/IO.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Point_set_3.h>

#include <cassert>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point;
typedef Kernel::FT FT;
typedef CGAL::Point_set_3<Point> Point_set;
typedef CGAL::Octree::Octree
        <Point_set, typename Point_set::Point_map>
        Octree;

void test_1_node() {

  // Define the dataset
  Point_set points;
  points.insert({-1, -1, -1});
  auto point_map = points.point_map();

  // Create the octree
  Octree octree(points, point_map);
  octree.refine(10, 1);

  // Compare the top (only) node
  assert(octree.bbox(octree.root()) == CGAL::Bbox_3(-1, -1, -1, -1, -1, -1));
}

void test_9_nodes() {

  // Define the dataset
  Point_set points;
  points.insert({-1, -1, -1});
  points.insert({1, 1, 1});
  auto point_map = points.point_map();

  // Create the octree
  Octree octree(points, point_map, 1.1);
  octree.refine(10, 1);

  // Compare the top node
  assert(octree.bbox(octree.root()) == CGAL::Bbox_3(-1.1, -1.1, -1.1, 1.1, 1.1, 1.1));
}

int main(void) {

  test_1_node();
  test_9_nodes();

  return EXIT_SUCCESS;
}
