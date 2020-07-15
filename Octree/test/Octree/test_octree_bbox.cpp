
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

  // Compare the child nodes
  assert(octree.bbox(octree.root()[0]) == CGAL::Bbox_3(-1.1, -1.1, -1.1, 0, 0, 0));
  assert(octree.bbox(octree.root()[1]) == CGAL::Bbox_3(0, -1.1, -1.1, 1.1, 0, 0));
  assert(octree.bbox(octree.root()[2]) == CGAL::Bbox_3(-1.1, 0, -1.1, 0, 1.1, 0));
  assert(octree.bbox(octree.root()[3]) == CGAL::Bbox_3(0, 0, -1.1, 1.1, 1.1, 0));
  assert(octree.bbox(octree.root()[4]) == CGAL::Bbox_3(-1.1, -1.1, 0, 0, 0, 1.1));
  assert(octree.bbox(octree.root()[5]) == CGAL::Bbox_3(0, -1.1, 0, 1.1, 0, 1.1));
  assert(octree.bbox(octree.root()[6]) == CGAL::Bbox_3(-1.1, 0, 0, 0, 1.1, 1.1));
  assert(octree.bbox(octree.root()[7]) == CGAL::Bbox_3(0, 0, 0, 1.1, 1.1, 1.1));
}

int main(void) {

  test_1_node();
  test_9_nodes();

  return EXIT_SUCCESS;
}
