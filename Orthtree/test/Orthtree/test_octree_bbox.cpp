
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

void test_1_node() {

  // Define the dataset
  Point_set points;
  points.insert({-1, -1, -1});

  // Create the octree
  Octree octree(points, points.point_map());
  octree.refine(10, 1);

  // Compare the top (only) node
  assert(octree.bbox(octree.root()) == CGAL::Bbox_3(-1, -1, -1, -1, -1, -1));
}

void test_9_nodes() {

  // Define the dataset
  Point_set points;
  points.insert({-1, -1, -1});
  points.insert({1, 1, 1});

  // Create the octree
  Octree octree(points, points.point_map(), 1.1);
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

void test_25_nodes() {

  // Define the dataset
  Point_set points;
  points.insert({-1, -1, -1});
  points.insert({1, 1, 1});
  points.insert({-1, -1, -0.5});
  points.insert({1, 0.5, 1});

  // Create the octree
  Octree octree(points, points.point_map(), 1.5);
  octree.refine(10, 1);

  // Compare the top node
  assert(octree.bbox(octree.root()) == CGAL::Bbox_3(-1.5, -1.5, -1.5, 1.5, 1.5, 1.5));

  // Compare the child nodes
  assert(octree.bbox(octree.root()[0]) == CGAL::Bbox_3(-1.5, -1.5, -1.5, 0, 0, 0));
  assert(octree.bbox(octree.root()[1]) == CGAL::Bbox_3(0, -1.5, -1.5, 1.5, 0, 0));
  assert(octree.bbox(octree.root()[2]) == CGAL::Bbox_3(-1.5, 0, -1.5, 0, 1.5, 0));
  assert(octree.bbox(octree.root()[3]) == CGAL::Bbox_3(0, 0, -1.5, 1.5, 1.5, 0));
  assert(octree.bbox(octree.root()[4]) == CGAL::Bbox_3(-1.5, -1.5, 0, 0, 0, 1.5));
  assert(octree.bbox(octree.root()[5]) == CGAL::Bbox_3(0, -1.5, 0, 1.5, 0, 1.5));
  assert(octree.bbox(octree.root()[6]) == CGAL::Bbox_3(-1.5, 0, 0, 0, 1.5, 1.5));
  assert(octree.bbox(octree.root()[7]) == CGAL::Bbox_3(0, 0, 0, 1.5, 1.5, 1.5));

  // Compare children of the first child
  assert(octree.bbox(octree.root()[0][0]) == CGAL::Bbox_3(-1.5, -1.5, -1.5, -0.75, -0.75, -0.75));
  assert(octree.bbox(octree.root()[0][1]) == CGAL::Bbox_3(-0.75, -1.5, -1.5, 0, -0.75, -0.75));
  assert(octree.bbox(octree.root()[0][2]) == CGAL::Bbox_3(-1.5, -0.75, -1.5, -0.75, 0, -0.75));
  assert(octree.bbox(octree.root()[0][3]) == CGAL::Bbox_3(-0.75, -0.75, -1.5, 0, 0, -0.75));
  assert(octree.bbox(octree.root()[0][4]) == CGAL::Bbox_3(-1.5, -1.5, -0.75, -0.75, -0.75, 0));
  assert(octree.bbox(octree.root()[0][5]) == CGAL::Bbox_3(-0.75, -1.5, -0.75, 0, -0.75, 0));
  assert(octree.bbox(octree.root()[0][6]) == CGAL::Bbox_3(-1.5, -0.75, -0.75, -0.75, 0, 0));
  assert(octree.bbox(octree.root()[0][7]) == CGAL::Bbox_3(-0.75, -0.75, -0.75, 0, 0, 0));

  // Compare children of the last child
  assert(octree.bbox(octree.root()[7][0]) == CGAL::Bbox_3(0, 0, 0, 0.75, 0.75, 0.75));
  assert(octree.bbox(octree.root()[7][1]) == CGAL::Bbox_3(0.75, 0, 0, 1.5, 0.75, 0.75));
  assert(octree.bbox(octree.root()[7][2]) == CGAL::Bbox_3(0, 0.75, 0, 0.75, 1.5, 0.75));
  assert(octree.bbox(octree.root()[7][3]) == CGAL::Bbox_3(0.75, 0.75, 0, 1.5, 1.5, 0.75));
  assert(octree.bbox(octree.root()[7][4]) == CGAL::Bbox_3(0, 0, 0.75, 0.75, 0.75, 1.5));
  assert(octree.bbox(octree.root()[7][5]) == CGAL::Bbox_3(0.75, 0, 0.75, 1.5, 0.75, 1.5));
  assert(octree.bbox(octree.root()[7][6]) == CGAL::Bbox_3(0, 0.75, 0.75, 0.75, 1.5, 1.5));
  assert(octree.bbox(octree.root()[7][7]) == CGAL::Bbox_3(0.75, 0.75, 0.75, 1.5, 1.5, 1.5));
}

int main(void) {

  test_1_node();
  test_9_nodes();
  test_25_nodes();

  return EXIT_SUCCESS;
}
