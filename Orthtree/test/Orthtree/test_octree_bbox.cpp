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

void test_1_node() {

  // Define the dataset
  Point_set points;
  points.insert({-1, -1, -1});

  // Create the octree
  Octree octree(points, points.point_map());
  octree.refine(10, 1);

  Octree::Bbox expected_bbox{-1, -1, -1, -1, -1, -1};

  // Compare the top (only) node
  assert(octree.bbox(octree.root()) == Octree::Bbox(-1, -1, -1, -1, -1, -1));
}

void test_9_nodes() {

  // Define the dataset
  Point_set points;
  points.insert({-1, -1, -1});
  points.insert({1, 1, 1});

  // Create the octree
  Octree octree(points, points.point_map());
  octree.refine(10, 1);

  // Compare the top node
  assert(octree.bbox(octree.root()) == Octree::Bbox(-1, -1, -1, 1, 1, 1));

  // Compare the child nodes
  assert(octree.bbox(octree.node(0)) == Octree::Bbox(-1, -1, -1, 0, 0, 0));
  assert(octree.bbox(octree.node(1)) == Octree::Bbox(0, -1, -1, 1, 0, 0));
  assert(octree.bbox(octree.node(2)) == Octree::Bbox(-1, 0, -1, 0, 1, 0));
  assert(octree.bbox(octree.node(3)) == Octree::Bbox(0, 0, -1, 1, 1, 0));
  assert(octree.bbox(octree.node(4)) == Octree::Bbox(-1, -1, 0, 0, 0, 1));
  assert(octree.bbox(octree.node(5)) == Octree::Bbox(0, -1, 0, 1, 0, 1));
  assert(octree.bbox(octree.node(6)) == Octree::Bbox(-1, 0, 0, 0, 1, 1));
  assert(octree.bbox(octree.node(7)) == Octree::Bbox(0, 0, 0, 1, 1, 1));
}

void test_25_nodes() {

  // Define the dataset
  Point_set points;
  points.insert({-1, -1, -1});
  points.insert({1, 1, 1});
  points.insert({-1, -1, -0.5});
  points.insert({1, 0.5, 1});

  // Create the octree
  Octree octree(points, points.point_map());
  octree.refine(10, 1);

  // Compare the top node
  assert(octree.bbox(octree.root()) == Octree::Bbox(-1, -1, -1, 1, 1, 1));

  // Compare the child nodes
  assert(octree.bbox(octree.node(0)) == Octree::Bbox(-1, -1, -1, 0, 0, 0));
  assert(octree.bbox(octree.node(1)) == Octree::Bbox(0, -1, -1, 1, 0, 0));
  assert(octree.bbox(octree.node(2)) == Octree::Bbox(-1, 0, -1, 0, 1, 0));
  assert(octree.bbox(octree.node(3)) == Octree::Bbox(0, 0, -1, 1, 1, 0));
  assert(octree.bbox(octree.node(4)) == Octree::Bbox(-1, -1, 0, 0, 0, 1));
  assert(octree.bbox(octree.node(5)) == Octree::Bbox(0, -1, 0, 1, 0, 1));
  assert(octree.bbox(octree.node(6)) == Octree::Bbox(-1, 0, 0, 0, 1, 1));
  assert(octree.bbox(octree.node(7)) == Octree::Bbox(0, 0, 0, 1, 1, 1));

  // Compare children of the first child
  assert(octree.bbox(octree.node(0, 0)) ==
         Octree::Bbox(-1, -1, -1, -0.5, -0.5, -0.5));
  assert(octree.bbox(octree.node(0, 1)) ==
         Octree::Bbox(-0.5, -1, -1, 0, -0.5, -0.5));
  assert(octree.bbox(octree.node(0, 2)) ==
         Octree::Bbox(-1, -0.5, -1, -0.5, 0, -0.5));
  assert(octree.bbox(octree.node(0, 3)) ==
         Octree::Bbox(-0.5, -0.5, -1, 0, 0, -0.5));
  assert(octree.bbox(octree.node(0, 4)) ==
         Octree::Bbox(-1, -1, -0.5, -0.5, -0.5, 0));
  assert(octree.bbox(octree.node(0, 5)) ==
         Octree::Bbox(-0.5, -1, -0.5, 0, -0.5, 0));
  assert(octree.bbox(octree.node(0, 6)) ==
         Octree::Bbox(-1, -0.5, -0.5, -0.5, 0, 0));
  assert(octree.bbox(octree.node(0, 7)) ==
         Octree::Bbox(-0.5, -0.5, -0.5, 0, 0, 0));

  // Compare children of the last child
  assert(octree.bbox(octree.node(7, 0)) ==
         Octree::Bbox(0, 0, 0, 0.5, 0.5, 0.5));
  assert(octree.bbox(octree.node(7, 1)) ==
         Octree::Bbox(0.5, 0, 0, 1, 0.5, 0.5));
  assert(octree.bbox(octree.node(7, 2)) ==
         Octree::Bbox(0, 0.5, 0, 0.5, 1, 0.5));
  assert(octree.bbox(octree.node(7, 3)) ==
         Octree::Bbox(0.5, 0.5, 0, 1, 1, 0.5));
  assert(octree.bbox(octree.node(7, 4)) ==
         Octree::Bbox(0, 0, 0.5, 0.5, 0.5, 1));
  assert(octree.bbox(octree.node(7, 5)) ==
         Octree::Bbox(0.5, 0, 0.5, 1, 0.5, 1));
  assert(octree.bbox(octree.node(7, 6)) ==
         Octree::Bbox(0, 0.5, 0.5, 0.5, 1, 1));
  assert(octree.bbox(octree.node(7, 7)) ==
         Octree::Bbox(0.5, 0.5, 0.5, 1, 1, 1));

  // All child nodes should share a vertex
  auto center_of_last_child = octree.bbox(octree.node(7, 7)).vertex(0);
  assert(octree.bbox(octree.node(7, 0)).vertex(7) == center_of_last_child);
  assert(octree.bbox(octree.node(7, 1)).vertex(4) == center_of_last_child);
  assert(octree.bbox(octree.node(7, 2)).vertex(6) == center_of_last_child);
  assert(octree.bbox(octree.node(7, 3)).vertex(5) == center_of_last_child);
  assert(octree.bbox(octree.node(7, 4)).vertex(2) == center_of_last_child);
  assert(octree.bbox(octree.node(7, 5)).vertex(3) == center_of_last_child);
  assert(octree.bbox(octree.node(7, 6)).vertex(1) == center_of_last_child);

  // Nodes of different sizes should also share vertices
  assert(octree.bbox(octree.node(7, 0)).vertex(0) == octree.bbox(octree.node(0, 7)).vertex(7));
}

int main(void) {

  test_1_node();
  test_9_nodes();
  test_25_nodes();

  return EXIT_SUCCESS;
}
