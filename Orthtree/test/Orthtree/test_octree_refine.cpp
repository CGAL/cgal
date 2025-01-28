
#define CGAL_TRACE_STREAM std::cerr

#include <CGAL/Octree.h>
#include <CGAL/Point_set_3.h>

#include <CGAL/Simple_cartesian.h>
#include <iostream>
#include <cassert>

using Kernel = CGAL::Simple_cartesian<double>;
using Point = Kernel::Point_3;
using Point_set = CGAL::Point_set_3<Point>;
using Octree = CGAL::Octree<Kernel, Point_set, typename Point_set::Point_map>;

class Split_nth_child_of_root {
  std::size_t m_n;
public:

  Split_nth_child_of_root(std::size_t n) : m_n(n) {}

  template <typename Node>
  bool operator()(const Node& node) const {
    return (node.depth() == 1 && node.local_coordinates().to_ulong() == m_n);
  }

  template<typename Node_index, typename Tree>
  bool operator()(Node_index i, const Tree &tree) const {
    return (tree.depth(i) == 1 && tree.local_coordinates(i).to_ulong() == m_n);
  }
};

void test_1_point() {

  // Define the dataset
  Point_set points;
  points.insert({-1, -1, -1});

  // Create the octree
  Octree octree(points, points.point_map());
  octree.refine(10, 1);

  // Check that the root node was never split
  assert(octree.is_leaf(octree.root()));
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
  Octree other(points, points.point_map());
  other.split(other.root());
  assert(Octree::is_topology_equal(other, octree));
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

  Octree other(points, points.point_map());
  other.split(other.root());
  other.split(other.node(3));
  other.split(other.node(7));
  assert(Octree::is_topology_equal(other, octree));
  assert(2 == octree.depth());

  // Applying another splitting criterion shouldn't reset the tree.
  octree.refine(Split_nth_child_of_root(2));
  other.split(other.node(2));
  assert(Octree::is_topology_equal(other, octree));

}

int main(void) {


  test_1_point();
  test_2_points();
  test_4_points();

  return EXIT_SUCCESS;
}
