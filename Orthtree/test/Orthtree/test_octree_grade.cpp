
#define CGAL_TRACE_STREAM std::cerr

#include <iostream>
#include <CGAL/Octree.h>
#include <CGAL/Orthtree/Traversals.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Point_set_3.h>

#include <cassert>
#include <CGAL/point_generators_3.h>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point;
typedef Kernel::FT FT;
typedef CGAL::Point_set_3<Point> Point_set;
typedef CGAL::Octree<Kernel, Point_set, typename Point_set::Point_map> Octree;
typedef CGAL::Orthtrees::Leaves_traversal Leaves_traversal;

std::size_t count_jumps(Octree &octree) {

  std::size_t jumps = 0;

  for (auto &node : octree.traverse(Leaves_traversal())) {

    for (int direction = 0; direction < 6; ++direction) {

      auto adjacent_node = node.adjacent_node(direction);

      if (adjacent_node.is_null())
        continue;

      if ((node.depth() - adjacent_node.depth()) > 1)
        jumps++;
    }
  }

  return jumps;
}

void test(std::size_t dataset_size) {

  // Create a dataset
  Point_set points;
  CGAL::Random_points_in_cube_3<Point> generator;
  points.reserve(dataset_size);
  for (std::size_t i = 0; i < dataset_size; ++i)
    points.insert(*(generator++));

  // Build an octree
  Octree octree(points, points.point_map());

  // Refine the octree
  octree.refine();

  // Count the jumps in depth
  auto jumps = count_jumps(octree);
  std::cout << "un-graded octree has " << jumps << " jumps" << std::endl;

  // Grade the octree
  octree.grade();

  // Count the jumps in depth
  jumps = count_jumps(octree);
  std::cout << "graded octree has " << jumps << " jumps" << std::endl;
  assert(jumps == 0);
}

int main(void) {

  test(100000);

  return EXIT_SUCCESS;
}
