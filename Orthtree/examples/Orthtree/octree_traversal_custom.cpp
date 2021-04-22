#include <fstream>
#include <iostream>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Octree.h>
#include <CGAL/Point_set_3.h>
#include <CGAL/Point_set_3/IO.h>

// Type Declarations
typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point;
typedef CGAL::Point_set_3<Point> Point_set;
typedef Point_set::Point_map Point_map;

typedef CGAL::Octree<Kernel, Point_set, Point_map> Octree;
typedef Octree::Node Node;

struct First_branch_traversal
{
  Node first (Node root) const
  {
    return root;
  }

  Node next (Node n) const
  {
    if (n.is_leaf())
      return Node();
    return n[0];
  }
};

int main(int argc, char **argv) {

  // Point set will be used to hold our points
  Point_set points;

  // Load points from a file.
  std::ifstream stream((argc > 1) ? argv[1] : "data/cube.pwn");
  stream >> points;
  if (0 == points.number_of_points()) {

    std::cerr << "Error: cannot read file" << std::endl;
    return EXIT_FAILURE;
  }
  std::cout << "loaded " << points.number_of_points() << " points\n" << std::endl;

  // Create an octree from the points
  Octree octree(points, points.point_map());

  // Build the octree
  octree.refine();

  // Print out the first branch using custom traversal
  for (auto &node : octree.traverse<First_branch_traversal>())
    std::cout << node << std::endl;

  return EXIT_SUCCESS;
}
