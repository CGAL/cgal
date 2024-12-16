#include <fstream>
#include <iostream>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Octree.h>
#include <CGAL/Point_set_3.h>
#include <CGAL/Point_set_3/IO.h>

// Type Declarations
using Kernel = CGAL::Simple_cartesian<double>;
using Point = Kernel::Point_3;
using Point_set = CGAL::Point_set_3<Point>;
using Point_map = Point_set::Point_map;
using Octree = CGAL::Octree<Kernel, Point_set, Point_map>;

int main(int argc, char** argv) {

  // Point set will be used to hold our points
  Point_set points;

  // Load points from a file.
  std::ifstream stream((argc > 1) ? argv[1] : CGAL::data_file_path("points_3/cube.pwn"));
  stream >> points;
  if (0 == points.number_of_points()) {

    std::cerr << "Error: cannot read file" << std::endl;
    return EXIT_FAILURE;
  }
  std::cout << "loaded " << points.number_of_points() << " points\n" << std::endl;

  // Create an octree from the points
  Octree octree(points, points.point_map());

  // Build the octree using the default arguments
  octree.refine();

  // Print out a few nodes
  std::cout << "Navigation relative to the root node" << std::endl;
  std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
  std::cout << "the root node: " << std::endl;
  std::cout << octree.to_string(octree.root()) << std::endl;
  std::cout << "the first child of the root node: " << std::endl;
  std::cout << octree.to_string(octree.child(octree.root(), 0)) << std::endl;
  std::cout << "the fifth child: " << std::endl;
  std::cout << octree.to_string(octree.child(octree.root(), 4)) << std::endl;
  std::cout << "the fifth child, accessed without going through root: " << std::endl;
  std::cout << octree.to_string(octree.node(4)) << std::endl;
  std::cout << "the second child of the fourth child: " << std::endl;
  std::cout << octree.to_string(octree.child(octree.child(octree.root(), 4), 1)) << std::endl;
  std::cout << "the second child of the fourth child, accessed without going through root: " << std::endl;
  std::cout << octree.to_string(octree.node(4, 1)) << std::endl;
  std::cout << std::endl;

  // Retrieve one of the deeper children
  Octree::Node_index cur = octree.node(4, 3);
  std::cout << "Navigation relative to a child node" << std::endl;
  std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
  std::cout << "the third child of the fourth child: " << std::endl;
  std::cout << octree.to_string(cur) << std::endl;
  std::cout << "the third child: " << std::endl;
  std::cout << octree.to_string(octree.parent(cur)) << std::endl;
  std::cout << "the next sibling of the third child of the fourth child: " << std::endl;
  std::cout << octree.to_string(octree.child(octree.parent(cur), octree.local_coordinates(cur).to_ulong() + 1))
            << std::endl;

  return EXIT_SUCCESS;
}
