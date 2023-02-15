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

int main(int argc, char **argv) {

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
  std::cout << octree.root() << std::endl;
  std::cout << "the first child of the root node: " << std::endl;
  std::cout << octree.root()[0] << std::endl;
  std::cout << "the fifth child: " << std::endl;
  std::cout << octree.root()[4] << std::endl;
  std::cout << "the fifth child, accessed without the root keyword: " << std::endl;
  std::cout << octree[4] << std::endl;
  std::cout << "the second child of the fourth child: " << std::endl;
  std::cout << octree.root()[4][1] << std::endl;
  std::cout << "the second child of the fourth child, accessed without the root keyword: " << std::endl;
  std::cout << octree[4][1] << std::endl;
  std::cout << std::endl;

  // Retrieve one of the deeper children
  Octree::Node cur = octree[3][2];
  std::cout << "Navigation relative to a child node" << std::endl;
  std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
  std::cout << "the third child of the fourth child: " << std::endl;
  std::cout << cur << std::endl;
  std::cout << "the third child: " << std::endl;
  std::cout << cur.parent() << std::endl;
  std::cout << "the next sibling of the third child of the fourth child: " << std::endl;
  std::cout << cur.parent()[cur.local_coordinates().to_ulong() + 1] << std::endl;

  return EXIT_SUCCESS;
}
