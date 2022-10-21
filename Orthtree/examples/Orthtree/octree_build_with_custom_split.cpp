#include <fstream>
#include <iostream>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Octree.h>
#include <CGAL/Point_set_3.h>
#include <CGAL/Point_set_3/IO.h>

// Type Declarations
typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point;
typedef Kernel::FT FT;
typedef CGAL::Point_set_3<Point> Point_set;
typedef Point_set::Point_map Point_map;

typedef CGAL::Octree<Kernel, Point_set, Point_map> Octree;

// Split Predicate
// The predicate is a functor which returns a boolean value, whether a node needs to be split or not
struct Split_by_ratio {

  std::size_t ratio;

  Split_by_ratio(std::size_t ratio) : ratio(ratio) {}

  template<class Node>
  bool operator()(const Node &n) const {
    return n.size() > (ratio * n.depth());
  }
};

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
  std::cout << "loaded " << points.number_of_points() << " points" << std::endl;

  // Create an octree from the points
  Octree octree(points, points.point_map());

  // Build the octree using our custom split predicate
  octree.refine(Split_by_ratio(2));

  // Print out the tree
  std::cout << octree;

  return EXIT_SUCCESS;
}
