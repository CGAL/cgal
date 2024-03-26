#include <fstream>
#include <iostream>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Octree.h>
#include <CGAL/Point_set_3.h>
#include <CGAL/Point_set_3/IO.h>

// Type Declarations
using Kernel = CGAL::Simple_cartesian<double>;
using Point = Kernel::Point_3;
using FT = Kernel::FT;
using Point_set = CGAL::Point_set_3<Point>;
using Point_map = Point_set::Point_map;
using Octree = CGAL::Octree<Kernel, Point_set, Point_map>;

// Split Predicate
// The predicate is a functor which returns a Boolean value, whether a node needs to be split or not
struct Split_by_ratio {

  std::size_t ratio;

  explicit Split_by_ratio(std::size_t ratio) : ratio(ratio) {}

  template<typename Node_index, typename Tree>
  bool operator()(Node_index i, const Tree &tree) const {
    std::size_t num_points = tree.data(i).size();
    std::size_t depth = tree.depth(i);
    return num_points > (ratio * depth);
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
