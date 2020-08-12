#include <fstream>
#include <iostream>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Octree.h>
#include <CGAL/Octree/IO.h>

// Type Declarations
typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point;
typedef std::vector<Point> Point_vector;
typedef CGAL::Octree::Octree<Point_vector> Octree;

int main(int argc, char **argv) {

  // Here, our point set is a vector
  Point_vector points;

  // Add a few points to the vector
  points.emplace_back(1, 1, 1);
  points.emplace_back(2, 1, -11);
  points.emplace_back(2, 1, 1);
  points.emplace_back(1, -2, 1);
  points.emplace_back(1, 1, 1);
  points.emplace_back(-1, 1, 1);

  // Create an octree from the points
  Octree octree(points);

  // Build the octree with a small bucket size, using a more verbose method
  octree.refine(CGAL::Octree::Split_criterion::Max_depth_or_bucket_size(10, 2));

  // Print out the tree
  std::cout << octree;

  return EXIT_SUCCESS;
}
