#define CGAL_TRACE_STREAM std::cerr

#include <iostream>
#include <CGAL/Octree.h>
#include <CGAL/Octree/IO.h>
#include <CGAL/Octree/Traversal.h>
#include <CGAL/Simple_cartesian.h>

#include <cassert>
#include <CGAL/point_generators_3.h>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point;
typedef std::vector<Point> Point_vector;
typedef CGAL::Octree::Octree<Point_vector> Octree;

int main(void) {

  // Fill a vector with points
  Point_vector points;
  points.emplace_back(1, 1, 1);
  points.emplace_back(-1, 1, 1);
  points.emplace_back(1, -1, 1);
  points.emplace_back(-1, -1, 1);
  points.emplace_back(1, 1, -1);
  points.emplace_back(-1, 1, -1);
  points.emplace_back(1, -1, -1);
  points.emplace_back(-1, -1, -1);
  points.emplace_back(0.9, -1, -1);
  points.emplace_back(0.9, -0.9, -1);
  points.emplace_back(0.9, -0.95, -1);
  points.emplace_back(0.9, -0.9, -0.9);
  points.emplace_back(0.9, -0.95, -0.9);
  points.emplace_back(0.9, -1, -1);
  points.emplace_back(-0.9, -1, -1);

  // Create an octree from the vector
  Octree octree(points);

  // Build the octree
  octree.refine(10, 2);

  // Intersection with a point (not particularly useful)
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  {
    std::vector<const Octree::Node *> nodes{};

  }

  // Intersection with a sphere
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  {
    std::vector<const Octree::Node *> nodes{};

  }

  // Intersection with a ray
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  {
    std::vector<const Octree::Node *> nodes{};

  }

  // Print out the tree
  std::cout << octree;

  return EXIT_SUCCESS;
}
