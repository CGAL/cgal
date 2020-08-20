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
    // Set the search point
    auto query = Point{1, 1, 1};

    // Get a list of nodes intersected
    std::vector<const Octree::Node *> nodes{};
    octree.intersecting_nodes(query, std::back_inserter(nodes));

    // A point should only intersect one node
    assert(1 == nodes.size());

    // That node should be the node leaf that contains the point
    assert(octree.locate(Point{1, 1, 1}) == *nodes[0]);
  }

  // Intersection with a sphere
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  {
    // Set the search point
    auto query = Kernel::Sphere_3(Point{1, 1, 1}, 4.0);

    // Get a list of nodes intersected
    std::vector<const Octree::Node *> nodes{};
    octree.intersecting_nodes(query, std::back_inserter(nodes));

  }

  // Intersection with a ray
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  {
    // Set the search point
    auto query = Kernel::Ray_3(Point{1, 1, 1}, Point{0, 0, 0});

    // Get a list of nodes intersected
    std::vector<const Octree::Node *> nodes{};
    octree.intersecting_nodes(query, std::back_inserter(nodes));

  }

  // Print out the tree
  std::cout << octree;

  return EXIT_SUCCESS;
}
