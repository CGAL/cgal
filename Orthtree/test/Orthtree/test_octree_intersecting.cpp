#define CGAL_TRACE_STREAM std::cerr

#include <CGAL/Octree.h>
#include <CGAL/Orthtree/Traversals.h>

#include <CGAL/point_generators_3.h>
#include <CGAL/Simple_cartesian.h>
#include <iostream>
#include <cassert>

using Kernel = CGAL::Simple_cartesian<double>;
using Point = Kernel::Point_3;
using Point_vector = std::vector<Point>;
using Octree = CGAL::Octree<Kernel, Point_vector>;

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
  Octree octree(Octree::Traits{points});

  // Build the octree
  octree.refine(10, 2);

  // Intersection with a point (not particularly useful)
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  {
    // Set the search point
    auto query = Point{1, 1, 1};

    // Get a list of nodes intersected
    std::vector<Octree::Node_index> nodes{};
    octree.intersected_nodes(query, std::back_inserter(nodes));

    // A point should only intersect one node
    assert(1 == nodes.size());

    // That node should be the node leaf that contains the point
    assert(octree.locate(Point(1, 1, 1)) == nodes[0]);
  }

  // Intersection with a sphere
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  {
    // Set the search point
    auto query = Kernel::Sphere_3(Point{1, 0.5, 1}, 1.0);

    // Get a list of nodes intersected
    std::vector<Octree::Node_index> nodes{};
    octree.intersected_nodes(query, std::back_inserter(nodes));

    // Check the results
    assert(4 == nodes.size());
    assert(octree.node(Octree::Traits::RIGHT_TOP_BACK) == nodes[0]);
    assert(octree.node(Octree::Traits::RIGHT_BOTTOM_FRONT) == nodes[1]);
    assert(octree.node(Octree::Traits::LEFT_TOP_FRONT) == nodes[2]);
    assert(octree.node(Octree::Traits::RIGHT_TOP_FRONT) == nodes[3]);
  }

  // Intersection with a ray
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  {
    // Set the search point
    auto query = Kernel::Ray_3(Point{1, 1, 1}, Point{0, 0, 0});

    // Get a list of nodes intersected
    std::vector<Octree::Node_index> nodes{};
    octree.intersected_nodes(query, std::back_inserter(nodes));

    // Check the results
    assert(8 == nodes.size());
    assert(octree.node(Octree::Traits::LEFT_BOTTOM_BACK) == nodes[0]);
    assert(
      octree.node(Octree::Traits::RIGHT_BOTTOM_BACK, Octree::Traits::LEFT_TOP_FRONT)
      == nodes[1]
    );
    assert(octree.node(Octree::Traits::LEFT_TOP_BACK) == nodes[2]);
    assert(octree.node(Octree::Traits::RIGHT_TOP_BACK) == nodes[3]);
    assert(octree.node(Octree::Traits::LEFT_BOTTOM_FRONT) == nodes[4]);
    assert(octree.node(Octree::Traits::RIGHT_BOTTOM_FRONT) == nodes[5]);
    assert(octree.node(Octree::Traits::LEFT_TOP_FRONT) == nodes[6]);
    assert(octree.node(Octree::Traits::RIGHT_TOP_FRONT) == nodes[7]);
  }

  return EXIT_SUCCESS;
}
