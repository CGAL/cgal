#include <fstream>
#include <iostream>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Octree.h>
#include <CGAL/Point_set_3.h>
#include <CGAL/Point_set_3/IO.h>

#include <boost/iterator/function_output_iterator.hpp>

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
  std::cout << "loaded " << points.number_of_points() << " points" << std::endl;

  // Create an octree from the points
  Octree octree(points, points.point_map());

  // Build the octree
  octree.refine(10, 20);

  // Find the nearest points to a few locations
  std::vector<Point> points_to_find = {
    {0,         0,         0},
    {1,         1,         1},
    {-1,        -1,        -1},
    {-0.46026,  -0.25353,  0.32051},
    {-0.460261, -0.253533, 0.320513}
  };

  for (const Point& p : points_to_find)
    octree.nearest_k_neighbors
    (p, 1, // k=1 to find the single closest point
      boost::make_function_output_iterator
      ([&](const Point_set::Index& nearest)
        {
          std::cout << "the nearest point to (" << p <<
            ") is (" << points.point(nearest) << ")" << std::endl;
        }));

  typename Octree::Sphere s(points_to_find[0], 0.0375);
  std::cout << std::endl << "Closest points within the sphere around " << s.center() << " with squared radius of " << s.squared_radius() << ":" << std::endl;
  octree.neighbors_within_radius
  (s,
    boost::make_function_output_iterator
    ([&](const Point_set::Index& nearest)
      {
        std::cout << points.point(nearest) << "    dist: " << (Point(0, 0, 0) - points.point(nearest)).squared_length() << std::endl;
      }));

  std::size_t k = 3;
  std::cout << std::endl << "The up to " << k << " closest points to(" << s.center() << ") within a squared radius of " << s.squared_radius() << " are:" << std::endl;
  octree.nearest_k_neighbors_within_radius
  (s, k,
    boost::make_function_output_iterator
    ([&](const Point_set::Index& nearest)
      {
        std::cout << "(" << points.point(nearest) << "    dist: " << (Point(0, 0, 0) - points.point(nearest)).squared_length() << std::endl;
      }));

  return EXIT_SUCCESS;
}
