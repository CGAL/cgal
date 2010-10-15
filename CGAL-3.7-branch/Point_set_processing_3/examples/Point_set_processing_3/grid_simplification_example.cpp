#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/grid_simplify_point_set.h>
#include <CGAL/IO/read_xyz_points.h>

#include <vector>
#include <fstream>

// types
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point;

int main(void)
{
  // Reads a .xyz point set file in points[].
  std::vector<Point> points;
  std::ifstream stream("data/oni.xyz");
  if (!stream ||
      !CGAL::read_xyz_points(stream, std::back_inserter(points)))
  {
    std::cerr << "Error: cannot read file data/oni.xyz" << std::endl;
    return EXIT_FAILURE;
  }

  // simplification by clustering using erase-remove idiom
  double cell_size = 0.001;
  points.erase(CGAL::grid_simplify_point_set(points.begin(), points.end(), cell_size),
               points.end());

  // Optional: after erase(), use Scott Meyer's "swap trick" to trim excess capacity
  std::vector<Point>(points).swap(points);

  return EXIT_SUCCESS;
}

