#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/property_map.h>
#include <CGAL/remove_outliers.h>
#include <CGAL/IO/read_xyz_points.h>

#include <vector>
#include <fstream>

// types
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point;

int main(void)
{
  // Reads a .xyz point set file in points[].
  // The Dereference_property_map property map can be omitted here as it is the default value.
  std::vector<Point> points;
  std::ifstream stream("data/oni.xyz");
  if (!stream ||
      !CGAL::read_xyz_points(stream, std::back_inserter(points),
                             CGAL::Dereference_property_map<Point>()))
  {
    std::cerr << "Error: cannot read file data/oni.xyz" << std::endl;
    return EXIT_FAILURE;
  }

  // Removes outliers using erase-remove idiom.
  // The Dereference_property_map property map can be omitted here as it is the default value.
  const double removed_percentage = 5.0; // percentage of points to remove
  const int nb_neighbors = 24; // considers 24 nearest neighbor points
  points.erase(CGAL::remove_outliers(points.begin(), points.end(),
                                     CGAL::Dereference_property_map<Point>(),
                                     nb_neighbors, removed_percentage), 
               points.end());

  // Optional: after erase(), use Scott Meyer's "swap trick" to trim excess capacity
  std::vector<Point>(points).swap(points);

  return EXIT_SUCCESS;
}

