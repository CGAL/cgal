#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/property_map.h>
#include <CGAL/compute_average_spacing.h>
#include <CGAL/remove_outliers.h>
#include <CGAL/IO/read_xyz_points.h>

#include <vector>
#include <fstream>
#include <iostream>

// types
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point;

int main(int argc, char*argv[])
{
  const char* fname = (argc>1)?argv[1]:"data/oni.xyz";
  // Reads a .xyz point set file in points[].
  // The Identity_property_map property map can be omitted here as it is the default value.
  std::vector<Point> points;
  std::ifstream stream(fname);
  if (!stream ||
      !CGAL::read_xyz_points(stream, std::back_inserter(points),
                             CGAL::parameters::point_map(CGAL::Identity_property_map<Point>())))
  {
    std::cerr << "Error: cannot read file " << fname << std::endl;
    return EXIT_FAILURE;
  }

  // Removes outliers using erase-remove idiom.
  // The Identity_property_map property map can be omitted here as it is the default value.
  const int nb_neighbors = 24; // considers 24 nearest neighbor points

  // Estimate scale of the point set with average spacing
  const double average_spacing = CGAL::compute_average_spacing<CGAL::Sequential_tag>
    (points, nb_neighbors);

  //////////////////
  // FIRST OPTION //
  // I don't know the ratio of outliers present in the point set
  std::vector<Point>::iterator first_to_remove
    = CGAL::remove_outliers<CGAL::Parallel_if_available_tag>
    (points,
     nb_neighbors,
     CGAL::parameters::threshold_percent (100.). // No limit on the number of outliers to remove
     threshold_distance (2. * average_spacing)); // Point with distance above 2*average_spacing are considered outliers

  std::cerr << (100. * std::distance(first_to_remove, points.end()) / (double)(points.size()))
            << "% of the points are considered outliers when using a distance threshold of "
            << 2. * average_spacing << std::endl;


  ///////////////////
  // SECOND OPTION //
  // I know the ratio of outliers present in the point set
  const double removed_percentage = 5.0; // percentage of points to remove

  points.erase(CGAL::remove_outliers<CGAL::Parallel_if_available_tag>
               (points,
                nb_neighbors,
                CGAL::parameters::threshold_percent(removed_percentage). // Minimum percentage to remove
                threshold_distance(0.)), // No distance threshold (can be omitted)
               points.end());

  // Optional: after erase(), use Scott Meyer's "swap trick" to trim excess capacity
  std::vector<Point>(points).swap(points);

  return EXIT_SUCCESS;
}
