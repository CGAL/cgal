#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/random_simplify_point_set.h>
#include <CGAL/IO/read_points.h>
#include <CGAL/IO/write_points.h>

#include <vector>
#include <fstream>
#include <iostream>

// types
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point;

int main(int argc, char*argv[])
{
  const std::string fname = (argc>1) ? argv[1] : CGAL::data_file_path("points_3/oni.pwn");

  // Reads a point set file in points[].
  std::vector<Point> points;
  if(!CGAL::IO::read_points(fname, std::back_inserter(points)))
  {
    std::cerr << "Error: cannot read file " << fname << std::endl;
    return EXIT_FAILURE;
  }

  // Randomly simplifies using erase-remove idiom
  const double removed_percentage = 97.0; // percentage of points to remove
  points.erase(CGAL::random_simplify_point_set(points, removed_percentage), points.end());

  // Optional: after erase(), use Scott Meyer's "swap trick" to trim excess capacity
  std::vector<Point>(points).swap(points);

  // Saves point set.
  const std::string output_filename = (argc>2) ? argv[2] : "Three_lady_copy.xyz";
  if(!CGAL::IO::write_points(output_filename, points, CGAL::parameters::stream_precision(17)))
    return EXIT_FAILURE;

  return EXIT_SUCCESS;
}
