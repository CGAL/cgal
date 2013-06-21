#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/regularize_and_simplify_point_set.h>
#include <CGAL/IO/read_xyz_points.h>
#include <CGAL/IO/write_xyz_points.h>

#include <vector>
#include <fstream>

// types
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point;

int main(void)
{
  // Reads a .xyz point set file in points[].
  std::vector<Point> points;
  std::ifstream stream("data/sphere_20k.xyz");
  if (!stream ||
      !CGAL::read_xyz_points(stream, std::back_inserter(points)))
  {
    std::cerr << "Error: cannot read file data/oni.xyz" << std::endl;
    return EXIT_FAILURE;
  }

  // Randomly simplifies using erase-remove idiom
  const double retain_percentage = 5.0; // percentage of points to remove

  std::vector<Point> points_sampled;
  points_sampled.assign(points.size() * (retain_percentage / 100.), Point());
  
  std::copy(CGAL::regularize_and_simplify_point_set(points.begin(), points.end(), retain_percentage, 450, 50), points.end(), points_sampled.begin());

  // Optional: after erase(), use Scott Meyer's "swap trick" to trim excess capacity
  std::vector<Point>(points).swap(points);

  // Saves point set.
  // Note: write_xyz_points_and_normals() requires an output iterator
  // over points as well as property maps to access each
  // point position and normal.
  std::ofstream out("data/sphere_20k_sampled.xyz");
  if (!out ||
	  !CGAL::write_xyz_points(
	  out, points_sampled.begin(), points_sampled.end()))
  {
	  return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}

