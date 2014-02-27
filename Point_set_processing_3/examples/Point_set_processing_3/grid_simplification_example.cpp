#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/grid_simplify_point_set.h>
#include <CGAL/random_simplify_point_set.h>
#include <CGAL/IO/read_xyz_points.h>
#include <CGAL/IO/write_xyz_points.h>
#include <CGAL/Timer.h>

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

  CGAL::Timer task_timer;
  task_timer.start();
  std::cout << "Run algorithm example: " << std::endl;

  // simplification by clustering using erase-remove idiom
  double cell_size = 0.12;
  points.erase(CGAL::grid_simplify_point_set(points.begin(), points.end(), cell_size),
               points.end());
  //points.erase(CGAL::random_simplify_point_set(points.begin(), points.end(), 90),
  //  points.end());

  std::cout << "total done: " << task_timer.time() << " seconds, "  << std::endl;

  // Optional: after erase(), use Scott Meyer's "swap trick" to trim excess capacity
  std::vector<Point>(points).swap(points);

  // Saves point set.
  // Note: write_xyz_points_and_normals() requires an output iterator
  // over points as well as property maps to access each
  // point position and normal.
  std::ofstream out("data/sphere_20k_grid_downsampled.xyz");
  if (!out ||
	  !CGAL::write_xyz_points(
	  out, points.begin(), points.end()))
  {
	  return EXIT_FAILURE;
  }

  system("Pause");

  return EXIT_SUCCESS;
}

