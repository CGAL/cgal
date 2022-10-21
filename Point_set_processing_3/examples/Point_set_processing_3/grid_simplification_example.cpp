#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/grid_simplify_point_set.h>
#include <CGAL/IO/read_points.h>

#include <vector>
#include <fstream>

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

  // simplification by clustering using erase-remove idiom
  double cell_size = 0.03;
  unsigned int min_points_per_cell = 3;

  auto iterator_to_first_to_remove
    = CGAL::grid_simplify_point_set
    (points, cell_size,
     CGAL::parameters::min_points_per_cell(min_points_per_cell)); // optional

  points.erase(iterator_to_first_to_remove, points.end());

  // Optional: after erase(), shrink_to_fit to trim excess capacity
  points.shrink_to_fit();

  return EXIT_SUCCESS;
}
