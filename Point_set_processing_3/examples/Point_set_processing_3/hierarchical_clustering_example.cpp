#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/hierarchical_clustering.h>
#include <CGAL/IO/read_xyz_points.h>

#include <vector>
#include <fstream>

// types
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point;

int main(int argc, char*argv[])
{
  // Reads a .xyz point set file in points[].
  std::vector<Point> points;
  const char* fname = (argc>1)?argv[1]:"data/oni.xyz";
  std::ifstream stream(fname);
  if (!stream ||
      !CGAL::read_xyz_points(stream, std::back_inserter(points)))
  {
    std::cerr << "Error: cannot read file " << fname << std::endl;
    return EXIT_FAILURE;
  }

  std::vector<Point> output; // Algorithm generate a new set of points
  CGAL::hierarchical_clustering (points.begin (), points.end (),
				 std::back_inserter (output));

  return EXIT_SUCCESS;
}

