#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/estimate_scale.h>
#include <CGAL/IO/read_xyz_points.h>

#include <vector>
#include <fstream>
#include <boost/tuple/tuple.hpp>

// Types
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point;

int main (int argc, char** argv)
{
  const char* fname = (argc>1)?argv[1]:"data/sphere_20k.xyz";
  // Reads a .xyz point set file in points.
  // As the point is the second element of the tuple (that is with index 1)
  // we use a property map that accesses the 1st element of the tuple.
    
  std::vector<Point> points;
  std::ifstream stream(fname);
  if (!stream ||
      !CGAL::read_xyz_points(
                             stream, std::back_inserter(points)))
    {
      std::cerr << "Error: cannot read file " << fname << std::endl;
      return EXIT_FAILURE;
    }

  // estimate global scale
  std::size_t scale = CGAL::estimate_global_k_neighbor_scale (points.begin(), points.end(),
                                                              CGAL::Identity_property_map<Point>(),
                                                              Kernel());
  std::cout << "Global K scale: " << scale << std::endl;

  return EXIT_SUCCESS;
}

