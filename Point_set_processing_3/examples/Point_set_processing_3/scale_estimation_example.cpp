#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/estimate_scale.h>
#include <CGAL/IO/read_xyz_points.h>

#include <vector>
#include <fstream>
#include <boost/tuple/tuple.hpp>

// Types
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point_3;
typedef Kernel::Point_2 Point_2;

int main (int argc, char** argv)
{

  // 3D CASE
  {
    const char* fname = (argc>1)?argv[1]:"data/sphere_20k.xyz";
    // Reads a .xyz point set file in points.
    // As the point is the second element of the tuple (that is with index 1)
    // we use a property map that accesses the 1st element of the tuple.
  
    std::vector<Point_3> points;
    std::ifstream stream(fname);
    if (!stream ||
        !CGAL::read_xyz_points(
                               stream, std::back_inserter(points)))
      {
        std::cerr << "Error: cannot read file " << fname << std::endl;
        return EXIT_FAILURE;
      }

    // estimate global scale
    std::size_t k_scale = CGAL::estimate_global_k_neighbor_scale (points.begin(), points.end(),
                                                                  CGAL::Identity_property_map<Point_3>(),
                                                                  Kernel());
    std::cout << "Global K scale: " << k_scale << std::endl;

    FT range_scale = CGAL::estimate_global_range_scale (points.begin(), points.end(),
                                                        CGAL::Identity_property_map<Point_3>(),
                                                        Kernel());
    std::cout << "Global range scale: " << range_scale << std::endl;
  }

  // 2D CASE
  {
    std::vector<Point_2> points;
    for (std::size_t i = 0; i < 5000; ++ i)
      {
        double angle = CGAL_PI * 2. * (rand() / (double)RAND_MAX);
        double radius = 1. + 0.1 * (rand() / (double)RAND_MAX);
        points.push_back (Point_2 (radius * std::cos (angle),
                                   radius * std::sin (angle)));
      }

        
    // estimate global scale
    std::size_t k_scale = CGAL::estimate_global_k_neighbor_scale (points.begin(), points.end(),
                                                                  CGAL::Identity_property_map<Point_2>(),
                                                                  Kernel());
    std::cout << "Global K scale: " << k_scale << std::endl;

    FT range_scale = CGAL::estimate_global_range_scale (points.begin(), points.end(),
                                                        CGAL::Identity_property_map<Point_2>(),
                                                        Kernel());
    std::cout << "Global range scale: " << range_scale << std::endl;

  }
  return EXIT_SUCCESS;
}

