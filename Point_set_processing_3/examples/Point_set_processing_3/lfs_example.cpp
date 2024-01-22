#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/estimate_lfs.h>
#include <CGAL/IO/read_points.h>

#include <vector>
#include <utility> // defines std::pair


// types
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;

// Point with normal vector stored in a std::pair.
typedef std::pair<Point, Vector> Point_with_normal;

// Concurrency
typedef CGAL::Parallel_if_available_tag Concurrency_tag;


int main(void)
{

  // read xyz
  const std::string filename = "./frog.xyz";

  std::vector<Point_with_normal> points;
  if(!CGAL::IO::read_points(filename,
                            std::back_inserter(points),
                            CGAL::parameters::point_map(CGAL::First_of_pair_property_map<Point_with_normal>())))
  {
    std::cerr << "Error: cannot read file " << filename<< std::endl;
    return EXIT_FAILURE;
  }
  
  unsigned int jet_k = 24;
  std::size_t N_rays = 60;
  FT apex_angle = 30;
  std::vector<FT> lfses = CGAL::estimate_local_feature_size<Concurrency_tag>(points, jet_k, N_rays, apex_angle,
    CGAL::parameters::point_map(CGAL::First_of_pair_property_map<Point_with_normal>())
                                                .normal_map(CGAL::Second_of_pair_property_map<Point_with_normal>()));

  for (const auto &lfs : lfses)
    std::cerr << lfs << "\n";


  return EXIT_SUCCESS;
}
