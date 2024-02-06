#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/estimate_lfs.h>
#include <CGAL/IO/read_points.h>

#include <vector>
#include <utility>
#include <tuple>

// types
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;

typedef std::tuple<Point, Vector, FT> Point_with_normal_and_lfs;

// Concurrency
typedef CGAL::Parallel_if_available_tag Concurrency_tag;

int main(void)
{

  // read xyz
  const std::string filename = CGAL::data_file_path("points_3/kitten.xyz");

  std::vector<Point_with_normal_and_lfs> points;
  if(!CGAL::IO::read_points(filename,
                            std::back_inserter(points),
                            CGAL::parameters::point_map(CGAL::Nth_of_tuple_property_map<0, Point_with_normal_and_lfs>())
                            .normal_map(CGAL::Nth_of_tuple_property_map<1, Point_with_normal_and_lfs>())))
  {
    std::cerr << "Error: cannot read file " << filename<< std::endl;
    return EXIT_FAILURE;
  }

  unsigned int jet_k = 24;
  std::size_t N_rays = 60;
  FT apex_angle = 30;
  
  auto lfs_map = CGAL::Nth_of_tuple_property_map<2, Point_with_normal_and_lfs>();
  CGAL::estimate_local_feature_size<Concurrency_tag>(points,
                                                     jet_k,
                                                     N_rays,
                                                     apex_angle,
                                                     lfs_map,
                                                     CGAL::parameters::point_map(CGAL::Nth_of_tuple_property_map<0, Point_with_normal_and_lfs>())
                                                                      .normal_map(CGAL::Nth_of_tuple_property_map<1, Point_with_normal_and_lfs>()));

  for (const auto &pts : points){
    std::cerr << std::get<2>(pts) << std::endl;
  }

  return EXIT_SUCCESS;
}
