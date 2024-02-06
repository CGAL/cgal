#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/estimate_lfs.h>
#include <CGAL/IO/read_points.h>
#include <CGAL/Point_set_3.h>
#include <CGAL/Point_set_3/IO.h>

#include <vector>



// types
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;

typedef CGAL::Point_set_3<Point, Vector> Point_set;
typedef Point_set::Property_map<Point> Point_map;
typedef Point_set::Property_map<FT> FT_map;

// types for K nearest neighbors search structure
typedef CGAL::Point_set_processing_3::internal::Neighbor_query<Kernel, Point_set&, Point_map> Neighbor_query;

// Concurrency
typedef CGAL::Parallel_if_available_tag Concurrency_tag;
//typedef CGAL::Parallel_tag Concurrency_tag;
typedef CGAL::Sequential_tag Concurrency_tag;

int main(void)
{
  // read xyz
  const std::string filename = CGAL::data_file_path("points_3/kitten.xyz");

  Point_set point_set;
  point_set.add_normal_map();

  FT_map lfs_map;
  boost::tie (lfs_map, boost::tuples::ignore) = point_set.add_property_map<FT> ("LFS", 0.);

  if (!CGAL::IO::read_point_set(filename, point_set))
  {
    std::cerr << "Error: cannot read file " << filename<< std::endl;
    return EXIT_FAILURE;
  }


  unsigned int jet_k = 24;
  std::size_t N_rays = 60;
  FT apex_angle = 30;
  CGAL::estimate_local_feature_size<Concurrency_tag>(point_set, jet_k, N_rays, apex_angle, lfs_map,
    CGAL::parameters::point_map(point_set.point_map())
                                                .normal_map(point_set.normal_map()));
  // print
  for (Point_set::iterator it = point_set.begin(); it != point_set.end(); ++ it)
   {
      Point p = point_set.point(*it);
      Vector n = point_set.normal(*it);
      FT lfs = lfs_map[*it];
      std::cerr << "point:" << p << " normal:" << n << " lfs:" << lfs << std::endl;
   }

  std::cerr << "this is a test!" << std::endl;

  return EXIT_SUCCESS;
}
