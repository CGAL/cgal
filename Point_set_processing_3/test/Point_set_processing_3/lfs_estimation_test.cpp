#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/estimate_local_feature_size.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Random.h>

#include <vector>
#include <utility>
#include <tuple>
#include <math.h>

// types
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;
typedef CGAL::Random_points_on_sphere_3<Point> Random_points_on_sphere_3;

typedef std::tuple<Point, Vector, FT> Point_with_normal_and_lfs;


// Concurrency
typedef CGAL::Parallel_if_available_tag Concurrency_tag;

void test (std::vector<Point>& input,
           std::ptrdiff_t result0 = 1, int result1 = 1, int result2 = 1, int result3 = 1, int result4 = 1)
{
  ;
}


int main(void)
{
  const int num_points = 1000; // number of sampled points
  FT r = 1.0; // unit sphere

  CGAL::Random rand = CGAL::Random(23);

  // std::vector<Point> point_coordinates;
  // point_coordinates.reserve(num_points);
  // std::copy_n(Random_points_on_sphere_3(r, rand), num_points, std::back_inserter(point_coordinates));

  std::vector<Point_with_normal_and_lfs> points;
  std::generate_n(std::back_inserter(points),
                    num_points,
                    [&]() -> Point_with_normal_and_lfs {
                        Point point = *CGAL::Random_points_on_sphere_3<Point>(r, rand);
                        Vector normal = point - CGAL::ORIGIN;
                        return std::make_tuple(point, normal, FT(0));
                    });


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

  // evalute the LFS on the sphere
  // the GT should be the radius
  FT mse = 0.0;
  for (const auto &pts : points)
  {
    FT lfs = std::get<2>(pts);
    mse += std::pow((lfs - r), 2);
  }
  mse /= num_points;

  std::cout << "Mean Squared Error (MSE): " << mse << std::endl;

  return EXIT_SUCCESS;
}

