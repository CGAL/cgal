#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Point_set_3.h>
#include <CGAL/grid_simplify_point_set.h>
#include <CGAL/point_generators_3.h>

#include <boost/unordered_set.hpp>

#include <fstream>
#include <limits>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;

typedef CGAL::Point_set_3<Point> Point_set;

int main() {

  std::ifstream stream_xyz{CGAL::data_file_path("points_3/cube.xyz")};
  Point_set points_xyz;
  stream_xyz >> points_xyz;
  assert(points_xyz.number_of_points() == 8);
  assert(!points_xyz.has_normal_map());

  std::ifstream stream_pwn{CGAL::data_file_path("points_3/cube.pwn")};
  Point_set points_pwn;
  stream_pwn >> points_pwn;
  assert(points_pwn.number_of_points() == 50000);
  assert(points_pwn.has_normal_map());

}
