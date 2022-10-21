#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Point_set_3.h>
#include <CGAL/IO/read_points.h>

#include <fstream>
#include <limits>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;

typedef CGAL::Point_set_3<Point> Point_set;

int main (int argc, char** argv)
{
  const std::string filename = argc > 1 ? argv[1] : CGAL::data_file_path("meshes/camel.off");

  Point_set point_set;
  point_set.add_normal_map();

  // Reading input in OFF format
  if(!CGAL::IO::read_points(filename, point_set.index_back_inserter(),
                        CGAL::parameters::point_map(point_set.point_push_map())
                                         .normal_map(point_set.normal_push_map())))
  {
    std::cerr << "Can't read input file " << std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
