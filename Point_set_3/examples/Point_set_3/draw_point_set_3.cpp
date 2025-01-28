#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Point_set_3.h>
#include <CGAL/draw_point_set_3.h>

#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;
typedef CGAL::Point_set_3<Point> Point_set;

int main (int argc, char** argv)
{
  const std::string filename = argc > 1 ? argv[1] : CGAL::data_file_path("points_3/oni.pwn");

  Point_set point_set;
  if(!CGAL::IO::read_point_set(filename, point_set))
  {
    std::cerr << "Can't read input file " << filename << std::endl;
    return EXIT_FAILURE;
  }

  CGAL::draw(point_set);

  return EXIT_SUCCESS;
}
