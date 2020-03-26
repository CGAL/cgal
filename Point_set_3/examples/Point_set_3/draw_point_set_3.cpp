#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Point_set_3.h>
#include <CGAL/Point_set_3/IO.h>
#include <CGAL/draw_point_set_3.h>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;
typedef CGAL::Point_set_3<Point> Point_set;

int main (int argc, char** argv)
{
  std::ifstream f (argc > 1 ? argv[1] : "data/oni.xyz");
  Point_set point_set;

  // Reading input in XYZ format
  if (!f || !CGAL::read_xyz_point_set (f, point_set))
  {
    std::cerr<<"Can't read input file "
             <<(argc > 1 ? argv[1] : "data/oni.xyz")<< std::endl;
    return EXIT_FAILURE;
  }

  CGAL::draw(point_set);

  return EXIT_SUCCESS;
}
