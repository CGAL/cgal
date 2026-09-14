
//#define CGAL_FORCE_IFORMAT_DOUBLE

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Timer.h>
#include <CGAL/Point_set_3.h>

#include <iostream>
#include <string>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;

typedef CGAL::Point_set_3<Point> Point_set;

int main (int argc, char** argv)
{
  const std::string fname = argc > 1 ? argv[1] : CGAL::data_file_path("points_3/oniAscii.ply");

  CGAL::Timer timer;
  timer.start();
  // Reading input
  Point_set point_set;
  std::string comment;
  if(!CGAL::IO::read_PLY(fname, point_set, comment, CGAL::parameters::use_binary_mode(false)))
  {
    std::cerr << "Can't read input file " << std::endl;
    return EXIT_FAILURE;
  }


  std::cout << "Read " << point_set.size() << " points in " << timer.time() << " seconds.\n";

  return 0;
}
