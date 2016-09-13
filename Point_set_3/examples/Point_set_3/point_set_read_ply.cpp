#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Point_set_3.h>
#include <CGAL/IO/read_ply_points.h>
#include <CGAL/IO/read_ply_point_set_3.h>

#include <fstream>
#include <limits>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;

typedef CGAL::Point_set_3<Kernel> Point_set;
typedef CGAL::Ply_interpreter_point_set_3<Kernel> Ply_interpreter;

int main (int argc, char** argv)
{
  std::ifstream f (argc > 1 ? argv[1] : "data/TODO.ply");

  Point_set point_set;
  Ply_interpreter interpreter(point_set);
  
  if (!f ||
      !CGAL::read_ply_custom_points (f, interpreter, Kernel()))
    {
      std::cerr << "Can't read input file " << std::endl;
    }

  std::cerr << point_set.info();
  
  return 0;
}
