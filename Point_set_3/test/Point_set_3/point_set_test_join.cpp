#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Point_set_3.h>
#include <CGAL/IO/read_xyz_points.h>
#include <CGAL/IO/write_xyz_points.h>
#include <CGAL/grid_simplify_point_set.h>

#include <fstream>
#include <limits>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;

typedef CGAL::Point_set_3<Point> Point_set;
typedef CGAL::cpp11::array<unsigned char, 3> Color;

std::size_t nb_test = 0;
std::size_t nb_success = 0;

void test (bool expr, const char* msg)
{
  ++ nb_test;
  if (!expr)
    std::cerr << "Error on test " << nb_test << ": " << msg << std::endl;
  else
    ++ nb_success;
}

void print_point_set (const Point_set& ps, const char* msg)
{
  std::cerr << msg << std::endl;
  if (ps.has_normal_map())
    for (Point_set::const_iterator it = ps.begin(); it != ps.end(); ++ it)
      std::cerr << *it << ": " << ps.point(*it) 
                << ", normal " << ps.normal(*it) << std::endl;
  else
    for (Point_set::const_iterator it = ps.begin(); it != ps.end(); ++ it)
      std::cerr << *it << ": " << ps.point(*it) << std::endl;
}


int main (int, char**)
{
  Point_set ps1, ps2;
  ps1.add_normal_map();

  for (std::size_t i = 0; i < 5; ++ i)
    ps1.insert (Point (i, i, i), Vector (i, i, i));

  ps1.remove (ps1.end() - 3);
  
  for (std::size_t i = 5; i < 10; ++ i)
    ps2.insert (Point (i, i, i));

  ps2.remove (ps2.end() - 3);

  print_point_set (ps1, "PS1 = ");
  print_point_set (ps2, "PS2 = ");

  ps1 += ps2;
  print_point_set (ps1, "JOINT PS1 = ");

  return EXIT_SUCCESS;
};
