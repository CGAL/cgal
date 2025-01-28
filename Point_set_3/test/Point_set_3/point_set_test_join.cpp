#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Point_set_3.h>
#include <CGAL/IO/write_xyz_points.h>
#include <CGAL/grid_simplify_point_set.h>

#include <fstream>
#include <limits>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;

typedef CGAL::Point_set_3<Point> Point_set;
typedef std::array<unsigned char, 3> Color;

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
  std::optional<Point_set::Property_map<int>> intensity
    = ps.property_map<int>("intensity");

  std::cerr << msg << std::endl;
  for (Point_set::const_iterator it = ps.begin(); it != ps.end(); ++ it)
  {
    std::cerr << *it << ": " << ps.point(*it);
    if (ps.has_normal_map())
      std::cerr << ", normal " << ps.normal(*it);
    if (intensity.has_value())
      std::cerr << ", intensity " << intensity.value()[*it];
    std::cerr << std::endl;
  }
}

int main (int, char**)
{
  Point_set ps1, ps2;
  ps1.add_normal_map();

  for (std::size_t i = 0; i < 5; ++ i)
    ps1.insert (Point (double(i), double(i), double(i)), Vector (double(i), double(i), double(i)));

  ps1.remove (ps1.end() - 3);

  for (std::size_t i = 5; i < 10; ++ i)
    ps2.insert (Point (double(i), double(i), double(i)));

  ps2.remove (ps2.end() - 3);

  print_point_set (ps1, "PS1 = ");
  print_point_set (ps2, "PS2 = ");

  ps1 += ps2;
  print_point_set (ps1, "JOINT PS1 = ");

  Point_set ps3;
  ps3.add_normal_map();

  Point_set::Property_map<int> intensity;
  bool okay;

  boost::tie (intensity, okay) = ps3.add_property_map<int>("intensity", 0);
  assert (okay);

  Point_set::iterator it = ps3.insert (Point (double(0), double(1), double(2)),
                                       Vector (double(3), double(4), double(5)));
  intensity[*it] = 42;

  print_point_set (ps3, "PS3 = ");
  ps1.copy_properties (ps3);

  print_point_set (ps1, "PS1 with PS3 properties = ");
  ps1.insert (ps3, *it);
  print_point_set (ps1, "PS1 with PS3 properties + PS3 item copied = ");

  return EXIT_SUCCESS;
}
