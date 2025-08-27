#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/property_map.h>
#include <CGAL/IO/read_las_points.h>
#include <CGAL/IO/write_las_points.h>

#include <utility>
#include <vector>
#include <fstream>

// types
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point;
typedef std::array<unsigned short, 4> Color;
typedef std::pair<Point, Color> PointWithColor;

int main(int argc, char*argv[])
{
  const char* fname = "colored_points.las";

  std::ofstream os(fname, std::ios::binary);

  std::vector<PointWithColor> points; // store points
  points.push_back(std::make_pair(Point(0, 0, 0), Color{ 65535, 0, 0, 0 }));
  points.push_back(std::make_pair(Point(1, 0, 0), Color{ 0, 65535, 0, 0 }));
  points.push_back(std::make_pair(Point(0, 1, 0), Color{ 0, 0, 65535, 0 }));
  points.push_back(std::make_pair(Point(1, 1, 0), Color{ 0, 65535, 65535, 0 }));
  points.push_back(std::make_pair(Point(1, 1, 1), Color{ 65535, 65535, 0, 0 }));

  // Writes a .las point set file with colors
  if(!CGAL::IO::write_LAS_with_properties(os, points,
                                         CGAL::IO::make_las_point_writer(CGAL::First_of_pair_property_map<PointWithColor>()),
                                         std::make_tuple(CGAL::Second_of_pair_property_map<PointWithColor>(),
                                                         CGAL::IO::LAS_property::R(),
                                                         CGAL::IO::LAS_property::G(),
                                                         CGAL::IO::LAS_property::B(),
                                                         CGAL::IO::LAS_property::I())))
  {
    std::cerr << "Error: cannot write file " << fname << std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
