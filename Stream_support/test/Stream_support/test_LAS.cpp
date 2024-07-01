#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/property_map.h>
#include <CGAL/IO/read_las_points.h>

#include <utility>
#include <vector>
#include <fstream>

// types
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point;
typedef std::array<unsigned short, 4> Color;
typedef std::pair<Point, Color> PointWithColor;

int main(int argc, char* argv[])
{
  std::string fname = (argc > 1) ? argv[1] : CGAL::data_file_path("points_3/colors.las");

  // Reads a .las point set file with normal vectors and colors
  std::vector<PointWithColor> points; // store points
  if (!CGAL::IO::read_LAS_with_properties(fname.c_str(), std::back_inserter(points),
    CGAL::IO::make_las_point_reader(CGAL::First_of_pair_property_map<PointWithColor>()),
    std::make_tuple(CGAL::Second_of_pair_property_map<PointWithColor>(),
      CGAL::Construct_array(),
      CGAL::IO::LAS_property::R(),
      CGAL::IO::LAS_property::G(),
      CGAL::IO::LAS_property::B(),
      CGAL::IO::LAS_property::I())))
  {
    std::cerr << "Error: cannot read file " << fname << std::endl;
    return EXIT_FAILURE;
  }

  CGAL_assertion(points.size() == 3);
  CGAL_assertion(points[0].second[0] == 255);
  CGAL_assertion(points[0].second[1] == 0);
  CGAL_assertion(points[0].second[2] == 0);

  CGAL_assertion(points[1].second[0] == 0);
  CGAL_assertion(points[1].second[1] == 255);
  CGAL_assertion(points[1].second[2] == 0);

  return EXIT_SUCCESS;
}
