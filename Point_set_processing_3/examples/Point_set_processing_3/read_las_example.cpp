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
typedef CGAL::cpp11::array<unsigned short, 4> Color;
typedef std::pair<Point, Color> PointWithColor;

int main(int argc, char*argv[])
{
  const char* fname = (argc>1) ? argv[1] : "data/pig_points.las";
  // Reads a .las point set file with normal vectors and colors

  std::vector<PointWithColor> points; // store points
  std::ifstream in(fname, std::ios_base::binary);

  if (!in ||
      !CGAL::read_las_points_with_properties
      (in,
       std::back_inserter (points),
       CGAL::make_las_point_reader (CGAL::First_of_pair_property_map<PointWithColor>()),
       std::make_tuple
       (CGAL::Second_of_pair_property_map<PointWithColor>(),
        CGAL::Construct_array(),
        CGAL::LAS_property::R(),
        CGAL::LAS_property::G(),
        CGAL::LAS_property::B(),
        CGAL::LAS_property::I())))
    {
      std::cerr << "Error: cannot read file " << fname << std::endl;
      return EXIT_FAILURE;
    }

  for (std::size_t i = 0; i < points.size(); ++ i)
    std::cout << points[i].first << std::endl;
  
  return EXIT_SUCCESS;
}
