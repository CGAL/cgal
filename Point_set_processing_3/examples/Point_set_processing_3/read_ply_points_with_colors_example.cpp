#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/property_map.h>
#include <CGAL/IO/read_ply_points.h>

#include <utility>
#include <vector>
#include <fstream>

// types
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;
typedef std::array<unsigned char, 3> Color;

// Point with normal, color and intensity
typedef std::tuple<Point, Vector, Color, int> PNCI;
typedef CGAL::Nth_of_tuple_property_map<0, PNCI> Point_map;
typedef CGAL::Nth_of_tuple_property_map<1, PNCI> Normal_map;
typedef CGAL::Nth_of_tuple_property_map<2, PNCI> Color_map;
typedef CGAL::Nth_of_tuple_property_map<3, PNCI> Intensity_map;

int main(int argc, char*argv[])
{
  const std::string fname = (argc>1) ? argv[1] : CGAL::data_file_path("points_3/colors.ply");

  // Reads a .ply point set file with normal vectors and colors
  std::vector<PNCI> points; // store points
  std::ifstream in(fname);
  if(!CGAL::IO::read_PLY_with_properties(in, std::back_inserter(points),
                                         CGAL::make_ply_point_reader(Point_map()),
                                         std::make_pair(Intensity_map(), CGAL::PLY_property<int>("intensity")),
                                         std::make_tuple(Color_map(),
                                                         CGAL::Construct_array(),
                                                         CGAL::IO::PLY_property<unsigned char>("red"),
                                                         CGAL::IO::PLY_property<unsigned char>("green"),
                                                         CGAL::IO::PLY_property<unsigned char>("blue")),
                                         CGAL::IO::make_ply_normal_reader(Normal_map())))
  {
    std::cerr << "Error: cannot read file " << fname << std::endl;
    return EXIT_FAILURE;
  }

  // Display points read
  for(std::size_t i = 0; i < points.size (); ++ i)
  {
    const Point& p = get<0>(points[i]);
    const Vector& n = get<1>(points[i]);
    const Color& c = get<2>(points[i]);
    int I = get<3>(points[i]);
    std::cerr << "Point (" << p
              << ") with normal (" << n
              << "), color (" << int(c[0]) << " " << int(c[1]) << " " << int(c[2])
              << ") and intensity " << I << std::endl;
  }

  return EXIT_SUCCESS;
}
