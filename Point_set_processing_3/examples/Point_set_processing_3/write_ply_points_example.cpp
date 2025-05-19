#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/property_map.h>
#include <CGAL/IO/write_ply_points.h>

#include <utility>
#include <vector>
#include <fstream>

// types
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;
typedef std::array<unsigned char, 4> Color;

// Point with normal, color and intensity
typedef std::tuple<Point, Color, int> PCI;
typedef CGAL::Nth_of_tuple_property_map<0, PCI> Point_map;
typedef CGAL::Nth_of_tuple_property_map<1, PCI> Color_map;
typedef CGAL::Nth_of_tuple_property_map<2, PCI> Intensity_map;

// Define how a color should be stored
namespace CGAL {

template< class F >
struct Output_rep< ::Color, F > {
  const ::Color& c;
  static const bool is_specialized = true;
  Output_rep (const ::Color& c) : c(c)
  { }
  std::ostream& operator() (std::ostream& out) const
  {
    if (IO::is_ascii(out))
      out << int(c[0]) << " " << int(c[1]) << " " << int(c[2]) << " " << int(c[3]);
    else
      out.write(reinterpret_cast<const char*>(&c), sizeof(c));
    return out;
  }
};

} // namespace CGAL

int main(int, char**)
{
  std::vector<PCI> points; // store points

  for (int i = 0; i < 10; ++ i)
    points.push_back (std::make_tuple (Point (i / 10., i / 20., i / 30.),
                                               CGAL::make_array (static_cast<unsigned char>(255 / (i + 1)),
                                                                 static_cast<unsigned char>(192 / (i + 1)),
                                                                 static_cast<unsigned char>(128 / (i + 1)),
                                                                 static_cast<unsigned char>(64 / (i + 1))),
                                               i));

  std::ofstream f("out.ply", std::ios::binary);
  CGAL::IO::set_binary_mode(f); // The PLY file will be written in the binary format

  CGAL::IO::write_PLY_with_properties(f, points,
                                      CGAL::make_ply_point_writer (Point_map()),
                                      std::make_tuple(Color_map(),
                                                      CGAL::IO::PLY_property<unsigned char>("red"),
                                                      CGAL::IO::PLY_property<unsigned char>("green"),
                                                      CGAL::IO::PLY_property<unsigned char>("blue"),
                                                      CGAL::IO::PLY_property<unsigned char>("alpha")),
                                  std::make_pair(Intensity_map(), CGAL::IO::PLY_property<int>("intensity")));

  return EXIT_SUCCESS;
}
