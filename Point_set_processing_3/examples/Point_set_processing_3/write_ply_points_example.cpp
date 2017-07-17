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
typedef CGAL::cpp11::array<unsigned char, 4> Color;

// Point with normal, color and intensity
typedef CGAL::cpp11::tuple<Point, Color, int> PCI;
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
      if (is_ascii(out))
        out << int(c[0]) << " " << int(c[1]) << " " << int(c[2]) << " " << int(c[3]);
      else
        out.write(reinterpret_cast<const char*>(&c), sizeof(c));
      return out;
    }
  }; 
}


int main(int, char**)
{
  std::vector<PCI> points; // store points

  for (int i = 0; i < 10; ++ i)
    points.push_back (CGAL::cpp11::make_tuple (Point (i / 10., i / 20., i / 30.),
                                               CGAL::make_array ((unsigned char)(255 / (i + 1)),
                                                                 (unsigned char)(192 / (i + 1)),
                                                                 (unsigned char)(128 / (i + 1)),
                                                                 (unsigned char)(64 / (i + 1))),
                                               i));

  std::ofstream f("out.ply");
  CGAL::set_binary_mode(f); // The PLY file will be written in the binary format
  
  CGAL::write_ply_points_with_properties
    (f, points.begin(), points.end(),
     CGAL::make_ply_point_writer (Point_map()),
     std::make_tuple(Color_map(),
                     CGAL::PLY_property<unsigned char>("red"),
                     CGAL::PLY_property<unsigned char>("green"),
                     CGAL::PLY_property<unsigned char>("blue"),
                     CGAL::PLY_property<unsigned char>("alpha")),
     std::make_pair (Intensity_map(), CGAL::PLY_property<int>("intensity")));
  
  return EXIT_SUCCESS;
}
