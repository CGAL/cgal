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


namespace std
{
  std::ostream& operator<< (std::ostream& stream, const Color& c)
  {
    stream << int(c[0]) << " " << int(c[1]) << " " << int(c[2]) << " " << int(c[3]);
    return stream;
  }
}


int main(int, char**)
{
  std::vector<PCI> points; // store points

  for (std::size_t i = 0; i < 10; ++ i)
    points.push_back (CGAL::cpp11::make_tuple (Point (i / 10., i / 20., i / 30.),
                                               CGAL::make_array ((unsigned char)(255 / (i + 1)),
                                                                 (unsigned char)(192 / (i + 1)),
                                                                 (unsigned char)(128 / (i + 1)),
                                                                 (unsigned char)(64 / (i + 1))),
                                               i));

  std::ofstream f("out.ply");
  if (CGAL::get_mode(f) == CGAL::IO::ASCII)
    std::cerr << "Okay!" << std::endl;
  
  CGAL::write_ply_points_with_properties
    (f, points.begin(), points.end(),
     CGAL::Ply::point_writer (Point_map()),
     CGAL::cpp11::make_tuple(Color_map(),
                             CGAL::Ply::Property<unsigned char>("red"),
                             CGAL::Ply::Property<unsigned char>("green"),
                             CGAL::Ply::Property<unsigned char>("blue"),
                             CGAL::Ply::Property<unsigned char>("alpha")),
     std::make_pair (Intensity_map(), CGAL::Ply::Property<int>("intensity")));
  
  return EXIT_SUCCESS;
}
