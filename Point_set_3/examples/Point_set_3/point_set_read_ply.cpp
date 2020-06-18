#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Point_set_3.h>
#include <CGAL/Point_set_3/IO.h>

#include <fstream>
#include <limits>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;

typedef CGAL::Point_set_3<Point> Point_set;

int main (int argc, char** argv)
{
  std::ifstream f (argc > 1 ? argv[1] : "data/example.ply",
                   std::ios_base::binary); // Mandatory on Windows if input is binary PLY

  Point_set point_set;

  if (!f || !CGAL::read_ply_point_set (f, point_set)) // same as `f >> point_set`
  {
    std::cerr << "Can't read input file " << std::endl;
  }

  // Shows which properties are defined
  std::vector<std::string> properties = point_set.properties();
  std::cerr << "Properties:" << std::endl;
  for (std::size_t i = 0; i < properties.size(); ++ i)
    std::cerr << " * " << properties[i] << std::endl;

  // Recover "label" property of type int
  Point_set::Property_map<boost::int32_t> label_prop;
  bool found = false;
  boost::tie (label_prop, found)  = point_set.property_map<boost::int32_t> ("label");

  if (found)
    {
      std::cerr << "Point set has an integer \"label\" property with values:" << std::endl;
      for (Point_set::iterator it = point_set.begin(); it != point_set.end(); ++ it)
        std::cerr << " * " << label_prop[*it] << std::endl;
    }

  if (argc > 2 && strcmp (argv[2], "-b") == 0) // Optional binary output
  {
    std::ofstream out ("out.ply",
                       std::ios_base::binary); // Mandatory on Windows
    CGAL::set_binary_mode (out); // Select binary mode (ASCII is default)
    out << point_set; // same as `CGAL::write_ply_point_set (out, point_set)`
  }
  else // ASCII output
  {
    std::ofstream out ("out.ply");
    out.precision(17); // Use sufficient precision in ASCII
    CGAL::write_ply_point_set (out, point_set); // same as `out << point_set`
  }

  return 0;
}
