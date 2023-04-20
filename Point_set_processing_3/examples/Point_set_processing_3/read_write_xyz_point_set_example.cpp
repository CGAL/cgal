#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/property_map.h>
#include <CGAL/IO/read_xyz_points.h>
#include <CGAL/IO/write_xyz_points.h>

#include <utility> // defines std::pair
#include <vector>
#include <fstream>
#include <iostream>

// types
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;

// Point with normal vector stored as a std::pair.
typedef std::pair<Point, Vector> Pwn;

int main(int argc, char*argv[])
{
  const std::string fname = (argc>1) ? argv[1] : CGAL::data_file_path("points_3/oni.pwn");

  // Reads a .xyz point set file in points[].
  // Note: read_points() requires an output iterator
  // over points and as well as property maps to access each
  // point position and normal.
  std::vector<Pwn> points;
  if(!CGAL::IO::read_XYZ(fname,
                         std::back_inserter(points),
                         CGAL::parameters::point_map(CGAL::First_of_pair_property_map<Pwn>())
                                          .normal_map(CGAL::Second_of_pair_property_map<Pwn>())))
  {
    std::cerr << "Error: cannot read file " << fname << std::endl;
    return EXIT_FAILURE;
  }

  // Saves point set.
  // Note: write_XYZ() requires property maps to access each
  // point position and normal.
  if(!CGAL::IO::write_XYZ("oni_copy.xyz", points,
                          CGAL::parameters::point_map(CGAL::First_of_pair_property_map<Pwn>())
                                           .normal_map(CGAL::Second_of_pair_property_map<Pwn>())
                                           .stream_precision(17)))
    return EXIT_FAILURE;

  return EXIT_SUCCESS;
}
