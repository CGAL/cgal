#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/IO/read_xyz_points.h>
#include <CGAL/Point_with_normal_3.h>
#include <CGAL/property_map.h>

#include <CGAL/Shape_detection_3.h>

#include <iostream>
#include <fstream>

// Type declarations
typedef CGAL::Exact_predicates_inexact_constructions_kernel  Kernel;
typedef std::pair<Kernel::Point_3, Kernel::Vector_3>         Point_with_normal;
typedef std::vector<Point_with_normal>                       Pwn_vector;
typedef CGAL::First_of_pair_property_map<Point_with_normal>  Point_map;
typedef CGAL::Second_of_pair_property_map<Point_with_normal> Normal_map;

// In Shape_detection_traits the basic types, i.e., Point and Vector types
// as well as iterator type and property maps, are defined.
typedef CGAL::Shape_detection_3::Shape_detection_traits
  <Kernel, Pwn_vector, Point_map, Normal_map>                Traits;
typedef CGAL::Shape_detection_3::Region_growing<Traits>      Region_growing;
typedef CGAL::Shape_detection_3::Plane<Traits>               Plane;

int main(int argc, char** argv)
{
  Pwn_vector points;
  std::ifstream stream(argc>1 ? argv[1] : "data/cube.pwn");

  if (!stream || 
    !CGAL::read_xyz_points_and_normals(stream,
      std::back_inserter(points),
      Point_map(),
      Normal_map()))
  {
      std::cerr << "Error: cannot read file cube.pwn" << std::endl;
      return EXIT_FAILURE;
  }

  Region_growing region_growing;
  region_growing.set_input(points);
  region_growing.add_shape_factory<Plane>();
  region_growing.detect();

  std::cout << region_growing.shapes().end() - region_growing.shapes().begin() << " shapes detected." << std::endl;

  return EXIT_SUCCESS;
}
