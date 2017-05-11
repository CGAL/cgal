#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/IO/read_xyz_points.h>
#include <CGAL/Point_with_normal_3.h>
#include <CGAL/property_map.h>

#include <CGAL/Point_set_3.h>
#include <CGAL/Point_set_3/Point_set_processing_3.h>
#include <CGAL/Point_set_3/IO.h>
#include <CGAL/Shape_detection_3.h>

#include <iostream>
#include <fstream>

// Type declarations
typedef CGAL::Exact_predicates_inexact_constructions_kernel  Kernel;
typedef CGAL::Point_set_3<typename Kernel::Point_3>          Point_set;
typedef typename Point_set::Point_map                        Point_map;
typedef typename Point_set::Vector_map                       Vector_map;

// In Shape_detection_traits the basic types, i.e., Point and Vector types
// as well as iterator type and property maps, are defined.
typedef CGAL::Shape_detection_3::Shape_detection_traits
  <Kernel, Point_set, Point_map, Vector_map>                 Traits;
typedef CGAL::Shape_detection_3::Region_growing<Traits>      Region_growing;
typedef CGAL::Shape_detection_3::Plane<Traits>               Plane;

int main(int argc, char** argv)
{
  Point_set points;

  std::ifstream stream(argc>1 ? argv[1] : "data/cube.pwn");
  stream >> points;

  if (points.empty())
  {
    std::cerr << "Error: cannot read file " << std::endl;
    return EXIT_FAILURE;
  }

  if (!(points.has_normal_map()))
  {
    std::cout << "Normal estimation" << std::endl;
    CGAL::jet_estimate_normals<CGAL::Sequential_tag> (points, 12);
  }

  std::cout << "Shape detection" << std::endl;
  Region_growing region_growing;
  region_growing.set_input(points, points.point_map(), points.normal_map());
  region_growing.add_shape_factory<Plane>();
  region_growing.detect();

  std::cout << region_growing.shapes().end() - region_growing.shapes().begin() << " shapes detected." << std::endl;

  return EXIT_SUCCESS;
}
