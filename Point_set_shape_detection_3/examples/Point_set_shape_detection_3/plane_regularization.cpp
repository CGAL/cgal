#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/IO/read_xyz_points.h>
#include <CGAL/Point_with_normal_3.h>
#include <CGAL/property_map.h>

#include <CGAL/Shape_detection_3.h>
#include <CGAL/regularize_planes.h>

#include <iostream>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel  Kernel;
typedef std::pair<Kernel::Point_3, Kernel::Vector_3>         Point_with_normal;
typedef std::vector<Point_with_normal>                       Pwn_vector;
typedef CGAL::First_of_pair_property_map<Point_with_normal>  Point_map;
typedef CGAL::Second_of_pair_property_map<Point_with_normal> Normal_map;

typedef CGAL::Shape_detection_3::Shape_detection_traits
  <Kernel, Pwn_vector, Point_map, Normal_map>                Traits;
typedef CGAL::Shape_detection_3::Region_growing<Traits>      Region_growing;
typedef CGAL::Shape_detection_3::Plane<Traits>               Plane;

int main(int argc, char** argv) 
{
  Pwn_vector points;
  std::ifstream stream(argc > 1 ? argv[1] : "data/cube.pwn");

  if (!stream || 
    !CGAL::read_xyz_points_and_normals(stream,
      std::back_inserter(points),
      Point_map(),
      Normal_map()))
  {
      std::cerr << "Error: cannot read file cube.pwn" << std::endl;
      return EXIT_FAILURE;
  }

  // Call RANSAC shape detection with planes
  Region_growing region_growing;
  region_growing.set_input(points);
  region_growing.add_shape_factory<Plane>();
  region_growing.detect();

  typename Region_growing::Plane_range planes = region_growing.planes();
  // Regularize detected planes
  CGAL::regularize_planes (points,
                           Point_map(),
                           planes,
                           CGAL::Shape_detection_3::Plane_map<Traits>(),
                           CGAL::Shape_detection_3::Point_to_shape_index_map<Traits>(points, planes),
                           true, // Regularize parallelism
                           true, // Regularize orthogonality
                           false, // Do not regularize coplanarity
                           true, // Regularize Z-symmetry (default)
                           10); // 10 degrees of tolerance for parallelism/orthogonality
  
  return EXIT_SUCCESS;
}
