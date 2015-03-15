#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/IO/read_xyz_points.h>
#include <CGAL/Point_with_normal_3.h>
#include <CGAL/property_map.h>

#include <CGAL/Efficient_RANSAC_3.h>
#include <CGAL/Plane_shape.h>

#include <iostream>
#include <fstream>

// Type declarations
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Point_with_normal_3<Kernel>                   Point_with_normal;
typedef std::vector<Point_with_normal>                      Pwn_vector;
typedef CGAL::Identity_property_map<Point_with_normal>      Point_pmap;
typedef CGAL::Normal_of_point_with_normal_pmap<Kernel>      Normal_pmap;

// In Efficient_RANSAC_traits_3 the basic types, i.e., Point and Vector types
// as well as iterator type and property maps, are defined.
typedef CGAL::Efficient_RANSAC_traits_3<Kernel,
  Pwn_vector::iterator, Point_pmap, Normal_pmap>            Traits;
typedef CGAL::Efficient_RANSAC_3<Traits>                     Efficient_RANSAC;

int main() 
{
  // Points with normals.
  Pwn_vector points;

  // Loads point set from a file. 
  // read_xyz_points_and_normals takes an OutputIterator for storing the points
  // and a property map to store the normal vector with each point.
  std::ifstream stream("cube.pwn");

  if (!stream || 
    !CGAL::read_xyz_points_and_normals(stream,
      std::back_inserter(points),
      Normal_pmap()))
  {
      std::cerr << "Error: cannot read file cube.pwn" << std::endl;
      return EXIT_FAILURE;
  }

  // Instantiates shape detection engine and provides input data.
  Efficient_RANSAC sd(points.begin(), points.end(),
	                 Point_pmap(), Normal_pmap());

  // Registers planar shapes via template method.
  sd.add_shape_factory<CGAL::Plane_shape<Traits> >();

  // Detects registered shapes with default parameters.
  sd.detect();

  // Prints number of detected shapes.
  std::cout << sd.number_of_shapes() << " shapes detected." << std::endl;

  return EXIT_SUCCESS;
}
