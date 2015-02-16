#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/IO/read_xyz_points.h>
#include <CGAL/Point_with_normal_3.h>
#include <CGAL/property_map.h>
#include <CGAL/Timer.h>
#include <CGAL/number_utils.h>

#include <CGAL/Shape_detection_3.h>
#include <CGAL/Plane_shape.h>
#include <CGAL/Cylinder_shape.h>
#include <CGAL/Cone_shape.h>
#include <CGAL/Sphere_shape.h>
#include <CGAL/Torus_shape.h>

#include <iostream>
#include <fstream>

// Type declarations
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT                                          FT;
typedef CGAL::Point_with_normal_3<Kernel>                   Point_with_normal;
typedef std::vector<Point_with_normal>                      Pwn_vector;
typedef CGAL::Identity_property_map<Point_with_normal>      Point_pmap;
typedef CGAL::Normal_of_point_with_normal_pmap<Kernel>      Normal_pmap;

// In Shape_detection_traits_3 the basic types, i.e., Point and Vector types
// as well as iterator type and property maps, are defined.
typedef CGAL::Shape_detection_traits_3<Kernel,
  Pwn_vector::iterator, Point_pmap, Normal_pmap>            Traits;
typedef CGAL::Shape_detection_3<Traits>                     Shape_detection;


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

  // Measures time before setting up the shape detection.
  CGAL::Timer time;
  time.start();

  // Instantiates shape detection engine and provides input data.
  Shape_detection sd(points.begin(),
    points.end(), Point_pmap(), Normal_pmap());

  // Registers detection of planes

  sd.add_shape_factory(new 
    CGAL::Shape_factory<CGAL::Plane_shape<Traits> >);

  // Detects shapes.
  sd.detect();

  // Measures time after detection.
  time.stop();

  // Prints number of assigned shapes and unsassigned points.
  std::cout << "time: " << time.time() * 1000 << "ms" << std::endl;
  std::cout << sd.number_of_shapes() << " primitives, "
    << sd.number_of_unassigned_points()
    << " unassigned points" << std::endl;

  // Shape_detection_3::shapes_begin() provides
  // an iterator to the detected shapes.
  Shape_detection::Shape_iterator it = sd.shapes_begin();
  while (it != sd.shapes_end()) {

    const Shape_detection::Shape *shape = *it;
    // Using Shape_base::info() for printing 
    // the parameters of the detected shape.
    std::cout << (*it)->info();

    // Sums distances of points to detected shapes.
    FT sum_distances = 0;

    // Iterates through point indices assigned to each detected shape.
    std::vector<size_t>::const_iterator
      index_it = (*it)->assigned_points().begin();

    while (index_it != (*it)->assigned_points().end()) {
      
      // Retrieves point
      const Point_with_normal &p = *(points.begin() + (*index_it));

      // Adds Euclidean distance between point and shape.
      sum_distances += CGAL::sqrt((*it)->squared_distance(p));

      // Proceeds with next point.
      index_it++;
    }

    // Computes and prints average distance.
    FT average_distance = sum_distances / shape->assigned_points().size();
    std::cout << " average distance: " << average_distance << std::endl;

    // Proceeds with next detected shape.
    it++;
  }

  return EXIT_SUCCESS;
}
