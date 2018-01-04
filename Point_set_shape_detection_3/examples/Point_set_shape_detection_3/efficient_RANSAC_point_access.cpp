#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/IO/read_xyz_points.h>
#include <CGAL/Point_with_normal_3.h>
#include <CGAL/property_map.h>
#include <CGAL/Timer.h>
#include <CGAL/number_utils.h>

#include <CGAL/Shape_detection_3.h>

#include <iostream>
#include <fstream>

// Type declarations
typedef CGAL::Exact_predicates_inexact_constructions_kernel  Kernel;
typedef Kernel::FT                                           FT;
typedef std::pair<Kernel::Point_3, Kernel::Vector_3>         Point_with_normal;
typedef std::vector<Point_with_normal>                       Pwn_vector;
typedef CGAL::First_of_pair_property_map<Point_with_normal>  Point_map;
typedef CGAL::Second_of_pair_property_map<Point_with_normal> Normal_map;

// In Efficient_RANSAC_traits the basic types, i.e., Point and Vector types
// as well as iterator type and property maps, are defined.
typedef CGAL::Shape_detection_3::Efficient_RANSAC_traits<Kernel,
  Pwn_vector, Point_map, Normal_map>            Traits;
typedef CGAL::Shape_detection_3::Efficient_RANSAC<Traits>   Efficient_ransac;
typedef CGAL::Shape_detection_3::Plane<Traits>              Plane;


int main()
{
  // Points with normals.
  Pwn_vector points;

  // Loads point set from a file. 
  // read_xyz_points_and_normals takes an OutputIterator for storing the points
  // and a property map to store the normal vector with each point.
  std::ifstream stream("data/cube.pwn");

  if (!stream ||
    !CGAL::read_xyz_points_and_normals(stream,
      std::back_inserter(points),
      Point_map(),
      Normal_map()))
  {
      std::cerr << "Error: cannot read file cube.pwn" << std::endl;
      return EXIT_FAILURE;
  }

  // Instantiates shape detection engine.
  Efficient_ransac ransac;

  // Provides the input data.
  ransac.set_input(points);

  // Registers detection of planes
  ransac.add_shape_factory<Plane>();

  // Measures time before setting up the shape detection.
  CGAL::Timer time;
  time.start();

  // Build internal data structures.
  ransac.preprocess();

  // Measures time after preprocessing.
  time.stop();

  std::cout << "preprocessing took: " << time.time() * 1000 << "ms" << std::endl;

  // Perform detection several times and choose result with highest coverage.
  Efficient_ransac::Shape_range shapes = ransac.shapes();
  FT best_coverage = 0;

  for (size_t i = 0;i<3;i++) {
    // Reset timer.
    time.reset();
    time.start();

    // Detects shapes.
    ransac.detect();

    // Measures time after detection.
    time.stop();

    // Compute coverage, i.e. ratio of the points assigned to a shape.
    FT coverage = FT(points.size() - ransac.number_of_unassigned_points())
                  / FT(points.size());

    // Prints number of assigned shapes and unsassigned points.
    std::cout << "time: " << time.time() * 1000 << "ms" << std::endl;
    std::cout << ransac.shapes().end() - ransac.shapes().begin() << " primitives, "
      << coverage << " coverage" << std::endl;
    
    // Choose result with highest coverage.
    if (coverage > best_coverage) {
      best_coverage = coverage;
      // Efficient_ransac::shapes() provides
      // an iterator range to the detected shapes. 
      shapes = ransac.shapes();
    }
  }

  Efficient_ransac::Shape_range::iterator it = shapes.begin();

  while (it != shapes.end()) {
    boost::shared_ptr<Efficient_ransac::Shape> shape = *it;
    // Using Shape_base::info() for printing 
    // the parameters of the detected shape.
    std::cout << (*it)->info();

    // Sums distances of points to detected shapes.
    FT sum_distances = 0;

    // Iterates through point indices assigned to each detected shape.
    std::vector<std::size_t>::const_iterator
      index_it = (*it)->indices_of_assigned_points().begin();

    while (index_it != (*it)->indices_of_assigned_points().end()) {
      
      // Retrieves point
      const Point_with_normal &p = *(points.begin() + (*index_it));

      // Adds Euclidean distance between point and shape.
      sum_distances += CGAL::sqrt((*it)->squared_distance(p.first));

      // Proceeds with next point.
      index_it++;
    }

    // Computes and prints average distance.
    FT average_distance = sum_distances / shape->indices_of_assigned_points().size();
    std::cout << " average distance: " << average_distance << std::endl;

    // Proceeds with next detected shape.
    it++;
  }

  return EXIT_SUCCESS;
}
