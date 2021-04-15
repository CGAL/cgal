#include <fstream>
#include <iostream>

#include <CGAL/Timer.h>
#include <CGAL/number_utils.h>
#include <CGAL/property_map.h>
#include <CGAL/IO/read_xyz_points.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Shape_detection/Efficient_RANSAC.h>

// Type declarations.
typedef CGAL::Exact_predicates_inexact_constructions_kernel  Kernel;
typedef Kernel::FT                                           FT;
typedef std::pair<Kernel::Point_3, Kernel::Vector_3>         Point_with_normal;
typedef std::vector<Point_with_normal>                       Pwn_vector;
typedef CGAL::First_of_pair_property_map<Point_with_normal>  Point_map;
typedef CGAL::Second_of_pair_property_map<Point_with_normal> Normal_map;

typedef CGAL::Shape_detection::Efficient_RANSAC_traits
<Kernel, Pwn_vector, Point_map, Normal_map>             Traits;
typedef CGAL::Shape_detection::Efficient_RANSAC<Traits> Efficient_ransac;
typedef CGAL::Shape_detection::Plane<Traits>            Plane;

int main(int argc, char** argv) {

  // Points with normals.
  Pwn_vector points;

  // Load point set from a file.
  std::ifstream stream((argc > 1) ? argv[1] : "data/cube.pwn");

  if (!stream ||
    !CGAL::read_xyz_points(
      stream,
      std::back_inserter(points),
      CGAL::parameters::point_map(Point_map()).
      normal_map(Normal_map()))) {

    std::cerr << "Error: cannot read file cube.pwn!" << std::endl;
    return EXIT_FAILURE;
  }

  // Instantiate shape detection engine.
  Efficient_ransac ransac;

  // Provide input data.
  ransac.set_input(points);

  // Register detection of planes.
  ransac.add_shape_factory<Plane>();

  // Measure time before setting up the shape detection.
  CGAL::Timer time;
  time.start();

  // Build internal data structures.
  ransac.preprocess();

  // Measure time after preprocessing.
  time.stop();

  std::cout << "preprocessing took: " << time.time() * 1000 << "ms" << std::endl;

  // Perform detection several times and choose result with the highest coverage.
  Efficient_ransac::Shape_range shapes = ransac.shapes();

  FT best_coverage = 0;
  for (std::size_t i = 0; i < 3; ++i) {

    // Reset timer.
    time.reset();
    time.start();

    // Detect shapes.
    ransac.detect();

    // Measure time after detection.
    time.stop();

    // Compute coverage, i.e. ratio of the points assigned to a shape.
    FT coverage =
    FT(points.size() - ransac.number_of_unassigned_points()) / FT(points.size());

    // Print number of assigned shapes and unassigned points.
    std::cout << "time: " << time.time() * 1000 << "ms" << std::endl;
    std::cout << ransac.shapes().end() - ransac.shapes().begin()
    << " primitives, " << coverage << " coverage" << std::endl;

    // Choose result with the highest coverage.
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

    // Use Shape_base::info() to print the parameters of the detected shape.
    std::cout << (*it)->info();

    // Sums distances of points to the detected shapes.
    FT sum_distances = 0;

    // Iterate through point indices assigned to each detected shape.
    std::vector<std::size_t>::const_iterator
    index_it = (*it)->indices_of_assigned_points().begin();

    while (index_it != (*it)->indices_of_assigned_points().end()) {

      // Retrieve point.
      const Point_with_normal& p = *(points.begin() + (*index_it));

      // Adds Euclidean distance between point and shape.
      sum_distances += CGAL::sqrt((*it)->squared_distance(p.first));

      // Proceed with the next point.
      index_it++;
    }

    // Compute and print the average distance.
    FT average_distance = sum_distances / shape->indices_of_assigned_points().size();
    std::cout << " average distance: " << average_distance << std::endl;

    // Proceed with the next detected shape.
    it++;
  }
  return EXIT_SUCCESS;
}
