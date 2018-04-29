#if defined (_MSC_VER) && !defined (_WIN64)
#pragma warning(disable:4244) // boost::number_distance::distance()
                              // converts 64 to 32 bits integers
#endif

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/IO/read_xyz_points.h>
#include <CGAL/Point_with_normal_3.h>
#include <CGAL/property_map.h>

#include <CGAL/Shape_detection_3.h>

#include <CGAL/Timer.h>

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
typedef CGAL::Shape_detection_3::Efficient_RANSAC<Traits>    Efficient_ransac;
typedef CGAL::Shape_detection_3::Region_growing<Traits>      Region_growing;
typedef CGAL::Shape_detection_3::Plane<Traits>               Plane;

struct Timeout_callback
{
  mutable int nb;
  mutable CGAL::Timer timer;
  const double limit;
  
  Timeout_callback(double limit) : nb(0), limit(limit)
  {
    timer.start();
  }
  
  bool operator()(double advancement) const
  {
    // Avoid calling time() at every single iteration, which could
    // impact performances very badly
    ++ nb;
    if (nb % 1000 != 0)
      return true;

    // If the limit is reach, interrupt the algorithm
    if (timer.time() > limit)
    {
      std::cerr << "Algorithm takes too long, exiting ("
                << 100. * advancement << "% done)" << std::endl;
      return false;
    }

    return true;
  }
};


// This program both works for RANSAC and Region Growing
template <typename ShapeDetection>
int run(const char* filename)
{
  Pwn_vector points;
  std::ifstream stream(filename);

  if (!stream || 
    !CGAL::read_xyz_points(stream,
      std::back_inserter(points),
      CGAL::parameters::point_map(Point_map()).
      normal_map(Normal_map())))
  {
      std::cerr << "Error: cannot read file cube.pwn" << std::endl;
      return EXIT_FAILURE;
  }

  ShapeDetection shape_detection;
  shape_detection.set_input(points);
  shape_detection.template add_shape_factory<Plane>();

  // Create callback that interrupts the algorithm if it takes more than half a second
  Timeout_callback timeout_callback(0.5);
  
  // Detects registered shapes with default parameters.
  shape_detection.detect(typename ShapeDetection::Parameters(),
                         timeout_callback);

  return EXIT_SUCCESS;
}


int main (int argc, char** argv)
{
  if (argc > 1 && std::string(argv[1]) == "-r")
  {
    std::cout << "Efficient RANSAC" << std::endl;
    return run<Efficient_ransac> ((argc > 2) ? argv[2] : "data/cube.pwn");
  }

  std::cout << "Region Growing" << std::endl;
  return run<Region_growing> ((argc > 1) ? argv[1] : "data/cube.pwn");
}
