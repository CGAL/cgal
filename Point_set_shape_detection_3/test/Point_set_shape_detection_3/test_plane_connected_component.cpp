#include "generators.h"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Simple_cartesian.h>

#include <CGAL/Shape_detection_3.h>
#include <CGAL/Point_with_normal_3.h>
#include <CGAL/property_map.h>

const int rounds = 50;

template <class K>
bool test_plane_connected_component() {
  typedef typename CGAL::Point_with_normal_3<K>               Pwn;
  typedef typename CGAL::Point_3<K>                           Point;
  typedef typename CGAL::Vector_3<K>                          Vector;
  typedef std::vector<Pwn>                                    Pwn_vector;
  typedef typename CGAL::Identity_property_map<Pwn>           Point_map;
  typedef typename CGAL::Normal_of_point_with_normal_pmap<K>  Normal_map;

  typedef typename CGAL::Shape_detection_3::Efficient_RANSAC_traits<K,
    Pwn_vector, Point_map, Normal_map>                        Traits;

  typedef typename CGAL::Shape_detection_3::Efficient_RANSAC<Traits>
    Efficient_ransac;

  typedef typename CGAL::Shape_detection_3::Plane<Traits>     Plane;

  std::size_t success = 0;

  for (std::size_t i = 0;i<rounds;i++) {
    Pwn_vector points;

    K::FT dist = 0;
    Vector normal;
    CGAL::Bbox_3 bbox(-10, -10, -10, 10, 10, 10);

    std::size_t index = 0;

    // Sample 4 rectangles with 0.05 spacing between points
    // and 0.2 spacing between rectangles.
    Vector offset[] = {Vector(0, 0, 0), Vector(1.2, 0, 0),
      Vector(0, 1.2, 0), Vector(1.2, 1.2, 0)};

    for (std::size_t j = 0;j<4;j++) {
      for (std::size_t x = 0;x<=20;x++)
        for (std::size_t y = 0;y<=20;y++)
          points.push_back(Pwn(Point(x * 0.05, y * 0.05, 0) + offset[j],
                                Vector(0, 0, 1)));
    }
        
    Efficient_ransac ransac;

    ransac.add_shape_factory<Plane>();
    
    ransac.set_input(points);

    // For the first half the rounds chose a high cluster_epsilon to find one
    // shape and for the second half choose a small cluster_epsilon to find
    // four separated shapes.
    
    Efficient_ransac::Parameters parameters;
    parameters.probability = 0.05f;
    parameters.min_points = 100;
    parameters.epsilon = 0.002f;
    parameters.normal_threshold = 0.9f;

    if (i <= rounds/2)
      parameters.cluster_epsilon = 0.201f;
    else
      parameters.cluster_epsilon = 0.051f;

    if (!ransac.detect(parameters)) {
      std::cout << " aborted" << std::endl;
      return false;
    }
    
    Efficient_ransac::Shape_range shapes = ransac.shapes();

    if (i <= rounds/2 && shapes.size() != 1)
      continue;

    if (i > rounds/2 && shapes.size() != 4)
      continue;

    success++;
  }

  if (success >= rounds * 0.8) {
    std::cout << " succeeded" << std::endl;
    return true;
  }
  else {
    std::cout << " failed" << std::endl;
    return false;
  }
}


int main() {
  bool success = true;

  std::cout << "test_plane_connected_component<CGAL::Simple_cartesian<float>> ";
  if (!test_plane_connected_component<CGAL::Simple_cartesian<float>>()) 
    success = false;

  std::cout << "test_plane_connected_component<CGAL::Simple_cartesian<double>> ";
  if (!test_plane_connected_component<CGAL::Simple_cartesian<double>>())
    success = false;

  std::cout << "test_plane_connected_component<CGAL::Exact_predicates_inexact_constructions_kernel> ";
  if (!test_plane_connected_component<CGAL::Exact_predicates_inexact_constructions_kernel>()) 
    success = false;

  return (success) ? EXIT_SUCCESS : EXIT_FAILURE;
}
