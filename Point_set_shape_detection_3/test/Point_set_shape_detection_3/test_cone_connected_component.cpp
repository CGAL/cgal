#include "generators.h"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Simple_cartesian.h>

#include <CGAL/Shape_detection_3.h>
#include <CGAL/Point_with_normal_3.h>
#include <CGAL/property_map.h>


template <class K>
bool test_cone_connected_component() {
  const int NB_ROUNDS = 10;
  const int NB_POINTS = 2000;

  typedef typename K::FT                                      FT;
  typedef CGAL::Point_with_normal_3<K>                        Pwn;
  typedef CGAL::Point_3<K>                                    Point;
  typedef CGAL::Vector_3<K>                                   Vector;
  typedef std::vector<Pwn>                                    Pwn_vector;
  typedef CGAL::Identity_property_map<Pwn>                    Point_map;
  typedef CGAL::Normal_of_point_with_normal_pmap<K>           Normal_map;

  typedef CGAL::Shape_detection_3::Shape_detection_traits<
    K, Pwn_vector, Point_map, Normal_map>                     Traits;

  typedef CGAL::Shape_detection_3::Efficient_RANSAC<Traits>   Efficient_ransac;
  typedef CGAL::Shape_detection_3::Cone<Traits>               Cone;

  std::size_t success = 0;

  for (int i = 0; i < NB_ROUNDS ; i++) {
    Pwn_vector points;

    // generate random points on a cone
    FT spacing;

    do {
      points.clear();
      CGAL::Bbox_3 bbox(-10, -10, -10, 10, 10, 10);
      Vector axis = random_normal<K>();
      Point apex = random_point_in<K>(bbox);
      FT angle = random_float((FT) 0.1, (FT) 1.4);
      FT start = random_float((FT) 0, (FT) 2.0);
      FT end = start + random_float((FT) 0.5, (FT) 2.0);
      FT mid = (start + end) / (FT)2.0;

      sample_cone(NB_POINTS, apex, axis, angle, start, end,
        std::back_inserter(points)); 

      CGAL::Vector_3<K> n = random_normal<K>();
      CGAL::Plane_3<K> pl(apex + mid * axis, n);

      spacing = FT(0.75) * (mid / cos(angle)) * sin(angle);

      filter_by_distance(pl, spacing * FT(0.5), points);
    } while (points.size() < NB_POINTS / 2.0);

    Efficient_ransac ransac;

    ransac.template add_shape_factory<Cone>();

    ransac.set_input(points);

    // Same parameters as for the parameters unit tests, besides
    // the cluster_epsilon.
    typename Efficient_ransac::Parameters parameters;
    parameters.probability = 0.05f;
    parameters.min_points = points.size()/4;
    parameters.epsilon = 0.002f;
    parameters.normal_threshold = 0.9f;

    // The first half of rounds choose a high cluster_epsilon to get only
    // a single shape and a lower cluster_epsilon for the second half
    // to get two separated shapes.
    if (i < NB_ROUNDS/2)
      parameters.cluster_epsilon = spacing * FT(1.1);
    else
      parameters.cluster_epsilon = spacing * FT(0.35);

    if (!ransac.detect(parameters)) {
      std::cout << " aborted" << std::endl;
      return false;
    }

    typename Efficient_ransac::Shape_range shapes = ransac.shapes();
    
    if (i < NB_ROUNDS/2 && shapes.size() != 1)
      continue;

    if (i >= NB_ROUNDS/2 && shapes.size() != 2)
      continue;

    success++;
  }

  if (success >= NB_ROUNDS * 0.8) {
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

  std::cout << "test_cone_connected_component<CGAL::Simple_cartesian<float>> ";
  if (!test_cone_connected_component<CGAL::Simple_cartesian<float> >()) 
    success = false;

  std::cout << "test_cone_connected_component<CGAL::Simple_cartesian<double>> ";
  if (!test_cone_connected_component<CGAL::Simple_cartesian<double> >())
    success = false;

  std::cout << "test_cone_connected_component<CGAL::Exact_predicates_inexact_constructions_kernel> ";
  if (!test_cone_connected_component<CGAL::Exact_predicates_inexact_constructions_kernel>()) 
    success = false;

  return (success) ? EXIT_SUCCESS : EXIT_FAILURE;
}
