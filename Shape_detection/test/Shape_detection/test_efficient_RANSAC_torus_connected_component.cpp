#include "include/test_efficient_RANSAC_generators.h"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Simple_cartesian.h>

#include <CGAL/Shape_detection/Efficient_RANSAC.h>
#include <CGAL/Point_with_normal_3.h>
#include <CGAL/property_map.h>

template <class K>
bool test_torus_connected_component() {
  const int NB_ROUNDS = 10;
  const int NB_POINTS = 2000;

  typedef typename K::FT                                      FT;
  typedef CGAL::Point_with_normal_3<K>                        Pwn;
  typedef CGAL::Point_3<K>                                    Point;
  typedef CGAL::Vector_3<K>                                   Vector;
  typedef std::vector<Pwn>                                    Pwn_vector;
  typedef CGAL::Identity_property_map<Pwn>                    Point_map;
  typedef CGAL::Normal_of_point_with_normal_map<K>            Normal_map;

  typedef CGAL::Shape_detection::Efficient_RANSAC_traits<K, Pwn_vector, Point_map, Normal_map> Traits;
  typedef CGAL::Shape_detection::Efficient_RANSAC<Traits>   Efficient_ransac;
  typedef CGAL::Shape_detection::Torus<Traits>              Torus;

  std::size_t success = 0;

  for (int i = 0;i<NB_ROUNDS;i++) {
    Pwn_vector points;

    // Generate random points on torus.
    CGAL::Bbox_3 bbox(-10, -10, -10, 10, 10, 10);
    FT minor_radius = (FT) 0.7;
    FT major_radius = (FT) 2.0;
    Vector axis = random_normal<K>();
    Point center = random_point_in<K>(bbox);

    sample_torus(NB_POINTS, center, axis,
      major_radius, minor_radius, std::back_inserter(points));

    CGAL::Vector_3<K> n = random_normal<K>();
    n = CGAL::cross_product(axis, n);
    n = n * (FT) 1.0 / (CGAL::sqrt(n.squared_length()));
    CGAL::Plane_3<K> pl(center, n);

    FT spacing = (FT) 1;

    filter_by_distance(pl, spacing * FT(0.5), points);

    Efficient_ransac ransac;

    ransac.template add_shape_factory<Torus>();

    ransac.set_input(points);

    // Same parameters as for the parameters unit tests, besides
    // the cluster_epsilon.
    typename Efficient_ransac::Parameters parameters;
    parameters.probability = 0.05f;
    parameters.min_points = points.size()/10;
    parameters.epsilon = 0.002f;
    parameters.normal_threshold = 0.9f;

    // The first half of rounds choose a high cluster_epsilon to get only
    // a single shape and a lower cluster_epsilon for the second half
    // to get two separated shapes.
    if (i < NB_ROUNDS/2)
      parameters.cluster_epsilon = spacing * (FT) 1.5;
    else
      parameters.cluster_epsilon = spacing * (FT) 0.9;

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

  std::cout << "test_torus_connected_component<CGAL::Simple_cartesian<float>> ";
  if (!test_torus_connected_component<CGAL::Simple_cartesian<float> >())
    success = false;

  std::cout << "test_torus_connected_component<CGAL::Simple_cartesian<double>> ";
  if (!test_torus_connected_component<CGAL::Simple_cartesian<double> >())
    success = false;

  std::cout << "test_torus_connected_component<CGAL::Exact_predicates_inexact_constructions_kernel> ";
  if (!test_torus_connected_component<CGAL::Exact_predicates_inexact_constructions_kernel>())
    success = false;

  return (success) ? EXIT_SUCCESS : EXIT_FAILURE;
}
