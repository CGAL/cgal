#include "include/test_efficient_RANSAC_generators.h"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Simple_cartesian.h>

#include <CGAL/Shape_detection/Efficient_RANSAC.h>
#include <CGAL/Point_with_normal_3.h>
#include <CGAL/property_map.h>

template <class K>
bool test_plane_connected_component() {
  const int NB_ROUNDS = 10;

  typedef typename K::FT                                      FT;
  typedef CGAL::Point_with_normal_3<K>                        Pwn;
  typedef CGAL::Point_3<K>                                    Point;
  typedef CGAL::Vector_3<K>                                   Vector;
  typedef std::vector<Pwn>                                    Pwn_vector;
  typedef CGAL::Identity_property_map<Pwn>                    Point_map;
  typedef CGAL::Normal_of_point_with_normal_map<K>            Normal_map;

  typedef typename CGAL::Shape_detection::Efficient_RANSAC_traits<K, Pwn_vector, Point_map, Normal_map> Traits;

  typedef typename CGAL::Shape_detection::Efficient_RANSAC<Traits> Efficient_ransac;
  typedef typename CGAL::Shape_detection::Plane<Traits> Plane;

  std::size_t success = 0;

  for (int i = 0 ; i < NB_ROUNDS ; i++) {
    Pwn_vector points;

    Vector normal;
    CGAL::Bbox_3 bbox(-10, -10, -10, 10, 10, 10);

    // Sample 4 rectangles with 0.05 spacing between points
    // and 0.2 spacing between rectangles.
    Vector offset[] = {Vector((FT) 0, (FT) 0, (FT) 0),
      Vector((FT) 1.2, (FT) 0, (FT) 0),
      Vector((FT) 0, (FT) 1.2, (FT) 0),
      Vector((FT) 1.2, (FT) 1.2, (FT) 0)};

    for (std::size_t j = 0;j<4;j++) {
      for (std::size_t x = 0;x<=20;x++)
        for (std::size_t y = 0;y<=20;y++)
          points.push_back(Pwn(Point(FT(x * 0.05), FT(y * 0.05), (FT) 1.0) + offset[j],
                                Vector((FT) 0, (FT) 0, (FT) 1)));
    }

    Efficient_ransac ransac;

    ransac.template add_shape_factory<Plane>();

    ransac.set_input(points);

    // Same parameters as for the parameters unit tests, besides
    // the cluster_epsilon.
    typename Efficient_ransac::Parameters parameters;
    parameters.probability = 0.05f;
    parameters.min_points = 100;
    parameters.epsilon = 0.002f;
    parameters.normal_threshold = 0.9f;

    // For the first half the rounds chose a high cluster_epsilon to find one
    // shape and for the second half choose a small cluster_epsilon to find
    // four separated shapes.
    if (i < NB_ROUNDS/2)
      parameters.cluster_epsilon = 0.201f;
    else
      parameters.cluster_epsilon = 0.051f;

    if (!ransac.detect(parameters)) {
      std::cout << " aborted" << std::endl;
      return false;
    }

    typename Efficient_ransac::Shape_range shapes = ransac.shapes();

    if (i < NB_ROUNDS/2 && shapes.size() != 1)
      continue;

    if (i >= NB_ROUNDS/2 && shapes.size() != 4)
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

  std::cout << "test_plane_connected_component<CGAL::Simple_cartesian<float>> ";
  if (!test_plane_connected_component<CGAL::Simple_cartesian<float> >())
    success = false;

  std::cout << "test_plane_connected_component<CGAL::Simple_cartesian<double>> ";
  if (!test_plane_connected_component<CGAL::Simple_cartesian<double> >())
    success = false;

  std::cout << "test_plane_connected_component<CGAL::Exact_predicates_inexact_constructions_kernel> ";
  if (!test_plane_connected_component<CGAL::Exact_predicates_inexact_constructions_kernel>())
    success = false;

  return (success) ? EXIT_SUCCESS : EXIT_FAILURE;
}
