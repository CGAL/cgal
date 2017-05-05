#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Point_with_normal_3.h>
#include <CGAL/property_map.h>

#include <CGAL/Shape_detection_3.h>
#include <CGAL/structure_point_set.h>

#include <CGAL/Random.h>

#include <iostream>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel  Kernel;
typedef Kernel::Point_3                                      Point;
typedef Kernel::Vector_3                                     Vector;
typedef Kernel::Plane_3                                      Plane;
typedef std::pair<Point, Vector>                             Point_with_normal;
typedef std::vector<Point_with_normal>                       Pwn_vector;
typedef CGAL::First_of_pair_property_map<Point_with_normal>  Point_map;
typedef CGAL::Second_of_pair_property_map<Point_with_normal> Normal_map;

typedef CGAL::Shape_detection_3::Shape_detection_traits
  <Kernel, Pwn_vector, Point_map, Normal_map>                Traits;
typedef CGAL::Shape_detection_3::Efficient_RANSAC<Traits>    Efficient_ransac;

typedef CGAL::Point_set_with_structure<Traits>               Points_with_structure;

template <typename OutputIterator>
void generate_random_points (const Point& origin, const Vector& base1, const Vector& base2,
                             std::size_t nb_pts, OutputIterator output)
{
  Vector normal = CGAL::cross_product (base1, base2);
  normal = normal / std::sqrt (normal * normal);

  for (std::size_t i = 0; i < nb_pts; ++ i)
    {
      Point point = origin
        + CGAL::get_default_random().get_double() * base1
        + CGAL::get_default_random().get_double() * base2;

      *(output ++) = std::make_pair (point, normal);
    }
}


int main() 
{
  Vector vx (1., 0., 0.),
         vy (0., 1., 0.),
         vz (0., 0., 1.);
  
  Efficient_ransac ransac;
  ransac.add_shape_factory<CGAL::Shape_detection_3::Plane<Traits> >();
  
  const std::size_t nb_pts = 1000;
  
  Efficient_ransac::Parameters op;
  op.probability = 0.05;
  op.min_points = nb_pts / 2;
  op.epsilon = 0.02;
  op.cluster_epsilon = 0.05;
  op.normal_threshold = 0.8;
  
  Pwn_vector points;

  generate_random_points (Point (0., 0., 0.), vx, vy,
                          5000, std::back_inserter (points));
  generate_random_points (Point (0., 0., 0.), vx, vz,
                          5000, std::back_inserter (points));
  generate_random_points (Point (0., 0., 0.), vy, vz,
                          5000, std::back_inserter (points));
  generate_random_points (Point (0., 0., 1.), vx, vy,
                          5000, std::back_inserter (points));
  generate_random_points (Point (0., 1., 0.), vx, vz,
                          5000, std::back_inserter (points));
  generate_random_points (Point (1., 0., 0.), vy, vz,
                          5000, std::back_inserter (points));


  ransac.set_input(points);
  ransac.detect(op);

  Points_with_structure pss (points.begin(), points.end(), ransac, op.cluster_epsilon);

  std::vector<Point> vertices;

  for (std::size_t i = 0; i < pss.size(); ++ i)
    if (pss.adjacency (i).size () == 3)
      vertices.push_back (pss.point (i));

  if (vertices.size () != 8)
    {
      std::cerr << "Error: 8 point should have been structural vertices." << std::endl;
      return EXIT_FAILURE;
    }

  std::vector<Point> ground_truth;
  ground_truth.push_back (Point (0., 0., 0.));
  ground_truth.push_back (Point (0., 0., 1.));
  ground_truth.push_back (Point (0., 1., 0.));
  ground_truth.push_back (Point (0., 1., 1.));
  ground_truth.push_back (Point (1., 0., 0.));
  ground_truth.push_back (Point (1., 0., 1.));
  ground_truth.push_back (Point (1., 1., 0.));
  ground_truth.push_back (Point (1., 1., 1.));
  std::vector<bool> found (ground_truth.size(), false);
  std::size_t nb_found = 0;
  
  for (std::size_t i = 0; i < vertices.size(); ++ i)
    for (std::size_t j = 0; j < ground_truth.size(); ++ j)
      {
        if (found[j])
          continue;

        if (CGAL::squared_distance (ground_truth[j], vertices[i]) < 1e-6)
          {
            found[j] = true;
            ++ nb_found;
            break;
          }
      }

  if (nb_found != ground_truth.size())
    {
      std::cerr << "Error: the following vert(ex/ices) was/were not found:" << std::endl;
      for (std::size_t i = 0; i < ground_truth.size(); ++ i)
        if (!(found[i]))
          std::cerr << " * " << ground_truth[i] << std::endl;
      return EXIT_FAILURE;
    }
      
  
  return EXIT_SUCCESS;
}
