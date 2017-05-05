#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Point_with_normal_3.h>
#include <CGAL/property_map.h>

#include <CGAL/Shape_detection_3.h>
#include <CGAL/regularize_planes.h>
#include <CGAL/Random.h>
#include <CGAL/Aff_transformation_3.h>

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

typedef CGAL::Aff_transformation_3<Kernel>                   Transform;


template <typename OutputIterator>
void generate_random_points (const Point& origin, const Vector& base1, const Vector& base2,
                             std::size_t nb_pts, OutputIterator output)
{
  Vector normal = CGAL::cross_product (base1, base2);
  normal = normal / std::sqrt (normal * normal);
  //  Vector noise = 0.00001 * normal;

  for (std::size_t i = 0; i < nb_pts; ++ i)
    {
      Point point = origin
        + CGAL::get_default_random().get_double() * base1
        + CGAL::get_default_random().get_double() * base2;
        //        + CGAL::get_default_random().get_double() * noise;
      *(output ++) = std::make_pair (point, normal);
    }
}

double to_rad (const double& deg) { return deg * CGAL_PI / 180.; }
double to_deg (const double& rad) { return rad * 180. / CGAL_PI; }

void rotate_vector (Vector& vector, int axis, double angle)
{
  double cos_a = std::cos (angle);
  double sin_a = std::sin (angle);

  Transform trans;
  if (axis == 0) // around x
    {
      trans = Transform (1.,  0.,     0.,
                         0.,  cos_a, -sin_a,
                         0.,  sin_a,  cos_a);
    }
  else if (axis == 1) // around y
    {
      trans = Transform ( cos_a, 0., sin_a,
                          0.,    1., 0.,
                         -sin_a, 0., cos_a);
    }
  else if (axis == 2) // around z
    {
      trans = Transform ( cos_a, -sin_a, 0.,
                          sin_a,  cos_a, 0.,
                          0.,     0.,    1.);
    }
  else
    abort();

  vector = trans (vector);
  
}

std::vector<Plane> get_ransac_planes (const Efficient_ransac& ransac)
{
  std::vector<Plane> out;
  BOOST_FOREACH(boost::shared_ptr<Efficient_ransac::Shape> shape, ransac.shapes())
    {
      out.push_back ((Plane)(*(dynamic_cast<CGAL::Shape_detection_3::Plane<Traits>*>(shape.get ()))));
    }
  return out;
}

bool planes_differ (const Plane& a, const Plane& b)
{
  return ((std::fabs (a.a() - b.a()) > 1e-6) ||
          (std::fabs (a.b() - b.b()) > 1e-6) ||
          (std::fabs (a.c() - b.c()) > 1e-6) ||
          (std::fabs (a.d() - b.d()) > 1e-6));
}

void check_ransac_size (const Efficient_ransac& ransac, std::size_t nb)
{
  if (ransac.shapes().size() != nb)
    {
      std::cerr << "Error: " << ransac.shapes().size() << " detected plane(s) instead of " << nb << "." << std::endl;
      abort();
    }
}

void check_planes_unchanged (const std::vector<Plane>& before, const std::vector<Plane>& after)
{
  for (std::size_t i = 0; i < before.size(); ++ i)
    if (planes_differ (before[i], after[i]))
      std::cerr << "Error: [" << before[i] << "] was altered as [" << after[i] << "]" << std::endl;
}

void check_planes_changed (const std::vector<Plane>& before, const std::vector<Plane>& after)
{
  for (std::size_t i = 0; i < before.size(); ++ i)
    if (planes_differ (before[i], after[i]))
      return;
  std::cerr << "Error: no plane has been altered by regularization while at least one should have." << std::endl;
}


int main() 
{
  Vector vx (1., 0., 0.),
         vy (0., 1., 0.),
         vz (0., 0., 1.);
  
  Vector dvx = vx;
  rotate_vector (dvx, 1, to_rad(10.));


  Efficient_ransac ransac;
  
  const std::size_t nb_pts = 5000;
  
  Efficient_ransac::Parameters op;
  op.probability = 0.05;
  op.min_points = nb_pts / 2;
  op.epsilon = 0.02;
  op.cluster_epsilon = 0.2;
  op.normal_threshold = 0.8;
  
  // Test parallelism
  std::cerr << "Testing parallelism..." << std::endl;
  {
    Pwn_vector points;

    generate_random_points (Point (0., 0., 0.), vx, vy,
                            5000, std::back_inserter (points));

      
    generate_random_points (Point (0., 0., 0.5), dvx, vy,
                            5000, std::back_inserter (points));

    ransac.set_input(points);
    ransac.add_shape_factory<CGAL::Shape_detection_3::Plane<Traits> >();
    ransac.detect(op);
    check_ransac_size (ransac, 2);

    std::vector<Plane> before = get_ransac_planes(ransac);
    CGAL::regularize_planes (ransac, true, false, false, false, 5.);
    std::vector<Plane> after = get_ransac_planes(ransac);

    std::cerr << " * Nothing should change now..." << std::endl;
    check_planes_unchanged (before, after);
    
    CGAL::regularize_planes (ransac, true, false, false, false, 15.);
    after = get_ransac_planes(ransac);

    std::cerr << " * Something should change now..." << std::endl;
    check_planes_changed (before, after);
    ransac.clear ();
  }

  // Test orthogonality
  std::cerr << "Testing orthogonality..." << std::endl;
  {
    Pwn_vector points;

    generate_random_points (Point (0., 0., 0.), vy, vz,
                            5000, std::back_inserter (points));
      
    generate_random_points (Point (0.5, 0., 0.), dvx, vy,
                            5000, std::back_inserter (points));

    
    ransac.set_input(points);
    ransac.add_shape_factory<CGAL::Shape_detection_3::Plane<Traits> >();
    ransac.detect(op);
    check_ransac_size (ransac, 2);

    std::vector<Plane> before = get_ransac_planes(ransac);
    CGAL::regularize_planes (ransac, false, true, false, true, 5.);
    std::vector<Plane> after = get_ransac_planes(ransac);

    std::cerr << " * Nothing should change now..." << std::endl;
    check_planes_unchanged (before, after);
    
    CGAL::regularize_planes (ransac, false, true, false, true, 15.);
    after = get_ransac_planes(ransac);

    std::cerr << " * Something should change now..." << std::endl;
    check_planes_changed (before, after);
    ransac.clear ();
  }

  // Test coplanarity
  std::cerr << "Testing coplanarity..." << std::endl;
  {
    Pwn_vector points;

    generate_random_points (Point (0., 0., 0.), vx, vy,
                            5000, std::back_inserter (points));
      
    generate_random_points (Point (1., 1., 0.2), vx, vy,
                            5000, std::back_inserter (points));

    ransac.set_input(points);
    ransac.add_shape_factory<CGAL::Shape_detection_3::Plane<Traits> >();
    ransac.detect(op);
    check_ransac_size (ransac, 2);

    std::vector<Plane> before = get_ransac_planes(ransac);
    CGAL::regularize_planes (ransac, true, false, true, false, 5., 0.1);
    std::vector<Plane> after = get_ransac_planes(ransac);

    std::cerr << " * Nothing should change now..." << std::endl;
    check_planes_unchanged (before, after);
    
    CGAL::regularize_planes (ransac, true, false, true, false, 5., 0.3);
    after = get_ransac_planes(ransac);

    std::cerr << " * Something should change now..." << std::endl;
    check_planes_changed (before, after);
    ransac.clear ();
  }

  // Test symmetry
  std::cerr << "Testing symmetry..." << std::endl;
  {
    Pwn_vector points;
    
    Vector dvx1 = vx;
    Vector dvx2 = vx;
    rotate_vector (dvx1, 1, to_rad(40.));
    rotate_vector (dvx2, 1, to_rad(-50.));

    generate_random_points (Point (0., 0., -0.5), dvx1, vy,
                            5000, std::back_inserter (points));
      
    generate_random_points (Point (0., 0., 0.5), dvx2, vy,
                            5000, std::back_inserter (points));

    ransac.set_input(points);
    ransac.add_shape_factory<CGAL::Shape_detection_3::Plane<Traits> >();
    ransac.detect(op);
    check_ransac_size (ransac, 2);

    std::vector<Plane> before = get_ransac_planes(ransac);
    CGAL::regularize_planes (ransac, false, false, false, true, 5., 0.01, Vector(1., 0., 0.));
    std::vector<Plane> after = get_ransac_planes(ransac);

    std::cerr << " * Nothing should change now..." << std::endl;
    check_planes_unchanged (before, after);
    
    CGAL::regularize_planes (ransac, false, false, false, true, 15., 0.01, Vector(1., 0., 0.));
    after = get_ransac_planes(ransac);

    std::cerr << " * Something should change now..." << std::endl;
    check_planes_changed (before, after);

    ransac.clear ();
  }

  
  
  return EXIT_SUCCESS;
}
