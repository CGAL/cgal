#include "include/utils.h"
#include "include/Saver.h"
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Shape_detection/Efficient_RANSAC.h>
#include <CGAL/Shape_regularization/regularize_planes.h>

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;

using Point  = typename Kernel::Point_3;
using Vector = typename Kernel::Vector_3;
using Plane  = typename Kernel::Plane_3;

using Point_with_normal = std::pair<Point, Vector>;
using Pwn_vector        = std::vector<Point_with_normal>;
using Point_map         = CGAL::First_of_pair_property_map<Point_with_normal>;
using Normal_map        = CGAL::Second_of_pair_property_map<Point_with_normal>;
using Transform         = CGAL::Aff_transformation_3<Kernel>;

using Traits = CGAL::Shape_detection::Efficient_RANSAC_traits<
  Kernel, Pwn_vector, Point_map, Normal_map>;
using Efficient_RANSAC = CGAL::Shape_detection::Efficient_RANSAC<Traits>;

template<typename OutputIterator>
void generate_random_points(
  const Point& origin,
  const Vector& base1,
  const Vector& base2,
  const std::size_t nb_pts,
  OutputIterator output) {

  Vector normal = CGAL::cross_product(base1, base2);
  normal /= std::sqrt(normal * normal);

  // Vector noise = 0.00001 * normal;
  for (std::size_t i = 0; i < nb_pts; ++i) {
    const Point point = origin
      + CGAL::get_default_random().get_double() * base1
      + CGAL::get_default_random().get_double() * base2;
      // + CGAL::get_default_random().get_double() * noise;
      *(output++) = std::make_pair(point, normal);
  }
}

double to_rad(const double deg) { return deg * CGAL_PI / 180.0; }
double to_deg(const double rad) { return rad * 180.0 / CGAL_PI; }

void rotate_vector(
  Vector& vector,
  const int axis,
  const double angle) {

  const double cos_a = std::cos(angle);
  const double sin_a = std::sin(angle);

  Transform transform;
  if (axis == 0) {      // around x
    transform = Transform(
      1.0,   0.0,    0.0,
      0.0, cos_a, -sin_a,
      0.0, sin_a,  cos_a);
  }
  else if (axis == 1) { // around y
    transform = Transform(
      cos_a, 0.0, sin_a,
        0.0, 1.0,   0.0,
     -sin_a, 0.0, cos_a);
  }
  else if (axis == 2) { // around z
    transform = Transform(
      cos_a, -sin_a, 0.0,
      sin_a,  cos_a, 0.0,
        0.0,    0.0, 1.0);
  }
  else abort();
  vector = transform(vector);
}

std::vector<Plane> get_ransac_planes(
  const Efficient_RANSAC& ransac) {

  std::vector<Plane> planes;
  for (const auto& shape : ransac.shapes()) {
    planes.push_back((Plane)(*(dynamic_cast<
    CGAL::Shape_detection::Plane<Traits>*>(shape.get()))));
  }
  return planes;
}

bool planes_difference(
  const Plane& a,
  const Plane& b) {

  return (
    (CGAL::abs(a.a() - b.a()) > 1e-6) ||
    (CGAL::abs(a.b() - b.b()) > 1e-6) ||
    (CGAL::abs(a.c() - b.c()) > 1e-6) ||
    (CGAL::abs(a.d() - b.d()) > 1e-6) );
}

void check_ransac_size(
  const Efficient_RANSAC& ransac,
  const std::size_t nb) {

  if (ransac.shapes().size() != nb) {
    std::cerr << "Error: " << ransac.shapes().size() <<
    " detected plane(s) instead of " << nb << "." << std::endl;
    abort();
  }
}

void check_planes_unchanged(
  const std::vector<Plane>& before,
  const std::vector<Plane>& after) {

  for (std::size_t i = 0; i < before.size(); ++i)
    if (planes_difference(before[i], after[i]))
      std::cerr << "Error: [" << before[i] <<
      "] was altered as [" << after[i] << "]" << std::endl;
}

void check_planes_changed(
  const std::vector<Plane>& before,
  const std::vector<Plane>& after) {

  for (std::size_t i = 0; i < before.size(); ++i) {
    if (planes_difference(before[i], after[i])) {
      return;
    }
  }
  std::cerr << "Error: no plane has been altered by regularization" <<
  " while at least one should have." << std::endl;
}

int main() {

  Vector vx(1.0, 0.0, 0.0),
         vy(0.0, 1.0, 0.0),
         vz(0.0, 0.0, 1.0);

  Vector dvx = vx;
  rotate_vector(dvx, 1, to_rad(10.0));

  Efficient_RANSAC ransac;
  const std::size_t nb_pts = 5000;

  Efficient_RANSAC::Parameters op;
  op.probability = 0.05;
  op.min_points = nb_pts / 2;
  op.epsilon = 0.02;
  op.cluster_epsilon = 0.2;
  op.normal_threshold = 0.8;

  // Test parallelism.
  std::cerr << "Testing parallelism..." << std::endl;
  {
    Pwn_vector points;
    generate_random_points(
      Point(0.0, 0.0, 0.0), vx, vy,
      5000, std::back_inserter(points));

    generate_random_points(
      Point(0.0, 0.0, 0.5), dvx, vy,
      5000, std::back_inserter(points));

    ransac.set_input(points);
    ransac.add_shape_factory< CGAL::Shape_detection::Plane<Traits> >();
    ransac.detect(op);
    check_ransac_size(ransac, 2);

    const std::vector<Plane> before = get_ransac_planes(ransac);

    // Test regularization.
    Efficient_RANSAC::Plane_range planes = ransac.planes();
    CGAL::Shape_regularization::Planes::regularize_planes(
      points,
      Point_map(),
      planes,
      CGAL::Shape_detection::Plane_map<Traits>(),
      CGAL::Shape_detection::Point_to_shape_index_map<Traits>(points, planes),
      true, false, false, false, 5.0);

    std::vector<Plane> after = get_ransac_planes(ransac);

    std::cerr << "* nothing should change now..." << std::endl;
    check_planes_unchanged(before, after);

    CGAL::Shape_regularization::Planes::regularize_planes(
      planes,
      points,
      CGAL::parameters::
      plane_map(CGAL::Shape_detection::Plane_map<Traits>()).
      point_map(Point_map()).
      plane_index_map(
        CGAL::Shape_detection::Point_to_shape_index_map<Traits>(points, planes)).
      regularize_parallelism(true).
      regularize_orthogonality(false).
      regularize_coplanarity(false).
      regularize_axis_symmetry(false).
      maximum_angle(15.0));

    after = get_ransac_planes(ransac);

    std::cerr << "* something should change now..." << std::endl;
    check_planes_changed(before, after);
    ransac.clear();
  }

  // Test orthogonality.
  std::cerr << "Testing orthogonality..." << std::endl;
  {
    Pwn_vector points;
    generate_random_points(
      Point(0.0, 0.0, 0.0), vy, vz,
      5000, std::back_inserter(points));

    generate_random_points(
      Point(0.5, 0.0, 0.0), dvx, vy,
      5000, std::back_inserter(points));

    ransac.set_input(points);
    ransac.add_shape_factory< CGAL::Shape_detection::Plane<Traits> >();
    ransac.detect(op);
    check_ransac_size(ransac, 2);

    const std::vector<Plane> before = get_ransac_planes(ransac);

    Efficient_RANSAC::Plane_range planes = ransac.planes();
    CGAL::Shape_regularization::Planes::regularize_planes(
      points,
      Point_map(),
      planes,
      CGAL::Shape_detection::Plane_map<Traits>(),
      CGAL::Shape_detection::Point_to_shape_index_map<Traits>(points, planes),
      false, true, false, true, 5.0);

    std::vector<Plane> after = get_ransac_planes(ransac);

    std::cerr << "* nothing should change now..." << std::endl;
    check_planes_unchanged(before, after);

    CGAL::Shape_regularization::Planes::regularize_planes(
      planes,
      CGAL::Shape_detection::Plane_map<Traits>(),
      points,
      Point_map(),
      CGAL::parameters::plane_index_map(
        CGAL::Shape_detection::Point_to_shape_index_map<Traits>(points, planes)).
      regularize_parallelism(false).
      regularize_orthogonality(true).
      regularize_coplanarity(false).
      regularize_axis_symmetry(true).
      maximum_angle(15.0));

    after = get_ransac_planes(ransac);

    std::cerr << "* something should change now..." << std::endl;
    check_planes_changed(before, after);
    ransac.clear();
  }

  // Test coplanarity.
  std::cerr << "Testing coplanarity..." << std::endl;
  {
    Pwn_vector points;
    generate_random_points(
      Point(0.0, 0.0, 0.0), vx, vy,
      5000, std::back_inserter(points));

    generate_random_points(
      Point(1.0, 1.0, 0.2), vx, vy,
      5000, std::back_inserter(points));

    ransac.set_input(points);
    ransac.add_shape_factory< CGAL::Shape_detection::Plane<Traits> >();
    ransac.detect(op);
    check_ransac_size(ransac, 2);

    const std::vector<Plane> before = get_ransac_planes(ransac);

    Efficient_RANSAC::Plane_range planes = ransac.planes();
    CGAL::Shape_regularization::Planes::regularize_planes(
      points,
      Point_map(),
      planes,
      CGAL::Shape_detection::Plane_map<Traits>(),
      CGAL::Shape_detection::Point_to_shape_index_map<Traits>(points, planes),
      true, false, true, false, 5.0, 0.1);

    std::vector<Plane> after = get_ransac_planes(ransac);

    std::cerr << "* nothing should change now..." << std::endl;
    check_planes_unchanged(before, after);

    CGAL::Shape_regularization::Planes::regularize_planes(
      planes,
      points,
      CGAL::parameters::
      point_map(Point_map()).
      plane_map(CGAL::Shape_detection::Plane_map<Traits>()).
      plane_index_map(
        CGAL::Shape_detection::Point_to_shape_index_map<Traits>(points, planes)).
      regularize_parallelism(true).
      regularize_orthogonality(false).
      regularize_coplanarity(true).
      regularize_axis_symmetry(false).
      maximum_angle(5.0).
      maximum_offset(0.3));

    after = get_ransac_planes(ransac);

    std::cerr << "* something should change now..." << std::endl;
    check_planes_changed(before, after);
    ransac.clear();
  }

  // Test symmetry.
  std::cerr << "Testing symmetry..." << std::endl;
  {
    Pwn_vector points;
    Vector dvx1 = vx;
    Vector dvx2 = vx;
    rotate_vector(dvx1, 1, to_rad( 40.0));
    rotate_vector(dvx2, 1, to_rad(-50.0));

    generate_random_points(
      Point(0.0, 0.0, -0.5), dvx1, vy,
      5000, std::back_inserter(points));

    generate_random_points(
      Point(0.0, 0.0, 0.5), dvx2, vy,
      5000, std::back_inserter(points));

    ransac.set_input(points);
    ransac.add_shape_factory< CGAL::Shape_detection::Plane<Traits> >();
    ransac.detect(op);
    check_ransac_size(ransac, 2);

    const std::vector<Plane> before = get_ransac_planes(ransac);

    Efficient_RANSAC::Plane_range planes = ransac.planes();
    CGAL::Shape_regularization::Planes::regularize_planes(
      points,
      Point_map(),
      planes,
      CGAL::Shape_detection::Plane_map<Traits>(),
      CGAL::Shape_detection::Point_to_shape_index_map<Traits>(points, planes),
      false, false, false, true, 5.0, 0.01, Vector(1.0, 0.0, 0.0));

    std::vector<Plane> after = get_ransac_planes(ransac);

    std::cerr << "* nothing should change now..." << std::endl;
    check_planes_unchanged(before, after);

    CGAL::Shape_regularization::Planes::regularize_planes(
      planes,
      CGAL::Shape_detection::Plane_map<Traits>(),
      points,
      Point_map(),
      CGAL::parameters::plane_index_map(
        CGAL::Shape_detection::Point_to_shape_index_map<Traits>(points, planes)).
      regularize_parallelism(false).
      regularize_orthogonality(false).
      regularize_coplanarity(false).
      regularize_axis_symmetry(true).
      maximum_angle(15.0).
      symmetry_direction(Vector(1.0, 0.0, 0.0)));

    after = get_ransac_planes(ransac);

    std::cerr << "* something should change now..." << std::endl;
    check_planes_changed(before, after);
    ransac.clear();
  }
  return EXIT_SUCCESS;
}
