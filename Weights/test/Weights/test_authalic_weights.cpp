#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Weights/authalic_weights.h>
#include <CGAL/Weights/utils.h>

// Typedefs.
using SCK   = CGAL::Simple_cartesian<double>;
using EPICK = CGAL::Exact_predicates_inexact_constructions_kernel;
using EPECK = CGAL::Exact_predicates_exact_constructions_kernel;

template<typename Kernel>
bool test_query(
  const typename Kernel::FT x,
  const typename Kernel::FT y) {

  using Point_2 = typename Kernel::Point_2;
  using Point_3 = typename Kernel::Point_3;

  // 2D configuration.
  const Point_2 t2 = Point_2(-1,  0);
  const Point_2 r2 = Point_2( 0, -1);
  const Point_2 p2 = Point_2( 1,  0);
  const Point_2 q2 = Point_2( x,  y);

  // 3D configuration.
  const Point_3 t3 = Point_3(-1,  0, 1);
  const Point_3 r3 = Point_3( 0, -1, 1);
  const Point_3 p3 = Point_3( 1,  0, 1);
  const Point_3 q3 = Point_3( x,  y, 1);

  const auto a2 =
    CGAL::Weights::authalic_weight(t2, r2, p2, q2);
  const auto b2 =
    CGAL::Weights::half_authalic_weight(
      CGAL::Weights::cotangent(t2, r2, q2),
      CGAL::Weights::squared_distance(q2, r2)) +
    CGAL::Weights::half_authalic_weight(
      CGAL::Weights::cotangent(q2, r2, p2),
      CGAL::Weights::squared_distance(q2, r2));
  assert(a2 == b2);

  const auto a3 =
    CGAL::Weights::authalic_weight(t3, r3, p3, q3);
  const auto b3 =
    CGAL::Weights::half_authalic_weight(
      CGAL::Weights::cotangent(t3, r3, q3),
      CGAL::Weights::squared_distance(q3, r3)) +
    CGAL::Weights::half_authalic_weight(
      CGAL::Weights::cotangent(q3, r3, p3),
      CGAL::Weights::squared_distance(q3, r3));
  assert(a3 == b3);

  assert(a2 == a3 && b2 == b3);
  return true;
}

template<typename Kernel>
bool test_authalic_weights() {

  assert(test_query<Kernel>(0, 0));
  return true;
}

int main() {

  assert(test_authalic_weights<SCK>());
  assert(test_authalic_weights<EPICK>());
  assert(test_authalic_weights<EPECK>());
  std::cout << "* test_authalic_weights: SUCCESS" << std::endl;
  return EXIT_SUCCESS;
}
