#include <CGAL/Simple_cartesian.h>
#include <CGAL/Kernel/global_functions_3.h>
#include <iostream>
#include <utility>
using K = CGAL::Simple_cartesian<double>;
using Point_3 = K::Point_3;

struct query {
  Point_3 p;
  double expected_angle;
};

void sign_test()
{
  K::Point_3 a(0,0,0), b(1,0,0), c(0,1, 0), d(0,0,1);

  assert( CGAL::approximate_dihedral_angle(a, b, c, d) > 0);
  assert( CGAL::approximate_dihedral_angle(c, a, b, d) > 0);
  assert( CGAL::approximate_dihedral_angle(a, d, b, c) > 0);
  assert( CGAL::approximate_dihedral_angle(c, b, d, a) > 0);
  assert( CGAL::approximate_dihedral_angle(d, b, a, c) > 0);
  assert( CGAL::approximate_dihedral_angle(d, c, b, a) > 0);

  assert( CGAL::approximate_dihedral_angle(a, b, d, c) < 0);
  assert( CGAL::approximate_dihedral_angle(c, a, d, b) < 0);
  assert( CGAL::approximate_dihedral_angle(a, d, c, b) < 0);
  assert( CGAL::approximate_dihedral_angle(c, b, a, d) < 0);
  assert( CGAL::approximate_dihedral_angle(d, b, c, a) < 0);
  assert( CGAL::approximate_dihedral_angle(d, c, a, b) < 0);
}

auto almost_equal_angle(double a, double b) {
  return (std::min)(std::abs(a - b), std::abs(a + 360 - b)) < 0.1;
}

void test_regular_tetrahedron()
{
  auto half_root_of_2 = std::sqrt(2) / 2;

  // Regular tetrahedron
  Point_3 a{ -1,  0, -half_root_of_2};
  Point_3 b{  1,  0, -half_root_of_2};
  Point_3 c{  0,  1,  half_root_of_2};
  Point_3 d{  0, -1,  half_root_of_2};
  assert(orientation(a, b, c, d) == CGAL::POSITIVE);
  assert(almost_equal_angle(CGAL::approximate_dihedral_angle(a, b, c, d), 70.5288));
}

int main() {
  std::cout.precision(17);
  sign_test();
  test_regular_tetrahedron();

  Point_3 a = {0, 0, 0};
  Point_3 b = {0, -1, 0}; // ab is oriented so that it sees the plan xz positively.

  const query queries[] = {
    { {  1, 0,  0},    0.},
    { {  1, 0,  1},   45.},
    { {  0, 0,  1},   90.},
    { { -1, 0,  1},  135.},
    { { -1, 0,  0},  180.},
    { { -1, 0, -1}, -135.},
    { {  0, 0, -1},  -90.},
    { {  1, 0, -1},  -45.},
  };

  auto cnt = 0u;
  for(double yc = -10; yc < 10; yc += 0.1) {
    // c can be any point in the half-plane xy, with x>0
    Point_3 c{1, yc, 0};
    // std::cout << "c = " << c << '\n';
    for(const auto& query : queries) {
      for(double yp = -10; yp < 10; yp += 0.3) {
        const auto& expected = query.expected_angle;
        const Point_3 p{query.p.x(), yp, query.p.z()};
        // std::cout << "p = " << p << '\n';
        auto approx = CGAL::approximate_dihedral_angle(a, b, c, p);
        // std::cout << approx << "  -- " << expected << '\n';
        if(!almost_equal_angle(approx, expected)) {
          std::cout << "ERROR:\n";
          std::cout << "CGAL::approximate_dihedral_angle(" << a << ", " << b << ", " << c << ", " << p << ") = " << approx << '\n';
          std::cout << "expected: " << expected << '\n';
          return 1;
        }
        ++cnt;
      }
    }
  }
  std::cout << "OK (" << cnt << " tests)\n";
  assert(cnt > 10000);
}
