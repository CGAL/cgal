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

int main() {
  sign_test();

  Point_3 a = {0, 0, 0};
  Point_3 b = {0, 1, 0};
  Point_3 c = {1, 0, 0};

  const query queries[] = {
    { {  1, 0,  0},    0.},
    { {  1, 0,  1},   -45.},
    { {  0, 0,  1},   -90.},
    { { -1, 0,  1},  -135.},
    { { -1, 0,  0},  -180.},
    { { -1, 0, -1}, 135.},
    { {  0, 0, -1},  90.},
    { {  1, 0, -1},  45.},
  };

  for(auto query: queries) {
    const auto& expected = query.expected_angle;
    const auto& p = query.p;
    auto approx = CGAL::approximate_dihedral_angle(a, b, c, p);
    std::cout << approx << "  -- " << expected << '\n';
    assert( std::abs(approx - expected) < 0.1 );
  }
}
