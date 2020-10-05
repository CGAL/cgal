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

int main() {
  Point_3 a = {0, 0, 0};
  Point_3 b = {0, 1, 0};
  Point_3 c = {1, 0, 0};

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

  for(auto query: queries) {
    const auto& expected = query.expected_angle;
    const auto& p = query.p;
    auto approx = CGAL::approximate_dihedral_angle(a, b, c, p);
    std::cout << approx << "  -- " << expected << '\n';
    assert( std::abs(approx - expected) < 0.1 );
  }
};
