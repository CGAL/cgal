#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/Linear_cell_complex/lcc_scaled_jacobian.h>
#include <array>
#include <cmath>
#include <cstdlib>
#include <iostream>

// Scaled jacobian of hexahedra: a cube, an inverted corner, a mirrored
// (inverted) hexahedron and a degenerate one.

using LCC=CGAL::Linear_cell_complex_for_combinatorial_map<3>;
using Point=LCC::Point;

void check(bool b, const char* what)
{
  if(!b)
  {
    std::cerr<<"lcc_scaled_jacobian_test FAILED: "<<what<<std::endl;
    std::exit(EXIT_FAILURE);
  }
}

double jacobian(const std::array<Point, 8>& p)
{
  LCC lcc;
  LCC::Dart_descriptor d=lcc.make_hexahedron(p[0], p[1], p[2], p[3],
                                             p[4], p[5], p[6], p[7]);
  return CGAL::scaled_jacobian_of_hexahedron(lcc, d);
}

int main()
{
  const std::array<Point, 8> cube={Point(0,0,0), Point(1,0,0), Point(1,1,0), Point(0,1,0),
                                   Point(0,1,1), Point(0,0,1), Point(1,0,1), Point(1,1,1)};
  const double j_cube=jacobian(cube);
  check(std::abs(j_cube-1.)<1e-12, "std::abs(j_cube-1.)<1e-12");

  std::array<Point, 8> pushed=cube;
  pushed[0]=Point(0.9,0.9,0.9); // first corner pushed inside
  check(jacobian(pushed)<0., "jacobian(pushed)<0.");

  std::array<Point, 8> mirrored=cube;
  for(Point& q: mirrored) { q=Point(-q.x(), q.y(), q.z()); }
  check(std::abs(jacobian(mirrored)+1.)<1e-12, "std::abs(jacobian(mirrored)+1.)<1e-12");

  std::array<Point, 8> degenerate=cube;
  degenerate[1]=degenerate[0]; // an edge of length 0
  const double j_degenerate=jacobian(degenerate);
  check(std::isfinite(j_degenerate) && j_degenerate<=0., "std::isfinite(j_degenerate) && j_degenerate<=0.");

  std::cout<<"lcc_scaled_jacobian_test: OK"<<std::endl;
  return EXIT_SUCCESS;
}
