#include <CGAL/Exact_rational.h>
#include <CGAL/Extended_cartesian.h>
#include <CGAL/Nef_polyhedron_2.h>
#include <cassert>

typedef CGAL::Exact_rational FT;
typedef CGAL::Extended_cartesian<FT> Extended_kernel;
typedef CGAL::Nef_polyhedron_2<Extended_kernel> Nef_polyhedron;
typedef Nef_polyhedron::Point Point;
typedef Nef_polyhedron::Line  Line;

int main() {

  Nef_polyhedron N1(Nef_polyhedron::COMPLETE);

  Line l(2.1,4.8,2.0); // l : 2.1x + 4.8y + 2 = 0
  Nef_polyhedron N2(l,Nef_polyhedron::INCLUDED);
  Nef_polyhedron N3 = N2.complement();
  assert(N1 == N2.join(N3));

  Point p1(0.1,0.), p2(10.8,10.25), p3(-20.18,15.14);
  Point triangle[3] = { p1, p2, p3 };
  Nef_polyhedron N4(triangle, triangle+3);
  Nef_polyhedron N5 = N2.intersection(N4);
  assert(N5 <= N2 && N5 <= N4);

  return 0;
}
