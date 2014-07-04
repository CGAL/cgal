#include <CGAL/Exact_integer.h>
#include <CGAL/Filtered_extended_homogeneous.h>
#include <CGAL/Nef_polyhedron_2.h>

typedef CGAL::Exact_integer RT;
typedef CGAL::Filtered_extended_homogeneous<RT> Extended_kernel;
typedef CGAL::Nef_polyhedron_2<Extended_kernel> Nef_polyhedron;
typedef Nef_polyhedron::Point Point;
typedef Nef_polyhedron::Line  Line;

int main() {

  Nef_polyhedron N1(Nef_polyhedron::COMPLETE);

  Line l(2,4,2); // l : 2x + 4y + 2 = 0
  Nef_polyhedron N2(l,Nef_polyhedron::INCLUDED);
  Nef_polyhedron N3 = N2.complement();
  CGAL_assertion(N1 == N2.join(N3));

  Point p1(0,0), p2(10,10), p3(-20,15);
  Point triangle[3] = { p1, p2, p3 };
  Nef_polyhedron N4(triangle, triangle+3);
  Nef_polyhedron N5 = N2.intersection(N4);
  CGAL_assertion(N5 <= N2 && N5 <= N4);

  return 0;
}
