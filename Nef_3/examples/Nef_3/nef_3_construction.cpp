#include <CGAL/Exact_integer.h>
#include <CGAL/Extended_homogeneous.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <cassert>

typedef CGAL::Extended_homogeneous<CGAL::Exact_integer>  Kernel;
typedef CGAL::Nef_polyhedron_3<Kernel>  Nef_polyhedron;
typedef Nef_polyhedron::Plane_3  Plane_3;

int main() {
  Nef_polyhedron N0;
  Nef_polyhedron N1(Nef_polyhedron::EMPTY);
  Nef_polyhedron N2(Nef_polyhedron::COMPLETE);
  Nef_polyhedron N3(Plane_3( 1, 2, 5,-1));
  Nef_polyhedron N4(Plane_3( 1, 2, 5,-1), Nef_polyhedron::INCLUDED);
  Nef_polyhedron N5(Plane_3( 1, 2, 5,-1), Nef_polyhedron::EXCLUDED);

  assert(N0 == N1);
  assert(N3 == N4);
  assert(N0 != N2);
  assert(N3 != N5);

  assert(N4 >= N5);
  assert(N5 <= N4);
  assert(N4 > N5);
  assert(N5 < N4);

  N5 = N5.closure();
  assert(N4 >= N5);
  assert(N4 <= N5);

  return 0;
}
