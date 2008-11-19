#include <CGAL/Gmpz.h>
#include <CGAL/Extended_homogeneous.h>
#include <CGAL/Nef_polyhedron_3.h>

typedef CGAL::Extended_homogeneous<CGAL::Gmpz>  Kernel;
typedef CGAL::Nef_polyhedron_3<Kernel>  Nef_polyhedron;
typedef Nef_polyhedron::Plane_3  Plane_3;

int main() {
  Nef_polyhedron N0;
  Nef_polyhedron N1(Nef_polyhedron::EMPTY);
  Nef_polyhedron N2(Nef_polyhedron::COMPLETE);
  Nef_polyhedron N3(Plane_3( 1, 2, 5,-1));
  Nef_polyhedron N4(Plane_3( 1, 2, 5,-1), Nef_polyhedron::INCLUDED);
  Nef_polyhedron N5(Plane_3( 1, 2, 5,-1), Nef_polyhedron::EXCLUDED);

  CGAL_assertion(N0 == N1);
  CGAL_assertion(N3 == N4);
  CGAL_assertion(N0 != N2);
  CGAL_assertion(N3 != N5);

  CGAL_assertion(N4 >= N5);
  CGAL_assertion(N5 <= N4);
  CGAL_assertion(N4 > N5);
  CGAL_assertion(N5 < N4);

  N5 = N5.closure();
  CGAL_assertion(N4 >= N5);
  CGAL_assertion(N4 <= N5);

  return 0;
}
