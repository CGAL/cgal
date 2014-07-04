#include <CGAL/Exact_integer.h>
#include <CGAL/Extended_homogeneous.h>
#include <CGAL/Nef_polyhedron_3.h>

typedef CGAL::Extended_homogeneous<CGAL::Exact_integer>  Kernel;
typedef CGAL::Nef_polyhedron_3<Kernel>  Nef_polyhedron;
typedef Kernel::Plane_3  Plane_3;

int main() {

  Nef_polyhedron N1(Plane_3( 1, 0, 0,-1));
  Nef_polyhedron N2(Plane_3(-1, 0, 0,-1));
  Nef_polyhedron N3(Plane_3( 0, 1, 0,-1));
  Nef_polyhedron N4(Plane_3( 0,-1, 0,-1));
  Nef_polyhedron N5(Plane_3( 0, 0, 1,-1));
  Nef_polyhedron N6(Plane_3( 0, 0,-1,-1));

  Nef_polyhedron I1(!N1 + !N2);  // open slice in yz-plane
  Nef_polyhedron I2(N3 - !N4);   // closed slice in xz-plane
  Nef_polyhedron I3(N5 ^ N6);    // open slice in yz-plane
  Nef_polyhedron Cube1(I2 * !I1);
  Cube1 *= !I3;
  Nef_polyhedron Cube2 = N1 * N2 * N3 * N4 * N5 * N6;

  CGAL_assertion(Cube1 == Cube2);  // both are closed cube
  CGAL_assertion(Cube1 == Cube1.closure());
  CGAL_assertion(Cube1 == Cube1.regularization());
  CGAL_assertion((N1 - N1.boundary()) == N1.interior());
  CGAL_assertion(I1.closure() == I1.complement().interior().complement());
  CGAL_assertion(I1.regularization() == I1.interior().closure());
  return 0;
}
