#include <CGAL/Exact_integer.h>
#include <CGAL/Extended_homogeneous.h>
#include <CGAL/Nef_polyhedron_3.h>

typedef CGAL::Exact_integer  NT;
typedef CGAL::Extended_homogeneous<NT>  Kernel;
typedef CGAL::Nef_polyhedron_3<Kernel>  Nef_polyhedron;
typedef Nef_polyhedron::Plane_3  Plane_3;

int main() {

  Nef_polyhedron N1(Plane_3( 1, 0, 0,-1),Nef_polyhedron::INCLUDED);
  Nef_polyhedron N2(Plane_3(-1, 0, 0,-1),Nef_polyhedron::INCLUDED);
  Nef_polyhedron N3(Plane_3( 0, 1, 0,-1),Nef_polyhedron::INCLUDED);
  Nef_polyhedron N4(Plane_3( 0,-1, 0,-1),Nef_polyhedron::INCLUDED);
  Nef_polyhedron N5(Plane_3( 0, 0, 1,-1),Nef_polyhedron::INCLUDED);
  Nef_polyhedron N6(Plane_3( 0, 0,-1,-1),Nef_polyhedron::INCLUDED);

  Nef_polyhedron I1(!N1+!N2);
  Nef_polyhedron I2(N3-!N4);
  Nef_polyhedron I3(N5^N6);
  Nef_polyhedron Cube1(I2 *!I1);
  Cube1 *= !I3;
  Nef_polyhedron Cube2 = N1 * N2 * N3 * N4 * N5 * N6;
  CGAL_assertion (Cube1 == Cube2);
  return 0;
}
