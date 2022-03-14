#include <CGAL/Exact_integer.h>
#include <CGAL/Extended_homogeneous.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <cassert>

typedef CGAL::Exact_integer  NT;
typedef CGAL::Extended_homogeneous<NT>  Kernel;
typedef CGAL::Nef_polyhedron_3<Kernel>  Nef_polyhedron;
//typedef Nef_polyhedron::Plane_3  Plane_3;
typedef Kernel::Plane_3 Plane_3;
int main() {

  Nef_polyhedron N1(Plane_3(2,5,7,11), Nef_polyhedron::INCLUDED);
  Nef_polyhedron N2(Plane_3(2,5,7,11), Nef_polyhedron::EXCLUDED);

  assert(N1 >= N2);
  assert(N2 <= N1);
  assert(N1 != N2);
  assert(N1 > N2);
  assert(N2 < N1);

  N2 = N2.closure();
  assert(N1==N2);
  assert(N1>=N2);
  assert(N1<=N2);

  return 0;
}
