#include <CGAL/Exact_integer.h>
#include <CGAL/Homogeneous.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <cassert>

typedef CGAL::Homogeneous<CGAL::Exact_integer>  Kernel;
typedef CGAL::Nef_polyhedron_3<Kernel> Nef_polyhedron;

int main() {
  Nef_polyhedron N0(Nef_polyhedron::EMPTY);
  Nef_polyhedron N1(Nef_polyhedron::COMPLETE);

  assert(N0 == N1.complement());
  assert(N0 != N1);
  return 0;
}
