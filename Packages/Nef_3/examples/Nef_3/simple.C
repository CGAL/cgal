// file: examples/Nef_3/simple.C

#include <CGAL/Gmpq.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Nef_polyhedron_3.h>

typedef CGAL::Cartesian<CGAL::Gmpq>  Kernel;
typedef CGAL::Nef_polyhedron_3<Kernel> Nef_polyhedron;

int main() {
  Nef_polyhedron N0(Nef_polyhedron::EMPTY);
  Nef_polyhedron N1(Nef_polyhedron::COMPLETE);

  CGAL_assertion (N0 == N1.complement());
  return 0;
}
