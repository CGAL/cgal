#include <CGAL/Gmpz.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Nef_polyhedron_3.h>

typedef CGAL::Gmpz  NT;
typedef CGAL::Cartesian<NT>  Kernel;
typedef CGAL::Nef_polyhedron_3<Kernel> Nef_polyhedron;

int main() {

  Nef_polyhedron N0(Nef_polyhedron::EMPTY);
  Nef_polyhedron N1(Nef_polyhedron::COMPLETE);

  N1 = N1.complement();
  CGAL_assertion (N0 == N1);
  return 0;
}
