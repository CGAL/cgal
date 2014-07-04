#include <CGAL/Exact_integer.h>
#include <CGAL/Homogeneous.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/IO/Nef_polyhedron_iostream_3.h>

typedef CGAL::Exact_integer  NT;
typedef CGAL::Homogeneous<NT>  Kernel;
typedef CGAL::Nef_polyhedron_3<Kernel>  Nef_polyhedron;

int main() {

  Nef_polyhedron N;
  std::cin >> N;

  CGAL_assertion((N - N.boundary()) == N.interior());
  CGAL_assertion(N.closure() == N.complement().interior().complement());
  CGAL_assertion(N.regularization() == N.interior().closure());

  return 0;
}
