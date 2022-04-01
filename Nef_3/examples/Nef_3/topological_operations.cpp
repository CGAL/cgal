#include <CGAL/Exact_integer.h>
#include <CGAL/Homogeneous.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/IO/Nef_polyhedron_iostream_3.h>
#include <cassert>

typedef CGAL::Exact_integer  NT;
typedef CGAL::Homogeneous<NT>  Kernel;
typedef CGAL::Nef_polyhedron_3<Kernel>  Nef_polyhedron;

int main() {

  Nef_polyhedron N;
  std::cin >> N;

  assert((N - N.boundary()) == N.interior());
  assert(N.closure() == N.complement().interior().complement());
  assert(N.regularization() == N.interior().closure());

  return 0;
}
