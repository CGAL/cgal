// examples/Nef_3/visualization_SNC.C
// -------------------------------------
#include <CGAL/Gmpz.h>
#include <CGAL/Homogeneous.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/Nef_3/Visualizor.h>
#include <CGAL/IO/Nef_polyhedron_iostream_3.h>

typedef CGAL::Gmpz RT;
typedef CGAL::Homogeneous<RT> Kernel;
typedef CGAL::Nef_polyhedron_3<Kernel> Nef_polyhedron_3;
typedef CGAL::Visualizor<Nef_polyhedron_3> Visualizor;

int main() {

  Nef_polyhedron_3 N;
  std::cin >> N;

  Visualizor V(N);
  V.draw();
  CGAL::OGL::start_viewer();

  return 0;
}
