// examples/Nef_3/visualization_SM.C
// -------------------------------------
#include <CGAL/Gmpz.h>
#include <CGAL/Homogeneous.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/IO/Nef_polyhedron_iostream_3.h>
#include <CGAL/IO/Nef_polyhedron_S2_OGLUT_stream.h>

typedef CGAL::Gmpz RT;
typedef CGAL::Homogeneous<RT> Kernel;
typedef CGAL::Nef_polyhedron_3<Kernel> Nef_polyhedron_3;
typedef Nef_polyhedron_3::Vertex_const_iterator Vertex_const_iterator;
typedef Nef_polyhedron_3::Nef_polyhedron_S2 Nef_polyhedron_S2;


int main() {

  Nef_polyhedron_3 N;
  std::cin >> N;

  Vertex_const_iterator v = N.vertices_begin();
  Nef_polyhedron_S2 S(N.get_sphere_map(v));
  CGAL::ogl << S;
  CGAL::ogl << "Visualization example";
  CGAL::ogl.display();;

  return 0;
}
