#include <CGAL/Gmpz.h>
#include <CGAL/Extended_homogeneous.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/IO/Nef_polyhedron_iostream_3.h>

typedef CGAL::Gmpz  NT;
//instead of
//typedef CGAL::Extended_homogeneous<NT>  Kernel;
// workaround for VC++
struct Kernel : public CGAL::Extended_homogeneous<NT> {};

typedef CGAL::Nef_polyhedron_3<Kernel>  Nef_polyhedron;
typedef Nef_polyhedron::RT  RT;
typedef Nef_polyhedron::Point_3  Point_3;
typedef Nef_polyhedron::Plane_3  Plane_3;
typedef Nef_polyhedron::Vertex_const_iterator  Vertex_const_iterator;

int main() {

  Nef_polyhedron N;
  std::cin >> N;

  Vertex_const_iterator v;
  for(v = N.vertices_begin(); v != N.vertices_end(); ++v) {
    Point_3 p(v->point());
    if(p.hx().degree() > 0 || p.hy().degree() > 0 || p.hz().degree() > 0)
      std::cout << "extended vertex at " << p << std::endl;
    else
      std::cout << "standard vertex at " << p << std::endl;

    if(p == Point_3(RT(0,1), RT(0,1), RT(0,1)))
       std::cout << "  found vertex (right,back,top) of the infimaximal box"
                 << std::endl;
  }

  return 0;
}
