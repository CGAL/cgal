#include <CGAL/Exact_integer.h>
#include <CGAL/Homogeneous.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/IO/Nef_polyhedron_iostream_3.h>

typedef CGAL::Homogeneous<CGAL::Exact_integer> Kernel;
typedef CGAL::Nef_polyhedron_3<Kernel> Nef_polyhedron_3;
  typedef Kernel::Point_3 Point_3;

int main() {
  //We've put the typedefs here as VC7 gives us an ICE if they are global typedefs
  typedef Nef_polyhedron_3::Vertex_const_handle Vertex_const_handle;
  typedef Nef_polyhedron_3::Halfedge_const_handle Halfedge_const_handle;
  typedef Nef_polyhedron_3::Halffacet_const_handle Halffacet_const_handle;
  typedef Nef_polyhedron_3::Volume_const_handle Volume_const_handle;
  typedef Nef_polyhedron_3::Object_handle Object_handle;


  Nef_polyhedron_3 N;
  std::cin >> N;

  Vertex_const_handle v;
  Halfedge_const_handle e;
  Halffacet_const_handle f;
  Volume_const_handle c;
  Object_handle o = N.locate(Point_3(0,0,0));
  if(CGAL::assign(v,o))
    std::cout << "Locating vertex" << std::endl;
  else if(CGAL::assign(e,o))
    std::cout << "Locating edge" << std::endl;
  else if(CGAL::assign(f,o))
    std::cout << "Locating facet" << std::endl;
  else if(CGAL::assign(c,o))
    std::cout << "Locating volume" << std::endl;
  //other cases can not occur

  return 0;
}
