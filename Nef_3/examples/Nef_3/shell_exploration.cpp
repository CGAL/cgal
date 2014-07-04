#include <CGAL/Exact_integer.h>
#include <CGAL/Homogeneous.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/IO/Nef_polyhedron_iostream_3.h>

typedef CGAL::Homogeneous<CGAL::Exact_integer> Kernel;
typedef CGAL::Nef_polyhedron_3<Kernel> Nef_polyhedron;
typedef Nef_polyhedron::Vertex_const_handle Vertex_const_handle;
typedef Nef_polyhedron::Halfedge_const_handle Halfedge_const_handle;
typedef Nef_polyhedron::Halffacet_const_handle Halffacet_const_handle;
typedef Nef_polyhedron::SHalfedge_const_handle SHalfedge_const_handle;
typedef Nef_polyhedron::SHalfloop_const_handle SHalfloop_const_handle;
typedef Nef_polyhedron::SFace_const_handle SFace_const_handle;
typedef Nef_polyhedron::Volume_const_iterator Volume_const_iterator;
typedef Nef_polyhedron::Shell_entry_const_iterator Shell_entry_const_iterator;
typedef Kernel::Point_3 Point_3;

class Shell_explorer {
  bool first;
  Vertex_const_handle v_min;

public:
  Shell_explorer()
    : first(true) {}

  void visit(Vertex_const_handle v) {
    if(first || CGAL::lexicographically_xyz_smaller(v->point(),v_min->point())) {
      v_min = v;
      first=false;
    }
  }

  void visit(Halfedge_const_handle ) {}
  void visit(Halffacet_const_handle ) {}
  void visit(SHalfedge_const_handle ) {}
  void visit(SHalfloop_const_handle ) {}
  void visit(SFace_const_handle ) {}

  Vertex_const_handle& minimal_vertex() { return v_min; }
  void reset_minimal_vertex() { first = true; }
};

int main() {
  Nef_polyhedron N;
  std::cin >> N;

  int ic = 0;
  Volume_const_iterator c;
  Shell_explorer SE;
  CGAL_forall_volumes(c,N) {
    std::cout << "Volume " << ic++ << std::endl;
    int is = 0;
    Shell_entry_const_iterator it;
    CGAL_forall_shells_of(it,c) {
      SE.reset_minimal_vertex();
      N.visit_shell_objects(SFace_const_handle(it),SE);
      Point_3 p(SE.minimal_vertex()->point());
      std::cout << "  minimal vertex of shell " << is++
                << " is at " << p << std::endl;
    }
  }
}
