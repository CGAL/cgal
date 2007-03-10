#include <CGAL/Gmpz.h>
#include <CGAL/Homogeneous.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/IO/Nef_polyhedron_iostream_3.h>

typedef CGAL::Homogeneous<CGAL::Gmpz> Kernel;
typedef CGAL::Nef_polyhedron_3<Kernel> Nef_polyhedron_3;


int main() {

  // We've put the typedefs here as VC7 gives us an ICE if they are global typedefs
  typedef Nef_polyhedron_3::Vertex_const_iterator Vertex_const_iterator;
  typedef Nef_polyhedron_3::Nef_polyhedron_S2 Nef_polyhedron_S2;
  typedef Nef_polyhedron_S2::SVertex_const_handle SVertex_const_handle;
  typedef Nef_polyhedron_S2::SHalfedge_const_handle SHalfedge_const_handle;
  typedef Nef_polyhedron_S2::SHalfloop_const_handle SHalfloop_const_handle;
  typedef Nef_polyhedron_S2::SFace_const_iterator SFace_const_iterator;
  typedef Nef_polyhedron_S2::SFace_cycle_const_iterator
    SFace_cycle_const_iterator;

  Nef_polyhedron_3 N;
  std::cin >> N;

  Vertex_const_iterator v = N.vertices_begin();
  Nef_polyhedron_S2 S(N.get_sphere_map(v));

  int i=0;
  SFace_const_iterator sf;
  for(sf = S.sfaces_begin(); sf != S.sfaces_end(); sf++) {
    SFace_cycle_const_iterator it;
    std::cout << "the sface cycles of sface " << i++ << " start with an\n";
    for(it = sf->sface_cycles_begin(); it != sf->sface_cycles_end(); it++) {
      if (it.is_svertex())
        std::cout << "  svertex at position "
                  << SVertex_const_handle(it)->point() << std::endl;
      else if (it.is_shalfedge())
        std::cout << "  shalfedge from "
                  << SHalfedge_const_handle(it)->source()->point() << " to "
                  << SHalfedge_const_handle(it)->target()->point() << std::endl;
      else if (it.is_shalfloop())
        std::cout << "  shalfloop lying in the plane "
                  << SHalfloop_const_handle(it)->circle() << std::endl;
      // other cases can not occur.
    }
  }

  return 0;
}
