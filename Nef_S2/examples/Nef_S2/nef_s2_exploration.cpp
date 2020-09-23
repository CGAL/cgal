#include <CGAL/Exact_rational.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Nef_polyhedron_S2.h>
#include <CGAL/Nef_S2/create_random_Nef_S2.h>

typedef CGAL::Exact_rational FT;
typedef CGAL::Cartesian<FT> Kernel;
typedef CGAL::Nef_polyhedron_S2<Kernel> Nef_polyhedron_S2;
typedef Nef_polyhedron_S2::SVertex_const_handle SVertex_const_handle;
typedef Nef_polyhedron_S2::SHalfedge_const_handle SHalfedge_const_handle;
typedef Nef_polyhedron_S2::SHalfloop_const_handle SHalfloop_const_handle;
typedef Nef_polyhedron_S2::SFace_const_iterator SFace_const_iterator;
typedef Nef_polyhedron_S2::SFace_cycle_const_iterator
                           SFace_cycle_const_iterator;

int main() {

  Nef_polyhedron_S2 S;
  CGAL::create_random_Nef_S2(S,5);

  int i=0;
  SFace_const_iterator sf;
  CGAL_forall_sfaces(sf,S) {
    SFace_cycle_const_iterator it;
    std::cout << "the sface cycles of sface " << i++;
    std::cout << " start with an " << std::endl;
    CGAL_forall_sface_cycles_of(it,sf) {
      if (it.is_svertex()) {
        std::cout << "  svertex at position ";
        std::cout << SVertex_const_handle(it)->point() << std::endl;
      }
      else if (it.is_shalfedge()) {
        std::cout << "  shalfedge from ";
        std::cout << SHalfedge_const_handle(it)->source()->point() << " to ";
        std::cout << SHalfedge_const_handle(it)->target()->point() << std::endl;
      }
      else if (it.is_shalfloop()) {
        std::cout << "  shalfloop lying in the plane ";
        std::cout << SHalfloop_const_handle(it)->circle() << std::endl;
      }
      else
        std::cout << "something is wrong" << std::endl;
    }
  }
  return 0;
}
