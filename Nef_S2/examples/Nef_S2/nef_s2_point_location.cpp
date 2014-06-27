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
typedef Nef_polyhedron_S2::SFace_const_handle SFace_const_handle;
typedef Nef_polyhedron_S2::Object_handle Object_handle;
typedef Nef_polyhedron_S2::Sphere_point Sphere_point;

int main() {

  Nef_polyhedron_S2 S;
  CGAL::create_random_Nef_S2(S,5);

  SVertex_const_handle sv;
  SHalfedge_const_handle se;
  SHalfloop_const_handle sl;
  SFace_const_handle sf;
  Object_handle o = S.locate(Sphere_point(1,0,0));
  if(CGAL::assign(sv,o))
    std::cout << "Locating svertex" << std::endl;
  else if(CGAL::assign(se,o))
    std::cout << "Locating shalfedge" << std::endl;
  else if(CGAL::assign(sl,o))
    std::cout << "Locating shalfloop" << std::endl;
  else if(CGAL::assign(sf,o))
    std::cout << "Locating sface" << std::endl;
  else {
    std::cout << "something wrong" << std::endl;
    return 1;
  }
  return 0;
}
