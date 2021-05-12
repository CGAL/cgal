#include <CGAL/leda_integer.h>
#include <CGAL/Homogeneous.h>
#include <CGAL/Nef_polyhedron_S2.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Nef_polyhedron_iostream_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Nef_S2/Gausian_map.h>
#include <CGAL/Nef_S2/gausian_map_to_polyhedron_3.h>
#include <CGAL/Nef_S2/SM_point_locator.h>

#include <CGAL/convexity_check_3.h>

#include <fstream>
#include <sstream>

#define CGAL_NEF3_SPHERE_SWEEP_OPTIMIZATION_OFF

typedef leda_integer NT;
typedef CGAL::Homogeneous<NT> Kernel;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef CGAL::Nef_polyhedron_3<Kernel> Nef_polyhedron_3;
typedef CGAL::Gausian_map<Kernel> Gausian_map;
typedef Gausian_map::SM_const_decorator SM_decorator;
typedef SM_decorator::Sphere_point Sphere_point;
typedef CGAL::SM_point_locator<SM_decorator> SM_point_locator;

struct Plane_equation {
  template <class Facet>
  typename Facet::Plane_3 operator()( Facet& f) {
    typename Facet::Halfedge_handle h = f.halfedge();
    typedef typename Facet::Plane_3  Plane;
    return Plane( h->vertex()->point(),
                  h->next()->vertex()->point(),
                  h->next()->next()->vertex()->point());
  }
};

int main(int argc, char* argv[]) {

  CGAL_assertion(argc==2);
  std::ifstream in1(argv[1]);

  /*
  Polyhedron P1;
  in1 >> P1;

  std::transform( P1.facets_begin(), P1.facets_end(), P1.planes_begin(),
                  Plane_equation());

  CGAL_assertion(is_strongly_convex_3(P1));
  */

  Nef_polyhedron_3 N1;
  in1 >> N1;
  Gausian_map G1(N1);

  CGAL::Timer tboth, ttop, tbottom, ptop, pbottom;
  tboth.start();
  G1.locate_top_and_bottom();
  G1.get_top();
  G1.get_bottom();
  tboth.stop();

  ttop.start();
  G1.locate_top();
  ttop.stop();

  tbottom.start();
  G1.locate_bottom();
  tbottom.stop();

  SM_point_locator pl(G1.sphere_map());
  ptop.start();
  pl.locate(Sphere_point(0,0,1));
  ptop.stop();

  pbottom.start();
  pl.locate(Sphere_point(0,0,-1));
  pbottom.stop();

  std::cerr << "New method - both at once: " << tboth.time() << std::endl;
  std::cerr << "New method - top only: " << ttop.time() << std::endl;
  std::cerr << "New method - bottom only: " << tbottom.time() << std::endl;
  std::cerr << "Trivial point location - top only: " << ptop.time() << std::endl;
  std::cerr << "Trivial point location - bottom only: " << pbottom.time() << std::endl;

  // G1.visualize();
}
