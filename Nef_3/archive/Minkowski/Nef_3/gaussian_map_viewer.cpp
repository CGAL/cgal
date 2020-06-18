#include <CGAL/leda_integer.h>
#include <CGAL/Homogeneous.h>
#include <CGAL/Nef_polyhedron_S2.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Nef_polyhedron_iostream_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Nef_S2/Gausian_map.h>
#include <CGAL/Nef_S2/gausian_map_to_polyhedron_3.h>

#include <CGAL/convexity_check_3.h>

#include <fstream>
#include <sstream>

#define CGAL_NEF3_SPHERE_SWEEP_OPTIMIZATION_OFF

typedef leda_integer NT;
typedef CGAL::Homogeneous<NT> Kernel;
typedef CGAL::Nef_polyhedron_S2<Kernel> Nef_polyhedron_S2;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef CGAL::Gausian_map<Kernel> Gausian_map;

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

  Polyhedron P1;
  in1 >> P1;

  std::transform( P1.facets_begin(), P1.facets_end(), P1.planes_begin(),
                  Plane_equation());

  CGAL_assertion(is_strongly_convex_3(P1));

  Gausian_map G1(P1);
  G1.visualize();
}
