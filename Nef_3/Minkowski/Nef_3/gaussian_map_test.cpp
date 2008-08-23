#include <CGAL/leda_integer.h>
#include <CGAL/Homogeneous.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Nef_polyhedron_iostream_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Nef_S2/Gausian_map.h>
#include <CGAL/Nef_S2/gausian_map_to_polyhedron_3.h>
#include <CGAL/Nef_S2/gausian_map_to_nef_3.h>

#include <CGAL/convexity_check_3.h>

#include <CGAL/IO/Qt_widget_Nef_3.h>
#include <qapplication.h>

#include <fstream>
#include <sstream>

#define CGAL_NEF3_SPHERE_SWEEP_OPTIMIZATION_OFF

typedef leda_integer NT;
typedef CGAL::Homogeneous<NT> Kernel;
typedef CGAL::Nef_polyhedron_3<Kernel> Nef_polyhedron;
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

  CGAL_assertion(argc==3);
  std::ifstream in1(argv[1]);
  std::ifstream in2(argv[2]);
  
  Polyhedron P1,P2;
  in1 >> P1;
  in2 >> P2;

  std::transform( P1.facets_begin(), P1.facets_end(), P1.planes_begin(),
		  Plane_equation());

  std::transform( P2.facets_begin(), P2.facets_end(), P2.planes_begin(),
		  Plane_equation());

  CGAL_assertion(is_strongly_convex_3(P1));
  CGAL_assertion(is_strongly_convex_3(P2));

  double t0,t1,t2;
  CGAL::Timer t;
  t.start();

  // create Gausian map of P1
  Gausian_map G1(P1);
  G1.dump();

  G1.dump();
  CGAL::gausian_map_to_nef_3<Kernel, Nef_polyhedron::Items, Nef_polyhedron::Mark> Converter1(G1);
  Nef_polyhedron N1;
  N1.delegate(Converter1,true);

  t0 = t.time();

  // create Gausian map of P2
  Gausian_map G2(P2);
  G2.dump();

  t1 = t.time();

  Gausian_map G;
  G.minkowski_sum(G1,G2);

  t2 = t.time();

  /*
  // conversion from gausian_map to Polyhedron_3
  gausian_map_to_polyhedron_3<Kernel, Polyhedron::HDS> Converter(G);
  Polyhedron P;
  P.delegate(Converter);
  */

  G.dump();

  // conversion from gausian_map to Nef_3
  CGAL::gausian_map_to_nef_3<Kernel, Nef_polyhedron::Items, Nef_polyhedron::Mark> Converter(G);
  Nef_polyhedron N;
  N.delegate(Converter,true);

  std::cerr << N;

  t.stop();
  std::cerr << "Runtime Minkowski Sum: " << t0 << std::endl;
  std::cerr << "Runtime Minkowski Sum: " << t1 << std::endl;
  std::cerr << "Runtime Minkowski Sum: " << t2 << std::endl;
  std::cerr << "Runtime Minkowski Sum: " << t.time() << std::endl;

  // Visualization of result
  //  Nef_polyhedron N(P);  
  QApplication a(argc, argv);
  CGAL::Qt_widget_Nef_3<Nef_polyhedron>* w = 
    new CGAL::Qt_widget_Nef_3<Nef_polyhedron>(N);
  a.setMainWidget(w);
  w->show();
  a.exec();
}
