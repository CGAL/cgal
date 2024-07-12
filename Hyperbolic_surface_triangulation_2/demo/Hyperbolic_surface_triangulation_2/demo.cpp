#include "window.h"

#include <CGAL/Gmpq.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Hyperbolic_Delaunay_triangulation_traits_2.h>
#include <CGAL/Hyperbolic_surface_traits_2.h>
#include <CGAL/Hyperbolic_fundamental_domain_2.h>
#include <CGAL/Hyperbolic_fundamental_domain_factory_2.h>
#include <CGAL/Hyperbolic_surface_triangulation_2.h>

using namespace CGAL;

typedef Cartesian<Gmpq>                                                 Kernel;
typedef Hyperbolic_Delaunay_triangulation_traits_2<Kernel>              ParentTraits;
typedef Hyperbolic_surface_traits_2<ParentTraits>                      Traits;
typedef Hyperbolic_fundamental_domain_2<Traits>                         Domain;
typedef Hyperbolic_fundamental_domain_factory_2<Traits>                 Factory;
typedef Hyperbolic_surface_triangulation_2<Traits>                      Triangulation;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv){
  // 1. Generate the triangulation
  Factory factory = Factory(time(NULL));
  Domain domain = factory.generate_domain_g2();
  Triangulation triangulation = Triangulation(domain);
  triangulation.make_delaunay();

  // 2. Draw the triangulation
  QApplication app(argc, argv);
  app.setApplicationName("Hyperbolic surfaces triangulation 2 Demo");

  DemoWindow window;
  window.item().draw_triangulation(triangulation);
  window.show();

  QStringList args = app.arguments();
  args.removeAt(0);
  return app.exec();
}
