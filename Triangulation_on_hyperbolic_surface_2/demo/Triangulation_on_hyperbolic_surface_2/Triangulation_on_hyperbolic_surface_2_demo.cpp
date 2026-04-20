#include "window.h"

#include <CGAL/Exact_rational.h>
#include <CGAL/Gmpq.h>
#include <CGAL/Simple_cartesian.h>

#include <CGAL/Hyperbolic_Delaunay_triangulation_CK_traits_2.h>
#include <CGAL/Hyperbolic_surface_traits_2.h>
#include <CGAL/Hyperbolic_fundamental_domain_factory_2.h>
#include <CGAL/Triangulation_on_hyperbolic_surface_2.h>

using namespace CGAL;

typedef Gmpq                                                      NumberType;
typedef Circular_kernel_2<Simple_cartesian<NumberType>, Algebraic_kernel_for_circles_2_2<NumberType>> Kernel;
typedef Hyperbolic_Delaunay_triangulation_CK_traits_2<Kernel>     ParentTraits;

// typedef Simple_cartesian<NumberType>                                Kernel;
// typedef Hyperbolic_Delaunay_triangulation_traits_2<Kernel>              ParentTraits;
typedef Hyperbolic_surface_traits_2<ParentTraits>                       Traits;
typedef Hyperbolic_fundamental_domain_2<Traits>                         Domain;
typedef Hyperbolic_fundamental_domain_factory_2<Traits>                 Factory;
typedef Triangulation_on_hyperbolic_surface_2<Traits>                   Triangulation;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv)
{
  // 1. Generate the triangulation
  Factory factory;
  Domain domain = factory.make_hyperbolic_fundamental_domain_g2(time(NULL));
  Triangulation triangulation = Triangulation(domain);
  triangulation.make_Delaunay();

  // 2. Draw the triangulation
  QApplication app(argc, argv);
  app.setApplicationName("Triangulation on hyperbolic surface 2 Demo");

  DemoWindow window;
  window.item().draw_triangulation(triangulation);
  window.show();

  QStringList args = app.arguments();
  args.removeAt(0);
  return app.exec();
}
