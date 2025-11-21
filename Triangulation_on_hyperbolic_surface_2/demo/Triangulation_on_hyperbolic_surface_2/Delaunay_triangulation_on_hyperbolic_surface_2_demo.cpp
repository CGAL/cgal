#include "window.h"

#include <CGAL/Gmpq.h>
#include <CGAL/Circular_kernel_2.h>
#include <CGAL/Algebraic_kernel_for_circles_2_2.h>
#include <CGAL/Hyperbolic_Delaunay_triangulation_CK_traits_2.h>
#include <CGAL/Hyperbolic_surface_traits_2.h>
#include <CGAL/Hyperbolic_fundamental_domain_factory_2.h>
#include <CGAL/Triangulation_on_hyperbolic_surface_2.h>

#include <CGAL/Timer.h>

using namespace CGAL;

typedef CGAL::Gmpq                                                      NumberType;
typedef CGAL::Circular_kernel_2<CGAL::Simple_cartesian<NumberType>,CGAL::Algebraic_kernel_for_circles_2_2<NumberType>> Kernel;
typedef CGAL::Hyperbolic_Delaunay_triangulation_CK_traits_2<Kernel>     ParentTraits;
typedef Hyperbolic_surface_traits_2<ParentTraits>                       Traits;

typedef CGAL::Hyperbolic_fundamental_domain_2<Traits>                                                           Domain;
typedef CGAL::Hyperbolic_isometry_2<Traits>                                                                     Isometry;
typedef CGAL::Hyperbolic_fundamental_domain_factory_2<Traits>                                                   Factory;
typedef CGAL::Delaunay_triangulation_on_hyperbolic_surface_2<Traits>                                            Delaunay_triangulation;
typedef typename Delaunay_triangulation::Anchor                                                                 Anchor;

/*
  HOW TO USE THIS DEMO
  ./Delaunay_triangulation_on_hyperbolic_surface_2_demo [epsilon] [surface seed] [precision]
  With no arguments, uses the default values defined below.
*/

// DEFAULT VALUES
double epsilon = 0.25;
int seed = time(NULL);
int p = 1;

int main(int argc, char *argv[])
{
  // 1. Parse args and generate the input
  if (argc > 1) {
    epsilon = std::stod(argv[1]);
  }

  Domain domain;
  if (argc <= 2) {
    std::cout << "Using random seed " << seed << std::endl;
  } else {
    seed = atoi(argv[2]);
  }
  Factory factory;
  std::cout << "Generating surface with seed " << seed << "..." << std::endl;
  domain = factory.make_hyperbolic_fundamental_domain_g2(seed);
  Delaunay_triangulation dt = Delaunay_triangulation(domain);

  if (argc > 3) {
    p = atoi(argv[3]);
  }

  // 2. Get a vertex
  // So that if you run the demo on a same surface but with different values of epsilon,
  // the drawing will be centered at the same vertex and it will look similar.
  Point v0 = dt.anchor().vertices[0];

  // 3. Compute epsilon-net and display useful info
  if constexpr(!std::is_same<NumberType, CGAL::Gmpq>::value) {
    std::cout << "WARNING: Not using the CGAL::Gmpq number type. Precision will be ignored and to_double approximation will be used instead." << std::endl;
  }
  std::cout << "Computing a " << epsilon << "-net with floating-point precision " << p*53 << "..." << std::endl;
  CGAL::Timer timer;
  timer.start();
  std::cout << "Is epsilon-net? " << dt.epsilon_net(epsilon, p) << std::endl;
  timer.stop();
  std::cout << "Done in " << timer.time() << " seconds." << std::endl;
  dt.combinatorial_map().display_characteristics(std::cout) << std::endl;

  // 4. SET THE FIRST ANCHOR OF THE DRAWING
  Anchor anchor = dt.locate(v0);
  int index = 0;
  for (int i = 0; i < 3; i++) {
    if (v0 == anchor.vertices[i]) {
      index = i;
    }
  }

  Anchor start = Anchor();
  start.dart = anchor.dart;
  for (int i = 0; i < 3; i++) {
    start.vertices[i] = anchor.vertices[(i + index) % 3];
    if (i < index) {
      start.dart = dt.Base::ccw(start.dart);
    }
  }

  // 5. DRAW
  QApplication app(argc, argv);
  app.setApplicationName("Delaunay triangulation on hyperbolic surface 2 Demo");
  DemoWindow window;
  window.item().draw_triangulation(dt, start);
  window.show();
  QStringList args = app.arguments();
  args.removeAt(0);
  return app.exec();
}
