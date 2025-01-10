#include <CGAL/Exact_rational.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Hyperbolic_Delaunay_triangulation_traits_2.h>
#include <CGAL/Hyperbolic_surface_traits_2.h>
#include <CGAL/Hyperbolic_fundamental_domain_factory_2.h>
#include <CGAL/Delaunay_Triangulation_on_hyperbolic_surface_2.h>
#include <time.h>

typedef CGAL::Exact_rational                                         Rational;
typedef CGAL::Simple_cartesian<Rational>                             Kernel;
typedef CGAL::Hyperbolic_Delaunay_triangulation_traits_2<Kernel>     ParentTraits;
typedef CGAL::Hyperbolic_surface_traits_2<ParentTraits>              Traits;
typedef CGAL::Hyperbolic_fundamental_domain_2<Traits>                Domain;
typedef CGAL::Hyperbolic_fundamental_domain_factory_2<Traits>        Factory;
typedef CGAL::Triangulation_on_hyperbolic_surface_2<Traits>          Triangulation;
typedef CGAL::Delaunay_triangulation_on_hyperbolic_surface_2<Traits> Delaunay_triangulation;

int main(){
  // Generates the domain:
  Factory factory = Factory();
  Domain domain = factory.make_hyperbolic_fundamental_domain_g2(time(NULL));

  // Triangulates the domain:
  Triangulation triangulation = Triangulation(domain);

  // Saves the triangulation:
  std::ofstream  output_file = std::ofstream ("OutputTriangulation.txt");
  output_file << triangulation;
  output_file.close();

  // Prints the triangulation:
  std::cout << triangulation << std::endl;

  // Generates a Delaunay triangulation
  Delaunay_triangulation dt = Delaunay_triangulation(domain);

  triangulation.has_anchor();
  dt.has_anchor();
  dt.make_Delaunay();
 // Prints the triangulation:
  std::cout << dt << std::endl;

  return 0;
}
