#include <CGAL/Exact_rational.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Hyperbolic_Delaunay_triangulation_traits_2.h>
#include <CGAL/Hyperbolic_surface_traits_2.h>
#include <CGAL/Hyperbolic_fundamental_domain_2.h>
#include <CGAL/Hyperbolic_fundamental_domain_factory_2.h>
#include <CGAL/Hyperbolic_surface_triangulation_2.h>
#include <time.h>

using namespace CGAL;
typedef Simple_cartesian<Exact_rational>                       Kernel;
typedef Hyperbolic_Delaunay_triangulation_traits_2<Kernel>           ParentTraits;
typedef Hyperbolic_surface_traits_2<ParentTraits>                    Traits;
typedef Hyperbolic_fundamental_domain_2<Traits>                      Domain;
typedef Hyperbolic_fundamental_domain_factory_2<Traits>              Factory;
typedef Hyperbolic_surface_triangulation_2<Traits>                   Triangulation;

int main(){
  // Generates the domain:
  Factory factory = Factory(time(NULL));
  Domain domain = factory.generate_domain_g2();

  // Triangulates the domain:
  Triangulation triangulation = Triangulation(domain);

  // Applies the Delaunay flip algorithm to the triangulation:
  triangulation.make_delaunay();

  // Saves the triangulation:
  std::ofstream  output_file = std::ofstream ("OutputTriangulation.txt");
  output_file << triangulation;
  output_file.close();

  // Prints the triangulation:
  std::cout << triangulation << std::endl;

  return 0;
}
