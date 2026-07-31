#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Hyperbolic_Delaunay_triangulation_traits_2.h>
#include <CGAL/Hyperbolic_surface_traits_2.h>
#include <CGAL/Hyperbolic_fundamental_domain_factory_2.h>
#include <CGAL/Triangulation_on_hyperbolic_surface_2.h>
#include <CGAL/Triangulation_on_hyperbolic_surface_2_IO.h>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel          Kernel;
typedef CGAL::Hyperbolic_Delaunay_triangulation_traits_2<Kernel>     ParentTraits;
typedef CGAL::Hyperbolic_surface_traits_2<ParentTraits>              Traits;
typedef CGAL::Hyperbolic_fundamental_domain_2<Traits>                Domain;
typedef CGAL::Hyperbolic_fundamental_domain_factory_2<Traits>        Factory;
typedef CGAL::Triangulation_on_hyperbolic_surface_2<Traits>          Triangulation;

int main() {

  // Generates the domain:
  Factory factory = Factory();
  Domain domain = factory.make_hyperbolic_fundamental_domain_g2(time(NULL)); // get a random

  // Triangulates the domain:
  Triangulation triangulation = Triangulation(domain);

  // Applies the Delaunay flip algorithm to the triangulation:
  triangulation.make_Delaunay();
  triangulation.combinatorial_map().display_characteristics(std::cout) << std::endl;

  // Saves the triangulation:
  std::ofstream  output_file = std::ofstream ("OutputTriangulation.txt");
  output_file << triangulation;
  output_file.close();

  // Prints the triangulation:
  std::cout << triangulation << std::endl;

  return EXIT_SUCCESS;
}
