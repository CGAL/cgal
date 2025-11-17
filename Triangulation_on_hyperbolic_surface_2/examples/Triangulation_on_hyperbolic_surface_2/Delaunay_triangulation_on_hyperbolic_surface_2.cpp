#include <CGAL/Exact_rational.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Circular_kernel_2.h>
#include <CGAL/Algebraic_kernel_for_circles_2_2.h>
#include <CGAL/Hyperbolic_Delaunay_triangulation_CK_traits_2.h>
#include <CGAL/Hyperbolic_surface_traits_2.h>

#include <CGAL/Delaunay_triangulation_on_hyperbolic_surface_2.h>
#include <CGAL/Triangulation_on_hyperbolic_surface_2_IO.h>
#include <CGAL/Hyperbolic_fundamental_domain_factory_2.h>

#include <time.h>

typedef CGAL::Exact_rational                                                                                            Rational;
typedef CGAL::Circular_kernel_2<CGAL::Simple_cartesian<Rational>,CGAL::Algebraic_kernel_for_circles_2_2<Rational>>		Kernel;
typedef CGAL::Hyperbolic_Delaunay_triangulation_CK_traits_2<Kernel>                                             		ParentTraits;
typedef CGAL::Hyperbolic_surface_traits_2<ParentTraits>                                                         		Traits;
typedef CGAL::Hyperbolic_fundamental_domain_2<Traits>                													Domain;
typedef CGAL::Hyperbolic_fundamental_domain_factory_2<Traits>        													Factory;
typedef CGAL::Delaunay_triangulation_on_hyperbolic_surface_2<Traits>          											Delaunay_triangulation;

int main() {
  // Generates the domain:
  Factory factory = Factory();
  Domain domain = factory.make_hyperbolic_fundamental_domain_g2(time(NULL)); // get a random seed with time(NULL)

  // Triangulates the domain:
  Delaunay_triangulation dt = Delaunay_triangulation(domain);

  // Saves the triangulation:
  std::ofstream output_file = std::ofstream ("OutputDelaunayTriangulation.txt");
  output_file << dt;
  output_file.close();

  // Prints the triangulation:
  std::cout << dt << std::endl;

  return 0;
}
