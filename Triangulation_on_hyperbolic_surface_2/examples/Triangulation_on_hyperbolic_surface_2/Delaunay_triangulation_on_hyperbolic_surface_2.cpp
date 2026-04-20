#include <CGAL/Exact_rational.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Circular_kernel_2.h>
#include <CGAL/Algebraic_kernel_for_circles_2_2.h>
#include <CGAL/Hyperbolic_Delaunay_triangulation_CK_traits_2.h>
#include <CGAL/Hyperbolic_surface_traits_2.h>
#include <CGAL/Delaunay_triangulation_on_hyperbolic_surface_2.h>
#include <CGAL/Triangulation_on_hyperbolic_surface_2_IO.h>
#include <CGAL/Hyperbolic_fundamental_domain_factory_2.h>

typedef CGAL::Gmpq    Rational;
typedef CGAL::Circular_kernel_2<CGAL::Simple_cartesian<Rational>,CGAL::Algebraic_kernel_for_circles_2_2<Rational>>    Kernel;
typedef CGAL::Hyperbolic_Delaunay_triangulation_CK_traits_2<Kernel>                                                 ParentTraits;
typedef CGAL::Hyperbolic_surface_traits_2<ParentTraits>                                                             Traits;
typedef CGAL::Hyperbolic_fundamental_domain_2<Traits>                                         Domain;
typedef CGAL::Hyperbolic_fundamental_domain_factory_2<Traits>                                 Factory;
typedef CGAL::Delaunay_triangulation_on_hyperbolic_surface_2<Traits>                                Delaunay_triangulation;

int main() {
  // Generates the domain:
  Factory factory = Factory();
  Domain domain = factory.make_hyperbolic_fundamental_domain_g2(time(NULL)); // get a random seed with time(NULL)

  // Triangulates the domain:
  Delaunay_triangulation dt = Delaunay_triangulation(domain);

  // Compute epsilon-net and display useful info
  if constexpr(!std::is_same<Rational, CGAL::Gmpq>::value) {
      std::cout << "WARNING: Not using the CGAL::Gmpq number type. Precision will be ignored and to_double approximation will be used instead." << std::endl;
    }
  std::cout << "Computing a " << 0.25 << "-net with floating-point precision " << 2*53 << "..." << std::endl;
  std::cout << "Is epsilon-net? " << dt.epsilon_net(0.25, 2) << std::endl;
  std::cout << "Done " << std::endl;
  dt.combinatorial_map().display_characteristics(std::cout) << std::endl;

  // Locate a point in the triangulation
  Delaunay_triangulation::Point query(1/2,1/2);
  Delaunay_triangulation::Locate_type lt;
  unsigned li=0;
  unsigned ld=0;
  dt.locate(query, lt, li, ld, dt.anchor());
  std::cout << "The number of triangles traversed by the walk used for the point location algorithm is " << ld << std::endl;

  // Saves the triangulation:
  std::ofstream output_file = std::ofstream ("OutputDelaunayTriangulation.txt");
  output_file << dt;
  output_file.close();

  return 0;
}
