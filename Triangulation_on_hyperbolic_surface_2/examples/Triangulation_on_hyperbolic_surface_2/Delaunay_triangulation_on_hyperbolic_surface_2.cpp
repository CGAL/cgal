//#include <CGAL/Simple_cartesian.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Circular_kernel_2.h>
#include <CGAL/Algebraic_kernel_for_circles_2_2.h>
#include <CGAL/Hyperbolic_Delaunay_triangulation_CK_traits_2.h>
#include <CGAL/Hyperbolic_surface_Delaunay_traits_2.h>
#include <CGAL/Delaunay_triangulation_on_hyperbolic_surface_2.h>
#include <CGAL/Hyperbolic_fundamental_domain_factory_2.h>
#include <CGAL/Hyperbolic_Dirichlet_domain_2.h>
#include <CGAL/Triangulation_on_hyperbolic_surface_2_IO.h>

using Rational = CGAL::Gmpq;
//using Kernel = CGAL::Circular_kernel_2<CGAL::Simple_cartesian<Rational>,CGAL::Algebraic_kernel_for_circles_2_2<Rational>>;
using Kernel = CGAL::Circular_kernel_2<CGAL::Cartesian<Rational>,CGAL::Algebraic_kernel_for_circles_2_2<Rational>>;
using ParentTraits = CGAL::Hyperbolic_Delaunay_triangulation_CK_traits_2<Kernel>;
using HSTraits = CGAL::Hyperbolic_surface_traits_2<ParentTraits>;
using Traits = CGAL::Hyperbolic_surface_Delaunay_traits_2<HSTraits>;
using Domain = CGAL::Hyperbolic_fundamental_domain_2<Traits>;
using Factory = CGAL::Hyperbolic_fundamental_domain_factory_2<Traits>;
using Delaunay_triangulation = CGAL::Delaunay_triangulation_on_hyperbolic_surface_2<Traits>;

int main() {
  // Generate the domain
  Factory factory = Factory();
  Domain domain = factory.make_hyperbolic_fundamental_domain_g2(time(NULL)); // get a random seed
  Traits gt = Traits();
  unsigned p;
  double eps;

  std::cout << "/////// Delaunay of a random genus 2 surface given by an octagon ///////\n";
  Delaunay_triangulation dt = Delaunay_triangulation(gt, domain);
  dt.combinatorial_map().display_characteristics(std::cout) << std::endl;
  std::cout << "packing_value() = " << dt.packing_value() << std::endl;
  std::cout << "covering_value() = " << dt.covering_value() << std::endl;

  std::cout << "/////// Epsilon net of this random genus 2 surface ///////\n";
  p = 1;
  eps = 0.5;
  dt.set_circumcenter_approximation_precision(p);
  std::cout << "Computing a " << eps << "-net with floating-point precision " << 53*p << "...\n";
  std::cout << "Is epsilon-net? " << dt.construct_epsilon_net(eps) << std::endl;
  std::cout << "packing_value() = " << dt.packing_value() << std::endl;
  std::cout << "covering_value() = " << dt.covering_value() << std::endl;
  dt.combinatorial_map().display_characteristics(std::cout) << std::endl;

  // Locate a point in the triangulation
  Delaunay_triangulation::Point query(1/2,1/2);
  Delaunay_triangulation::Locate_type lt;
  unsigned li=0;
  unsigned ld=0;
  dt.locate(query, lt, li, ld, dt.anchor());
  std::cout << "The number of triangles traversed by the walk used for the point location algorithm is " << ld << std::endl;

  std::cout << "/////// Epsilon net of a genus 5 surface (this fails with double precision)///////\n";
  std::ifstream input_file("../FM-genus-5.txt", std::ios::in);
  if (!input_file.is_open()) {
    std::cerr << "Could not open input file ./FM-genus-5.txt!\n";
    return EXIT_FAILURE;
  }
  input_file >> domain;
  Delaunay_triangulation dt_g5 = Delaunay_triangulation(gt, domain);
  dt_g5.combinatorial_map().display_characteristics(std::cout) << std::endl;
  eps = 1;
  p = 1;
  dt_g5.set_circumcenter_approximation_precision(p);
  std::cout << "packing_value() = " << dt_g5.packing_value() << std::endl;
  std::cout << "Computing a " << eps << "-net with floating-point precision " << 53*p << "...\n";
  std::cout << "Is epsilon-net? " << dt_g5.construct_epsilon_net(eps) << std::endl;
  std::cout << "packing_value() = " << dt_g5.packing_value() << std::endl;
  std::cout << "covering_value() = " << dt_g5.covering_value() << std::endl;
  dt_g5.combinatorial_map().display_characteristics(std::cout) << std::endl;

  std::cout << "/////// Epsilon net of a genus 5 surface higher precision (this works with quad precision)///////\n";
  dt_g5.~Delaunay_triangulation();
  new (&dt_g5) Delaunay_triangulation(gt, domain);
  p = 2;
  dt_g5.set_circumcenter_approximation_precision(p);
  dt_g5.set_circumcenter_approximation_precision(p);
  std::cout << "Computing a " << eps << "-net with floating-point precision " << 53*p << "...\n";
  std::cout << "Is epsilon-net? " << dt_g5.construct_epsilon_net(eps) << std::endl;
  std::cout << "packing_value() = " << dt_g5.packing_value() << std::endl;
  std::cout << "covering_value() = " << dt_g5.covering_value() << std::endl;
  dt_g5.combinatorial_map().display_characteristics(std::cout) << std::endl;

  return EXIT_SUCCESS;
}
