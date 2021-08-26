#include <CGAL/Real_timer.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Barycentric_coordinates_2/Delaunay_domain_2.h>
#include <CGAL/Barycentric_coordinates_2/Harmonic_coordinates_2.h>

using Kernel   = CGAL::Exact_predicates_inexact_constructions_kernel;
using FT       = typename Kernel::FT;
using Point_2  = typename Kernel::Point_2;
using Timer    = CGAL::Real_timer;
using Vertices = std::vector<Point_2>;

using Domain = CGAL::Barycentric_coordinates::Delaunay_domain_2<Vertices, Kernel>;
using HMC2   = CGAL::Barycentric_coordinates::Harmonic_coordinates_2<Vertices, Domain, Kernel>;

void generate_regular_polygon(
  const std::size_t n,
  const double radius,
  Vertices& vertices) {

  vertices.clear();
  vertices.reserve(n);
  for (std::size_t i = 0; i < n; ++i) {
    vertices.push_back(Point_2(
      static_cast<FT>(+radius * sin((CGAL_PI / n) + ((i * 2.0 * CGAL_PI) / n))),
      static_cast<FT>(-radius * cos((CGAL_PI / n) + ((i * 2.0 * CGAL_PI) / n)))));
  }
}

int main() {

  Timer timer;
  std::cout.precision(10);
  const double radius = 1.0;
  const std::size_t number_of_runs = 1;

  const std::vector<std::size_t> ns = {
    5, 10, 25, 50, 100, 500, 1000
  };
  const std::vector<FT> scales = {
    0.00718, 0.0086, 0.00885, 0.00886, 0.00887, 0.00888, 0.00888
  };

  std::list<Point_2> seeds;
  seeds.push_back(Point_2(0.0, 0.0));
  for (std::size_t i = 0; i < ns.size(); ++i) {

    Vertices vertices;
    generate_regular_polygon(ns[i], radius, vertices);

    Domain domain(vertices);
    domain.create(scales[i], seeds);
    HMC2 harmonic_coordinates_2(vertices, domain);

    double setup = 0.0, factorize = 0.0, solve = 0.0;
    for (std::size_t k = 0; k < number_of_runs; ++k) {
      harmonic_coordinates_2.clear();

      timer.reset(); timer.start();
      harmonic_coordinates_2.setup();
      timer.stop();
      setup += timer.time();

      timer.reset(); timer.start();
      harmonic_coordinates_2.factorize();
      timer.stop();
      factorize += timer.time();

      timer.reset(); timer.start();
      harmonic_coordinates_2.solve();
      timer.stop();
      solve += timer.time();
    }
    setup /= static_cast<double>(number_of_runs);
    factorize /= static_cast<double>(number_of_runs);
    solve /= static_cast<double>(number_of_runs);

    const double total = setup + factorize + solve;
    std::cout << "benchmark_hm_n_vertices, num vertices/num queries: " <<
      ns[i] << "/" << domain.number_of_vertices() << std::endl;
    std::cout <<
      "benchmark_hm_n_vertices, compute (CPU time setup/factorize/solve/total): " <<
    setup << "/" << factorize << "/" << solve << "/" << total << " seconds" << std::endl;
  }
  return EXIT_SUCCESS;
}
