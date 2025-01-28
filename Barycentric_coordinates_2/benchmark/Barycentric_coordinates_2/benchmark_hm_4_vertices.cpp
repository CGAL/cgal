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

int save_raw_data() {

  Timer timer;
  std::cout.precision(10);

  const std::vector<Point_2> vertices = {
    Point_2(0, 0), Point_2(1, 0),
    Point_2(1, 1), Point_2(0, 1)
  };

  const FT min_scale = 0.00158;
  const FT max_scale = 0.17500;
  FT step = (max_scale - min_scale) / FT(100);

  std::list<Point_2> seeds;
  seeds.push_back(Point_2(0.5, 0.5));

  Domain domain(vertices);
  HMC2 harmonic_coordinates_2(vertices, domain);

  std::size_t count = 0;
  for (FT scale = max_scale; scale >= min_scale; scale -= step) {
    if (scale < min_scale) break;
    if (count == 80) step /= FT(10);
    ++count;

    domain.clear();
    domain.create(scale, seeds);
    harmonic_coordinates_2.clear();

    timer.reset(); timer.start();
    harmonic_coordinates_2.setup();
    timer.stop();
    const double setup = timer.time();

    timer.reset(); timer.start();
    harmonic_coordinates_2.factorize();
    timer.stop();
    const double factorize = timer.time();

    timer.reset(); timer.start();
    harmonic_coordinates_2.solve();
    timer.stop();
    const double solve = timer.time();

    std::cout << domain.number_of_vertices() << " " <<
    setup << " " << factorize << " " << solve << std::endl;
  }
  return EXIT_SUCCESS;
}

int main() {

  // return save_raw_data();

  Timer timer;
  std::cout.precision(10);
  const std::size_t number_of_runs = 1;

  const std::vector<Point_2> vertices = {
    Point_2(0, 0), Point_2(1, 0),
    Point_2(1, 1), Point_2(0, 1)
  };

  const std::vector<FT> scales = {
    0.175, 0.072, 0.05, 0.01, 0.0071, 0.00502, 0.00224, 0.00158
  };

  std::list<Point_2> seeds;
  seeds.push_back(Point_2(0.5, 0.5));

  Domain domain(vertices);
  for (const FT scale : scales) {

    domain.clear();
    domain.create(scale, seeds);
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
    std::cout << "benchmark_hm_4_vertices, num queries: " <<
      domain.number_of_vertices() << std::endl;
    std::cout <<
      "benchmark_hm_4_vertices, compute (CPU time setup/factorize/solve/total): " <<
    setup << "/" << factorize << "/" << solve << "/" << total << " seconds" << std::endl;
  }
  return EXIT_SUCCESS;
}
