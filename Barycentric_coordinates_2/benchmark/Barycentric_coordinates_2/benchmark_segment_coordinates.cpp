#include <CGAL/Real_timer.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Barycentric_coordinates_2/segment_coordinates_2.h>

using Kernel  = CGAL::Exact_predicates_inexact_constructions_kernel;
using FT      = typename Kernel::FT;
using Point_2 = typename Kernel::Point_2;
using Timer   = CGAL::Real_timer;

int main() {

  const std::size_t number_of_iterations = 100000;
  const std::size_t number_of_points     = 10000;
  const std::size_t number_of_runs       = 1;

  const FT zero = FT(0);
  const FT one  = FT(1);
  const FT step = one / static_cast<FT>(number_of_points);

  const Point_2 p0 = Point_2(zero - step, zero);
  const Point_2 p1 = Point_2(one  + step, zero);

  Timer timer;
  std::vector<FT> coordinates(2);
  auto it = coordinates.begin();

  double time = 0.0;
  for (std::size_t i = 0; i < number_of_runs; ++i) {
    timer.start();
    for (std::size_t j = 0; j < number_of_iterations; ++j) {
      for (FT x = zero; x <= one; x += step) {
        const Point_2 query(x, zero);
        CGAL::Barycentric_coordinates::segment_coordinates_2(
          p0, p1, query, it);
      }
    }
    timer.stop();
    time += timer.time();
    timer.reset();
  }

  const double mean_time =
    time / static_cast<double>(number_of_runs);
  std::cout.precision(10);
  std::cout << "benchmark_segment_coordinates (CPU time): " <<
    mean_time << " seconds" << std::endl;
  return EXIT_SUCCESS;
}
