#include <CGAL/Real_timer.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Barycentric_coordinates_2/triangle_coordinates_2.h>

using Kernel  = CGAL::Exact_predicates_inexact_constructions_kernel;
using FT      = typename Kernel::FT;
using Point_2 = typename Kernel::Point_2;
using Timer   = CGAL::Real_timer;

int main() {

  const std::size_t number_of_x_coordinates = 100000;
  const std::size_t number_of_y_coordinates = 10000;
  const std::size_t number_of_runs          = 1;

  const FT zero = FT(0);
  const FT one  = FT(1);
  const FT two  = FT(2);
  const FT x_step = one / static_cast<FT>(number_of_x_coordinates);
  const FT y_step = one / static_cast<FT>(number_of_y_coordinates);

  const Point_2 p0 = Point_2(zero - x_step, zero - x_step);
  const Point_2 p1 = Point_2(two  + y_step, zero - x_step);
  const Point_2 p2 = Point_2(zero - x_step, two  + y_step);

  Timer timer;
  std::vector<FT> coordinates(3);
  auto it = coordinates.begin();

  double time = 0.0;
  for (std::size_t i = 0; i < number_of_runs; ++i) {
    timer.start();
    for (FT x = zero; x <= one; x += x_step) {
      for (FT y = zero; y <= one; y += y_step) {
        const Point_2 query(x, y);
        CGAL::Barycentric_coordinates::triangle_coordinates_2(
          p0, p1, p2, query, it);
      }
    }
    timer.stop();
    time += timer.time();
    timer.reset();
  }

  const double mean_time =
    time / static_cast<double>(number_of_runs);
  std::cout.precision(10);
  std::cout << "benchmark_triangle_coordinates (CPU time): " <<
    mean_time << " seconds" << std::endl;
  return EXIT_SUCCESS;
}
