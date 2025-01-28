#include <CGAL/Real_timer.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Barycentric_coordinates_2/Mean_value_coordinates_2.h>

using Kernel   = CGAL::Exact_predicates_inexact_constructions_kernel;
using FT       = typename Kernel::FT;
using Point_2  = typename Kernel::Point_2;
using Timer    = CGAL::Real_timer;
using Vertices = std::vector<Point_2>;

using MVC2 = CGAL::Barycentric_coordinates::Mean_value_coordinates_2<Vertices, Kernel>;

int main() {

  const std::size_t number_of_x_coordinates = 100000;
  const std::size_t number_of_y_coordinates = 1000;
  const std::size_t number_of_runs          = 1;

  const FT zero = FT(0);
  const FT one  = FT(1);
  const FT x_step = one / static_cast<FT>(number_of_x_coordinates);
  const FT y_step = one / static_cast<FT>(number_of_y_coordinates);

  const std::vector<Point_2> vertices = {
    Point_2(zero, zero - y_step), Point_2(one, zero - y_step),
    Point_2(FT(3) / FT(2), -FT(1) / FT(2)), Point_2(2, -FT(1) / FT(2)),
    Point_2(FT(5) / FT(2), 0), Point_2(2, FT(1) / FT(2)),
    Point_2(FT(5) / FT(2), 1), Point_2(3, FT(3) / FT(4)),
    Point_2(3, FT(5) / FT(4)), Point_2(FT(5) / FT(2), FT(7) / FT(4)),
    Point_2(3, FT(5) / FT(2)), Point_2(FT(5) / FT(2), FT(5) / FT(2)),
    Point_2(FT(9) / FT(4), 2), Point_2(FT(7) / FT(4), 2),
    Point_2(2, FT(5) / FT(2)), Point_2(FT(3) / FT(2), FT(5) / FT(2)),
    Point_2(FT(5) / FT(4), 2), Point_2(FT(3) / FT(4), 2),
    Point_2(1, FT(5) / FT(2)), Point_2(FT(1) / FT(2), FT(5) / FT(2)),
    Point_2(FT(1) / FT(4), 2), Point_2(-FT(1) / FT(4), 2),
    Point_2(zero, FT(5) / FT(2)), Point_2(-FT(1) / FT(2), FT(5) / FT(2)),
    Point_2(-FT(3) / FT(4), 2), Point_2(-FT(1) / FT(2), FT(3) / FT(2)),
    Point_2(-FT(5) / FT(4), FT(3) / FT(2)), Point_2(-FT(1) / FT(2), FT(3) / FT(4)),
    Point_2(-one, FT(1) / FT(2)), Point_2(-one, zero),
    Point_2(-FT(3) / FT(2), zero), Point_2(-FT(3) / FT(2), -FT(1) / FT(2)),
    Point_2(-FT(1) / FT(2), -FT(1) / FT(2)), Point_2(-FT(1) / FT(2), zero - y_step)
  };

  MVC2 mean_value_coordinates_2(
    vertices, CGAL::Barycentric_coordinates::Computation_policy_2::FAST);

  Timer timer;
  std::vector<FT> coordinates(34);
  auto it = coordinates.begin();

  double time = 0.0;
  for (std::size_t i = 0; i < number_of_runs; ++i) {
    timer.start();
    for (FT x = zero; x <= one; x += x_step) {
      for (FT y = zero; y <= one; y += y_step) {
        const Point_2 query(x, y);

           mean_value_coordinates_2(query, it);
        // mean_value_coordinates_2.weights(query, it);
      }
    }
    timer.stop();
    time += timer.time();
    timer.reset();
  }

  const double mean_time =
    time / static_cast<double>(number_of_runs);
  std::cout.precision(10);
  std::cout << "benchmark_mv_34_vertices (CPU time): " <<
    mean_time << " seconds" << std::endl;
  return EXIT_SUCCESS;
}
