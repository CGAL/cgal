#include <CGAL/Real_timer.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Barycentric_coordinates_3/tetrahedron_coordinates.h>

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using FT = typename Kernel::FT;
using Point_3 = typename Kernel::Point_3;
using Timer = CGAL::Real_timer;

int main() {

  const std::size_t number_of_x_coordinates = 100;
  const std::size_t number_of_y_coordinates = 100;
  const std::size_t number_of_z_coordinates = 100;

  const std::size_t number_of_runs = 10;

  const FT zero = FT(0);
  const FT one = FT(1);
  const FT four = FT(4);

  const FT x_step = one / static_cast<FT>(number_of_x_coordinates);
  const FT y_step = one / static_cast<FT>(number_of_y_coordinates);
  const FT z_step = one / static_cast<FT>(number_of_z_coordinates);

  const Point_3 p0 = Point_3(zero - x_step, zero - y_step, zero - z_step);
  const Point_3 p1 = Point_3(four + x_step, zero - y_step, zero - z_step);
  const Point_3 p2 = Point_3(zero - x_step, four + y_step, zero - z_step);
  const Point_3 p3 = Point_3(zero - x_step, zero - y_step, four + z_step);

  Timer timer;
  std::vector<FT> coordinates(4);
  double time = 0.0;

  for(std::size_t i = 0; i < number_of_runs; i++){

    timer.start();
    for (FT x = zero; x <= one; x += x_step){
      for (FT y = zero; y <= one; y += y_step){
        for (FT z = zero; z <= one; z += z_step){

          const Point_3 query(x, y, z);
          CGAL::Barycentric_coordinates::tetrahedron_coordinates(
            p0, p1, p2, p3, query, coordinates.begin());
        }
      }
    }

    timer.stop();
    time += timer.time();
    timer.reset();
  }

  const double mean_time =
    time / static_cast<double>(number_of_runs);

  std::cout.precision(10);
  std::cout << "benchmark_tretrahedron_coordinates (CPU time): " <<
    mean_time << " seconds" << std::endl;
  return EXIT_SUCCESS;
}
