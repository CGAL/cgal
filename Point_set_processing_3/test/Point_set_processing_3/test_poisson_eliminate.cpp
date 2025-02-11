#include <vector>
#include <iostream>
#include <fstream>
#include <filesystem>
#include <string>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/poisson_eliminate.h>
#include <CGAL/compute_average_spacing.h>
#include <CGAL/IO/read_points.h>
#include <CGAL/IO/write_points.h>
#include <CGAL/Real_timer.h>

template<class K>
void check(const std::vector<typename K::Point_3> &input, const std::vector<typename K::Point_3>& output, bool progressive) {
  using FT = typename K::FT;
  FT average_distance_input = CGAL::compute_average_spacing<CGAL::Parallel_if_available_tag>(input, 3);
  std::cout << " " << average_distance_input << " average point distance" << std::endl;

  if (output.size() < input.size()) {
    FT average_distance_after = CGAL::compute_average_spacing<CGAL::Parallel_if_available_tag>(output, 3);
    std::cout << " " << average_distance_after << " average point distance after elimination" << std::endl;
    assert(average_distance_after > average_distance_input);

    if (progressive) {
      FT average_distance_progressive = CGAL::compute_average_spacing<CGAL::Parallel_if_available_tag>(CGAL::make_range(output.begin(), output.begin() + output.size() / 2), 3);
      std::cout << " " << average_distance_progressive << " average point distance in the first half after progressive elimination" << std::endl;
      assert(average_distance_progressive > average_distance_after);
    }
  }
}

void no_check(const std::vector<CGAL::Exact_predicates_exact_constructions_kernel::Point_3>& input, const std::vector<CGAL::Exact_predicates_exact_constructions_kernel::Point_3>& output, bool progressive) {
}

template<class K, class Check_distance_functor>
void sampling(const std::string& filename, Check_distance_functor check_average_distance) {
  using Point_3 = typename K::Point_3;
  std::vector<Point_3> points;

  if (!CGAL::IO::read_points(
    filename,
    std::back_inserter(points))) {

    std::cerr << "Error: cannot read file!" << std::endl;
    return;
  }

  std::cout << " " << points.size() << " input points" << std::endl;

  std::size_t target_size = points.size() / 3;

  bool progressive = true;
  bool weight_limiting = false;
  bool tiling = false;
  std::vector<Point_3> output;

  CGAL::Real_timer timer;

  output.reserve(target_size);
  timer.start();
  CGAL::poisson_eliminate(points, target_size, std::back_inserter(output), CGAL::parameters::dimension(2).progressive(progressive).weight_limiting(weight_limiting).tiling(tiling));
  timer.stop();
  std::cout << timer.time() << std::endl;
  std::cout << " " << output.size() << " points after elimination" << std::endl;

  CGAL::IO::write_points("radar" + std::to_string(target_size) + "_2d.xyz", output, CGAL::parameters::stream_precision(17));

  assert(output.size() == target_size);

  check_average_distance(points, output, progressive);

  std::cout << " done" << std::endl;
}


int main(int argc, char* argv[])
{
  std::cout << "testing Simple_cartesian<double>" << std::endl;
  sampling<CGAL::Simple_cartesian<double>>(CGAL::data_file_path("points_3/radar.xyz"), check< CGAL::Simple_cartesian<double>>);

   std::cout << std::endl << "testing Exact_predicates_inexact_constructions_kernel" << std::endl;
   sampling<CGAL::Exact_predicates_inexact_constructions_kernel>(CGAL::data_file_path("points_3/radar.xyz"), check< CGAL::Exact_predicates_inexact_constructions_kernel>);

   std::cout << std::endl << "testing Exact_predicates_exact_constructions_kernel" << std::endl;
   sampling<CGAL::Exact_predicates_exact_constructions_kernel>(CGAL::data_file_path("points_3/radar.xyz"), no_check);

  return 0;
}
