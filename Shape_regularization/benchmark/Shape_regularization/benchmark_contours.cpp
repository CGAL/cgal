#include <CGAL/Random.h>
#include <CGAL/Real_timer.h>
#include <CGAL/random_polygon_2.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Shape_regularization.h>

#include "../../test/Shape_regularization/include/Saver.h"

using Kernel    = CGAL::Exact_predicates_inexact_constructions_kernel;
using FT        = typename Kernel::FT;
using Point_2   = typename Kernel::Point_2;
using Segment_2 = typename Kernel::Segment_2;
using Contour   = std::vector<Point_2>;
using Timer     = CGAL::Real_timer;

using Creator   = CGAL::Creator_uniform_2<int, Point_2>;
using Generator = CGAL::Random_points_in_square_2<Point_2, Creator>;

using CD = CGAL::Shape_regularization::Contours::Longest_direction_2<Kernel, Contour>;

void create_pattern(
  const std::size_t n,
  std::vector<Point_2>& contour,
  const FT d) {

  contour.clear();
  contour.reserve(n);

  FT step = FT(1);
  FT x = FT(0), y = FT(0);

  for (std::size_t i = 0; i < n / 2; ++i) {
    x = step * i;
    if (i % 2 == 0) y = FT(0);
    else y = d;
    const Point_2 vertex = Point_2(x, y);
    contour.push_back(vertex);
  }
  contour.push_back(Point_2(x, FT(2)));

  for (std::size_t i = 0; i < n / 2 - 1; ++i) {
    x -= step;
    if (i % 2 == 0) y = FT(2) - d;
    else y = FT(2);
    const Point_2 vertex = Point_2(x, y);
    contour.push_back(vertex);
  }
  assert(contour.size() == n);

  // CGAL::Shape_regularization::Tests::Saver<Kernel> saver;
  // saver.export_eps_closed_contour(contour, "pattern_closed", 100.0);
  // saver.export_eps_open_contour(contour, "pattern_open", 100.0);
  // exit(EXIT_SUCCESS);
}

void create_pseudo_random(
  const std::size_t n,
  std::vector<Point_2>& contour) {

  const FT d = FT(0);
  create_pattern(n, contour, d);

  // Perturb edges.
  CGAL::Random rand;
  for (std::size_t i = 0; i < contour.size(); ++i) {
    const std::size_t ip = (i + 1) % contour.size();
    Segment_2 edge(contour[i], contour[ip]);
    const FT angle = static_cast<FT>(rand.get_int(-15, 15));
    CGAL::Shape_regularization::internal::
      rotate_segment_2(angle, FT(0), edge);
    contour[i] = edge.source();
    contour[ip] = edge.target();
  }

  // CGAL::Shape_regularization::Tests::Saver<Kernel> saver;
  // saver.export_eps_closed_contour(contour, "pseudo_random_closed", 100.0);
  // saver.export_eps_open_contour(contour, "pseudo_random_open", 100.0);
  // exit(EXIT_SUCCESS);
}

void create_random(
  const std::size_t n,
  std::vector<Point_2>& contour) {

  contour.clear();
  contour.reserve(n);

  const FT radius = 100.0;
  std::vector<Point_2> point_set;
  point_set.reserve(n);

  CGAL::copy_n_unique(
    Generator(radius), n, std::back_inserter(point_set));
  CGAL::random_polygon_2(
    point_set.size(), std::back_inserter(contour), point_set.begin());

  // The above algorithm does not guarantee to have exactly n vertices!
  // assert(contour.size() == n);

  // CGAL::Shape_regularization::Tests::Saver<Kernel> saver;
  // saver.export_eps_closed_contour(contour, "random_closed", 1.0);
  // saver.export_eps_open_contour(contour, "random_open", 1.0);
  // exit(EXIT_SUCCESS);
}

void benchmark_contours(
  const std::size_t n,
  const bool simple_output,
  const std::size_t num_iters) {

  Timer timer;
  std::vector<Point_2> contour, regularized;
  const FT d = FT(1) / FT(10);

  create_pattern(n, contour, d);

  // create_pseudo_random(n, contour);
  // create_random(n, contour);

  const FT max_offset_2 = FT(1) / FT(5);
  regularized.clear();

  timer.start();
  CD closed_directions(contour, true);
  timer.stop();
  // const double longest_closed = timer.time();
  timer.reset();

  double closed_time = 0.0;
  for (std::size_t i = 0; i < num_iters; ++i) {
    regularized.clear();
    timer.start();
    CGAL::Shape_regularization::Contours::regularize_closed_contour(
      contour, closed_directions, std::back_inserter(regularized),
      CGAL::parameters::maximum_offset(max_offset_2));
    timer.stop();
    closed_time += timer.time();
    timer.reset();
  }
  closed_time /= static_cast<double>(num_iters);

  regularized.clear();
  timer.start();
  CD open_directions(contour, false);
  timer.stop();
  // const double longest_open = timer.time();
  timer.reset();

  double open_time = 0.0;
  for (std::size_t i = 0; i < num_iters; ++i) {
    regularized.clear();
    timer.start();
    CGAL::Shape_regularization::Contours::regularize_open_contour(
      contour, open_directions, std::back_inserter(regularized),
      CGAL::parameters::maximum_offset(max_offset_2));
    timer.stop();
    open_time += timer.time();
    timer.reset();
  }
  open_time /= static_cast<double>(num_iters);

  std::cout.precision(10);
  if (simple_output) {
    std::cout << contour.size() << " " << closed_time << " " << open_time << std::endl;
  } else {
    std::cout << "benchmark_contours " << contour.size() << " (CPU time " <<
    "closed/open): " << closed_time << "/" << open_time << " seconds" << std::endl;
  }
}

int main() {

  const std::size_t num_iters = 1;
  const std::vector<std::size_t> ns = {
    10, 100, 1000, 10000, 100000, 1000000, 10000000
  };
  for (const std::size_t n : ns) {
    benchmark_contours(n, false, num_iters);
  }

  // Dense results for plotting.
  // std::vector<std::size_t> ns;
  // for (std::size_t i = 10; i <= 10000; i += 10)
  //   ns.push_back(i);
  // for (const std::size_t n : ns)
  //   benchmark_contours(n, true, num_iters);

  return EXIT_SUCCESS;
}
