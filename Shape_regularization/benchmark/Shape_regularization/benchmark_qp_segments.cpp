#include <CGAL/Random.h>
#include <CGAL/Real_timer.h>
#include <CGAL/function_objects.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/Join_input_iterator.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Shape_regularization.h>

// #include "../../test/Shape_regularization/include/Saver.h"

using Kernel    = CGAL::Exact_predicates_inexact_constructions_kernel;
using FT        = typename Kernel::FT;
using Point_2   = typename Kernel::Point_2;
using Segment_2 = typename Kernel::Segment_2;
using Segments  = std::vector<Segment_2>;
using Indices   = std::vector<std::size_t>;
using Timer     = CGAL::Real_timer;

using PCreator = CGAL::Creator_uniform_2<FT, Point_2>;
using RS       = CGAL::Random_points_on_segment_2<Point_2, PCreator>;
using RC       = CGAL::Random_points_on_circle_2<Point_2, PCreator>;
using SCreator = CGAL::Creator_uniform_2<Point_2, Segment_2>;
using RSRC     = CGAL::Join_input_iterator_2<RS, RC, SCreator>;

using NQ = CGAL::Shape_regularization::Segments::Delaunay_neighbor_query_2<Kernel, Segments>;
using AR = CGAL::Shape_regularization::Segments::Angle_regularization_2<Kernel, Segments>;
using OR = CGAL::Shape_regularization::Segments::Offset_regularization_2<Kernel, Segments>;

void create_pattern(
  const std::size_t n,
  std::vector<Segment_2>& segments) {

  segments.clear();
  segments.reserve(n);

  FT step = FT(1), d = FT(1) / FT(10);
  FT x1 = FT(0), y1 = FT(0), x2 = FT(0), y2 = FT(0);

  for (std::size_t i = 0; i < n / 2; ++i) {
    x1 = step * i;
    y1 = FT(0);
    if (i % 2 == 0) x2 = step * i;
    else x2 = step * i - d;
    y2 = FT(1);
    const Point_2 source = Point_2(x1, y1);
    const Point_2 target = Point_2(x2, y2);
    segments.push_back(Segment_2(source, target));
  }

  for (std::size_t i = 0; i < n / 2; ++i) {
    x1 = step * i + d;
    y1 = FT(2);
    if (i % 2 == 0) x2 = step * i + d;
    else x2 = step * i;
    y2 = FT(4);
    const Point_2 source = Point_2(x1, y1);
    const Point_2 target = Point_2(x2, y2);
    segments.push_back(Segment_2(source, target));
  }
  assert(segments.size() == n);

  // CGAL::Shape_regularization::Tests::Saver<Kernel> saver;
  // saver.export_eps_segments(segments, "pattern", 100.0);
  // exit(EXIT_SUCCESS);
}

void create_random_in_square(
  const std::size_t n,
  std::vector<Segment_2>& segments) {

  segments.clear();
  segments.reserve(n);

  CGAL::Random rand;
  for (std::size_t i = 0; i < n / 2; ++i) {
    const FT x = static_cast<FT>(rand.get_int(-250, 250));
    const FT y = static_cast<FT>(rand.get_int(-250, 250));
    segments.push_back(Segment_2(
      Point_2(x, y), Point_2(x, -y)));
    segments.push_back(Segment_2(
      Point_2(-x, y), Point_2(x, y)));
  }
  assert(segments.size() == n);

  // Perturb segments.
  for (auto& segment : segments) {
    const FT angle = static_cast<FT>(rand.get_int(-15, 15));
    CGAL::Shape_regularization::internal::
      rotate_segment_2(angle, FT(0), segment);
  }

  // CGAL::Shape_regularization::Tests::Saver<Kernel> saver;
  // saver.export_eps_segments(segments, "pseudo_random", 1.0);
  // exit(EXIT_SUCCESS);
}

void create_random_in_circle(
  const std::size_t n,
  std::vector<Segment_2>& segments) {

  segments.clear();
  segments.reserve(n);

  const Point_2 source(0, -100);
  const Point_2 target(0, +100);
  const FT radius = 250.0;
  RS p1(source, target);
  RC p2(radius);

  RSRC generator(p1, p2);
  std::copy_n(generator, n, std::back_inserter(segments));
  assert(segments.size() == n);

  // CGAL::Shape_regularization::Tests::Saver<Kernel> saver;
  // saver.export_eps_segments(segments, "random", 1.0);
  // exit(EXIT_SUCCESS);
}

void benchmark_qp_segments(
  const std::size_t n,
  const bool regroup,
  const bool simple_output,
  const std::size_t num_iters) {

  const std::size_t m = 10; // number of segments in a group
  std::vector<Segment_2> segments;

  create_pattern(n, segments);

  // create_random_in_square(n, segments);
  // create_random_in_circle(n, segments);

  Timer timer;
  timer.start();
  NQ neighbor_query(segments);
  if (regroup) {
    Indices group;
    group.reserve(m);
    for (std::size_t i = 0; i < n;) {
      group.clear();
      for (std::size_t j = 0; j < m; ++j) {
        group.push_back(i + j);
      }
      neighbor_query.add_group(group);
      i += m;
    }
  }
  timer.stop();
  // const double delaunay_time = timer.time();
  timer.reset();

  const FT max_angle_2 = FT(10);
  timer.start();
  AR angle_regularization(
    segments, CGAL::parameters::maximum_angle(max_angle_2));
  if (regroup) {
    Indices group;
    group.reserve(m);
    for (std::size_t i = 0; i < n;) {
      group.clear();
      for (std::size_t j = 0; j < m; ++j) {
        group.push_back(i + j);
      }
      angle_regularization.add_group(group);
      i += m;
    }
  }
  timer.stop();
  // const double setup_angle_time = timer.time();
  timer.reset();

  double angle_time = 0.0;
  for (std::size_t i = 0; i < num_iters; ++i) {
    auto copied = segments;
    timer.start();
    CGAL::Shape_regularization::Segments::regularize_segments(
      copied, neighbor_query, angle_regularization);
    timer.stop();
    angle_time += timer.time();
    timer.reset();
  }
  angle_time /= static_cast<double>(num_iters);

  const FT max_offset_2 = FT(1) / FT(5);
  timer.start();
  OR offset_regularization(
    segments, CGAL::parameters::maximum_offset(max_offset_2));
  timer.stop();
  // const double setup_offset_time = timer.time();
  timer.reset();

  timer.start();
  std::vector<Indices> pgroups;
  angle_regularization.parallel_groups(
    std::back_inserter(pgroups));
  timer.stop();
  // const double init_group_time = timer.time();
  timer.reset();

  timer.start();
  neighbor_query.clear();
  for (const auto& pgroup : pgroups) {
    neighbor_query.add_group(pgroup);
    offset_regularization.add_group(pgroup);
  }
  timer.stop();
  // const double add_group_time = timer.time();
  timer.reset();

  double offset_time = 0.0;
  for (std::size_t i = 0; i < num_iters; ++i) {
    auto copied = segments;
    timer.start();
    CGAL::Shape_regularization::Segments::regularize_segments(
      copied, neighbor_query, offset_regularization);
    timer.stop();
    offset_time += timer.time();
    timer.reset();
  }
  offset_time /= static_cast<double>(num_iters);

  std::cout.precision(10);
  if (regroup && !simple_output) {
    std::cout << "grouped: " ;
  }

  // std::cout << "benchmark_qp_segments " << segments.size() << " (CPU time " <<
  // "delaunay/setup_angles/angles/setup_offsets/offsets): " <<
  //   delaunay_time << "/" <<
  //   setup_angle_time << "/" << angle_time << "/" <<
  //   setup_offset_time << "/" << offset_time <<
  // " seconds" << std::endl;

  if (!simple_output) {
    std::cout << "benchmark_qp_segments " << segments.size() << " (CPU time " <<
    "angles/offsets): " << angle_time << "/" << offset_time << " seconds" << std::endl;
  } else {
    if (!regroup) {
      std::cout << segments.size() << " " << angle_time << " " << offset_time << " ";
    } else {
      std::cout << angle_time << " " << offset_time << std::endl;
    }
  }
}

int main() {

  const std::size_t num_iters = 1;
  const std::vector<std::size_t> ns = {
    10, 50, 100, 500, 1000, 5000, 10000, 15000
  };
  std::cout << std::endl;
  for (const std::size_t n : ns) {
    benchmark_qp_segments(n, false, false, num_iters);
    benchmark_qp_segments(n, true , false, num_iters);
    std::cout << std::endl;
  }

  // Dense results for plotting.
  // std::vector<std::size_t> ns;
  // for (std::size_t i = 10; i <= 10000; i += 10)
  //   ns.push_back(i);
  // for (const std::size_t n : ns) {
  //   benchmark_qp_segments(n, false, true, num_iters);
  //   benchmark_qp_segments(n, true , true, num_iters);
  // }

  return EXIT_SUCCESS;
}
