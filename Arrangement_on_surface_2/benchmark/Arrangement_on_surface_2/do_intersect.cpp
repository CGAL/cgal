#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <chrono>
#include <iterator>

#include <boost/program_options.hpp>
#include <boost/iterator/function_output_iterator.hpp>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arr_polyline_traits_2.h>
#include <CGAL/Arr_non_caching_segment_basic_traits_2.h>
#include <CGAL/Arr_linear_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/draw_arrangement_2.h>
#include <CGAL/intersections.h>
#include <CGAL/Surface_sweep_2_algorithms.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/Join_input_iterator.h>
#include <CGAL/algorithm.h>

namespace po = boost::program_options;
namespace fs = std::filesystem;

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
// using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;
using Point = Kernel::Point_2;
using Segment = Kernel::Segment_2;

// The base traits handling segments
using Segment_traits = CGAL::Arr_segment_traits_2<Kernel>;
using Non_caching_segment_traits = CGAL::Arr_non_caching_segment_traits_2<Kernel>;

// The polyline traits wrapping the segment traits
using Polyline_traits = CGAL::Arr_polyline_traits_2<Segment_traits>;
using Polyline = Polyline_traits::Curve_2;
using X_monotone_polyline = Polyline_traits::X_monotone_curve_2;
using Polyline_point = Polyline_traits::Point_2;
using Polyline_multiplicity = Polyline_traits::Multiplicity;

// using Traits = CGAL::Arr_non_caching_segment_basic_traits_2<Kernel>;

using X_monotone_segment = Segment_traits::X_monotone_curve_2;

/*! Generates n random polylines within [0,1]x[0,1].
 * Strategy:
 * - Maintain a pool of existing vertices.
 * - To encourage endpoint intersections: Occasionally snap start/end points to the pool.
 * - To encourage interior intersections: Randomly place intermediate points.
 */
std::vector<Polyline> generate_polylines(std::size_t n, std::size_t seed) {
  std::vector<Polyline> polylines;
  std::vector<Point> vertex_pool;

  // Random Number Generation Setup
  std::mt19937 gen(seed);

  // Distribution for coordinates in [0, 1]
  std::uniform_real_distribution<double> coord_dist(0.0, 1.0);
  // Distribution for number of segments: 2, 3, or 4
  std::uniform_int_distribution<int> seg_count_dist(2, 4);
  // Coin flip for reusing existing vertices (50% chance)
  std::bernoulli_distribution reuse_dist(0.5);

  for (std::size_t i = 0; i < n; ++i) {
    auto num_segments = seg_count_dist(gen);
    auto num_points = num_segments + 1;

    std::vector<Point> points;
    points.reserve(num_points);

    for (std::size_t k = 0; k < num_points; ++k) {
      Point p;
      bool picked_existing = false;

      // Attempt to reuse an existing vertex for the Start or End of the polyline
      // This forces "Intersection at endpoints"
      bool is_endpoint = (k == 0 || k == num_points - 1);

      if (is_endpoint && !vertex_pool.empty() && reuse_dist(gen)) {
        // Pick a random index from the pool
        std::uniform_int_distribution<int> pool_idx_dist(0, vertex_pool.size() - 1);
        p = vertex_pool[pool_idx_dist(gen)];
        picked_existing = true;
      }
      else {
        // Generate a new random point in [0, 1]
        p = Point(coord_dist(gen), coord_dist(gen));
      }

      // Validation: CGAL Polylines cannot have zero-length segments.
      // If the new point is identical to the previous one, regenerate it.
      if (!points.empty() && p == points.back()) {
        // Simple retry strategy: perturb slightly or force new random
        p = Point(coord_dist(gen), coord_dist(gen));
      }

      points.push_back(p);
    }

    // Create the polyline using the traits constructor
    // Note: Arr_polyline_traits_2 constructs from a range of points
    Polyline polyline(points.begin(), points.end());
    polylines.push_back(polyline);

    // Add these points to the pool for future curves to potentially connect to
    vertex_pool.insert(vertex_pool.end(), points.begin(), points.end());
  }

  return polylines;
}

//!
template <typename Traits>
struct Test_functor {
  using X_monotone_curve = typename Traits::X_monotone_curve_2;

  bool& m_res;
  const X_monotone_curve& m_xcv1;
  const X_monotone_curve& m_xcv2;
  const Traits& m_traits;
  bool m_closed = true;
  Test_functor(bool& res, const X_monotone_curve& xcv1, const X_monotone_curve& xcv2, const Traits& traits,
               bool closed = true) :
    m_res(res),
    m_xcv1(xcv1),
    m_xcv2(xcv2),
    m_traits(traits),
    m_closed(closed)
  {}

  void operator()(std::pair<Polyline_point, Polyline_multiplicity> xres) {
    if (m_closed) {
      m_res = true;
      return;
    }

    if (m_res) return;

    auto eq = m_traits.equal_2_object();
    auto ctr_min_vertex = m_traits.construct_min_vertex_2_object();
    auto ctr_max_vertex = m_traits.construct_max_vertex_2_object();
    const auto& l1 = ctr_min_vertex(m_xcv1);
    const auto& r1 = ctr_max_vertex(m_xcv1);
    const auto& l2 = ctr_min_vertex(m_xcv2);
    const auto& r2 = ctr_max_vertex(m_xcv2);
    // std::cout << "YYYYYYYY\n";
    // std::cout << eq(xres.first, l1) << std::endl;
    // std::cout << eq(xres.first, r1) << std::endl;
    // std::cout << eq(xres.first, l2) << std::endl;
    // std::cout << eq(xres.first, r2) << std::endl;
    m_res = ((! eq(xres.first, l1) && ! eq(xres.first, r1)) || (! eq(xres.first, l2) && ! eq(xres.first, r2)));
    // std::cout << "m_res 2: " << m_res << std::endl;
  }
  void operator()(const X_monotone_polyline& poly) { m_res = true; }
};

// void print_exact_polyline(const X_monotone_polyline& xcv) {
//   std::cout << xcv.number_of_subcurves();
//   for (auto it = xcv.points_begin(); it != xcv.points_end(); ++it) {
//     const auto& p = *it;
//     std::cout << " " << p.x().exact() << " " << p.y().exact();
//   }
// }

//!
template <typename XCurves, typename Traits>
std::pair<bool, bool> verify_intersections(XCurves& xcurves, const Traits& traits, bool closed = true) {
  std::size_t cnt_do_intersect = 0;
  std::size_t cnt_intersect = 0;
  auto my_do_intersect = traits.do_intersect_2_object();
  auto my_intersect = traits.intersect_2_object();
  std::size_t counter = 0;
  std::size_t i = 0;

  for (const auto& s : xcurves) {
    ++i;
    std::size_t j = 0;
    for (const auto& t : xcurves) {
      ++j;
      if (j >= i) break;
      auto do_res = my_do_intersect(s, t, closed);
      if (do_res) ++cnt_do_intersect;

      bool res = false;
      my_intersect(s, t, boost::make_function_output_iterator(Test_functor(res, s, t, traits, closed)));
      if (res) ++cnt_intersect;
      if (do_res != res) {
#if 0
        print_exact_polyline(s);
        std::cout << std::endl;
        print_exact_polyline(t);
        std::cout << std::endl;
#endif
        std::cout << "s: " << s << std::endl;
        std::cout << "t: " << t << std::endl;
        std::cout << "do: " << do_res << ", cnt: " << res << std::endl;
        return std::make_pair(false, false);
      }
    }
  }

  // std::cout << "Do instersect counter: " << cnt_do_intersect << "\n";
  // std::cout << "Instersect counter:    " << cnt_intersect << "\n";
  return std::make_pair(true, cnt_intersect != 0);
}

//!
template <typename XCurves, typename Traits>
std::size_t count_intersections(XCurves& xcurves, const Traits& traits) {
  auto my_do_intersect = traits.do_intersect_2_object();
  std::size_t counter = 0;
  for (const auto& s : xcurves) {
    for (const auto& t : xcurves) if (my_do_intersect(s, t)) ++counter;
  }
  return counter;
}

/*! generates the curves.
 */
std::vector<X_monotone_polyline> generate(std::size_t n, std::size_t seed, const Polyline_traits& traits) {
  auto curves = generate_polylines(n, seed);
  auto make_x_monotone = traits.make_x_monotone_2_object();

  using Make_x_monotone_result = std::variant<Point, X_monotone_polyline>;

  std::vector<Make_x_monotone_result> xm_results;
  for (const auto& cv : curves) make_x_monotone(cv, std::back_inserter(xm_results));
  std::vector<X_monotone_polyline> xcurves;
  for (const auto& xm_res : xm_results) {
    auto* xcv_p = std::get_if<X_monotone_polyline>(&xm_res);
    xcurves.push_back(*xcv_p);
  }
  return xcurves;
}

/*! generates the curves.
 */
std::vector<Segment> generate_segments(std::size_t n, std::size_t seed) {
  std::vector<Segment> segments;
  segments.reserve(n);
  CGAL::Random rng(seed);
  using Point_creator = CGAL::Creator_uniform_2<double, Point>;
  CGAL::Random_points_in_square_2<Point, Point_creator> p1(100.0, rng);
  CGAL::Random_points_in_square_2<Point, Point_creator> p2(100.0, rng);
  using Seg_creator = CGAL::Creator_uniform_2<Point, Segment>;
  using Segment_generator = CGAL::Join_input_iterator_2<decltype(p1), decltype(p2), Seg_creator>;
  Segment_generator gen(p1, p2);
  std::copy_n(gen, n, std::back_inserter(segments));
  return segments;
}

/*! draws the arrangement induced by the curves.
 */
void draw(const std::vector<X_monotone_polyline>& xcurves) {
  using Arrangement = CGAL::Arrangement_2<Polyline_traits>;
  Arrangement arr;
  CGAL::insert(arr, xcurves.begin(), xcurves.end());

  CGAL::Graphics_scene_options<Arrangement, typename Arrangement::Vertex_const_handle,
                                 typename Arrangement::Halfedge_const_handle, typename Arrangement::Face_const_handle>
    gso;
  gso.colored_face = [](const Arrangement&, Arrangement::Face_const_handle) -> bool { return false; };
  gso.colored_edge = [](const Arrangement&, Arrangement::Halfedge_const_handle) -> bool { return true; };
  gso.edge_color = [](const Arrangement&, Arrangement::Halfedge_const_handle fh) -> CGAL::IO::Color {
    static int i = 0;
    unsigned char r = 0;
    unsigned char g = 0;
    unsigned char b = 0;
    if (i % 3 == 0) {
      r = 0xff;
      g = 0;
      b = 0;
    }
    else if (i % 3 == 1) {
      r = 0;
      g = 0xff;
      b = 0;
    }
    else {
      r = 0;
      g = 0;
      b = 0xff;
    }
    ++i;
    return CGAL::IO::Color(r, g, b);
  };
  CGAL::draw(arr, gso, "Intersections");
}

/*!
 */
int main(int argc, char* argv[]) {
  std::size_t n = 10;
  bool do_draw = false;
  std::size_t verbose = 0;
  std::random_device rd;
  std::size_t seed = rd();
  bool closed = false;
  bool random = true;
  bool verify = false;
  std::size_t repetitions = 100000;

  try {
    // Define options
    po::options_description desc("Allowed options");
    desc.add_options()
      ("closed,c", po::value<bool>(&verify)->default_value(false), "closed")
      ("help,h", "produce help message")
      ("draw,d", po::value<bool>(&do_draw)->implicit_value(true), "draw the meshes")
      ("number,n", po::value<std::size_t>(&n)->default_value(10), "number of polylines")
      ("seed", po::value<std::size_t>(&seed), "seed")
      ("random", po::value<bool>(&random)->default_value(true), "random")
      ("repetitions,r", po::value<std::size_t>(&repetitions)->default_value(100000), "repetitions")
      ("verbose,v", po::value<std::size_t>(&verbose)->default_value(0), "set verbosity level (0 = quiet")
      ("verify", po::value<bool>(&verify)->default_value(false), "verify")
      ;

    // Parse options
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::store(po::command_line_parser(argc, argv).options(desc).run(), vm);
    po::notify(vm);

    // Help
    if (vm.count("help")) {
      std::cout << desc << "\n";
      return 0;
    }
  }
  catch (std::exception& e) {
    std::cout << "Exception: " << e.what() << "\n";
    return 1;
  }

#if 0
  {
  // Verify Single
  X_monotone_polyline xcv0, xcv1, xcv2, xcv3, xcv4, xcv5;
  traits.push_back_2_object()(xcv0, X_monotone_segment(Point(2, 1), Point(0, 2)));
  traits.push_back_2_object()(xcv1, X_monotone_segment(Point(0, 1), Point(1, 1)));
  traits.push_back_2_object()(xcv2, X_monotone_segment(Point(1, 1), Point(3, 0)));
  traits.push_back_2_object()(xcv3, X_monotone_segment(Point(3, 0), Point(1, 1)));
  traits.push_back_2_object()(xcv4, X_monotone_segment(Point(2, 1), Point(3, 1)));
  std::vector<X_monotone_polyline> xcurves = {xcv1, xcv2, xcv3, xcv4};
  if (do_draw) draw(xcurves);
  auto [match, res1] = verify_intersections(xcurves, traits, closed);
  if (! match) {
    std::cerr << "Error: do_intersect() and intersect() do not match; seed (" << seed << ")\n";
    return -1;
  }
  auto res2 = CGAL::Surface_sweep_2::do_intersect(xcurves.begin(), xcurves.end(), closed, traits);
  if (res1 != res2) {
    std::cerr << "Error: Surface_sweep_2::do_intersect() and do_intersect() do not match; seed (" << seed << ")\n";
    return -1;
  }
  return 0;
  }
#endif

  if (verify) {
    Polyline_traits traits;
    std::size_t intersected = 0;
    for (std::size_t i = 0; i < repetitions; ++i) {
      if (i > 0) seed = rd();
      auto xcurves = generate(n, seed, traits);
      if (verbose > 1) std::cout << "Generates " << n << " random polylines, seed (" << seed << ")\n";
      if (verbose > 2) for (const auto& xcv : xcurves) std::cout << xcv << "\n";
      auto xsize = xcurves.size();
      if (verbose > 1) std::cout << "Operating on " << xsize << " x-monotone polylines\n";
      if (do_draw) draw(xcurves);
      auto [match, res1] = verify_intersections(xcurves, traits, closed);
      if (res1) ++intersected;
      if (verbose > 1) std::cout << "The x-monotonoe curves do" << ((res1) ? " " : " not ") << "intersect\n";
      if (! match) {
        std::cerr << "Error: do_intersect() and intersect() do not match; seed (" << seed << ")\n";
        return -1;
      }
      auto res2 = CGAL::Surface_sweep_2::do_intersect(xcurves.begin(), xcurves.end(), closed, traits);
      if (res1 != res2) {
        std::cerr << "Error: Surface_sweep_2::do_intersect() and do_intersect() do not match; seed (" << seed << ")\n";
        return -1;
      }
    }
    if (verbose > 0) std::cout << "Out of " << repetitions << " trials, " << intersected << " intersected\n";
    return 0;
  }

  // Bench
  auto segments = generate_segments(n, seed);
  if (verbose > 1) std::cout << "Generates " << n << " random segments\n";
  auto xsize = segments.size();
  std::cout << "Test intersections in " << xsize << " x-monotone curves" << std::endl;

  std::size_t intersection_count = 0;
  Non_caching_segment_traits traits;
  auto start = std::chrono::high_resolution_clock::now();
  for (std::size_t i = 0; i < repetitions; ++i) intersection_count = count_intersections(segments, traits);
  auto diff = std::chrono::high_resolution_clock::now() - start;
  auto duration = std::chrono::duration_cast<std::chrono::microseconds>(diff);
  std::cout << "Traits Do intersection time:   " << duration.count() << " microseconds\n";
  std::cout << "Intersections found: " << intersection_count << "\n";
  std::cout << "Not intersected: " << (xsize*xsize - intersection_count) << "\n";

  Kernel kernel;
  intersection_count = 0;
  start = std::chrono::high_resolution_clock::now();
  for (std::size_t i = 0; i < repetitions; ++i) intersection_count = count_intersections(segments, kernel);
  diff = std::chrono::high_resolution_clock::now() - start;
  duration = std::chrono::duration_cast<std::chrono::microseconds>(diff);
  std::cout << "Kernel Do intersection time:   " << duration.count() << " microseconds\n";
  std::cout << "Intersections found: " << intersection_count << "\n";
  std::cout << "Not intersected: " << (xsize*xsize - intersection_count) << "\n";

  start = std::chrono::high_resolution_clock::now();
  for (std::size_t i = 0; i < repetitions; ++i) {
    intersection_count = 0;
    for (const auto& s : segments) {
      for (const auto& t : segments) {
        if (CGAL::do_intersect(s, t)) ++intersection_count;
      }
    }
  }
  diff = std::chrono::high_resolution_clock::now() - start;
  duration = std::chrono::duration_cast<std::chrono::microseconds>(diff);
  std::cout << "Kernel Do intersection time:   " << duration.count() << " microseconds\n";
  std::cout << "Intersections found: " << intersection_count << "\n";
  std::cout << "Not intersected: " << (xsize*xsize - intersection_count) << "\n";

  return 0;
}
