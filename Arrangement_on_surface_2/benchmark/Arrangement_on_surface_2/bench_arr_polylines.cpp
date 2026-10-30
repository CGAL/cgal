// This program benchmarks the polyline traits and the construction of
// arrangements of polylines.
//
// The polylines are either generated at random or read from a file.  The
// program can
// 1. construct the arrangement of the polylines (aggregated insertion via the
//    surface sweep or incremental insertion with a choice of point location),
// 2. count the number of calls to each traits functor issued during the
//    construction (using Arr_counting_traits_2), and
// 3. measure individual traits functors applied to inputs derived from the
//    polylines (e.g., pairs of x-monotone polylines with overlapping bounding
//    boxes for Intersect_2).
//
// Run with --help for the list of options.
//
// The input (and output) file format is:
// <n>                                     // number of polylines
// <m_1> <x_1> <y_1> ... <x_m_1> <y_m_1>   // polyline #1
//   :
// <m_n> <x_1> <y_1> ... <x_m_n> <y_m_n>   // polyline #n
//
// The kernel is selected at compile time via CGAL_BENCH_KERNEL:
//   1 -- Exact_predicates_exact_constructions_kernel (default)
//   2 -- Simple_cartesian<Exact_rational>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <random>
#include <set>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

#include <boost/program_options.hpp>

#ifndef CGAL_BENCH_KERNEL
#define CGAL_BENCH_KERNEL 1
#endif

#if CGAL_BENCH_KERNEL == 1
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#elif CGAL_BENCH_KERNEL == 2
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_rational.h>
#else
#error "Unknown kernel"
#endif

#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arr_polyline_traits_2.h>
#include <CGAL/Arr_counting_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_naive_point_location.h>
#include <CGAL/Arr_walk_along_line_point_location.h>
#include <CGAL/Arr_landmarks_point_location.h>
#include <CGAL/Arr_trapezoid_ric_point_location.h>
#include <CGAL/box_intersection_d.h>

#if CGAL_BENCH_KERNEL == 1
using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;
const char* kernel_name = "Epeck";
#else
using Kernel = CGAL::Simple_cartesian<CGAL::Exact_rational>;
const char* kernel_name = "Simple_cartesian<Exact_rational>";
#endif

using Segment_traits = CGAL::Arr_segment_traits_2<Kernel>;
using Traits = CGAL::Arr_polyline_traits_2<Segment_traits>;
using Point = Traits::Point_2;
using Polyline = Traits::Curve_2;
using X_monotone_polyline = Traits::X_monotone_curve_2;
using Multiplicity = Traits::Multiplicity;
using Intersection_result = std::variant<std::pair<Point, Multiplicity>, X_monotone_polyline>;
using Make_x_monotone_result = std::variant<Point, X_monotone_polyline>;

using Counting_traits = CGAL::Arr_counting_traits_2<Traits>;

namespace po = boost::program_options;
using Clock = std::chrono::steady_clock;

//
struct Options {
  // Input
  std::string input_filename;
  std::string output_filename;
  std::string type = "walk";
  std::size_t num_polylines = 100;
  std::size_t num_points = 20;
  double size = 1000;
  double step = 50;
  double snap = 0;
  unsigned int seed = 0;

  // Construction
  std::string method = "sweep";
  bool count = false;
  bool check = false;
  std::size_t insert_points = 0;

  // Functors
  std::vector<std::string> functors;
  std::size_t max_pairs = 100000;
  std::size_t queries = 4;
  std::string split_at = "random";
  bool fresh_outputs = false;
  std::size_t repeats = 3;
  bool verbose = false;
};

//
template <typename T>
std::string join(const T& names) {
  std::string res;
  for (const auto& name : names) {
    if (! res.empty()) res += ", ";
    res += name;
  }
  return res;
}

const std::vector<std::string> generator_types = { "walk", "random", "monotone" };
const std::vector<std::string> methods = { "sweep", "sweep2", "naive", "walk", "landmarks", "trapezoid", "none" };
const std::vector<std::string> split_positions = { "random", "left", "middle", "right" };
const std::vector<std::string> functor_names = {
  "make_x_monotone", "min_vertex", "max_vertex", "is_vertical", "compare_endpoints_xy", "construct_opposite",
  "compare_y_at_x", "compare_y_at_x_left", "compare_y_at_x_right", "equal", "split_joint", "split_interior",
  "are_mergeable", "merge", "do_intersect", "intersect"
};

//! Parses the command line; returns false if the program should exit.
bool parse_options(int argc, char* argv[], Options& opts, int& status) {
  po::options_description input_opts("Input options");
  input_opts.add_options()
    ("input,i", po::value<std::string>(&opts.input_filename),
     "read the polylines from the given file instead of generating them")
    ("output,o", po::value<std::string>(&opts.output_filename), "write the polylines to the given file")
    ("type,t", po::value<std::string>(&opts.type)->default_value(opts.type),
     ("generator type; one of: " + join(generator_types) + "\n"
      "  walk     -- random walk with steps of random direction and length in (0, step]\n"
      "  random   -- vertices uniformly distributed in the box (long edges)\n"
      "  monotone -- x-monotone random walk").c_str())
    ("polylines,n", po::value<std::size_t>(&opts.num_polylines)->default_value(opts.num_polylines),
     "number of polylines")
    ("points,p", po::value<std::size_t>(&opts.num_points)->default_value(opts.num_points),
     "number of points per polyline")
    ("size", po::value<double>(&opts.size)->default_value(opts.size), "the vertices lie in the box [-size, size]^2")
    ("step", po::value<double>(&opts.step)->default_value(opts.step),
     "maximum step length of the walk and monotone generators")
    ("snap", po::value<double>(&opts.snap)->default_value(opts.snap),
     "if positive, round the coordinates to multiples of the given value (induces degeneracies)")
    ("seed,s", po::value<unsigned int>(&opts.seed)->default_value(opts.seed), "random seed");

  po::options_description arr_opts("Arrangement options");
  arr_opts.add_options()
    ("method,m", po::value<std::string>(&opts.method)->default_value(opts.method),
     ("construction method; one of: " + join(methods) + "\n"
      "  sweep  -- aggregated insertion\n"
      "  sweep2 -- aggregated insertion of the first half of the polylines, and then of the second half (into a non-empty "
      "arrangement)\n"
      "  the others are incremental insertion with the respective point location").c_str())
    ("count,c", po::bool_switch(&opts.count), "count the traits-functor calls during the construction")
    ("check", po::bool_switch(&opts.check), "check the validity of the arrangement")
    ("insert-points", po::value<std::size_t>(&opts.insert_points)->default_value(opts.insert_points),
     "insert points into the interiors of (up to) the given number of edges of the arrangement constructed by the "
     "sweep, one point in the interior of a random subcurve of every edge; every insertion splits the curve of the "
     "edge into new output polylines");

  po::options_description functor_opts("Functor options");
  functor_opts.add_options()
    ("functor,f", po::value<std::vector<std::string>>(&opts.functors)->composing(),
     ("benchmark the given functor (can be repeated); \"all\" for all; one of:\n  " + join(functor_names)).c_str())
    ("max-pairs", po::value<std::size_t>(&opts.max_pairs)->default_value(opts.max_pairs),
     "maximum number of x-monotone polyline pairs used by pair functors")
    ("queries", po::value<std::size_t>(&opts.queries)->default_value(opts.queries),
     "number of query points per x-monotone polyline for compare_y_at_x")
    ("split-at", po::value<std::string>(&opts.split_at)->default_value(opts.split_at),
     ("position of the split point along the x-monotone polyline (by rank of the subcurve, from left to right); "
      "one of: " + join(split_positions)).c_str())
    ("fresh-outputs", po::bool_switch(&opts.fresh_outputs),
     "split into new output polylines in every call (by default, the output polylines are reused, retaining their "
     "capacity, like the surface sweep does)")
    ("repeats,r", po::value<std::size_t>(&opts.repeats)->default_value(opts.repeats),
     "number of repetitions of each timed operation; the minimum time is reported");

  po::options_description general_opts("General options");
  general_opts.add_options()
    ("help,h", "print this help message")
    ("verbose,v", po::bool_switch(&opts.verbose), "print additional information");

  po::options_description all_opts("Usage: bench_arr_polylines [options]\nOptions");
  all_opts.add(general_opts).add(input_opts).add(arr_opts).add(functor_opts);

  po::variables_map vm;
  try {
    po::store(po::parse_command_line(argc, argv, all_opts), vm);
    po::notify(vm);
  }
  catch (const std::exception& e) {
    std::cerr << e.what() << std::endl;
    status = 1;
    return false;
  }

  if (vm.count("help")) {
    std::cout << all_opts << std::endl;
    status = 0;
    return false;
  }

  auto check = [&](const std::string& value, const std::vector<std::string>& valid, const char* what) {
    if (std::find(valid.begin(), valid.end(), value) != valid.end()) return true;
    std::cerr << "Invalid " << what << ": " << value << " (valid: " << join(valid) << ")" << std::endl;
    return false;
  };
  if (! check(opts.type, generator_types, "generator type") || ! check(opts.method, methods, "method") ||
      ! check(opts.split_at, split_positions, "split position")) {
    status = 1;
    return false;
  }
  std::vector<std::string> functors;
  for (const auto& f : opts.functors) {
    if (f == "all") functors.insert(functors.end(), functor_names.begin(), functor_names.end());
    else if (! check(f, functor_names, "functor")) {
      status = 1;
      return false;
    }
    else functors.push_back(f);
  }
  opts.functors = functors;
  if (opts.repeats == 0) opts.repeats = 1;
  if (opts.num_points < 2) {
    std::cerr << "The number of points per polyline must be at least 2" << std::endl;
    status = 1;
    return false;
  }
  return true;
}

//! Generates random polylines.
std::vector<Polyline> generate_polylines(const Options& opts, const Traits& traits) {
  std::mt19937 gen(opts.seed);
  std::uniform_real_distribution<double> coord(-opts.size, opts.size);
  std::uniform_real_distribution<double> angle(-CGAL_PI, CGAL_PI);
  std::uniform_real_distribution<double> unit(0, 1);

  auto clamp = [&](double v) { return std::clamp(v, -opts.size, opts.size); };
  auto snap = [&](double v) { return (opts.snap > 0) ? std::round(v / opts.snap) * opts.snap : v; };

  auto ctr_curve = traits.construct_curve_2_object();
  std::vector<Polyline> polylines;
  polylines.reserve(opts.num_polylines);
  std::vector<Point> pts;
  while (polylines.size() < opts.num_polylines) {
    pts.clear();
    double x = coord(gen);
    double y = coord(gen);
    for (std::size_t i = 0; i < opts.num_points; ++i) {
      if (i > 0) {
        if (opts.type == "random") {
          x = coord(gen);
          y = coord(gen);
        }
        else if (opts.type == "walk") {
          double a = angle(gen);
          double l = opts.step * (1 - unit(gen));
          x = clamp(x + l * std::cos(a));
          y = clamp(y + l * std::sin(a));
        }
        else {
          // monotone; when reaching the right side of the box, the polyline ends.
          x += opts.step * (1 - unit(gen));
          if (x > opts.size) break;
          y = clamp(y + opts.step * (2 * unit(gen) - 1));
        }
      }
      Point p(snap(x), snap(y));
      // Skip vertices that coincide with the previous vertex (possible after snapping or clamping).
      if (! pts.empty() && (pts.back() == p)) continue;
      pts.push_back(p);
    }
    if (pts.size() < 2) continue;
    polylines.push_back(ctr_curve(pts.begin(), pts.end()));
  }
  return polylines;
}

//! Reads polylines from a file.
bool read_polylines(const std::string& filename, const Traits& traits, std::vector<Polyline>& polylines) {
  std::ifstream is(filename);
  if (! is.is_open()) {
    std::cerr << "Failed to open " << filename << std::endl;
    return false;
  }
  auto ctr_curve = traits.construct_curve_2_object();
  std::size_t n;
  is >> n;
  std::vector<Point> pts;
  for (std::size_t i = 0; i < n; ++i) {
    std::size_t m;
    is >> m;
    pts.clear();
    for (std::size_t j = 0; j < m; ++j) {
      double x, y;
      is >> x >> y;
      pts.emplace_back(x, y);
    }
    if (! is) {
      std::cerr << "Failed to read polyline #" << i << " from " << filename << std::endl;
      return false;
    }
    polylines.push_back(ctr_curve(pts.begin(), pts.end()));
  }
  return true;
}

//! Writes polylines to a file.
bool write_polylines(const std::string& filename, const std::vector<Polyline>& polylines) {
  std::ofstream os(filename);
  if (! os.is_open()) {
    std::cerr << "Failed to open " << filename << std::endl;
    return false;
  }
  os << std::setprecision(17) << polylines.size() << std::endl;
  for (const auto& pl : polylines) {
    os << std::distance(pl.points_begin(), pl.points_end());
    for (auto it = pl.points_begin(); it != pl.points_end(); ++it)
      os << " " << CGAL::to_double(it->x()) << " " << CGAL::to_double(it->y());
    os << std::endl;
  }
  return true;
}

//! Runs `f` `repeats` times and returns the minimum duration in seconds.
template <typename F>
double time_min(std::size_t repeats, F f) {
  double best = std::numeric_limits<double>::max();
  for (std::size_t i = 0; i < repeats; ++i) {
    auto start = Clock::now();
    f();
    std::chrono::duration<double> d = Clock::now() - start;
    best = std::min(best, d.count());
  }
  return best;
}

//! Constructs the arrangement of the polylines with the given method.
template <typename Arrangement>
void construct(Arrangement& arr, const std::vector<Polyline>& polylines, const std::string& method) {
  if (method == "sweep") {
    CGAL::insert(arr, polylines.begin(), polylines.end());
    return;
  }
  if (method == "sweep2") {
    auto mid = polylines.begin() + polylines.size() / 2;
    CGAL::insert(arr, polylines.begin(), mid);
    CGAL::insert(arr, mid, polylines.end());
    return;
  }
  auto insert_all = [&](auto& pl) { for (const auto& cv : polylines) CGAL::insert(arr, cv, pl); };
  if (method == "naive") {
    CGAL::Arr_naive_point_location<Arrangement> pl(arr);
    insert_all(pl);
  }
  else if (method == "walk") {
    CGAL::Arr_walk_along_line_point_location<Arrangement> pl(arr);
    insert_all(pl);
  }
  else if (method == "landmarks") {
    CGAL::Arr_landmarks_point_location<Arrangement> pl(arr);
    insert_all(pl);
  }
  else if (method == "trapezoid") {
    // Td_traits copies the traits, but Arr_counting_traits_2 is not copyable.
    if constexpr (std::is_copy_constructible_v<typename Arrangement::Geometry_traits_2>) {
      CGAL::Arr_trapezoid_ric_point_location<Arrangement> pl(arr);
      insert_all(pl);
    }
    else std::cerr << "The trapezoid point location does not support the counting traits" << std::endl;
  }
}

//! Benchmarks the construction of the arrangement.
void bench_construction(const Options& opts, const Traits& traits, const std::vector<Polyline>& polylines) {
  using Arrangement = CGAL::Arrangement_2<Traits>;
  std::size_t nv = 0, ne = 0, nf = 0;
  bool valid = true;
  auto secs = time_min(opts.repeats, [&]() {
    Arrangement arr(&traits);
    construct(arr, polylines, opts.method);
    nv = arr.number_of_vertices();
    ne = arr.number_of_edges();
    nf = arr.number_of_faces();
    if (opts.check) valid = arr.is_valid();
  });
  std::cout << "Construction (" << opts.method << "): " << secs << " s" << std::endl
            << "  V = " << nv << ", E = " << ne << ", F = " << nf << std::endl;
  if (opts.check) std::cout << "  Validity: " << (valid ? "OK" : "INVALID") << std::endl;

  if (opts.count && (opts.method == "trapezoid"))
    std::cerr << "Counting is not supported with the trapezoid point location" << std::endl;
  else if (opts.count) {
    using Counting_arrangement = CGAL::Arrangement_2<Counting_traits>;
    Counting_traits cnt_traits;
    Counting_arrangement arr(&cnt_traits);
    construct(arr, polylines, opts.method);
    std::cout << "Traits-functor calls (" << opts.method << "):" << std::endl << cnt_traits;
  }
}

//! A point location that reports a given halfedge; it isolates the cost of
//! splitting edges by insert_point() from the cost of point location.
template <typename Arrangement>
class Given_halfedge_point_location {
public:
  using Halfedge_const_handle = typename Arrangement::Halfedge_const_handle;
  using Result_type = typename CGAL::Arr_point_location_result<Arrangement>::Type;

  Halfedge_const_handle he;

  Result_type locate(const Point&) const { return Result_type(he); }
};

//! Benchmarks the insertion of points into the interiors of edges.
void bench_point_insertion(const Options& opts, const Traits& traits, const std::vector<Polyline>& polylines) {
  using Arrangement = CGAL::Arrangement_2<Traits>;
  double best = std::numeric_limits<double>::max();
  std::size_t num_points = 0, num_subcurves = 0, nv = 0;
  bool valid = true;
  for (std::size_t r = 0; r < opts.repeats; ++r) {
    Arrangement arr(&traits);
    construct(arr, polylines, "sweep");

    // Select the edges and the points; the same ones in every repetition.
    std::mt19937 gen(opts.seed + 2);
    std::vector<std::pair<typename Arrangement::Halfedge_const_handle, Point>> targets;
    for (auto eit = arr.edges_begin(); eit != arr.edges_end(); ++eit) {
      if (targets.size() == opts.insert_points) break;
      const auto& xcv = eit->curve();
      const auto& seg = xcv[std::uniform_int_distribution<std::size_t>(0, xcv.number_of_subcurves() - 1)(gen)];
      targets.emplace_back(eit, CGAL::midpoint(seg.source(), seg.target()));
    }
    num_points = targets.size();
    num_subcurves = 0;
    for (const auto& [he, p] : targets) num_subcurves += he->curve().number_of_subcurves();

    Given_halfedge_point_location<Arrangement> pl;
    auto start = Clock::now();
    for (const auto& [he, p] : targets) {
      pl.he = he;
      CGAL::insert_point(arr, p, pl);
    }
    std::chrono::duration<double> d = Clock::now() - start;
    best = std::min(best, d.count());
    nv = arr.number_of_vertices();
    if (opts.check) valid = arr.is_valid();
  }
  std::cout << "Point insertion: " << best << " s" << std::endl
            << "  # points: " << num_points << ", avg. # subcurves of a split edge: "
            << (num_points ? double(num_subcurves) / double(num_points) : 0.0) << ", V = " << nv << std::endl;
  if (opts.check) std::cout << "  Validity: " << (valid ? "OK" : "INVALID") << std::endl;
}

//! The data derived from the polylines used by the functor benchmarks.
struct Functor_data {
  std::vector<X_monotone_polyline> xcvs;
  std::vector<std::pair<std::size_t, std::size_t>> pairs;               // indices into xcvs with overlapping bboxes
  std::vector<std::pair<Point, std::size_t>> queries;                   // query point and index into xcvs
  std::vector<std::pair<Point, std::size_t>> joints;                    // interior joint and index into xcvs
  std::vector<std::pair<Point, std::size_t>> interiors;                 // subcurve interior point and index into xcvs
  std::vector<std::pair<X_monotone_polyline, X_monotone_polyline>> mergeable;
  std::vector<std::tuple<Point, std::size_t, std::size_t>> left_queries;  // common point, curves defined to its left
  std::vector<std::tuple<Point, std::size_t, std::size_t>> right_queries; // common point, curves defined to its right
};

//! Prepares the data for the functor benchmarks.
void prepare(const Options& opts, const Traits& traits, const std::vector<Polyline>& polylines, Functor_data& data) {
  std::mt19937 gen(opts.seed + 1);
  std::uniform_real_distribution<double> unit(0, 1);

  // X-monotone polylines.
  auto make_x_monotone = traits.make_x_monotone_2_object();
  std::vector<Make_x_monotone_result> objs;
  for (const auto& cv : polylines) make_x_monotone(cv, std::back_inserter(objs));
  for (const auto& obj : objs)
    if (const auto* xcv = std::get_if<X_monotone_polyline>(&obj)) data.xcvs.push_back(*xcv);

  // Pairs of x-monotone polylines with overlapping bounding boxes; at most `max_pairs` (reservoir sampling).
  using Box = CGAL::Box_intersection_d::Box_with_info_d<double, 2, std::size_t>;
  std::vector<Box> boxes;
  boxes.reserve(data.xcvs.size());
  for (std::size_t i = 0; i < data.xcvs.size(); ++i) boxes.emplace_back(data.xcvs[i].bbox(), i);
  std::size_t num_overlapping = 0;
  auto callback = [&](const Box& b1, const Box& b2) {
    ++num_overlapping;
    std::pair<std::size_t, std::size_t> p(b1.info(), b2.info());
    if (data.pairs.size() < opts.max_pairs) {
      data.pairs.push_back(p);
      return;
    }
    std::uniform_int_distribution<std::size_t> idx(0, num_overlapping - 1);
    auto k = idx(gen);
    if (k < opts.max_pairs) data.pairs[k] = p;
  };
  CGAL::box_self_intersection_d<CGAL::Parallel_if_available_tag>(boxes.begin(), boxes.end(), callback);
  std::sort(data.pairs.begin(), data.pairs.end());

  // Query points in the x-range of each x-monotone polyline.
  for (std::size_t i = 0; i < data.xcvs.size(); ++i) {
    const auto& xcv = data.xcvs[i];
    auto bbox = xcv.bbox();
    double xmin = CGAL::to_double(xcv.points_begin()->x());
    double xmax = CGAL::to_double(std::prev(xcv.points_end())->x());
    if (xmin > xmax) std::swap(xmin, xmax);
    for (std::size_t j = 0; j < opts.queries; ++j) {
      double x = xmin + (xmax - xmin) * unit(gen);
      double y = bbox.ymin() + (bbox.ymax() - bbox.ymin()) * unit(gen);
      data.queries.emplace_back(Point(x, y), i);
    }
  }

  // Split points: interior joints and subcurve interior points. The subcurve
  // is selected by its rank in left-to-right order, in [first, n-1].
  auto rank = [&](std::size_t first, std::size_t n) -> std::size_t {
    if (opts.split_at == "left") return first;
    if (opts.split_at == "right") return n - 1;
    if (opts.split_at == "middle") return std::max(first, n / 2);
    std::uniform_int_distribution<std::size_t> idx(first, n - 1);
    return idx(gen);
  };
  for (std::size_t i = 0; i < data.xcvs.size(); ++i) {
    const auto& xcv = data.xcvs[i];
    auto n = xcv.number_of_subcurves();
    auto index = [&](std::size_t r) { return xcv[0].is_directed_right() ? r : n - 1 - r; };
    // The left endpoint of the subcurve of rank r > 0 is an interior joint.
    if (n > 1) data.joints.emplace_back(xcv[index(rank(1, n))].left(), i);
    const auto& seg = xcv[index(rank(0, n))];
    data.interiors.emplace_back(CGAL::midpoint(seg.source(), seg.target()), i);
  }

  // Mergeable pairs: the results of splitting.
  auto split = traits.split_2_object();
  for (const auto& [p, i] : data.joints) {
    X_monotone_polyline xcv1, xcv2;
    split(data.xcvs[i], p, xcv1, xcv2);
    data.mergeable.emplace_back(std::move(xcv1), std::move(xcv2));
  }

  // Compare_y_at_x_left/right queries: transversal intersection points of the pairs.
  auto intersect = traits.intersect_2_object();
  auto min_vertex = traits.construct_min_vertex_2_object();
  auto max_vertex = traits.construct_max_vertex_2_object();
  auto equal = traits.equal_2_object();
  std::vector<Intersection_result> res;
  for (const auto& [i, j] : data.pairs) {
    res.clear();
    intersect(data.xcvs[i], data.xcvs[j], std::back_inserter(res));
    for (const auto& r : res) {
      const auto* ip = std::get_if<std::pair<Point, Multiplicity>>(&r);
      if (! ip) continue;
      const auto& p = ip->first;
      if (! equal(p, min_vertex(data.xcvs[i])) && ! equal(p, min_vertex(data.xcvs[j])))
        data.left_queries.emplace_back(p, i, j);
      if (! equal(p, max_vertex(data.xcvs[i])) && ! equal(p, max_vertex(data.xcvs[j])))
        data.right_queries.emplace_back(p, i, j);
    }
  }

  if (opts.verbose) {
    std::size_t num_subcurves = 0;
    for (const auto& xcv : data.xcvs) num_subcurves += xcv.number_of_subcurves();
    std::cout << "Functor data:" << std::endl
              << "  # x-monotone polylines: " << data.xcvs.size() << " (" << num_subcurves << " subcurves)" << std::endl
              << "  # overlapping bbox pairs: " << num_overlapping << " (" << data.pairs.size() << " used)" << std::endl
              << "  # compare_y_at_x_left/right queries: " << data.left_queries.size() << "/"
              << data.right_queries.size() << std::endl;
  }
}

//! Benchmarks the requested functors.
void bench_functors(const Options& opts, const Traits& traits, const std::vector<Polyline>& polylines) {
  Functor_data data;
  prepare(opts, traits, polylines, data);

  std::cout << std::left << std::setw(24) << "Functor" << std::right << std::setw(12) << "calls"
            << std::setw(14) << "total [ms]" << std::setw(14) << "per call [ns]" << std::endl;

  // `f` performs all calls and returns a checksum, which prevents the compiler from eliminating them.
  std::size_t checksum = 0;
  const void* volatile sink = nullptr;
  auto run = [&](const std::string& name, std::size_t calls, auto f) {
    auto secs = time_min(opts.repeats, [&]() { checksum += f(); });
    std::cout << std::left << std::setw(24) << name << std::right << std::setw(12) << calls
              << std::setw(14) << std::fixed << std::setprecision(3) << secs * 1e3
              << std::setw(14) << std::setprecision(1) << ((calls > 0) ? secs * 1e9 / calls : 0.0)
              << std::defaultfloat << std::setprecision(6) << std::endl;
  };

  const auto& xcvs = data.xcvs;
  for (const auto& name : opts.functors) {
    if (name == "make_x_monotone") {
      auto make_x_monotone = traits.make_x_monotone_2_object();
      run(name, polylines.size(), [&]() {
        std::vector<Make_x_monotone_result> objs;
        for (const auto& cv : polylines) make_x_monotone(cv, std::back_inserter(objs));
        return objs.size();
      });
    }
    else if (name == "min_vertex") {
      auto min_vertex = traits.construct_min_vertex_2_object();
      run(name, xcvs.size(), [&]() {
        std::size_t s = 0;
        for (const auto& xcv : xcvs) sink = &min_vertex(xcv);
        return s;
      });
    }
    else if (name == "max_vertex") {
      auto max_vertex = traits.construct_max_vertex_2_object();
      run(name, xcvs.size(), [&]() {
        std::size_t s = 0;
        for (const auto& xcv : xcvs) sink = &max_vertex(xcv);
        return s;
      });
    }
    else if (name == "is_vertical") {
      auto is_vertical = traits.is_vertical_2_object();
      run(name, xcvs.size(), [&]() {
        std::size_t s = 0;
        for (const auto& xcv : xcvs) s += is_vertical(xcv);
        return s;
      });
    }
    else if (name == "compare_endpoints_xy") {
      auto cmp_endpoints = traits.compare_endpoints_xy_2_object();
      run(name, xcvs.size(), [&]() {
        std::size_t s = 0;
        for (const auto& xcv : xcvs) s += (cmp_endpoints(xcv) == CGAL::SMALLER);
        return s;
      });
    }
    else if (name == "construct_opposite") {
      auto opposite = traits.construct_opposite_2_object();
      run(name, xcvs.size(), [&]() {
        std::size_t s = 0;
        for (const auto& xcv : xcvs) s += opposite(xcv).number_of_subcurves();
        return s;
      });
    }
    else if (name == "compare_y_at_x") {
      auto cmp_y_at_x = traits.compare_y_at_x_2_object();
      run(name, data.queries.size(), [&]() {
        std::size_t s = 0;
        for (const auto& [p, i] : data.queries) s += (cmp_y_at_x(p, xcvs[i]) == CGAL::SMALLER);
        return s;
      });
    }
    else if (name == "compare_y_at_x_left") {
      auto cmp_y_at_x_left = traits.compare_y_at_x_left_2_object();
      run(name, data.left_queries.size(), [&]() {
        std::size_t s = 0;
        for (const auto& [p, i, j] : data.left_queries) s += (cmp_y_at_x_left(xcvs[i], xcvs[j], p) == CGAL::SMALLER);
        return s;
      });
    }
    else if (name == "compare_y_at_x_right") {
      auto cmp_y_at_x_right = traits.compare_y_at_x_right_2_object();
      run(name, data.right_queries.size(), [&]() {
        std::size_t s = 0;
        for (const auto& [p, i, j] : data.right_queries) s += (cmp_y_at_x_right(xcvs[i], xcvs[j], p) == CGAL::SMALLER);
        return s;
      });
    }
    else if (name == "equal") {
      // Compare each x-monotone polyline with a (deep) copy of itself, which is the expensive case.
      std::vector<X_monotone_polyline> copies;
      copies.reserve(xcvs.size());
      for (const auto& xcv : xcvs) copies.emplace_back(xcv.subcurves_begin(), xcv.subcurves_end());
      auto equal = traits.equal_2_object();
      run(name, xcvs.size(), [&]() {
        std::size_t s = 0;
        for (std::size_t i = 0; i < xcvs.size(); ++i) s += equal(xcvs[i], copies[i]);
        return s;
      });
    }
    else if ((name == "split_joint") || (name == "split_interior")) {
      const auto& queries = (name == "split_joint") ? data.joints : data.interiors;
      auto split = traits.split_2_object();
      run(name, queries.size(), [&]() {
        std::size_t s = 0;
        X_monotone_polyline xcv1, xcv2;
        for (const auto& [p, i] : queries) {
          if (opts.fresh_outputs) {
            X_monotone_polyline fresh1, fresh2;
            split(xcvs[i], p, fresh1, fresh2);
            s += fresh1.number_of_subcurves();
            continue;
          }
          split(xcvs[i], p, xcv1, xcv2);
          s += xcv1.number_of_subcurves();
        }
        return s;
      });
    }
    else if (name == "are_mergeable") {
      auto are_mergeable = traits.are_mergeable_2_object();
      run(name, data.mergeable.size(), [&]() {
        std::size_t s = 0;
        for (const auto& [xcv1, xcv2] : data.mergeable) s += are_mergeable(xcv1, xcv2);
        return s;
      });
    }
    else if (name == "merge") {
      auto merge = traits.merge_2_object();
      run(name, data.mergeable.size(), [&]() {
        std::size_t s = 0;
        X_monotone_polyline xcv;
        for (const auto& [xcv1, xcv2] : data.mergeable) {
          merge(xcv1, xcv2, xcv);
          s += xcv.number_of_subcurves();
        }
        return s;
      });
    }
    else if (name == "do_intersect") {
      auto do_intersect = traits.do_intersect_2_object();
      run(name, data.pairs.size(), [&]() {
        std::size_t s = 0;
        for (const auto& [i, j] : data.pairs) s += do_intersect(xcvs[i], xcvs[j]);
        return s;
      });
    }
    else if (name == "intersect") {
      auto intersect = traits.intersect_2_object();
      run(name, data.pairs.size(), [&]() {
        std::vector<Intersection_result> res;
        std::size_t s = 0;
        for (const auto& [i, j] : data.pairs) {
          res.clear();
          intersect(xcvs[i], xcvs[j], std::back_inserter(res));
          s += res.size();
        }
        return s;
      });
    }
  }
  if (opts.verbose) std::cout << "Checksum: " << checksum << std::endl;
}

int main(int argc, char* argv[]) {
  Options opts;
  int status;
  if (! parse_options(argc, argv, opts, status)) return status;

  Traits traits;
  std::vector<Polyline> polylines;
  if (! opts.input_filename.empty()) {
    if (! read_polylines(opts.input_filename, traits, polylines)) return 1;
  }
  else polylines = generate_polylines(opts, traits);
  if (! opts.output_filename.empty() && ! write_polylines(opts.output_filename, polylines)) return 1;

  std::size_t num_subcurves = 0;
  for (const auto& pl : polylines) num_subcurves += pl.number_of_subcurves();
  std::cout << "Kernel: " << kernel_name << std::endl;
  if (opts.input_filename.empty()) std::cout << "Generator: " << opts.type << ", seed: " << opts.seed << std::endl;
  else std::cout << "Input: " << opts.input_filename << std::endl;
  std::cout << "# polylines: " << polylines.size() << ", # segments: " << num_subcurves << std::endl;

  if (opts.method != "none") bench_construction(opts, traits, polylines);
  if (opts.insert_points > 0) bench_point_insertion(opts, traits, polylines);
  if (! opts.functors.empty()) bench_functors(opts, traits, polylines);
  return 0;
}
