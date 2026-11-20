// Copyright (c) 2026 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Benchmark for the 2D Arrangements package.
//
// The program measures the operations whose code paths go through the
// geometry traits held by an arrangement: construction, copying and
// assignment, insertion (aggregate and incremental), removal, point location,
// zone computation, overlay, and the arrangement-with-history operations.
// It is written so that the very same source compiles against a CGAL tree
// with or without the shared-pointer geometry-traits API; benchmarks that
// require the new API are skipped when it is absent. Build it once against
// each tree, run both executables with the same options, and compare the
// CSV files with compare_benchmarks.py.
//
// Every benchmark returns a checksum derived from the result of the timed
// work (e.g., the numbers of vertices, edges and faces). Two builds that
// perform identical work produce identical checksums, which guards against
// comparing timings of different computations.

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <memory>
#include <numeric>
#include <sstream>
#include <string>
#include <type_traits>
#include <variant>
#include <vector>

#include <CGAL/config.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Random.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arrangement_with_history_2.h>
#include <CGAL/Arr_overlay_2.h>
#include <CGAL/Arr_naive_point_location.h>
#include <CGAL/Arr_walk_along_line_point_location.h>
#include <CGAL/Arr_landmarks_point_location.h>
#include <CGAL/Arr_trapezoid_ric_point_location.h>

using Kernel     = CGAL::Exact_predicates_exact_constructions_kernel;
using Traits     = CGAL::Arr_segment_traits_2<Kernel>;
using Point      = Traits::Point_2;
using Segment    = Traits::Curve_2;
using Arr        = CGAL::Arrangement_2<Traits>;
using Aos        = Arr::Base;                   // Arrangement_on_surface_2
using Harr       = CGAL::Arrangement_with_history_2<Traits>;
using Shared_tr  = std::shared_ptr<const Traits>;

// Whether the arrangement type A offers the shared-pointer constructor.
template <typename A>
constexpr bool has_shared_ctor = std::is_constructible_v<A, Shared_tr>;

//-----------------------------------------------------------------------------
// Options
//-----------------------------------------------------------------------------
struct Options {
  unsigned int seed = 1;
  std::size_t  repetitions = 5;
  std::size_t  warmup = 1;
  std::size_t  segments = 1000;
  std::size_t  queries = 2000;
  std::size_t  small = 100000;    // number of objects in construction benches
  double       length = 0.3;     // maximal segment length (square is [-1,1]^2)
  std::string  filter;           // run only benchmarks whose name contains it
  std::string  csv;              // optional CSV output file
  bool         list = false;
};

static void usage(const char* prog) {
  std::cout <<
    "Usage: " << prog << " [options]\n"
    "  -s, --seed <n>          random seed (default 1)\n"
    "  -r, --repetitions <n>   timed repetitions per benchmark (default 5)\n"
    "  -w, --warmup <n>        untimed warm-up runs per benchmark (default 1)\n"
    "  -n, --segments <n>      number of input segments (default 1000)\n"
    "  -q, --queries <n>       number of point-location queries (default 2000;\n"
    "                          the naive strategy uses a tenth of them)\n"
    "  -m, --small <n>         objects created in construction benches (default 100000)\n"
    "  -l, --length <x>        maximal segment length; the input lies in\n"
    "                          [-1,1]^2, so this controls the number of\n"
    "                          intersections (default 0.3)\n"
    "  -f, --filter <str>      run only benchmarks whose name contains <str>\n"
    "  -o, --csv <file>        also write the results as CSV to <file>\n"
    "      --list              list the benchmarks and exit\n"
    "  -h, --help              print this message\n";
}

static bool parse(int argc, char* argv[], Options& opt) {
  auto value = [&](int& i) -> const char* {
    if (i + 1 >= argc) {
      std::cerr << "Missing value for " << argv[i] << "\n";
      std::exit(EXIT_FAILURE);
    }
    return argv[++i];
  };
  for (int i = 1; i < argc; ++i) {
    std::string a = argv[i];
    if (a == "-h" || a == "--help") { usage(argv[0]); std::exit(EXIT_SUCCESS); }
    else if (a == "-s" || a == "--seed") opt.seed = std::stoul(value(i));
    else if (a == "-r" || a == "--repetitions") opt.repetitions = std::stoul(value(i));
    else if (a == "-w" || a == "--warmup") opt.warmup = std::stoul(value(i));
    else if (a == "-n" || a == "--segments") opt.segments = std::stoul(value(i));
    else if (a == "-q" || a == "--queries") opt.queries = std::stoul(value(i));
    else if (a == "-m" || a == "--small") opt.small = std::stoul(value(i));
    else if (a == "-l" || a == "--length") opt.length = std::stod(value(i));
    else if (a == "-f" || a == "--filter") opt.filter = value(i);
    else if (a == "-o" || a == "--csv") opt.csv = value(i);
    else if (a == "--list") opt.list = true;
    else {
      std::cerr << "Unknown option " << a << "\n";
      usage(argv[0]);
      return false;
    }
  }
  if (opt.repetitions == 0 || opt.segments < 4 || opt.length <= 0.0) {
    std::cerr << "Invalid option values\n";
    return false;
  }
  return true;
}

//-----------------------------------------------------------------------------
// Input generation (Generator package)
//-----------------------------------------------------------------------------
// Segments are generated by drawing a source point uniformly in [-1,1]^2 and
// a displacement uniformly in a disc of radius `length`. Short segments keep
// the number of intersections, and hence the arrangement size, under control.
static std::vector<Segment>
random_segments(std::size_t n, double length, CGAL::Random& rnd) {
  using Kpoint = Kernel::Point_2;
  CGAL::Random_points_in_square_2<Kpoint> sources(1.0, rnd);
  CGAL::Random_points_in_disc_2<Kpoint>   offsets(length, rnd);
  std::vector<Segment> segs;
  segs.reserve(n);
  while (segs.size() < n) {
    const Kpoint p = *sources++;
    const Kpoint d = *offsets++;
    const Kpoint q(p.x() + d.x(), p.y() + d.y());
    if (p == q) continue;                     // skip degenerate segments
    segs.emplace_back(p, q);
  }
  return segs;
}

static std::vector<Point> random_points(std::size_t n, CGAL::Random& rnd) {
  CGAL::Random_points_in_square_2<Point> gen(1.0, rnd);
  std::vector<Point> pts;
  pts.reserve(n);
  std::copy_n(gen, n, std::back_inserter(pts));
  return pts;
}

//-----------------------------------------------------------------------------
// Timing infrastructure
//-----------------------------------------------------------------------------
// A benchmark brackets the work it wants measured with start()/stop(), so
// that per-repetition setup (e.g., copying an input arrangement) is excluded.
class Stopwatch {
  using Clock = std::chrono::steady_clock;
  Clock::time_point m_start;
  double m_elapsed = 0.0;   // milliseconds
public:
  void start() { m_start = Clock::now(); }
  void stop() {
    m_elapsed +=
      std::chrono::duration<double, std::milli>(Clock::now() - m_start).count();
  }
  double elapsed() const { return m_elapsed; }
};

using Bench_fn = std::function<std::uint64_t(Stopwatch&)>;

struct Benchmark {
  std::string name;
  std::string description;
  Bench_fn    fn;
};

struct Result {
  std::string         name;
  std::vector<double> times;   // milliseconds
  std::uint64_t       checksum = 0;
  bool                consistent = true;  // same checksum in every run
};

static double median(std::vector<double> v) {
  std::sort(v.begin(), v.end());
  const std::size_t n = v.size();
  return (n % 2) ? v[n / 2] : 0.5 * (v[n / 2 - 1] + v[n / 2]);
}

static double mean(const std::vector<double>& v)
{ return std::accumulate(v.begin(), v.end(), 0.0) / v.size(); }

static double stddev(const std::vector<double>& v) {
  if (v.size() < 2) return 0.0;
  const double m = mean(v);
  double s = 0.0;
  for (double x : v) s += (x - m) * (x - m);
  return std::sqrt(s / (v.size() - 1));
}

static Result run(const Benchmark& b, const Options& opt) {
  Result r;
  r.name = b.name;
  for (std::size_t i = 0; i < opt.warmup + opt.repetitions; ++i) {
    Stopwatch sw;
    const std::uint64_t cs = b.fn(sw);
    if (i == 0) r.checksum = cs;
    else if (cs != r.checksum) r.consistent = false;
    if (i >= opt.warmup) r.times.push_back(sw.elapsed());
  }
  return r;
}

//-----------------------------------------------------------------------------
// Helpers
//-----------------------------------------------------------------------------
// Combine counts into a checksum (order-dependent, FNV-1a style).
static std::uint64_t mix(std::uint64_t h, std::uint64_t v) {
  h ^= v + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
  return h;
}

template <typename A>
static std::uint64_t size_checksum(const A& arr) {
  std::uint64_t h = 0;
  h = mix(h, arr.number_of_vertices());
  h = mix(h, arr.number_of_edges());
  h = mix(h, arr.number_of_faces());
  return h;
}

// Prevent the optimizer from discarding results.
static volatile std::uint64_t g_sink = 0;

//-----------------------------------------------------------------------------
// Benchmarks
//-----------------------------------------------------------------------------
struct Input {
  Traits               traits;     // caller-owned traits (legacy path)
  std::vector<Segment> segments;
  std::vector<Segment> query_segments;
  std::vector<Point>   queries;
  Arr                  arr;        // arrangement of all segments
  Harr                 harr;       // arrangement with history of all segments
  Arr                  small_arr;  // arrangement of a few segments

  Input() : arr(&traits), harr(&traits), small_arr(&traits) {}
};

// Construction benchmarks, instantiated for several arrangement types.
template <typename A>
static void add_construction_benches(std::vector<Benchmark>& benches,
                                     const std::string& tag,
                                     Input& in, const Options& opt) {
  benches.push_back({"construct_default_" + tag,
    "construct and destroy default arrangements (traits owned)",
    [&opt](Stopwatch& sw) -> std::uint64_t {
      std::uint64_t h = 0;
      sw.start();
      for (std::size_t i = 0; i < opt.small; ++i) {
        A arr;
        h += arr.number_of_faces();
      }
      sw.stop();
      return h;
    }});

  benches.push_back({"construct_raw_" + tag,
    "construct and destroy arrangements from a raw traits pointer",
    [&opt, &in](Stopwatch& sw) -> std::uint64_t {
      std::uint64_t h = 0;
      sw.start();
      for (std::size_t i = 0; i < opt.small; ++i) {
        A arr(&in.traits);
        h += arr.number_of_faces();
      }
      sw.stop();
      return h;
    }});

  if constexpr (has_shared_ctor<A>) {
    benches.push_back({"construct_shared_" + tag,
      "construct and destroy arrangements from a shared traits pointer",
      [&opt](Stopwatch& sw) -> std::uint64_t {
        auto tr = std::make_shared<const Traits>();
        std::uint64_t h = 0;
        sw.start();
        for (std::size_t i = 0; i < opt.small; ++i) {
          A arr(tr);
          h += arr.number_of_faces();
        }
        sw.stop();
        return h;
      }});
  }

  benches.push_back({"copy_empty_owning_" + tag,
    "copy-construct empty arrangements that own their traits",
    [&opt](Stopwatch& sw) -> std::uint64_t {
      A src;
      std::uint64_t h = 0;
      sw.start();
      for (std::size_t i = 0; i < opt.small; ++i) {
        A arr(src);
        h += arr.number_of_faces();
      }
      sw.stop();
      return h;
    }});
}

static std::vector<Benchmark> make_benchmarks(Input& in, const Options& opt) {
  std::vector<Benchmark> b;

  // --- Construction and copying -------------------------------------------
  add_construction_benches<Aos>(b, "aos", in, opt);
  add_construction_benches<Arr>(b, "arr", in, opt);
  add_construction_benches<Harr>(b, "history", in, opt);

  b.push_back({"copy_small_owning",
    "copy-construct a small arrangement that owns its traits",
    [&](Stopwatch& sw) -> std::uint64_t {
      Arr src;
      CGAL::insert(src, in.segments.begin(), in.segments.begin() + 10);
      std::uint64_t h = 0;
      const std::size_t n = std::max<std::size_t>(1, opt.small / 10);
      sw.start();
      for (std::size_t i = 0; i < n; ++i) {
        Arr arr(src);
        h += arr.number_of_edges();
      }
      sw.stop();
      return h;
    }});

  b.push_back({"copy_large",
    "copy-construct the full arrangement",
    [&](Stopwatch& sw) -> std::uint64_t {
      sw.start();
      Arr arr(in.arr);
      sw.stop();
      return size_checksum(arr);
    }});

  b.push_back({"assign_large",
    "assign the full arrangement to a non-empty arrangement",
    [&](Stopwatch& sw) -> std::uint64_t {
      Arr arr(in.small_arr);
      sw.start();
      arr = in.arr;
      sw.stop();
      return size_checksum(arr);
    }});

  // --- Insertion ----------------------------------------------------------
  b.push_back({"insert_aggregate",
    "aggregate insertion of all segments (sweep line)",
    [&](Stopwatch& sw) -> std::uint64_t {
      Arr arr(&in.traits);
      sw.start();
      CGAL::insert(arr, in.segments.begin(), in.segments.end());
      sw.stop();
      return size_checksum(arr);
    }});

  b.push_back({"insert_aggregate_owning",
    "aggregate insertion into an arrangement that owns its traits",
    [&](Stopwatch& sw) -> std::uint64_t {
      Arr arr;
      sw.start();
      CGAL::insert(arr, in.segments.begin(), in.segments.end());
      sw.stop();
      return size_checksum(arr);
    }});

  b.push_back({"insert_incremental_walk",
    "incremental insertion of all segments (walk point location)",
    [&](Stopwatch& sw) -> std::uint64_t {
      Arr arr(&in.traits);
      CGAL::Arr_walk_along_line_point_location<Arr> pl(arr);
      sw.start();
      for (const auto& s : in.segments) CGAL::insert(arr, s, pl);
      sw.stop();
      return size_checksum(arr);
    }});

  b.push_back({"insert_history_aggregate",
    "aggregate insertion into an arrangement with history",
    [&](Stopwatch& sw) -> std::uint64_t {
      Harr arr(&in.traits);
      sw.start();
      CGAL::insert(arr, in.segments.begin(), in.segments.end());
      sw.stop();
      return mix(size_checksum(arr), arr.number_of_curves());
    }});

  // --- Removal ------------------------------------------------------------
  b.push_back({"remove_edges",
    "remove every other edge of the full arrangement",
    [&](Stopwatch& sw) -> std::uint64_t {
      Arr arr(in.arr);
      std::vector<Arr::Halfedge_handle> edges;
      std::size_t i = 0;
      for (auto eit = arr.edges_begin(); eit != arr.edges_end(); ++eit, ++i)
        if (i % 2 == 0) edges.push_back(eit);
      sw.start();
      for (auto e : edges) arr.remove_edge(e);
      sw.stop();
      return size_checksum(arr);
    }});

  b.push_back({"remove_curves_history",
    "remove every other curve of the arrangement with history",
    [&](Stopwatch& sw) -> std::uint64_t {
      Harr arr(in.harr);
      std::vector<Harr::Curve_handle> curves;
      std::size_t i = 0;
      for (auto cit = arr.curves_begin(); cit != arr.curves_end(); ++cit, ++i)
        if (i % 2 == 0) curves.push_back(cit);
      std::uint64_t removed = 0;
      sw.start();
      for (auto c : curves) removed += CGAL::remove_curve(arr, c);
      sw.stop();
      return mix(size_checksum(arr), removed);
    }});

  b.push_back({"split_merge_history",
    "split every edge at its midpoint and merge it back (with history)",
    [&](Stopwatch& sw) -> std::uint64_t {
      Harr arr(in.harr);
      std::vector<Harr::Halfedge_handle> edges;
      for (auto eit = arr.edges_begin(); eit != arr.edges_end(); ++eit)
        edges.push_back(eit);
      std::vector<Point> mids;
      mids.reserve(edges.size());
      for (auto e : edges)
        mids.push_back(CGAL::midpoint(e->source()->point(),
                                      e->target()->point()));
      std::uint64_t merged = 0;
      sw.start();
      for (std::size_t k = 0; k < edges.size(); ++k) {
        auto e1 = arr.split_edge(edges[k], mids[k]);
        auto e2 = e1->next();
        if (arr.are_mergeable(e1, e2)) {
          arr.merge_edge(e1, e2);
          ++merged;
        }
      }
      sw.stop();
      return mix(size_checksum(arr), merged);
    }});

  // --- Point location -----------------------------------------------------
  auto pl_checksum = [](const auto& obj, std::uint64_t h) {
    if (std::get_if<Arr::Face_const_handle>(&obj)) return mix(h, 1);
    if (std::get_if<Arr::Halfedge_const_handle>(&obj)) return mix(h, 2);
    return mix(h, 3);
  };

  auto add_pl = [&](const std::string& name, auto make_pl, bool time_build,
                    std::size_t divisor = 1) {
    b.push_back({"pl_" + name + "_query",
      "point-location queries (" + name +
        (divisor > 1 ? ", 1/" + std::to_string(divisor) + " of the queries)"
                     : ")"),
      [&in, make_pl, pl_checksum, divisor](Stopwatch& sw) -> std::uint64_t {
        auto pl = make_pl(in.arr);
        const std::size_t nq = std::max<std::size_t>(1,
                                                     in.queries.size() / divisor);
        std::uint64_t h = 0;
        sw.start();
        for (std::size_t k = 0; k < nq; ++k)
          h = pl_checksum(pl->locate(in.queries[k]), h);
        sw.stop();
        return h;
      }});
    if (time_build) {
      b.push_back({"pl_" + name + "_build",
        "construction of the point-location structure (" + name + ")",
        [&in, make_pl](Stopwatch& sw) -> std::uint64_t {
          sw.start();
          auto pl = make_pl(in.arr);
          sw.stop();
          g_sink = g_sink + reinterpret_cast<std::uintptr_t>(pl.get()) % 2;
          return 0;
        }});
    }
  };

  add_pl("naive", [](const Arr& a) {
    return std::make_unique<CGAL::Arr_naive_point_location<Arr>>(a); },
    false, 10);
  add_pl("walk", [](const Arr& a) {
    return std::make_unique<CGAL::Arr_walk_along_line_point_location<Arr>>(a); },
    false);
  add_pl("landmarks", [](const Arr& a) {
    return std::make_unique<CGAL::Arr_landmarks_point_location<Arr>>(a); }, true);
  add_pl("trapezoid", [](const Arr& a) {
    return std::make_unique<CGAL::Arr_trapezoid_ric_point_location<Arr>>(a); },
    true);

  // --- Zone ---------------------------------------------------------------
  b.push_back({"zone_do_intersect",
    "test query segments for intersection with the arrangement (zone)",
    [&](Stopwatch& sw) -> std::uint64_t {
      CGAL::Arr_walk_along_line_point_location<Arr> pl(in.arr);
      std::uint64_t h = 0;
      sw.start();
      for (const auto& s : in.query_segments)
        h = mix(h, CGAL::do_intersect(in.arr, s, pl) ? 1 : 0);
      sw.stop();
      return h;
    }});

  // --- Overlay ------------------------------------------------------------
  b.push_back({"overlay",
    "overlay the arrangements of the two halves of the input",
    [&](Stopwatch& sw) -> std::uint64_t {
      const auto mid = in.segments.begin() + in.segments.size() / 2;
      Arr a1(&in.traits), a2(&in.traits), res(&in.traits);
      CGAL::insert(a1, in.segments.begin(), mid);
      CGAL::insert(a2, mid, in.segments.end());
      sw.start();
      CGAL::overlay(a1, a2, res);
      sw.stop();
      return size_checksum(res);
    }});

  // --- Traits access ------------------------------------------------------
  // Mimics the internal pattern of fetching the traits from the arrangement
  // and invoking a predicate, so that any overhead of the accessor shows up.
  b.push_back({"traits_access",
    "fetch the traits from the arrangement and compare vertices",
    [&](Stopwatch& sw) -> std::uint64_t {
      std::vector<Arr::Vertex_const_handle> vs;
      for (auto vit = in.arr.vertices_begin(); vit != in.arr.vertices_end();
           ++vit)
        vs.push_back(vit);
      std::uint64_t h = 0;
      sw.start();
      for (int round = 0; round < 200; ++round)
        for (std::size_t k = 1; k < vs.size(); ++k) {
          const auto* tr = in.arr.traits_adaptor();
          h += static_cast<int>(tr->compare_xy_2_object()(vs[k - 1]->point(),
                                                         vs[k]->point())) + 1;
        }
      sw.stop();
      return h;
    }});

  return b;
}

//-----------------------------------------------------------------------------
// Main
//-----------------------------------------------------------------------------
int main(int argc, char* argv[]) {
  Options opt;
  if (!parse(argc, argv, opt)) return EXIT_FAILURE;

  // Generate the input. Points and segments come from independent streams
  // derived from the seed, so changing one size does not perturb the other.
  CGAL::Random seg_rnd(opt.seed);
  CGAL::Random qry_rnd(opt.seed + 1);
  Input in;
  in.segments = random_segments(opt.segments, opt.length, seg_rnd);
  in.query_segments = random_segments(opt.queries, opt.length, qry_rnd);
  in.queries = random_points(opt.queries, qry_rnd);
  CGAL::insert(in.arr, in.segments.begin(), in.segments.end());
  CGAL::insert(in.harr, in.segments.begin(), in.segments.end());
  CGAL::insert(in.small_arr, in.segments.begin(), in.segments.begin() + 3);

  std::vector<Benchmark> benches = make_benchmarks(in, opt);

  if (opt.list) {
    for (const auto& bm : benches)
      std::cout << std::left << std::setw(30) << bm.name << bm.description
                << "\n";
    return EXIT_SUCCESS;
  }

  std::cout << "CGAL " << CGAL_VERSION_STR
            << "  shared-pointer traits API: "
            << (has_shared_ctor<Aos> ? "yes" : "no") << "\n"
            << "seed " << opt.seed << ", repetitions " << opt.repetitions
            << ", warm-up " << opt.warmup << ", segments " << opt.segments
            << ", length " << opt.length << ", queries " << opt.queries
            << ", small " << opt.small << "\n"
            << "arrangement: " << in.arr.number_of_vertices() << " vertices, "
            << in.arr.number_of_edges() << " edges, "
            << in.arr.number_of_faces() << " faces\n\n";

  std::cout << std::left << std::setw(30) << "benchmark" << std::right
            << std::setw(12) << "min ms" << std::setw(12) << "median ms"
            << std::setw(12) << "mean ms" << std::setw(10) << "cv %"
            << "  checksum\n";

  std::ofstream csv;
  if (!opt.csv.empty()) {
    csv.open(opt.csv);
    if (!csv) {
      std::cerr << "Cannot open " << opt.csv << "\n";
      return EXIT_FAILURE;
    }
    csv << "# cgal=" << CGAL_VERSION_STR
        << " shared_api=" << (has_shared_ctor<Aos> ? 1 : 0)
        << " seed=" << opt.seed << " repetitions=" << opt.repetitions
        << " segments=" << opt.segments << " length=" << opt.length
        << " queries=" << opt.queries << " small=" << opt.small << "\n"
        << "name,repetitions,min_ms,median_ms,mean_ms,stddev_ms,checksum,"
           "times_ms\n";
  }

  bool all_consistent = true;
  for (const auto& bm : benches) {
    if (!opt.filter.empty() && bm.name.find(opt.filter) == std::string::npos)
      continue;
    const Result r = run(bm, opt);
    all_consistent = all_consistent && r.consistent;
    const double mn = *std::min_element(r.times.begin(), r.times.end());
    const double md = median(r.times);
    const double me = mean(r.times);
    const double sd = stddev(r.times);
    std::cout << std::left << std::setw(30) << r.name << std::right
              << std::fixed << std::setprecision(3)
              << std::setw(12) << mn << std::setw(12) << md
              << std::setw(12) << me << std::setprecision(1)
              << std::setw(10) << (me > 0 ? 100.0 * sd / me : 0.0)
              << "  " << std::hex << r.checksum << std::dec
              << (r.consistent ? "" : "  (INCONSISTENT)") << "\n";
    if (csv) {
      csv << r.name << ',' << r.times.size() << ',' << std::setprecision(6)
          << mn << ',' << md << ',' << me << ',' << sd << ','
          << std::hex << r.checksum << std::dec << ',';
      for (std::size_t k = 0; k < r.times.size(); ++k)
        csv << (k ? ";" : "") << r.times[k];
      csv << "\n";
    }
  }

  if (!all_consistent) {
    std::cerr << "\nWarning: some benchmarks produced different checksums "
                 "across repetitions.\n";
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
