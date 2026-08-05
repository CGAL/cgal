//! \file bench_lower_envelope.cpp
// Benchmark: measures the running time of lower_envelope_x_monotone_2()
// over randomly generated x-monotone line segments, for a range of input
// sizes.  The result is a table suitable for plotting a scaling curve.
//
// Usage:
//   bench_lower_envelope [--seed <uint>] [--reps <uint>]
//
// Output (tab-separated, one header line then one row per input size):
//   n   reps   mean_ms   stddev_ms   mean_vertices   mean_edges

#include <CGAL/Exact_rational.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arr_curve_data_traits_2.h>
#include <CGAL/Envelope_diagram_1.h>
#include <CGAL/envelope_2.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <random>
#include <string>
#include <iomanip>
#include <vector>

// ---------------------------------------------------------------------------
// Type setup.
// Arr_curve_data_traits_2 wraps Arr_segment_traits_2 and is required for
// compatibility with Envelope_diagram_1 and Env_divide_and_conquer_2: those
// templates type their internal curve-pointer containers from the diagram's
// X_monotone_curve_2, which must match the wrapped type used as input.
// We attach a dummy char label (unused) solely to satisfy the wrapper's
// interface — this matches the pattern used in the reference example.
// ---------------------------------------------------------------------------
using NT = CGAL::Exact_rational;
using Kernel = CGAL::Cartesian<NT>;
using Segment_traits  = CGAL::Arr_segment_traits_2<Kernel>;
using Traits_2 = CGAL::Arr_curve_data_traits_2<Segment_traits, char>;
using Point_2 = Traits_2::Point_2;
using Segment_2 = Traits_2::X_monotone_curve_2;  // the wrapped curve type
using Base_segment_2 = Segment_traits::X_monotone_curve_2;
using Diagram_1 = CGAL::Envelope_diagram_1<Traits_2>;

// ---------------------------------------------------------------------------
// Generate `n` random x-monotone segments.
// Each endpoint has integer coordinates in [0, coord_range].
// Segments with equal left and right x are discarded and regenerated so that
// every segment is strictly x-monotone (no vertical segments).
// ---------------------------------------------------------------------------
static std::vector<Segment_2>
generate_segments(std::size_t n, int coord_range, std::mt19937& rng) {
  std::uniform_int_distribution<int> dist(0, coord_range);
  std::vector<Segment_2> segs;
  segs.reserve(n);
  while (segs.size() < n) {
    int x1 = dist(rng), y1 = dist(rng);
    int x2 = dist(rng), y2 = dist(rng);
    if (x1 == x2) continue;          // skip vertical / degenerate segments
    if (x1 > x2) { std::swap(x1, x2); std::swap(y1, y2); } // ensure left < right
    // Construct the wrapped X_monotone_curve_2 from a base segment and a
    // dummy label.  The label is unused; it exists only to satisfy the
    // Arr_curve_data_traits_2 interface.
    Base_segment_2 base{Point_2{NT(x1), NT(y1)}, Point_2{NT(x2), NT(y2)}};
    segs.push_back(Segment_2{base, '\0'});
  }
  return segs;
}

// ---------------------------------------------------------------------------
// Run one lower-envelope computation and return {elapsed_ms, #vertices, #edges}.
// ---------------------------------------------------------------------------
struct Run_result {
  double ms;
  std::size_t vertices;
  std::size_t edges;
};

static Run_result run_once(std::vector<Segment_2>& segs, const Traits_2& traits) {
  Diagram_1 diag;
  auto t0 = std::chrono::steady_clock::now();
  CGAL::lower_envelope_x_monotone_2(segs.begin(), segs.end(), diag, traits);
  auto t1 = std::chrono::steady_clock::now();

  double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

  // Count diagram vertices and edges by walking left-to-right.
  // There is always one more edge than there are interior vertices.
  std::size_t nv = 0, ne = 1; // start with 1 for the leftmost edge
  {
    auto ec = diag.leftmost();
    while (ec != diag.rightmost()) {
      ++nv; // the vertex to the right of ec
      ++ne; // the edge to the right of that vertex
      ec = diag.right_edge(diag.right_vertex(ec));
    }
  }

  return {ms, nv, ne};
}

// ---------------------------------------------------------------------------
// Statistics helpers
// ---------------------------------------------------------------------------
static double mean(const std::vector<double>& v) {
  double s = 0;
  for (double x : v) s += x;
  return s / v.size();
}

static double stddev(const std::vector<double>& v, double m) {
  double s = 0;
  for (double x : v) s += (x-m)*(x-m);
  return std::sqrt(s / v.size());
}

// ---------------------------------------------------------------------------
// main
// ---------------------------------------------------------------------------
int main(int argc, char* argv[]) {
  // -- Parse optional arguments ---------------------------------------------
  unsigned seed = 42;
  unsigned reps = 5;          // repetitions per input size
  for (int i = 1; i < argc; ++i) {
    if (std::strcmp(argv[i], "--seed") == 0 && i+1 < argc) seed = static_cast<unsigned>(std::stoul(argv[++i]));
    else if (std::strcmp(argv[i], "--reps") == 0 && i+1 < argc) reps = static_cast<unsigned>(std::stoul(argv[++i]));
    else {
      std::cerr << "Usage: " << argv[0] << " [--seed <uint>] [--reps <uint>]\n";
      return 1;
    }
  }

  // -- Benchmark parameters -------------------------------------------------
  // Input sizes to test.  Geometric progression: 100, 200, 400, ..., 25600.
  const std::vector<std::size_t> sizes = {
    100, 200, 400, 800, 1600, 3200, 6400, 12800, 25600
  };

  // Coordinate range: sqrt(n) * 4 gives moderate intersection density without
  // degenerating to a fully overlapping clump.  We reuse a fixed range
  // (large enough for the biggest n) so results across sizes are comparable.
  const int coord_range = 10000;

  // -- Print header ---------------------------------------------------------
  std::cout << "# Lower envelope benchmark — Arr_segment_traits_2\n"
            << "# seed=" << seed << "  reps=" << reps << "\n"
            << "#\n"
            << std::left
            << std::setw(8) << "n"
            << std::setw(8) << "reps"
            << std::setw(14) << "mean_ms"
            << std::setw(14) << "stddev_ms"
            << std::setw(14) << "mean_verts"
            << std::setw(14) << "mean_edges"
            << "\n";

  // Shared traits object: reusing it improves speed if the traits caches data.
  Traits_2 traits;

  // -- Run benchmarks -------------------------------------------------------
  for (std::size_t n : sizes) {
    std::vector<double> times;
    times.reserve(reps);
    double sum_verts = 0, sum_edges = 0;

    for (unsigned r = 0; r < reps; ++r) {
      // Each repetition gets a distinct seed derived from the base seed, n,
      // and the repetition index, so the datasets are independent.
      std::mt19937 rng(seed ^ static_cast<unsigned>(n) ^ (r * 2654435761u));
      auto segs = generate_segments(n, coord_range, rng);

      auto res = run_once(segs, traits);
      times.push_back(res.ms);
      sum_verts += static_cast<double>(res.vertices);
      sum_edges += static_cast<double>(res.edges);
    }

    double m   = mean(times);
    double sd  = stddev(times, m);
    double mv  = sum_verts / reps;
    double me  = sum_edges / reps;

    std::cout << std::left
              << std::setw(8) << n
              << std::setw(8) << reps
              << std::setw(14) << std::fixed << std::setprecision(3) << m
              << std::setw(14) << sd
              << std::setw(14) << std::setprecision(1) << mv
              << std::setw(14) << me
              << "\n";
    std::cout.flush();
  }

  return 0;
}
