// Benchmark of hot-pixel snap rounding (CGAL::snap_rounding_2() with
// Hot_pixel_snap_rounding_traits_2).
//
// The input is either read from a file or generated:
//  - strip:  a strip of thin quads (rails, rungs, and one diagonal per quad)
//            rotated by a given angle; no two segments cross.
//  - random: random segments of a fixed length and random directions.
//
// Run with --help for the list of options.

#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <list>
#include <random>
#include <string>
#include <vector>

#include <boost/program_options.hpp>

#include <CGAL/Installation/internal/disable_deprecation_warnings_and_errors.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Hot_pixel_snap_rounding_traits_2.h>
#include <CGAL/Snap_rounding_2.h>
#include <CGAL/version.h>

namespace po = boost::program_options;

using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;
using Traits = CGAL::Hot_pixel_snap_rounding_traits_2<Kernel>;
using FT = Kernel::FT;
using Point_2 = Kernel::Point_2;
using Segment_2 = Kernel::Segment_2;
using Polyline_2 = std::list<Point_2>;
using Polylines = std::list<Polyline_2>;

struct Strip_params {
  double angle;
  int rails;
  int points;
  double step;
  double spacing;
  double jitter;
  double offset;
};

// Parallel rails with staggered vertices, connected by rungs and quad diagonals.
std::vector<Segment_2> generate_strip(const Strip_params& sp, unsigned seed) {
  const double a = sp.angle * CGAL_PI / 180.0;
  const double ux = std::cos(a), uy = std::sin(a);
  const double nx = -uy, ny = ux;
  std::mt19937 gen(seed);
  std::uniform_real_distribution<double> jitter(0.0, sp.jitter);

  std::vector<std::vector<Point_2>> v(sp.rails, std::vector<Point_2>(sp.points));
  for (int r = 0; r < sp.rails; ++r) {
    const double phase = std::fmod(r * 0.618034 * sp.step, sp.step);
    for (int p = 0; p < sp.points; ++p) {
      const double along = p * sp.step + phase;
      const double across = r * sp.spacing;
      v[r][p] = Point_2(sp.offset + along * ux + across * nx + jitter(gen),
                        sp.offset + along * uy + across * ny + jitter(gen));
    }
  }

  std::vector<Segment_2> segs;
  for (int r = 0; r < sp.rails; ++r) {
    for (int p = 0; p < sp.points; ++p) {
      if (p + 1 < sp.points) segs.emplace_back(v[r][p], v[r][p + 1]);
      if (r + 1 < sp.rails) segs.emplace_back(v[r][p], v[r + 1][p]);
      if ((r + 1 < sp.rails) && (p + 1 < sp.points)) segs.emplace_back(v[r][p], v[r + 1][p + 1]);
    }
  }
  return segs;
}

std::vector<Segment_2> generate_random(std::size_t n, double length, double extent, double offset, unsigned seed) {
  std::mt19937 gen(seed);
  std::uniform_real_distribution<double> pos(0.0, extent);
  std::uniform_real_distribution<double> dir(0.0, 2.0 * CGAL_PI);
  std::vector<Segment_2> segs;
  segs.reserve(n);
  for (std::size_t i = 0; i < n; ++i) {
    const double x = offset + pos(gen), y = offset + pos(gen), a = dir(gen);
    segs.emplace_back(Point_2(x, y), Point_2(x + length * std::cos(a), y + length * std::sin(a)));
  }
  return segs;
}

// File format: the number of segments followed by "x1 y1 x2 y2" per segment.
bool read_segments(const std::string& filename, std::vector<Segment_2>& segs) {
  std::ifstream is(filename);
  if (! is.is_open()) return false;
  std::size_t n;
  if (! (is >> n)) return false;
  segs.reserve(n);
  for (std::size_t i = 0; i < n; ++i) {
    FT x1, y1, x2, y2;
    if (! (is >> x1 >> y1 >> x2 >> y2)) return false;
    segs.emplace_back(Point_2(x1, y1), Point_2(x2, y2));
  }
  return true;
}

// An order-dependent hash of the output used to compare runs.
std::uint64_t hash_output(const Polylines& polylines) {
  std::uint64_t h = 1469598103934665603ull;
  auto mix = [&h](double d) {
    std::uint64_t bits;
    std::memcpy(&bits, &d, sizeof(bits));
    h = (h ^ bits) * 1099511628211ull;
  };
  for (const auto& pl : polylines) {
    mix(static_cast<double>(pl.size()));
    for (const auto& p : pl) {
      mix(CGAL::to_double(p.x()));
      mix(CGAL::to_double(p.y()));
    }
  }
  return h;
}

int main(int argc, char* argv[]) {
  po::options_description desc("Options");
  std::string input_file, generator;
  Strip_params sp;
  std::size_t random_n;
  double random_length, random_extent, pixel_size;
  unsigned int kd_trees, repetitions, seed;
  bool do_isr, grid_coords;
  std::vector<double> angles;

  desc.add_options()
    ("help,h", "print help message")
    ("input-file,i", po::value<std::string>(&input_file), "read the segments from a file")
    ("generator,g", po::value<std::string>(&generator)->default_value("strip"), "input generator: strip or random")
    ("angle,a", po::value<std::vector<double>>(&angles)->multitoken(),
     "strip angles in degrees (default: 35 0)")
    ("rails", po::value<int>(&sp.rails)->default_value(80), "strip: number of rails")
    ("points", po::value<int>(&sp.points)->default_value(160), "strip: number of vertices per rail")
    ("step", po::value<double>(&sp.step)->default_value(8000), "strip: distance between consecutive rail vertices")
    ("spacing", po::value<double>(&sp.spacing)->default_value(25), "strip: distance between rails")
    ("jitter", po::value<double>(&sp.jitter)->default_value(3), "strip: maximum vertex perturbation")
    ("offset", po::value<double>(&sp.offset)->default_value(4e8), "translation of the generated input")
    ("segments,n", po::value<std::size_t>(&random_n)->default_value(10000), "random: number of segments")
    ("length", po::value<double>(&random_length)->default_value(100), "random: segment length")
    ("extent", po::value<double>(&random_extent)->default_value(10000), "random: side of the square domain")
    ("seed", po::value<unsigned int>(&seed)->default_value(12345), "random seed")
    ("pixel-size,p", po::value<double>(&pixel_size)->default_value(1), "pixel size")
    ("kd-trees,k", po::value<unsigned int>(&kd_trees)->default_value(1), "number of kd-trees")
    ("isr", po::value<bool>(&do_isr)->default_value(false), "apply iterated snap rounding")
    ("grid-coordinates", po::value<bool>(&grid_coords)->default_value(false), "output integer grid coordinates")
    ("repetitions,r", po::value<unsigned int>(&repetitions)->default_value(1), "number of repetitions");

  po::variables_map vm;
  try {
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);
  }
  catch (const std::exception& e) {
    std::cerr << e.what() << "\n" << desc << std::endl;
    return 1;
  }
  if (vm.count("help")) {
    std::cout << desc << std::endl;
    return 0;
  }
  if (angles.empty()) angles = {35, 0};
  if (repetitions == 0) repetitions = 1;

  struct Input {
    std::string name;
    std::vector<Segment_2> segments;
  };
  std::vector<Input> inputs;
  if (vm.count("input-file")) {
    Input in{input_file, {}};
    if (! read_segments(input_file, in.segments)) {
      std::cerr << "Failed to read " << input_file << std::endl;
      return 1;
    }
    inputs.push_back(std::move(in));
  }
  else if (generator == "strip") {
    for (double angle : angles) {
      sp.angle = angle;
      inputs.push_back({"strip " + std::to_string(angle).substr(0, 5), generate_strip(sp, seed)});
    }
  }
  else if (generator == "random") {
    inputs.push_back({"random", generate_random(random_n, random_length, random_extent, sp.offset, seed)});
  }
  else {
    std::cerr << "Unknown generator " << generator << std::endl;
    return 1;
  }

  std::cout << "CGAL " << CGAL_VERSION_STR << ", pixel size " << pixel_size << ", kd-trees " << kd_trees
            << ", " << (do_isr ? "ISR" : "SR") << ", repetitions " << repetitions << std::endl;
  std::cout << std::fixed;
  for (const auto& in : inputs) {
    double total = 0, best = std::numeric_limits<double>::max();
    Polylines out;
    for (unsigned int i = 0; i < repetitions; ++i) {
      out.clear();
      auto start = std::chrono::steady_clock::now();
      CGAL::snap_rounding_2<Traits>(in.segments.begin(), in.segments.end(), out, FT(pixel_size),
                                    do_isr, grid_coords, kd_trees);
      std::chrono::duration<double> d = std::chrono::steady_clock::now() - start;
      total += d.count();
      best = (std::min)(best, d.count());
    }
    std::size_t num_points = 0;
    for (const auto& pl : out) num_points += pl.size();
    std::cout << std::left << std::setw(14) << in.name << std::right
              << " segments " << std::setw(8) << in.segments.size()
              << " | time best " << std::setprecision(3) << best << " s, mean " << total / repetitions << " s"
              << " | output points " << num_points
              << " | hash " << std::hex << hash_output(out) << std::dec << std::endl;
  }
  return 0;
}
