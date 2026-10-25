/*! \file arr_on_circle.cpp
 * Benchmarks arrangements on a circle, that is, `Arrangement_on_curve_1` with
 * `Circle_traits_1` and `Identified_topology_traits`. It times:
 * - the comparison of two directions (`Compare_x_1`), separately for distinct
 *   and for equal directions,
 * - sorting directions with `Compare_x_1`,
 * - incremental insertion of directions with the free function `insert()`,
 * - the overlay of two arrangements.
 * The arrangements store their vertices and edges either in lists (the
 * default) or in vectors; in the latter case, the arrangements built from
 * sorted directions reserve their storage in advance.
 * The directions are either given directly by their exact coordinates, or
 * constructed as the differences of pairs of points (as, e.g., the directions
 * of the edges of a polygon), in which case the kernel (Epeck) evaluates them
 * lazily. A fraction of the directions of the second arrangement of the
 * overlay duplicate (represented differently) directions of the first.
 */

#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include <boost/program_options.hpp>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Random.h>
#include <CGAL/Arrangement_on_curve_1/Circle_traits_1.h>
#include <CGAL/Arrangement_on_curve_1/Identified_topology_traits.h>
#include <CGAL/Arrangement_on_curve_1/Arrangement_on_curve_1.h>
#include <CGAL/Arrangement_on_curve_1/insert.h>
#include <CGAL/Arrangement_on_curve_1/overlay.h>

namespace po = boost::program_options;
namespace Aoc = CGAL::Arrangement_on_curve_1;

using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;
using Point_2 = Kernel::Point_2;
using Vector_2 = Kernel::Vector_2;
using Direction_2 = Kernel::Direction_2;
using Traits_1 = Aoc::Circle_traits_1<Kernel>;
template <bool UseVector>
using Arrangement =
  Aoc::Arrangement_on_curve_1<Traits_1, Aoc::Identified_topology_traits<Traits_1, void, void, UseVector>>;
using Clock = std::chrono::steady_clock;

// Generates directions that differ from the reference direction.
class Direction_generator {
public:
  Direction_generator(bool constructed, const Traits_1& traits, CGAL::Random& rnd) :
    m_constructed(constructed), m_traits(traits), m_rnd(rnd) {}

  // generates a random direction.
  Direction_2 operator()() {
    while (true) {
      Direction_2 d = m_constructed ? Direction_2(random_point() - random_point()) :
        Direction_2(m_rnd.get_double(-1, 1), m_rnd.get_double(-1, 1));
      if (! m_traits.is_on_x_identification_1_object()(d)) return d;
    }
  }

  // generates a direction equal to `d`, represented differently.
  Direction_2 duplicate(const Direction_2& d) {
    const double s = m_rnd.get_double(0.5, 4);
    if (m_constructed) return Direction_2(Vector_2(d.dx(), d.dy()) * s);
    return Direction_2(d.dx() * s, d.dy() * s);
  }

private:
  Point_2 random_point() { return Point_2(m_rnd.get_double(-10, 10), m_rnd.get_double(-10, 10)); }

  bool m_constructed;
  const Traits_1& m_traits;
  CGAL::Random& m_rnd;
};

// Times `f` `iterations` times (after a warm-up run); returns the median in ms.
template <typename F>
double median_ms(int iterations, F f) {
  f();
  std::vector<double> ts;
  for (int i = 0; i < iterations; ++i) {
    auto t0 = Clock::now();
    f();
    ts.push_back(std::chrono::duration<double, std::milli>(Clock::now() - t0).count());
  }
  std::sort(ts.begin(), ts.end());
  return ts[ts.size() / 2];
}

// Builds an arrangement from directions sorted by `Compare_x_1`, without
// duplicates, in linear time.
template <typename Arrangement>
void build_sorted(Arrangement& arr, const std::vector<Direction_2>& sorted) {
  if (sorted.empty()) return;
  arr.reserve(sorted.size());
  auto v = arr.insert_empty(sorted.front());
  for (auto it = std::next(sorted.begin()); it != sorted.end(); ++it) v = arr.insert_after(v, *it);
}

void print_row(const std::string& name, double ms, std::size_t num_compares) {
  std::cout << "  " << std::left << std::setw(36) << name << std::right << std::fixed
            << std::setw(12) << std::setprecision(3) << ms;
  if (num_compares > 0)
    std::cout << std::setw(12) << num_compares << std::setw(14) << std::setprecision(1)
              << 1e6 * ms / double(num_compares);
  std::cout << "\n";
}

// Times the operations on arrangements in list mode or in vector mode.
template <bool UseVector>
void run(std::shared_ptr<const Traits_1> traits, const std::vector<Direction_2>& dirs1,
         const std::vector<Direction_2>& dirs2, Direction_generator& gen, int iterations) {
  using Arr = Arrangement<UseVector>;
  const auto cmp = traits->compare_x_1_object();
  const std::size_t n = dirs1.size();

  // 1. Compare_x_1 on distinct and on equal directions.
  std::vector<std::pair<Direction_2, Direction_2>> distinct, equal;
  for (std::size_t i = 0; i < n; ++i) {
    distinct.emplace_back(dirs1[i], dirs2[(i * 7 + 1) % n]);
    equal.emplace_back(dirs1[i], gen.duplicate(dirs1[i]));
  }
  volatile int sink = 0;
  const int rounds = 100;
  auto compare_all = [&](const auto& pairs) {
    for (int r = 0; r < rounds; ++r)
      for (const auto& [d1, d2] : pairs) sink = sink + int(cmp(d1, d2));
  };
  print_row("Compare_x_1 (distinct directions)", median_ms(iterations, [&] { compare_all(distinct); }),
            rounds * distinct.size());
  print_row("Compare_x_1 (equal directions)", median_ms(iterations, [&] { compare_all(equal); }),
            rounds * equal.size());

  // 2. Sorting.
  auto less = [&](const Direction_2& a, const Direction_2& b) { return cmp(a, b) == CGAL::SMALLER; };
  std::size_t num_sort_compares = 0;
  {
    auto v = dirs1;
    std::sort(v.begin(), v.end(), [&](const Direction_2& a, const Direction_2& b) {
      ++num_sort_compares;
      return less(a, b);
    });
  }
  print_row("sort", median_ms(iterations, [&] { auto v = dirs1; std::sort(v.begin(), v.end(), less); }),
            num_sort_compares);

  // 3. Incremental insertion (each insertion locates the point by a linear
  // walk, which makes about n^2/4 comparisons in total).
  print_row("incremental insert (dirs1)", median_ms(iterations, [&] {
    Arr arr(traits);
    for (const auto& d : dirs1) Aoc::insert(arr, d);
  }), 0);

  // 4. Overlay of two arrangements built from the sorted directions.
  auto sorted_unique = [&](std::vector<Direction_2> v) {
    std::sort(v.begin(), v.end(), less);
    v.erase(std::unique(v.begin(), v.end(), [&](const Direction_2& a, const Direction_2& b)
                        { return cmp(a, b) == CGAL::EQUAL; }), v.end());
    return v;
  };
  const auto sorted1 = sorted_unique(dirs1);
  const auto sorted2 = sorted_unique(dirs2);
  Arr arr1(traits), arr2(traits);
  build_sorted(arr1, sorted1);
  build_sorted(arr2, sorted2);
  std::size_t num_res_vertices = 0;
  print_row("overlay", median_ms(iterations, [&] {
    Arr res(traits);
    Aoc::overlay(arr1, arr2, res);
    num_res_vertices = res.number_of_vertices();
  }), 0);
  print_row("build from sorted (both)", median_ms(iterations, [&] {
    Arr a1(traits), a2(traits);
    build_sorted(a1, sorted1);
    build_sorted(a2, sorted2);
  }), 0);

  std::cout << "\n  overlay: " << arr1.number_of_vertices() << " + " << arr2.number_of_vertices() << " vertices -> "
            << num_res_vertices << " vertices\n";
  (void)sink;
}

int main(int argc, char* argv[]) {
  po::options_description desc("Usage: arr_on_circle [options]\n"
                               "Benchmarks arrangements on a circle (Circle_traits_1).\nOptions");
  desc.add_options()
    ("help,h", "print this help message")
    ("seed,s", po::value<unsigned int>()->default_value(1337), "the random seed")
    ("iterations,i", po::value<int>()->default_value(5), "the number of timed runs (after a warm-up run)")
    ("directions,n", po::value<std::size_t>()->default_value(1000), "the number of directions per arrangement")
    ("duplicates", po::value<double>()->default_value(0.1),
     "the fraction of the directions of the second arrangement that duplicate directions of the first")
    ("reference", po::value<std::vector<double>>()->multitoken(),
     "the reference direction DX DY, that is, the identification point (default: 1 0)")
    ("exact", "use directions given by their exact coordinates; by default, the directions are constructed as "
     "differences of pairs of points")
    ("vector", "store the vertices and edges of the arrangements in vectors; by default, they are stored in lists");

  po::variables_map vm;
  try {
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);
    if (vm.count("help")) {
      std::cout << desc << "\n";
      return EXIT_SUCCESS;
    }
    if (vm["iterations"].as<int>() < 1) throw po::error("--iterations must be positive");
    if (vm["directions"].as<std::size_t>() < 2) throw po::error("--directions must be at least 2");
    const double f = vm["duplicates"].as<double>();
    if ((f < 0) || (f > 1)) throw po::error("--duplicates must be in [0, 1]");
    if (vm.count("reference") && (vm["reference"].as<std::vector<double>>().size() != 2))
      throw po::error("--reference expects two values");
  }
  catch (const std::exception& e) {
    std::cerr << "Error: " << e.what() << "\nRun with --help for the list of options.\n";
    return EXIT_FAILURE;
  }

  const unsigned int seed = vm["seed"].as<unsigned int>();
  const int iterations = vm["iterations"].as<int>();
  const std::size_t n = vm["directions"].as<std::size_t>();
  const double dup_fraction = vm["duplicates"].as<double>();
  const bool constructed = (vm.count("exact") == 0);
  Direction_2 ref(1, 0);
  if (vm.count("reference")) {
    const auto& r = vm["reference"].as<std::vector<double>>();
    ref = Direction_2(r[0], r[1]);
  }
  auto traits = std::make_shared<const Traits_1>(ref);

  // Generate the directions.
  CGAL::Random rnd(seed);
  Direction_generator gen(constructed, *traits, rnd);
  std::vector<Direction_2> dirs1, dirs2;
  for (std::size_t i = 0; i < n; ++i) dirs1.push_back(gen());
  for (std::size_t i = 0; i < n; ++i)
    dirs2.push_back((rnd.get_double() < dup_fraction) ? gen.duplicate(dirs1[rnd.get_int(0, int(n))]) : gen());

  std::cout << "Arrangement on a circle benchmark\n"
            << "  seed " << seed << ", " << iterations << " iterations, " << n << " directions per arrangement, "
            << (constructed ? "constructed" : "exact") << " directions, " << dup_fraction << " duplicates, "
            << "reference (" << ref.dx() << ", " << ref.dy() << "), "
            << (vm.count("vector") ? "vector" : "list") << " storage\n\n";
  std::cout << "  " << std::left << std::setw(36) << "operation" << std::right << std::setw(12) << "ms"
            << std::setw(12) << "compares" << std::setw(14) << "ns/compare" << "\n";

  if (vm.count("vector")) run<true>(traits, dirs1, dirs2, gen, iterations);
  else run<false>(traits, dirs1, dirs2, gen, iterations);
  return EXIT_SUCCESS;
}
