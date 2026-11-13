// Tests that Default_topology_traits selects the topology traits that matches
// the side categories of the geometry traits, that arrangements with the
// default topology traits work, that binary-search mode works, and that
// reserve() works.

#include <cstdlib>
#include <type_traits>

#include "test_utils.h"
#include <CGAL/Arrangement_on_curve_1/overlay.h>

// reserves storage for the given points, inserts them, and checks the
// result; in vector mode, also checks that the storage was not reallocated.
// Finally, reserves storage for twice as many vertices and checks that the
// arrangement is unchanged.
template <typename Arr>
bool test_reserve(Arr& arr, const std::vector<typename Arr::Point_1>& pts, std::string& why) {
  arr.reserve(pts.size());
  const auto& tt = arr.topology_traits();
  auto storage = [&]() {
    if constexpr (Arr::Topology_traits::use_vector)
      return std::make_pair(tt.raw_vertices().data(), tt.raw_edges().data());
    else return 0;
  };
  const auto before = storage();
  for (const auto& p : pts) Aoc::insert(arr, p);
  if (! is_valid(arr, pts, why)) return false;
  if (storage() != before) {
    why = "the storage was reallocated";
    return false;
  }
  arr.reserve(2 * pts.size());
  return is_valid(arr, pts, why);
}

template <bool UseVector>
void test_reserve_all(const std::string& mode, std::size_t& failures) {
  std::string why;
  {
    Test_arrangement<Closed_traits, UseVector> arr(std::make_shared<const Closed_traits>(FT(0), FT(10)));
    bool ok = test_reserve(arr, {FT(0), FT(3), FT(5), FT(7), FT(10)}, why);
    std::cout << (report(ok, mode + " reserve closed", why, failures) ? "PASS " : "FAIL ") << mode
              << " reserve closed" << std::endl;
  }
  {
    Test_arrangement<Identified_traits, UseVector> arr;
    bool ok = test_reserve(arr, {Direction_2(0, 1), Direction_2(-1, 1), Direction_2(-1, 0), Direction_2(0, -1),
                                 Direction_2(1, 0)}, why);
    std::cout << (report(ok, mode + " reserve identified", why, failures) ? "PASS " : "FAIL ") << mode
              << " reserve identified" << std::endl;
  }
}

template <typename Gt, typename Tt>
constexpr bool is_default_v = std::is_same_v<typename Aoc::Arrangement_on_curve_1<Gt>::Topology_traits, Tt>;

static_assert(is_default_v<Open_traits, Aoc::Open_topology_traits<Open_traits>>);
static_assert(is_default_v<Closed_traits, Aoc::Closed_topology_traits<Closed_traits>>);
static_assert(is_default_v<Left_closed_traits, Aoc::Left_closed_topology_traits<Left_closed_traits>>);
static_assert(is_default_v<Right_closed_traits, Aoc::Right_closed_topology_traits<Right_closed_traits>>);
static_assert(is_default_v<Identified_traits, Aoc::Identified_topology_traits<Identified_traits>>);
static_assert(std::is_same_v<Aoc::Default_topology_traits<Identified_traits, int, double, true>::type,
                             Aoc::Identified_topology_traits<Identified_traits, int, double, true>>);
static_assert(std::is_same_v<Aoc::Default_topology_traits<Closed_traits, void, int>::type,
                             Aoc::Closed_topology_traits<Closed_traits, void, int>>);

int main() {
  std::size_t failures = 0;
  std::string why;

  // Default topology traits of closed geometry traits, including overlay.
  {
    using Arr = Aoc::Arrangement_on_curve_1<Closed_traits>;
    auto tr = std::make_shared<const Closed_traits>(FT(0), FT(10));
    Arr arr_a(tr), arr_b(tr), arr_res(tr);
    for (int x : {0, 5}) Aoc::insert(arr_a, FT(x));
    for (int x : {10, 3}) Aoc::insert(arr_b, FT(x));
    Aoc::overlay(arr_a, arr_b, arr_res);
    bool ok = is_valid(arr_res, {FT(0), FT(3), FT(5), FT(10)}, why);
    std::cout << (report(ok, "default closed", why, failures) ? "PASS " : "FAIL ") << "default closed" << std::endl;
  }

  // Default topology traits of identified geometry traits.
  {
    using Arr = Aoc::Arrangement_on_curve_1<Identified_traits>;
    Arr arr;
    for (auto d : {Direction_2(0, 1), Direction_2(1, 0), Direction_2(-1, 0)}) Aoc::insert(arr, d);
    bool ok = is_valid(arr, {Direction_2(0, 1), Direction_2(1, 0), Direction_2(-1, 0)}, why);
    std::cout << (report(ok, "default identified", why, failures) ? "PASS " : "FAIL ") << "default identified"
              << std::endl;
  }

  // Binary-search mode; the vertices are inserted in increasing order.
  {
    using Arr =
      Aoc::Arrangement_on_curve_1<Closed_traits, Aoc::Closed_topology_traits<Closed_traits, void, void, true>, true>;
    Arr arr(std::make_shared<const Closed_traits>(FT(0), FT(10)));
    for (int x : {0, 4, 7, 10}) Aoc::insert(arr, FT(x));
    bool ok = is_valid(arr, {FT(0), FT(4), FT(7), FT(10)}, why);
    if (ok && (Aoc::locate(arr, FT(5)).index() != 1)) {
      ok = false;
      why = "an absent point is located at a vertex";
    }
    std::cout << (report(ok, "binary search", why, failures) ? "PASS " : "FAIL ") << "binary search" << std::endl;
  }

  test_reserve_all<false>("list", failures);
  test_reserve_all<true>("vector", failures);

  return (failures == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
