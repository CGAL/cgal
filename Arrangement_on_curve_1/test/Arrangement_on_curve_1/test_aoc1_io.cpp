// Tests the input and output operators for every topology, in list and vector
// mode: every arrangement constructed by inserting a permutation of a subset
// of a set of points is written, and read back into a non-empty arrangement,
// which must be replaced.  Also tests that malformed input is rejected.

#include <cstdlib>
#include <sstream>

#include "test_utils.h"
#include <CGAL/Arrangement_on_curve_1/IO/iostream.h>

template <typename Arr>
void test(const std::string& name, std::shared_ptr<const typename Arr::Geometry_traits_1> tr,
          const std::vector<typename Arr::Point_1>& pts, const typename Arr::Point_1& extra, std::size_t& failures) {
  std::size_t local_failures = 0;
  for_each_sequence(pts.size(), [&](const std::vector<int>& seq) {
    // In binary-search mode, vertices must be inserted in increasing order.
    if constexpr (Arr::binary_search_enabled) {
      if (! std::is_sorted(seq.begin(), seq.end())) return;
    }
    Arr arr(tr);
    std::vector<typename Arr::Point_1> inserted;
    for (int i : seq) {
      Aoc::insert(arr, pts[i]);
      inserted.push_back(pts[i]);
    }
    std::stringstream ss;
    ss << arr;
    const std::string text = ss.str();

    Arr read(tr);
    if constexpr (! Arr::binary_search_enabled) Aoc::insert(read, extra);
    ss >> read;
    std::stringstream ss2;
    ss2 << read;
    std::string why;
    bool ok = ! ss.fail() && is_valid(read, inserted, why);
    if (ok && (ss2.str() != text)) {
      ok = false;
      why = "the output of the read arrangement differs";
    }
    report(ok, name, why, local_failures);
  });
  failures += local_failures;
  std::cout << (local_failures ? "FAIL " : "PASS ") << name << " round trip" << std::endl;
}

template <typename Arr>
void test_rejection(const std::string& name, std::shared_ptr<const typename Arr::Geometry_traits_1> tr,
                    const std::string& text, std::size_t& failures) {
  Arr arr(tr);
  std::stringstream ss(text);
  ss >> arr;
  bool ok = ss.fail() && arr.empty() && (arr.number_of_edges() == 1);
  std::cout << (report(ok, name, "malformed input is accepted", failures) ? "PASS " : "FAIL ") << "rejects " << name
            << std::endl;
}

template <bool UseVector>
void test_all(const std::string& mode, std::size_t& failures) {
  using Closed_arr = Test_arrangement<Closed_traits, UseVector>;
  using Identified_arr = Test_arrangement<Identified_traits, UseVector>;
  auto closed_tr = std::make_shared<const Closed_traits>(FT(0), FT(10));
  auto identified_tr = std::make_shared<const Identified_traits>();

  test<Test_arrangement<Open_traits, UseVector>>(mode + " open", std::make_shared<const Open_traits>(),
                                                 {FT(-3), FT(2), FT(5), FT(7)}, FT(3), failures);
  test<Closed_arr>(mode + " closed", closed_tr, {FT(0), FT(2), FT(5), FT(7), FT(10)}, FT(3), failures);
  test<Test_arrangement<Left_closed_traits, UseVector>>(mode + " left-closed",
                                                        std::make_shared<const Left_closed_traits>(FT(0)),
                                                        {FT(0), FT(2), FT(5), FT(7)}, FT(3), failures);
  test<Test_arrangement<Right_closed_traits, UseVector>>(mode + " right-closed",
                                                         std::make_shared<const Right_closed_traits>(FT(10)),
                                                         {FT(2), FT(5), FT(7), FT(10)}, FT(3), failures);
  test<Identified_arr>(
    mode + " identified", identified_tr,
    {Direction_2(1, 0), Direction_2(0, 1), Direction_2(-1, 0), Direction_2(0, -1), Direction_2(1, 1)},
    Direction_2(-1, -1), failures);
  if constexpr (UseVector) {
    using Binary_search_arr =
      Aoc::Arrangement_on_curve_1<Closed_traits, Aoc::Closed_topology_traits<Closed_traits, void, void, true>, true>;
    test<Binary_search_arr>(mode + " closed binary-search", closed_tr, {FT(0), FT(2), FT(5), FT(7), FT(10)}, FT(3),
                            failures);
  }

  test_rejection<Closed_arr>(mode + " closed: out of order", closed_tr, "3\n2\n5\n4\n", failures);
  test_rejection<Closed_arr>(mode + " closed: duplicate", closed_tr, "2\n2\n2\n", failures);
  test_rejection<Closed_arr>(mode + " closed: left boundary not first", closed_tr, "2\n2\n0\n", failures);
  test_rejection<Closed_arr>(mode + " closed: right boundary not last", closed_tr, "2\n10\n5\n", failures);
  test_rejection<Closed_arr>(mode + " closed: truncated", closed_tr, "3\n2\n5\n", failures);
  test_rejection<Identified_arr>(mode + " identified: identification point not last", identified_tr,
                                 "2\n1 0\n0 1\n", failures);
  test_rejection<Identified_arr>(mode + " identified: out of order", identified_tr, "2\n0 -1\n0 1\n", failures);
}

int main() {
  std::size_t failures = 0;
  test_all<false>("list", failures);
  test_all<true>("vector", failures);
  return (failures == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
