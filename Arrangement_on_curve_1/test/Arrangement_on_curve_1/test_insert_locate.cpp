// Tests insert() and locate() for every topology, in list and vector mode:
// every permutation of every subset of a set of points is inserted; then every
// point is inserted again, which must not change the arrangement.

#include <cstdlib>

#include "test_utils.h"

template <typename Arr>
void test(const std::string& name, std::shared_ptr<const typename Arr::Geometry_traits_1> tr,
          const std::vector<typename Arr::Point_1>& pts, const typename Arr::Point_1& absent, std::size_t& failures) {
  std::size_t local_failures = 0;
  for_each_sequence(pts.size(), [&](const std::vector<int>& seq) {
    Arr arr(tr);
    std::vector<typename Arr::Point_1> inserted;
    for (int i : seq) {
      Aoc::insert(arr, pts[i]);
      inserted.push_back(pts[i]);
    }
    for (int i : seq) Aoc::insert(arr, pts[i]);
    std::string why;
    bool ok = is_valid(arr, inserted, why);
    if (ok && (Aoc::locate(arr, absent).index() != 1)) {
      ok = false;
      why = "an absent point is located at a vertex";
    }
    report(ok, name, why, local_failures);
  });
  failures += local_failures;
  std::cout << (local_failures ? "FAIL " : "PASS ") << name << std::endl;
}

template <bool UseVector>
void test_all(const std::string& mode, std::size_t& failures) {
  test<Test_arrangement<Open_traits, UseVector>>(mode + " open", std::make_shared<const Open_traits>(),
                                                 {FT(-3), FT(2), FT(5), FT(7)}, FT(1), failures);
  test<Test_arrangement<Closed_traits, UseVector>>(mode + " closed",
                                                   std::make_shared<const Closed_traits>(FT(0), FT(10)),
                                                   {FT(0), FT(2), FT(5), FT(7), FT(10)}, FT(1), failures);
  test<Test_arrangement<Left_closed_traits, UseVector>>(mode + " left-closed",
                                                        std::make_shared<const Left_closed_traits>(FT(0)),
                                                        {FT(0), FT(2), FT(5), FT(7)}, FT(1), failures);
  test<Test_arrangement<Right_closed_traits, UseVector>>(mode + " right-closed",
                                                         std::make_shared<const Right_closed_traits>(FT(10)),
                                                         {FT(2), FT(5), FT(7), FT(10)}, FT(1), failures);
  test<Test_arrangement<Identified_traits, UseVector>>(
    mode + " identified", std::make_shared<const Identified_traits>(),
    {Direction_2(1, 0), Direction_2(0, 1), Direction_2(-1, 0), Direction_2(0, -1), Direction_2(1, 1)},
    Direction_2(-3, -1), failures);
}

int main() {
  std::size_t failures = 0;
  test_all<false>("list", failures);
  test_all<true>("vector", failures);
  return (failures == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
