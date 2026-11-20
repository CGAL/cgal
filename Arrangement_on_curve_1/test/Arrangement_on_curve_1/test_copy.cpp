// Tests the copy constructor, the copy assignment operator, the move
// constructor, clear(), and swap() of the arrangement (and thus of the topology
// traits) for every topology, in list and vector mode.

#include <cstdlib>
#include <utility>

#include "test_utils.h"

template <typename Arr>
void test(const std::string& name, std::shared_ptr<const typename Arr::Geometry_traits_1> tr,
          const std::vector<typename Arr::Point_1>& pts, std::size_t& failures) {
  auto make = [&]() {
    Arr arr(tr);
    for (const auto& p : pts) Aoc::insert(arr, p);
    return arr;
  };
  std::string why;
  bool ok = true;

  Arr arr = make();
  Arr copy(arr);
  ok = ok && report(is_valid(copy, pts, why), name + " (copy constructor)", why, failures);

  Arr assigned = make();
  assigned.clear();
  assigned = arr;
  ok = ok && report(is_valid(assigned, pts, why), name + " (copy assignment)", why, failures);

  Arr moved(std::move(copy));
  ok = ok && report(is_valid(moved, pts, why), name + " (move constructor)", why, failures);

  Arr cleared = make();
  cleared.clear();
  ok = ok && report(is_valid(cleared, {}, why), name + " (clear)", why, failures);

  cleared.topology_traits().swap(assigned.topology_traits());
  ok = ok && report(is_valid(cleared, pts, why) && is_valid(assigned, {}, why), name + " (swap)", why, failures);

  // The source must be intact.
  ok = ok && report(is_valid(arr, pts, why), name + " (source)", why, failures);
  std::cout << (ok ? "PASS " : "FAIL ") << name << std::endl;
}

template <bool UseVector>
void test_all(const std::string& mode, std::size_t& failures) {
  test<Test_arrangement<Open_traits, UseVector>>(mode + " open", std::make_shared<const Open_traits>(),
                                                 {FT(4), FT(-3), FT(2)}, failures);
  test<Test_arrangement<Closed_traits, UseVector>>(mode + " closed",
                                                   std::make_shared<const Closed_traits>(FT(0), FT(10)),
                                                   {FT(4), FT(0), FT(10), FT(2)}, failures);
  test<Test_arrangement<Left_closed_traits, UseVector>>(mode + " left-closed",
                                                        std::make_shared<const Left_closed_traits>(FT(0)),
                                                        {FT(4), FT(0), FT(2)}, failures);
  test<Test_arrangement<Right_closed_traits, UseVector>>(mode + " right-closed",
                                                         std::make_shared<const Right_closed_traits>(FT(10)),
                                                         {FT(4), FT(10), FT(2)}, failures);
  test<Test_arrangement<Identified_traits, UseVector>>(mode + " identified",
                                                       std::make_shared<const Identified_traits>(),
                                                       {Direction_2(0, 1), Direction_2(1, 0), Direction_2(-1, 0)},
                                                       failures);
}

int main() {
  std::size_t failures = 0;
  test_all<false>("list", failures);
  test_all<true>("vector", failures);
  return (failures == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
