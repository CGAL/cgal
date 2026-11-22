// Tests remove() for every topology, in list and vector mode: every vertex of
// every arrangement constructed by inserting a permutation of a subset of a set
// of points is removed on its own, and all the vertices are removed one by
// one.  Also tests insert_before() and insert_after() next to the boundary
// vertices and the identification vertex.

#include <cstdlib>

#include "test_utils.h"

template <typename Arr>
void test(const std::string& name, std::shared_ptr<const typename Arr::Geometry_traits_1> tr,
          const std::vector<typename Arr::Point_1>& pts, std::size_t& failures) {
  std::size_t local_failures = 0;
  for_each_sequence(pts.size(), [&](const std::vector<int>& seq) {
    Arr arr(tr);
    std::vector<typename Arr::Point_1> inserted;
    for (int i : seq) {
      Aoc::insert(arr, pts[i]);
      inserted.push_back(pts[i]);
    }
    std::string why;

    // Remove every vertex on its own.
    for (std::size_t k = 0; k < inserted.size(); ++k) {
      Arr copy(arr);
      copy.remove(find_vertex(copy, inserted[k]));
      auto rest = inserted;
      rest.erase(rest.begin() + k);
      report(is_valid(copy, rest, why), name + " (single)", why, local_failures);
    }

    // Remove all the vertices, in insertion order.
    auto rest = inserted;
    while (! rest.empty()) {
      arr.remove(find_vertex(arr, rest.front()));
      rest.erase(rest.begin());
      if (! report(is_valid(arr, rest, why), name + " (sequential)", why, local_failures)) break;
    }
  });
  failures += local_failures;
  std::cout << (local_failures ? "FAIL " : "PASS ") << name << std::endl;
}

template <bool UseVector>
void test_neighbors(const std::string& mode, std::size_t& failures) {
  std::string why;
  {
    using Arr = Test_arrangement<Closed_traits, UseVector>;
    Arr arr(std::make_shared<const Closed_traits>(FT(0), FT(10)));
    auto vl = Aoc::insert(arr, FT(0));
    auto vr = Aoc::insert(arr, FT(10));
    auto v5 = arr.insert_after(vl, FT(5));      // next to the left boundary vertex
    arr.insert_before(vr, FT(7));               // next to the right boundary vertex
    arr.insert_after(v5, FT(6));
    arr.insert_before(v5, FT(3));
    bool ok = is_valid(arr, {FT(0), FT(3), FT(5), FT(6), FT(7), FT(10)}, why);
    if (ok) {
      Aoc::insert(arr, FT(9));
      Aoc::insert(arr, FT(1));
      arr.remove(find_vertex(arr, FT(10)));
      arr.remove(find_vertex(arr, FT(0)));
      ok = is_valid(arr, {FT(1), FT(3), FT(5), FT(6), FT(7), FT(9)}, why);
    }
    std::string name = mode + " closed insert_before/insert_after";
    std::cout << (report(ok, name, why, failures) ? "PASS " : "FAIL ") << name << std::endl;
  }
  {
    using Arr = Test_arrangement<Identified_traits, UseVector>;
    auto tr = std::make_shared<const Identified_traits>();
    Arr arr1(tr);                               // a single vertex, which is a self-loop
    auto v = Aoc::insert(arr1, Direction_2(0, 1));
    arr1.insert_before(v, Direction_2(1, 1));
    Arr arr2(tr);
    auto w = Aoc::insert(arr2, Direction_2(0, 1));
    arr2.insert_after(w, Direction_2(-1, 0));
    Arr arr3(tr);                               // next to the identification vertex
    Aoc::insert(arr3, Direction_2(0, 1));
    auto v_id = Aoc::insert(arr3, Direction_2(1, 0));
    arr3.insert_after(v_id, Direction_2(1, 1));
    arr3.insert_before(v_id, Direction_2(0, -1));
    Arr arr4(tr);                               // the identification vertex alone, which is a self-loop
    Aoc::insert(arr4, Direction_2(1, 0));
    arr4.remove(Aoc::insert(arr4, Direction_2(0, 1)));
    arr4.insert_after(find_vertex(arr4, Direction_2(1, 0)), Direction_2(-1, 0));
    bool ok = is_valid(arr1, {Direction_2(1, 1), Direction_2(0, 1)}, why) &&
      is_valid(arr2, {Direction_2(0, 1), Direction_2(-1, 0)}, why) &&
      is_valid(arr3, {Direction_2(1, 1), Direction_2(0, 1), Direction_2(0, -1), Direction_2(1, 0)}, why) &&
      is_valid(arr4, {Direction_2(-1, 0), Direction_2(1, 0)}, why);
    std::string name = mode + " identified insert_before/insert_after";
    std::cout << (report(ok, name, why, failures) ? "PASS " : "FAIL ") << name << std::endl;
  }
}

template <bool UseVector>
void test_all(const std::string& mode, std::size_t& failures) {
  test<Test_arrangement<Open_traits, UseVector>>(mode + " open", std::make_shared<const Open_traits>(),
                                                 {FT(-3), FT(2), FT(5), FT(7)}, failures);
  test<Test_arrangement<Closed_traits, UseVector>>(mode + " closed",
                                                   std::make_shared<const Closed_traits>(FT(0), FT(10)),
                                                   {FT(0), FT(2), FT(5), FT(7), FT(10)}, failures);
  test<Test_arrangement<Left_closed_traits, UseVector>>(mode + " left-closed",
                                                        std::make_shared<const Left_closed_traits>(FT(0)),
                                                        {FT(0), FT(2), FT(5), FT(7)}, failures);
  test<Test_arrangement<Right_closed_traits, UseVector>>(mode + " right-closed",
                                                         std::make_shared<const Right_closed_traits>(FT(10)),
                                                         {FT(2), FT(5), FT(7), FT(10)}, failures);
  test<Test_arrangement<Identified_traits, UseVector>>(
    mode + " identified", std::make_shared<const Identified_traits>(),
    {Direction_2(1, 0), Direction_2(0, 1), Direction_2(-1, 0), Direction_2(0, -1), Direction_2(1, 1)}, failures);
  test_neighbors<UseVector>(mode, failures);
}

int main() {
  std::size_t failures = 0;
  test_all<false>("list", failures);
  test_all<true>("vector", failures);
  return (failures == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
