// Tests that Copy_overlay_visitor copies the user data of the input cells to
// the resulting cells, that the allocator passed to the constructor is used,
// and that overlay works in binary-search mode.

#include <cstdlib>
#include <map>
#include <memory>
#include <string>
#include <utility>

#include "test_utils.h"
#include <CGAL/Arrangement_on_curve_1/Open_topology_traits.h>
#include <CGAL/Arrangement_on_curve_1/Closed_topology_traits.h>
#include <CGAL/Arrangement_on_curve_1/overlay.h>
#include <CGAL/Arrangement_on_curve_1/Copy_overlay_visitor.h>

// A stateful allocator that carries a tag, which identifies its origin.
template <typename T>
struct Tagged_allocator {
  using value_type = T;

  int m_tag = 0;

  Tagged_allocator() = default;
  explicit Tagged_allocator(int tag) : m_tag(tag) {}
  template <typename U> Tagged_allocator(const Tagged_allocator<U>& other) : m_tag(other.m_tag) {}

  T* allocate(std::size_t n) { return std::allocator<T>().allocate(n); }
  void deallocate(T* p, std::size_t n) { std::allocator<T>().deallocate(p, n); }

  template <typename U> bool operator==(const Tagged_allocator<U>& other) const { return m_tag == other.m_tag; }
  template <typename U> bool operator!=(const Tagged_allocator<U>& other) const { return m_tag != other.m_tag; }
};

// Overlays two arrangements on a line with vertex and edge data using
// Copy_overlay_visitor.  The data of a resulting vertex must be the data of the
// input vertex at the same point, where the first arrangement takes
// precedence, and the data of a resulting edge must be the data of the edge of
// the first arrangement that contains it.
void test_copy_visitor(std::size_t& failures) {
  using Topology_traits = Aoc::Open_topology_traits<Open_traits, std::string, int>;
  using Arr = Aoc::Arrangement_on_curve_1<Open_traits, Topology_traits>;
  auto tr = std::make_shared<const Open_traits>();
  Arr arr_a(tr), arr_b(tr), arr_res(tr);
  auto vm_a = arr_a.vertex_data_map();
  auto vm_b = arr_b.vertex_data_map();
  for (int x : {1, 3, 5}) put(vm_a, Aoc::insert(arr_a, FT(x)), "a" + std::to_string(x));
  for (int x : {3, 4}) put(vm_b, Aoc::insert(arr_b, FT(x)), "b" + std::to_string(x));
  {
    // Tag the edges of the first arrangement by their position along the line.
    int tag = 0;
    auto em_a = arr_a.edge_data_map();
    auto e = arr_a.left_boundary_edge();
    for (std::size_t i = 0; i < arr_a.number_of_edges(); ++i) {
      if (i > 0) e = arr_a.right_edge(arr_a.right_vertex(e));
      put(em_a, e, tag++);
    }
  }

  Aoc::Copy_overlay_visitor<Arr, Arr, Arr> visitor(arr_a, arr_b, arr_res);
  Aoc::overlay(arr_a, arr_b, arr_res, visitor);

  std::string why;
  bool ok = is_valid(arr_res, {FT(1), FT(3), FT(4), FT(5)}, why);
  auto vm_res = arr_res.vertex_data_map();
  const std::map<int, std::string> expected_labels{{1, "a1"}, {3, "a3"}, {4, "b4"}, {5, "a5"}};
  for (const auto& [x, label] : expected_labels) {
    if (ok && (get(vm_res, find_vertex(arr_res, FT(x))) != label)) {
      ok = false;
      why = "the data of the vertex at " + std::to_string(x);
    }
  }
  // A point in the interior of every resulting edge locates the edge of the
  // first arrangement that contains it.
  for (auto e : arr_res.edges()) {
    if (! ok) break;
    const bool has_l = arr_res.has_left_vertex(e), has_r = arr_res.has_right_vertex(e);
    FT l = has_l ? get(arr_res.vertex_point_map(), arr_res.left_vertex(e)) : FT(0);
    FT r = has_r ? get(arr_res.vertex_point_map(), arr_res.right_vertex(e)) : FT(0);
    FT q = (has_l && has_r) ? (l + r) / 2 : (has_l ? l + 1 : r - 1);
    auto loc = Aoc::locate(std::as_const(arr_a), q);
    if ((loc.index() != 1) || (get(arr_res.edge_data_map(), e) != get(arr_a.edge_data_map(), std::get<1>(loc)))) {
      ok = false;
      why = "the data of an edge";
    }
  }
  std::cout << (report(ok, "Copy_overlay_visitor", why, failures) ? "PASS " : "FAIL ") << "Copy_overlay_visitor"
            << std::endl;
}

// Constructs an arrangement with a tagged allocator, which must be used by its
// topology traits and survive copying.
void test_allocator(std::size_t& failures) {
  using Topology_traits = Aoc::Open_topology_traits<Open_traits, void, void, false, Tagged_allocator<char>>;
  using Arr = Aoc::Arrangement_on_curve_1<Open_traits, Topology_traits>;
  Arr arr(std::make_shared<const Open_traits>(), Tagged_allocator<char>(7));
  for (int x : {4, -1, 2}) Aoc::insert(arr, FT(x));
  Arr copy(arr);
  std::string why;
  bool ok = is_valid(arr, {FT(4), FT(-1), FT(2)}, why) && is_valid(copy, {FT(4), FT(-1), FT(2)}, why);
  if (ok && ((arr.get_allocator().m_tag != 7) || (copy.get_allocator().m_tag != 7))) {
    ok = false;
    why = "the allocator is not used";
  }
  std::cout << (report(ok, "allocator", why, failures) ? "PASS " : "FAIL ") << "allocator" << std::endl;
}

// Overlays two closed arrangements in binary-search mode.
void test_binary_search_overlay(std::size_t& failures) {
  using Arr =
    Aoc::Arrangement_on_curve_1<Closed_traits, Aoc::Closed_topology_traits<Closed_traits, void, void, true>, true>;
  auto tr = std::make_shared<const Closed_traits>(FT(0), FT(10));
  Arr arr_a(tr), arr_b(tr), arr_res(tr);
  for (int x : {0, 4, 10}) Aoc::insert(arr_a, FT(x));
  for (int x : {2, 4, 7}) Aoc::insert(arr_b, FT(x));
  Aoc::overlay(arr_a, arr_b, arr_res);
  std::string why;
  bool ok = is_valid(arr_res, {FT(0), FT(2), FT(4), FT(7), FT(10)}, why);
  std::cout << (report(ok, "binary-search overlay", why, failures) ? "PASS " : "FAIL ") << "binary-search overlay"
            << std::endl;
}

int main() {
  std::size_t failures = 0;
  test_copy_visitor(failures);
  test_allocator(failures);
  test_binary_search_overlay(failures);
  return (failures == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
