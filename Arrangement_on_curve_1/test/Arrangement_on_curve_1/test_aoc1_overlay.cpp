// Tests overlay() for every topology, in list and vector mode: every pair of
// arrangements constructed from subsets of a set of points is overlaid.  The
// result must be a valid arrangement of the union of the points, and the
// visitor must report every vertex and every edge of the result exactly once.

#include <cstdlib>
#include <map>

#include "test_utils.h"
#include <CGAL/Arrangement_on_curve_1/overlay.h>
#include <CGAL/Arrangement_on_curve_1/Direction_geom_traits_2_adaptor_1.h>
#include <CGAL/Arrangement_on_curve_1/Identified_geom_traits_2_adaptor_1.h>
#include <CGAL/Arr_geodesic_arc_on_sphere_traits_2.h>

// An overlay visitor that counts the reports of every result vertex and edge.
template <typename ArrangementR>
class Counting_visitor {
public:
  std::map<const void*, std::size_t> m_vertex_reports;
  std::map<const void*, std::size_t> m_edge_reports;

  template <typename A, typename B>
  void create_vertex_vv(A, B, typename ArrangementR::Vertex_descriptor v) { ++m_vertex_reports[key(v)]; }

  template <typename A, typename B>
  void create_vertex_ve(A, B, typename ArrangementR::Vertex_descriptor v) { ++m_vertex_reports[key(v)]; }

  template <typename A, typename B>
  void create_vertex_ev(A, B, typename ArrangementR::Vertex_descriptor v) { ++m_vertex_reports[key(v)]; }

  template <typename A, typename B>
  void create_edge_ee(A, B, typename ArrangementR::Edge_descriptor e) { ++m_edge_reports[key(e)]; }

  // determines whether every vertex and every edge of `arr` is reported exactly once.
  bool is_valid(const ArrangementR& arr) const {
    if ((m_vertex_reports.size() != arr.number_of_vertices()) || (m_edge_reports.size() != arr.number_of_edges()))
      return false;
    for (const auto& [k, n] : m_vertex_reports) if (n != 1) return false;
    for (const auto& [k, n] : m_edge_reports) if (n != 1) return false;
    return true;
  }

private:
  // A vector-mode descriptor is an index; a list-mode descriptor is an iterator.
  template <typename Descriptor>
  static const void* key(Descriptor d) {
    if constexpr (std::is_integral_v<Descriptor>) return reinterpret_cast<const void*>(d + 1);
    else return static_cast<const void*>(&*d);
  }
};

template <typename Arr>
void test(const std::string& name, std::shared_ptr<const typename Arr::Geometry_traits_1> tr,
          const std::vector<typename Arr::Point_1>& pts, std::size_t& failures) {
  std::size_t local_failures = 0;
  const std::size_t n = pts.size();
  for (unsigned mask_a = 0; mask_a < (1u << n); ++mask_a) {
    for (unsigned mask_b = 0; mask_b < (1u << n); ++mask_b) {
      Arr arr_a(tr), arr_b(tr), arr_res(tr);
      std::vector<typename Arr::Point_1> expected;
      for (unsigned i = 0; i < n; ++i) {
        if (mask_a & (1u << i)) Aoc::insert(arr_a, pts[i]);
        if (mask_b & (1u << i)) Aoc::insert(arr_b, pts[i]);
        if ((mask_a | mask_b) & (1u << i)) expected.push_back(pts[i]);
      }
      Counting_visitor<Arr> visitor;
      Aoc::overlay(arr_a, arr_b, arr_res, visitor);
      std::string why;
      bool ok = is_valid(arr_res, expected, why);
      if (ok && ! visitor.is_valid(arr_res)) {
        ok = false;
        why = "a cell is not reported exactly once";
      }
      report(ok, name, why, local_failures);
    }
  }
  failures += local_failures;
  std::cout << (local_failures ? "FAIL " : "PASS ") << name << std::endl;
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
  using Direction_traits = Aoc::Direction_geom_traits_2_adaptor_1<Kernel>;
  test<Test_arrangement<Direction_traits, UseVector>>(
    mode + " identified (directions)", std::make_shared<const Direction_traits>(Kernel(), Direction_2(1, 2)),
    {Direction_2(1, 2), Direction_2(0, 1), Direction_2(-1, -2), Direction_2(0, -1), Direction_2(1, 1)}, failures);
  // Points of the equator of the sphere; (-1, 0, 0) is the identification point.
  using Sphere_traits_2 = CGAL::Arr_geodesic_arc_on_sphere_traits_2<Kernel>;
  using Sphere_traits_1 = Aoc::Identified_geom_traits_2_adaptor_1<Sphere_traits_2>;
  auto sphere_tr = std::make_shared<const Sphere_traits_2>();
  auto ctr = sphere_tr->construct_point_2_object();
  test<Test_arrangement<Sphere_traits_1, UseVector>>(
    mode + " identified (equator)", std::make_shared<const Sphere_traits_1>(sphere_tr),
    {ctr(-1, 0, 0), ctr(1, 0, 0), ctr(0, 1, 0), ctr(0, -1, 0), ctr(1, 1, 0)}, failures);
}

int main() {
  std::size_t failures = 0;
  test_all<false>("list", failures);
  test_all<true>("vector", failures);
  return (failures == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
