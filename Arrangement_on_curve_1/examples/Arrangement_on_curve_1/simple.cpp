#include <iostream>
#include <string>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Arrangement_on_curve_1/Arrangement_on_curve_1.h>
#include <CGAL/Arrangement_on_curve_1/Unbounded_topology_traits.h>
#include <CGAL/Arrangement_on_curve_1/insert.h>

// Minimal geometry traits model utilizing direct coordinate FT parameters
template <typename Kernel>
class Aoc_ft_traits_1 {
public:
  using Point_1 = typename Kernel::FT;

  class Compare_x_1 {
  public:
    CGAL::Comparison_result operator()(const Point_1& p1, const Point_1& p2) const {
      if (p1 < p2) return CGAL::SMALLER;
      if (p1 > p2) return CGAL::LARGER;
      return CGAL::EQUAL;
    }
  };

  Compare_x_1 compare_x_1_object() const { return Compare_x_1(); }
};

using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;
using Geom_traits = Aoc_ft_traits_1<Kernel>;
using Point = Geom_traits::Point_1;
using Topo_traits = CGAL::Arrangement_on_curve_1::Unbounded_topology_traits<Point, std::string, void>;
using Arrangement = CGAL::Arrangement_on_curve_1::Arrangement_on_curve_1<Geom_traits, Topo_traits>;

int main() {
  auto traits_ptr = std::make_shared<const Geom_traits>();
  Arrangement arr(traits_ptr);

  std::cout << "Inserting raw coordinate fields into 1D track...\n";
  auto v1 = CGAL::Arrangement_on_curve_1::insert(arr, Kernel::FT(10.5));
  auto v2 = CGAL::Arrangement_on_curve_1::insert(arr, Kernel::FT(2.0));
  auto v3 = CGAL::Arrangement_on_curve_1::insert(arr, Kernel::FT(5.25));

  // Modify user-extended data fields on vertices
  v1->data() = "Right Node";
  v2->data() = "Left Node";
  v3->data() = "Center Node";

  // Traverse and inspect the arrangement
  std::cout << "\nResulting Sorted Sequence along the 2D Line:\n";
  auto v_pnt_map = arr.vertex_point_map();
  auto e = arr.unbounded_edge();
  while (arr.has_right_vertex(e)) {
    auto v = arr.right_vertex(e);
    std::cout << "  Vertex Point: (" << v->point() << ") | Extension ID: " << v->data() << "\n";
    e = arr.right_edge(v);
  }

  return 0;
}
