#include <iostream>
#include <string>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Arrangement_on_curve_1/Arrangement_on_curve_1.h>
#include <CGAL/Arrangement_on_curve_1/Ft_traits_1.h>
#include <CGAL/Arrangement_on_curve_1/insert.h>
#include <CGAL/Arrangement_on_curve_1/Unbounded_topology_traits.h>

using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;
using Geom_traits = CGAL::Arrangement_on_curve_1::Ft_traits_1<Kernel>;
using Point = Geom_traits::Point_1;
using Topo_traits = CGAL::Arrangement_on_curve_1::Unbounded_topology_traits<Point, std::string, void, true>;
using Arrangement = CGAL::Arrangement_on_curve_1::Arrangement_on_curve_1<Geom_traits, Topo_traits, true>;

int main() {
  auto traits_ptr = std::make_shared<const Geom_traits>();
  Arrangement arr(traits_ptr);

  std::cout << "Inserting raw coordinate fields into 1D track...\n";
  auto v2 = CGAL::Arrangement_on_curve_1::insert(arr, Kernel::FT(2.0));
  auto v3 = CGAL::Arrangement_on_curve_1::insert(arr, Kernel::FT(5.25));
  auto v1 = CGAL::Arrangement_on_curve_1::insert(arr, Kernel::FT(10.5));

  // Modify user-extended data fields on vertices
  auto v_data_map = arr.vertex_data_map();
  put(v_data_map, v1, "Right Node");
  put(v_data_map, v2, "Left Node");
  put(v_data_map, v3, "Center Node");

  // Traverse and inspect the arrangement
  std::cout << "\nResulting Sorted Sequence along the 1D Line:\n";
  auto v_pnt_map = arr.vertex_point_map();
  auto e = arr.unbounded_left_edge();
  while (arr.has_right_vertex(e)) {
    auto v = arr.right_vertex(e);
    std::cout << "  Vertex Point: (" << get(v_pnt_map, v) << ") | Extension ID: " << get(v_data_map, v) << "\n";
    e = arr.right_edge(v);
  }

  return 0;
}
