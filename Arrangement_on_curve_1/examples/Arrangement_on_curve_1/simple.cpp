#include <iostream>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Arrangement_on_curve_1/Arrangement_on_curve_1.h>
#include <CGAL/Arrangement_on_curve_1/insert.h>
#include <CGAL/Arrangement_on_curve_1/Unbounded_topology_traits.h>
#include <CGAL/Arrangement_on_curve_1/Line_2_traits_1.h>

using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;
using Traits = CGAL::Arrangement_on_curve_1::Line_2_traits_1<Kernel>;
using Topo = CGAL::Arrangement_on_curve_1::Unbounded_topology_traits<typename Traits::Point_1, int, void>;
using Arrangement = CGAL::Arrangement_on_curve_1::Arrangement_on_curve_1<Traits, Topo>;

int main() {
  // Define an infinite 2D line track: y = 0.5x
  Kernel::Point_2 origin(0, 0);
  Kernel::Point_2 direction_pt(2, 1);
  Kernel::Line_2 master_line(origin, direction_pt);

  Traits traits(master_line);
  Topo topo;
  Arrangement arr(traits, topo);

  // Collinear points along our line track
  Kernel::Point_2 p_middle(4, 2);
  Kernel::Point_2 p_far(8, 4);
  Kernel::Point_2 p_near(2, 1);

  std::cout << "Inserting collinear 2D points into the 1D arrangement...\n";
  auto v_mid = CGAL::Arrangement_on_curve_1::insert(arr, p_middle);
  auto v_far = CGAL::Arrangement_on_curve_1::insert(arr, p_far);
  auto v_near = CGAL::Arrangement_on_curve_1::insert(arr, p_near);

  // Fetch the data map and modify user properties via property map interfaces
  auto v_data_map = arr.vertex_data_map();
  put(v_data_map, v_mid, 42);
  put(v_data_map, v_far, 84);
  put(v_data_map, v_near, 21);

  std::cout << "\nResulting Sorted Sequence along the 2D Line:\n";
  auto v_pnt_map = arr.vertex_point_map();
  auto e = arr.unbounded_edge();
  while (arr.has_right_vertex(e)) {
    auto v = arr.right_vertex(e);
    std::cout << "  Vertex Point: (" << get(v_pnt_map, v) << ") | Extension ID: " << get(v_data_map, v) << "\n";
    e = arr.right_edge(v);
  }

  return 0;
}
