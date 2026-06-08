#include <iostream>
#include <string>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Arrangement_on_curve_1.h>
#include <CGAL/Arrangement_on_curve_1_functions.h>
#include <CGAL/Unbounded_topology_traits.h>
#include <CGAL/Line_3_traits_1.h>

using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;
using Traits = CGAL::Arrangement_on_curve_1::Line_3_traits_1<Kernel>;
using Topo = CGAL::Arrangement_on_curve_1::Unbounded_topology_traits<typename Traits::Point_1, std::string, void>;
using Arrangement = CGAL::Arrangement_on_curve_1::Arrangement_on_curve_1<Traits, Topo>;

int main() {
  // Define an infinite 3D line passing through (0,0,0) with direction (1,1,1)
  Kernel::Point_3 origin(0, 0, 0);
  Kernel::Point_3 direction_pt(1, 1, 1);
  Kernel::Line_3 master_line(origin, direction_pt);

  Traits traits(master_line);
  Topo topo;
  Arrangement arr(traits, topo);

  // Define 3D points that rest strictly on our 3D line trajectory
  Kernel::Point_3 p_middle(3, 3, 3);
  Kernel::Point_3 p_far(7, 7, 7);
  Kernel::Point_3 p_near(1, 1, 1);

  std::cout << "Inserting collinear 3D points into the 1D arrangement...\n";
  auto v_mid = CGAL::Arrangement_on_curve_1::insert(arr, p_middle);
  auto v_far = CGAL::Arrangement_on_curve_1::insert(arr, p_far);
  auto v_near = CGAL::Arrangement_on_curve_1::insert(arr, p_near);

  // Update string metadata using property maps
  auto v_data_map = arr.vertex_data_map();
  put(v_data_map, v_mid, std::string("Station B"));
  put(v_data_map, v_far, std::string("Station C"));
  put(v_data_map, v_near, std::string("Station A"));

  std::cout << "\nResulting Sorted Sequence along the 3D Line:\n";
  auto v_pnt_map = arr.vertex_point_map();
  auto e = arr.unbounded_edge();
  while (arr.has_right_vertex(e)) {
    auto v = arr.right_vertex(e);
    std::cout << "  Vertex Point: (" << get(v_pnt_map, v) << ") | Extension ID: " << get(v_data_map, v) << "\n";
    e = arr.right_edge(v);
  }

  return 0;
}
