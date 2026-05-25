#include <iostream>
#include <string>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Arrangement_on_curve_1.h>
#include <CGAL/Arrangement_on_curve_1_functions.h>
#include <CGAL/Line_3_traits_1.h>

using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;
using Traits = CGAL::Arrangement_on_curve_1::Line_3_traits_1<Kernel>;
using Arrangement = CGAL::Arrangement_on_curve_1::Arrangement_on_curve_1<Traits, std::string, void>;

int main() {
  // Define an infinite 3D line passing through (0,0,0) with direction (1,1,1)
  Kernel::Point_3 origin(0, 0, 0);
  Kernel::Point_3 direction_pt(1, 1, 1);
  Kernel::Line_3 master_line(origin, direction_pt);

  // Initialize our 1D topology space mapped to the 3D line geometry
  Traits traits(master_line);
  Arrangement arr(traits);

  // Define 3D points that rest strictly on our 3D line trajectory
  Kernel::Point_3 p_middle(3, 3, 3);
  Kernel::Point_3 p_far(7, 7, 7);
  Kernel::Point_3 p_near(1, 1, 1);

  std::cout << "Inserting collinear 3D points into the 1D arrangement...\n";
  auto v_mid  = CGAL::Arrangement_on_curve_1::insert(arr, p_middle);
  auto v_far  = CGAL::Arrangement_on_curve_1::insert(arr, p_far);
  auto v_near = CGAL::Arrangement_on_curve_1::insert(arr, p_near);

  // Attach metadata strings to the extended vertices
  v_mid->data()  = "Station B";
  v_far->data()  = "Station C";
  v_near->data() = "Station A";

  std::cout << "\nResulting Sorted Sequence along the 3D Line:\n";
  for (auto vit = arr.vertices_begin(); vit != arr.vertices_end(); ++vit) {
    std::cout << "  Vertex 3D Coordinate: (" << vit->point() << ") | Descriptor: " << vit->data() << "\n";
  }

  return 0;
}
