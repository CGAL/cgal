#include <iostream>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Arrangement_on_curve_1.h>
#include <CGAL/Arrangement_on_curve_1_functions.h>
#include <CGAL/Line_2_traits_1.h>

using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;
using Traits = CGAL::Arrangement_on_curve_1::Line_2_traits_1<Kernel>;
using Arrangement = CGAL::Arrangement_on_curve_1::Arrangement_on_curve_1<Traits, int, void>;

int main() {
  // Define an infinite 2D line track: y = 0.5x
  Kernel::Point_2 origin(0, 0);
  Kernel::Point_2 direction_pt(2, 1);
  Kernel::Line_2 master_line(origin, direction_pt);

  // Initialize traits using the capitalized class and file matching layout
  Traits traits(master_line);
  Arrangement arr(traits);

  // Collinear points along our line track
  Kernel::Point_2 p_middle(4, 2);
  Kernel::Point_2 p_far(8, 4);
  Kernel::Point_2 p_near(2, 1);

  std::cout << "Inserting collinear 2D points into the 1D arrangement...\n";
  auto v_mid  = CGAL::Arrangement_on_curve_1::insert(arr, p_middle);
  auto v_far  = CGAL::Arrangement_on_curve_1::insert(arr, p_far);
  auto v_near = CGAL::Arrangement_on_curve_1::insert(arr, p_near);

  v_mid->data()  = 42;
  v_far->data()  = 84;
  v_near->data() = 21;

  std::cout << "\nResulting Sorted Sequence along the 2D Line:\n";
  for (auto vit = arr.vertices_begin(); vit != arr.vertices_end(); ++vit) {
    std::cout << "  Vertex Point: (" << vit->point() << ") | Extension ID: " << vit->data() << "\n";
  }

  return 0;
}
