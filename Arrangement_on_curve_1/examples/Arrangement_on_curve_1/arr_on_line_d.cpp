#include <iostream>
#include <vector>

#include <CGAL/Cartesian_d.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>
#include <CGAL/Arrangement_on_curve_1.h>
#include <CGAL/Arrangement_on_curve_1_functions.h>
#include <CGAL/Line_d_traits_1.h>

// Use Quotient<MP_Float> to provide full Field operations (division) for Kernel_d
using NT          = CGAL::Quotient<CGAL::MP_Float>;
using Kernel      = CGAL::Cartesian_d<NT>;
using FT          = typename Kernel::FT;
using Traits      = CGAL::Arrangement_on_curve_1::Line_d_traits_1<Kernel>;
using Arrangement = CGAL::Arrangement_on_curve_1::Arrangement_on_curve_1<Traits, int, void>;

int main() {
  // Construct coordinates using the field-capable number type NT
  std::vector<NT> orig_coords = {NT(0), NT(0), NT(0), NT(0)};
  std::vector<NT> dir_coords = {NT(1), NT(2), NT(3), NT(4)};

  Kernel::Point_d origin(4, orig_coords.begin(), orig_coords.end());
  Kernel::Direction_d direction(Kernel::Vector_d(4, dir_coords.begin(), dir_coords.end()));
  Kernel::Line_d master_line(origin, direction);

  // Initialize traits and arrangement track
  Traits traits(master_line);
  Arrangement arr(traits);

  // Define our d-dimensional points using valid field coordinate sets
  std::vector<NT> c_mid = {NT(3), NT(6), NT(9), NT(12)};
  std::vector<NT> c_far = {NT(5), NT(10), NT(15), NT(20)};
  std::vector<NT> c_near = {NT(1), NT(2), NT(3), NT(4)};

  Kernel::Point_d p_middle(4, c_mid.begin(), c_mid.end());
  Kernel::Point_d p_far(4, c_far.begin(), c_far.end());
  Kernel::Point_d p_near(4, c_near.begin(), c_near.end());

  std::cout << "Inserting collinear 4D points into the 1D arrangement...\n";
  auto v_mid = CGAL::Arrangement_on_curve_1::insert(arr, p_middle);
  auto v_far = CGAL::Arrangement_on_curve_1::insert(arr, p_far);
  auto v_near = CGAL::Arrangement_on_curve_1::insert(arr, p_near);

  v_mid->data()  = 2;
  v_far->data()  = 3;
  v_near->data() = 1;

  std::cout << "\nResulting Sorted Sequence along the 4D Line:\n";
  for (auto vit = arr.vertices_begin(); vit != arr.vertices_end(); ++vit) {
    std::cout << "  Vertex 4D Coordinates: (";
    for (auto cit = vit->point().cartesian_begin(); cit != vit->point().cartesian_end(); ++cit)
      std::cout << " " << *cit;
    std::cout << " ) | Sort Priority ID: " << vit->data() << "\n";
  }

  return 0;
}
