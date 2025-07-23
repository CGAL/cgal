
#include <vector>

#include "arr_polylines.h"

int main() {
  Traits traits;
  auto cst_x_curve = traits.construct_x_monotone_curve_2_object();
  std::vector<Point> segs{Point(5, 5), Point(4, 4), Point(3, 3), Point(2, 2), Point(1, 1), Point(0, 0)};
  auto polyline = cst_x_curve(segs.begin(), segs.end());

  auto approx = traits.approximate_2_object();
  std::vector<Traits::Approximate_point_2> approx_points;
  // Specify l2r = true here.
  approx(polyline, 1, std::back_inserter(approx_points), true);

  std::cout << "Approximate points:\n";
  for(const auto& pt : approx_points) {
    std::cout << pt << "\n";
  }

  // Expected output:
  // Approximate points:
  // 0 0
  // 1 1
  // 2 2
  // 3 3
  // 4 4
  // 5 5

  // Got:
  // Approximate points:
  // 5 5
  // 4 4
  // 3 3
  // 2 2
  // 1 1
  // 0 0
}