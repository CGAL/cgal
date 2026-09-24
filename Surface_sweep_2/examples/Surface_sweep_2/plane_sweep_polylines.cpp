//! \file examples/Surface_sweep_2/plane_sweep.cpp
// Computing intersection points among curves using the surface-sweep alg.

#include <vector>
#include <chrono>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Surface_sweep_2_algorithms.h>
#include <CGAL/Arr_polyline_traits_2.h>

namespace ch = std::chrono;

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point = Kernel::Point_2;
using Segment_traits = CGAL::Arr_segment_traits_2<Kernel>;
using Traits = CGAL::Arr_polyline_traits_2<Segment_traits>;
using Curve = Traits::Curve_2;

int main() {
  Traits traits;
  auto polyline_construct = traits.construct_curve_2_object();

  // Construct the input curves.
  std::vector<Curve> curves;

  std::vector<Point> points1 = { Point(1, 1), Point(3,3), Point(0,6) };
  auto cv1 = polyline_construct(points1.begin(), points1.end());
  for (std::size_t i = 0; i < 2000000; ++i) curves.push_back(cv1);

  std::vector<Point> points2 = { Point(4, 0), Point(2,3), Point(5,6) };
  auto cv2 = polyline_construct(points2.begin(), points2.end());
  for (std::size_t i = 0; i < 2000000; ++i) curves.push_back(cv2);

  auto start = ch::high_resolution_clock::now();
  auto rc = CGAL::Surface_sweep_2::do_intersect(curves.begin(), curves.end());
  auto finish = ch::high_resolution_clock::now();
  auto duration = ch::duration_cast<ch::microseconds>(finish - start);
  std::cout << "time: " << duration.count() << std::endl;

  if (! rc) {
    std::cerr << "Failed to detect intersection\n";
    return -1;
  }

  return 0;
}
