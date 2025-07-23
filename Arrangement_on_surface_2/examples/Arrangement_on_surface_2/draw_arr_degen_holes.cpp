#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/draw_arrangement_2.h>

using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;
using Traits = CGAL::Arr_segment_traits_2<Kernel>;
using Point = Traits::Point_2;
using Arrangement = CGAL::Arrangement_2<Traits>;
using X_monotone_curve = Traits::X_monotone_curve_2;

int main() {
  Arrangement arr;
  auto traits = arr.traits();
  auto cst_x_curve = traits->construct_x_monotone_curve_2_object();

  // make a hexagon centered at the origin with radius 10
  double radius = 10.0;
  std::array<Point, 6> hexagon_points = {{
      {radius * cos(0 * CGAL_PI / 3), radius * sin(0 * CGAL_PI / 3)}, // 0°
      {radius * cos(1 * CGAL_PI / 3), radius * sin(1 * CGAL_PI / 3)}, // 60°
      {radius * cos(2 * CGAL_PI / 3), radius * sin(2 * CGAL_PI / 3)}, // 120°
      {radius * cos(3 * CGAL_PI / 3), radius * sin(3 * CGAL_PI / 3)}, // 180°
      {radius * cos(4 * CGAL_PI / 3), radius * sin(4 * CGAL_PI / 3)}, // 240°
      {radius * cos(5 * CGAL_PI / 3), radius * sin(5 * CGAL_PI / 3)}, // 300°
  }};
  std::array<X_monotone_curve, 6> hexagon;
  for(size_t i = 0; i < hexagon_points.size(); ++i) {
    size_t next_i = (i + 1) % hexagon_points.size();
    hexagon[i] = cst_x_curve(hexagon_points[i], hexagon_points[next_i]);
  }
  // rect hole
  auto hole_rectangle = {cst_x_curve({-2, -2}, {2, -2}), cst_x_curve({2, -2}, {2, 2}), cst_x_curve({2, 2}, {-2, 2}),
                         cst_x_curve({-2, 2}, {-2, -2})};
  // iso vertex inside rect hole
  auto iso_vertex_inside_hole = Point{0.5, 0.5};
  // degenerate segment below the rect hole
  auto degenerate_segment = cst_x_curve({0, -3}, {1, -3});

  CGAL::insert_point(arr, iso_vertex_inside_hole);
  CGAL::insert(arr, hexagon.begin(), hexagon.end());
  CGAL::insert(arr, hole_rectangle.begin(), hole_rectangle.end());
  CGAL::insert(arr, degenerate_segment);
  CGAL::draw(arr);
}