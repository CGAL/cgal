#include <vector>

#include <CGAL/draw_arrangement_2.h>

#include "arr_geodesic.h"
#include "arr_print.h"

void draw_face_crossing_boundary() {
  Arrangement arr;
  const auto& traits = *arr.geometry_traits();
  auto ctr_p = traits.construct_point_2_object();
  auto ctr_cv = traits.construct_curve_2_object();

  auto p1 = ctr_p(-0.95, 0.32, 0), p2 = ctr_p(-0.87, 0.02, 0.49), p3 = ctr_p(-0.93, -0.36, 0),
       p4 = ctr_p(-0.81, -0.03, -0.59);
  auto arcs = {ctr_cv(p1, p2), ctr_cv(p2, p3), ctr_cv(p3, p4), ctr_cv(p4, p1)};
  CGAL::insert(arr, arcs.begin(), arcs.end());

  print_arrangement_size(arr);
  CGAL::draw(arr, "face crossing boundary");
}

void draw_lakes() {
  Arrangement arr;
  const auto& traits = *arr.geometry_traits();
  auto ctr_p = traits.construct_point_2_object();
  auto ctr_cv = traits.construct_curve_2_object();

  auto poly1 = {
      ctr_p(-0.27, -0.053, -0.96), ctr_p(-0.76, -0.15, -0.63), ctr_p(-0.98, -0.19, -0.063), ctr_p(-0.98, -0.098, 0.2),
      ctr_p(-0.44, -0.18, 0.88),   ctr_p(0.39, -0.0049, 0.92), ctr_p(-0.01, 0.39, 0.92),    ctr_p(-0.54, 0.66, 0.53),
      ctr_p(-0.83, 0.56, 0.025),   ctr_p(-0.57, 0.32, -0.75),  ctr_p(-0.087, 0.048, -1),    ctr_p(-0.048, 0.088, -1),
      ctr_p(0.12, -0.14, -0.98),   ctr_p(-0.12, -0.14, -0.98),
  };
  auto poly2 = {ctr_p(-0.24, -0.53, -0.81), ctr_p(-0.47, -0.54, -0.69), ctr_p(-0.68, -0.65, -0.32),
                ctr_p(-0.71, -0.68, 0.2),   ctr_p(-0.54, -0.52, 0.67),  ctr_p(-0.18, -0.72, 0.67),
                ctr_p(0.31, -0.68, 0.67),   ctr_p(0.71, -0.69, 0.11),   ctr_p(0.6, -0.58, -0.56),
                ctr_p(0.21, -0.62, -0.75)};
  auto poly3 = {ctr_p(0.44, 0.27, -0.86), ctr_p(0.58, -0.063, -0.81), ctr_p(0.87, -0.094, -0.48),
                ctr_p(0.97, -0.1, 0.2),   ctr_p(0.46, 0.77, 0.45),    ctr_p(-0.023, 0.89, 0.45),
                ctr_p(-0.3, 0.95, 0.11),  ctr_p(-0.22, 0.69, -0.69),  ctr_p(-0.076, 0.35, -0.93)};
  auto poly4 = {
      ctr_p(0.4, 0.67, -0.63),  ctr_p(0.78, 0.39, -0.48),  ctr_p(0.92, 0.35, -0.15),
      ctr_p(0.52, 0.86, 0.025), ctr_p(0.068, 0.99, -0.15), ctr_p(0.22, 0.85, -0.48),
  };
  std::vector<Curve> arcs;
  std::vector<std::vector<Point>> polygons{poly1, poly2, poly3, poly4};
  for(const auto& poly : polygons) {
    for(size_t i = 0; i < poly.size(); ++i) {
      size_t next = (i + 1) % poly.size();
      arcs.push_back(ctr_cv(poly[i], poly[next]));
    }
  }

  CGAL::insert(arr, arcs.begin(), arcs.end());
  print_arrangement_size(arr);
  CGAL::draw(arr, "lakes");
}

void draw_guassian_map() {
  Arrangement arr;
  const auto& traits = *arr.geometry_traits();
  auto ctr_p = traits.construct_point_2_object();
  auto ctr_cv = traits.construct_curve_2_object();

  auto p1 = ctr_p(1, 1, 1), p2 = ctr_p(-1, -1, 1), p3 = ctr_p(-1, 1, -1), p4 = ctr_p(1, -1, -1);
  auto arcs = {ctr_cv(p1, p2), ctr_cv(p2, p3), ctr_cv(p3, p1), ctr_cv(p1, p4), ctr_cv(p2, p4), ctr_cv(p3, p4)};
  CGAL::insert(arr, arcs.begin(), arcs.end());
  print_arrangement_size(arr);
  CGAL::draw(arr, "guassian map of a tetrahedron");
}

void draw_random_arcs(int n) {
  Arrangement arr;
  const auto& traits = *arr.geometry_traits();
  auto ctr_p = traits.construct_point_2_object();
  auto ctr_cv = traits.construct_curve_2_object();

  CGAL::Random random;
  std::vector<Point> points;
  for(int i = 0; i < n; ++i) {
    double x = random.get_double(-1.0, 1.0);
    double y = random.get_double(-1.0, 1.0);
    double z = random.get_double(-1.0, 1.0);
    points.push_back(ctr_p(x, y, z));
  }
  std::vector<Curve> curves;
  for(int i = 0; i < n; ++i) {
    int j = random.get_int(0, n - 1);
    if(i == j) j = (j + 1) % n;
    curves.push_back(ctr_cv(points[i], points[j]));
  }

  CGAL::insert(arr, curves.begin(), curves.end());
  print_arrangement_size(arr);
  CGAL::draw(arr, (std::to_string(n) + " random arcs").c_str());
}

int main() {
  draw_face_crossing_boundary();
  draw_lakes();
  draw_guassian_map();
  draw_random_arcs(100);
  return 0;
}
