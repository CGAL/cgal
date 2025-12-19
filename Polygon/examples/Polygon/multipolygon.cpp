#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Multipolygon_with_holes_2.h>

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_2 = Kernel::Point_2;
using Polygon_2 = CGAL::Polygon_2<Kernel>;
using Polygon_with_holes_2 = CGAL::Polygon_with_holes_2<Kernel>;
using Multipolygon_with_holes_2 = CGAL::Multipolygon_with_holes_2<Kernel>;

int main() {

  Point_2 p1_outer[] = {Point_2(0,0), Point_2(1,0), Point_2(1,1), Point_2(0,1)};
  Point_2 p1_inner[] = {Point_2(0.2,0.2), Point_2(0.8,0.2), Point_2(0.8,0.8), Point_2(0.2,0.8)};

  Polygon_with_holes_2 p1(Polygon_2(p1_outer, p1_outer+4));
  Polygon_2 h(p1_inner, p1_inner+4);
  p1.add_hole(h);

  Point_2 p2_outer[] = {Point_2(0.4,0.4), Point_2(0.6,0.4), Point_2(0.6,0.6), Point_2(0.4,0.6)};
  Polygon_with_holes_2 p2(Polygon_2(p2_outer, p2_outer+4));

  Multipolygon_with_holes_2 mp;
  mp.add_polygon_with_holes(p1);
  mp.add_polygon_with_holes(p2);

  for (auto const& p: mp.polygons_with_holes()) {
    std::cout << p << std::endl;
  }

  return 0;
}
