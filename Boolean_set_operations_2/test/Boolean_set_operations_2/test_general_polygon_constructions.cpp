#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Boolean_set_operations_2/Polygon_conversions.h>
#include <CGAL/General_polygon_2.h>
#include <CGAL/Arr_polyline_traits_2.h>
#include <CGAL/Arr_segment_traits_2.h>

using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;
using Point_2 = Kernel::Point_2;
using Segment_2 = Kernel::Segment_2;
using Segment_traits = CGAL::Arr_segment_traits_2<Kernel>;
using Polyline_traits = CGAL::Arr_polyline_traits_2<Segment_traits>;

int main()
{
  std::vector<Point_2> points = {
    Point_2(0, 0),
    Point_2(1, 0),
    Point_2(2, 1),
    Point_2(1, 1),
    Point_2(0, 0.5)
  };

  auto point_to_segment = [&](const Point_2& p) -> Segment_2 {
    const Point_2* it = &p;
    if (it == &points.back()) return Segment_2(p, points[0]);
    return Segment_2(p, *(++it));
  };

  // Case 1: Segment-based GPS from Segment_2 range with transform iterator
  CGAL::General_polygon_2<Segment_traits>gp1
    (boost::make_transform_iterator(points.begin(), point_to_segment),
     boost::make_transform_iterator(points.end(), point_to_segment));

  // Case 5: Polyline-based GPS from Point range with internal
  // converter
  Polyline_traits::Curve_2 curve(boost::make_transform_iterator
                                 (points.begin(), point_to_segment),
                                 boost::make_transform_iterator
                                 (points.end(), point_to_segment));
  CGAL::General_polygon_2<Polyline_traits>gp2(curve);

  Polyline_traits traits;
  CGAL::General_polygon_2<Polyline_traits>gp3(traits.make_curve_2(points));

  CGAL::General_polygon_2<Polyline_traits>gp4(traits.make_curve_2(points, true));

  auto ctr = traits.construct_curve_2_object();
  auto gp5 = ctr(points.begin(), points.end());

  CGAL::Polygon_2<Kernel> gp6(points.begin(), points.end());
  auto gp7 = test_conversion(gp6);
  if (gp7 != gp6) {
    std::cerr << "Conversion failed\n";
    return 1;
  }

  return EXIT_SUCCESS;
}
