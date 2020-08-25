#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
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
  std::vector<Point_2> points
    = { Point_2 (0, 0),
        Point_2 (1, 0),
        Point_2 (2, 1),
        Point_2 (1, 1),
        Point_2 (0, 0.5) };

  auto point_to_segment
    = [&](const Point_2& p) -> Segment_2
  {
    const Point_2* it = &p;
    if (it == &points.back())
      return Segment_2 (p, points[0]);
    return Segment_2 (p, *(++ it));
  };

  // Case 1: Segment-based GPS from Segment_2 range with transform
  // iterator
  CGAL::General_polygon_2<Segment_traits>
    gps_1 (boost::make_transform_iterator
           (points.begin(), point_to_segment),
           boost::make_transform_iterator
           (points.end(), point_to_segment));

  // Case 2: Segment-based GPS from
  // Arr_segment_traits_2::X_monotone_curve_2 range with internal
  // converter
  CGAL::General_polygon_2<Segment_traits>
    gps_2 (CGAL::points_to_x_monotone_curves_begin<Segment_traits>(points),
           CGAL::points_to_x_monotone_curves_end<Segment_traits>(points));

  // Case 3: Segment-based GPS from Point range
  CGAL::General_polygon_2<Segment_traits>
    gps_3 (points.cbegin(), points.cend());

  // Case 4: Polyline-based GPS from
  // Arr_polyline_traits_2::X_monotone_curve_2 range with internal
  // converter
  CGAL::General_polygon_2<Polyline_traits>
    gps_4 (CGAL::points_to_x_monotone_curves_begin<Polyline_traits>(points),
           CGAL::points_to_x_monotone_curves_end<Polyline_traits>(points));

  // Case 5: Polyline-based GPS from Point range with internal
  // converter
  CGAL::General_polygon_2<Polyline_traits>
    gps_5 (points.cbegin(), points.cend());


  return EXIT_SUCCESS;
}
