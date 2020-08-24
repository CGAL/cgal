#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/General_polygon_2.h>
#include <CGAL/Arr_polyline_traits_2.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Boolean_set_operations_2/Points_to_x_monotone_curves_iterator.h>

using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;
using Point_2 = Kernel::Point_2;
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

  CGAL::General_polygon_2<Polyline_traits>
    gps (CGAL::points_to_x_monotone_curves_begin<Polyline_traits>(points),
         CGAL::points_to_x_monotone_curves_end<Polyline_traits>(points));


  return EXIT_SUCCESS;
}
