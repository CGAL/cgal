//! \file examples/Arrangement_2/ex_sweep_line.C
// Computing intersection points among curves using the sweep line.

#include <CGAL/Cartesian.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Sweep_line_2_algorithms.h>
#include <list>

typedef CGAL::Quotient<CGAL::MP_Float>                  NT;
typedef CGAL::Cartesian<NT>                             Kernel;
typedef CGAL::Arr_segment_traits_2<Kernel>              Traits_2;
typedef Traits_2::Point_2                               Point_2;
typedef Traits_2::Curve_2                               Segment_2;

int main()
{
  // Construct the input segments.
  std::list<Segment_2>   segments;
  
  segments.push_back (Segment_2 (Point_2 (1, 5), Point_2 (8, 5)));
  segments.push_back (Segment_2 (Point_2 (1, 1), Point_2 (8, 8)));
  segments.push_back (Segment_2 (Point_2 (3, 1), Point_2 (3, 8)));
  segments.push_back (Segment_2 (Point_2 (8, 5), Point_2 (8, 8)));

  // Compute all intersection points.
  Traits_2               traits;
  std::list<Point_2>     pts;

  CGAL::get_intersection_points (segments.begin(), segments.end(),
                                 std::back_inserter (pts), traits);
  
  // Print the result.
  std::copy (pts.begin(), pts.end(), std::ostream_iterator<Point_2>(std::cout, "\n"));

  return (0);
}
