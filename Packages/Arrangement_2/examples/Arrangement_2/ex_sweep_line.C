//! \file examples/Arrangement_2/ex_sweep_line.C
// Computing intersection points among curves using the sweep line.

#include <CGAL/Cartesian.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>
#include <CGAL/Sweep_line_2_algorithms.h>
#include <list>

typedef CGAL::Quotient<CGAL::MP_Float>                  NT;
typedef CGAL::Cartesian<NT>                             Kernel;
typedef Kernel::Point_2                                 Point_2;
typedef Kernel::Segment_2                               Segment_2;

int main()
{
  // Construct the input segments.
  Segment_2 segments[] = {Segment_2 (Point_2 (1, 5), Point_2 (8, 5)),
                          Segment_2 (Point_2 (1, 1), Point_2 (8, 8)),
                          Segment_2 (Point_2 (3, 1), Point_2 (3, 8)),
                          Segment_2 (Point_2 (8, 5), Point_2 (8, 8))};
  
  // Compute all intersection points.
  std::list<Point_2>     pts;

  CGAL::get_intersection_points (segments, segments + 4,
                                 std::back_inserter (pts));
  
  // Print the result.
  std::copy (pts.begin(), pts.end(),
             std::ostream_iterator<Point_2>(std::cout, "\n"));

  return (0);
}
