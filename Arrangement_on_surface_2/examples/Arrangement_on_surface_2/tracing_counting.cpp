//! \file examples/Arrangement_on_surface_2/tracing_counting.cpp
// Trace all traits operations and count them.

#include <CGAL/Cartesian.h>
#include <CGAL/Quotient.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_tracing_traits_2.h>
#include <CGAL/Arr_counting_traits_2.h>
#include <list>

typedef CGAL::Quotient<int>                           Number_type;
typedef CGAL::Cartesian<Number_type>                  Kernel;
typedef CGAL::Arr_segment_traits_2<Kernel>            Segment_traits_2;
typedef CGAL::Arr_tracing_traits_2<Segment_traits_2>  Tracing_traits_2;
typedef CGAL::Arr_counting_traits_2<Tracing_traits_2> Geom_traits_2;
typedef Geom_traits_2::Point_2                        Point_2;
typedef Geom_traits_2::X_monotone_curve_2             Segment_2;
typedef CGAL::Arrangement_2<Geom_traits_2>            Arrangement_2;

int main ()
{
  // Construct the arrangement of five intersecting segments.
  Geom_traits_2 traits;
  traits.disable_all_traces();
  traits.enable_trace(Tracing_traits_2::INTERSECT_OP);
  Arrangement_2           arr(&traits);;
  std::list<Segment_2>    segments;

  segments.push_back (Segment_2 (Point_2(1, 0), Point_2(2, 4)));
  segments.push_back (Segment_2 (Point_2(5, 0), Point_2(5, 5)));
  segments.push_back (Segment_2 (Point_2(1, 0), Point_2(5, 3)));  
  segments.push_back (Segment_2 (Point_2(0, 2), Point_2(6, 0)));
  segments.push_back (Segment_2 (Point_2(3, 0), Point_2(5, 5)));

  insert (arr, segments.begin(), segments.end());
  std::cout << traits << std::endl;

  return 0;
}
