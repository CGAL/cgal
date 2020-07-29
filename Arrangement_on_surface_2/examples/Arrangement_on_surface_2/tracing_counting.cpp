//! \file examples/Arrangement_on_surface_2/tracing_counting.cpp
// Trace all traits operations and count them.

#include <CGAL/Cartesian.h>
#include <CGAL/Quotient.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_tracing_traits_2.h>
#include <CGAL/Arr_counting_traits_2.h>
#include <list>

typedef CGAL::Quotient<int>                             Number_type;
typedef CGAL::Cartesian<Number_type>                    Kernel;
typedef CGAL::Arr_segment_traits_2<Kernel>              Segment_traits_2;
typedef CGAL::Arr_tracing_traits_2<Segment_traits_2>    Tracing_traits_2;
typedef CGAL::Arr_counting_traits_2<Tracing_traits_2>   Geom_traits_2;
typedef Geom_traits_2::Point_2                          Point_2;
typedef Geom_traits_2::X_monotone_curve_2               Segment_2;
typedef CGAL::Arrangement_2<Geom_traits_2>              Arrangement_2;

int main ()
{
  const Segment_2 s1(Point_2(0, 0), Point_2(2, 2));
  const Segment_2 s2(Point_2(2, 0), Point_2(0, 2));
  std::list<Segment_2> segments;
  segments.push_back(s1);
  segments.push_back(s2);

  Geom_traits_2 traits;
  traits.disable_all_traces();
  traits.enable_trace(Tracing_traits_2::INTERSECT_OP);

  // Construct an arrangement using aggregated insertion:
  Arrangement_2 arr1(&traits);
  insert(arr1, segments.begin(), segments.end());
  std::cout << traits << std::endl;
  traits.clear_counters();

  // Construct an arrangement using incremental insertion:
  Arrangement_2 arr2(&traits);
  insert(arr2, s1);
  insert(arr2, s2);
  std::cout << traits << std::endl;

  return 0;
}
