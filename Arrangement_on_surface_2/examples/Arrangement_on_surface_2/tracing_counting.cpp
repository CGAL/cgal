//! \file examples/Arrangement_on_surface_2/tracing_counting.cpp
// Trace all traits operations and count them.

#include <list>

#include <CGAL/basic.h>
#include <CGAL/Arr_tracing_traits_2.h>
#include <CGAL/Arr_counting_traits_2.h>

#include "arr_exact_construction_segments.h"

typedef CGAL::Arr_tracing_traits_2<Traits>              Tracing_traits;
typedef CGAL::Arr_counting_traits_2<Tracing_traits>     Geom_traits;
typedef CGAL::Arrangement_2<Geom_traits>                My_arrangement;

int main() {
  const Segment s1(Point(0, 0), Point(2, 2));
  const Segment s2(Point(2, 0), Point(0, 2));
  std::list<Segment> segments;
  segments.push_back(s1);
  segments.push_back(s2);

  Geom_traits traits;
  traits.disable_all_traces();
  traits.enable_trace(Tracing_traits::INTERSECT_OP);

  // Construct an arrangement using aggregated insertion:
  My_arrangement arr1(&traits);
  insert(arr1, segments.begin(), segments.end());
  std::cout << traits << std::endl;
  traits.clear_counters();

  // Construct an arrangement using incremental insertion:
  My_arrangement arr2(&traits);
  insert(arr2, s1);
  insert(arr2, s2);
  std::cout << traits << std::endl;

  return 0;
}
