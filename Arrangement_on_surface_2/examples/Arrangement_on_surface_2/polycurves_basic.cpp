//! \file examples/Arrangement_on_surface_2/polylines.cpp
// Constructing an arrangement of polylines.

#include <vector>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Arr_directional_non_caching_segment_basic_traits_2.h>
#include <CGAL/Arr_polycurve_basic_traits_2.h>
#include <CGAL/Arrangement_2.h>

#include "arr_print.h"

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef CGAL::Arr_directional_non_caching_segment_basic_traits_2<Kernel>
                                                            Subcurve_traits;
typedef CGAL::Arr_polycurve_basic_traits_2<Subcurve_traits> Geom_traits;
typedef Geom_traits::Point_2                                Point;
typedef Subcurve_traits::X_monotone_curve_2                 X_monotone_subcurve;
typedef Geom_traits::X_monotone_curve_2                     X_monotone_curve;
typedef CGAL::Arrangement_2<Geom_traits>                    Arrangement;

int main() {
  Geom_traits traits;
  Arrangement arr(&traits);

  auto ctr = traits.construct_x_monotone_curve_2_object();

  std::vector<X_monotone_subcurve> segs1;
  segs1.push_back(X_monotone_subcurve(Point(0, 0), Point(1, 1)));
  segs1.push_back(X_monotone_subcurve(Point(1, 1), Point(2, 2)));
  segs1.push_back(X_monotone_subcurve(Point(2, 2), Point(3, 1)));
  segs1.push_back(X_monotone_subcurve(Point(3, 1), Point(4, 0)));
  X_monotone_curve pc1 = ctr(segs1.begin(), segs1.end());

  std::vector<X_monotone_subcurve> segs2;
  segs2.push_back(X_monotone_subcurve(Point(0, 0), Point(1, 1)));
  segs2.push_back(X_monotone_subcurve(Point(1, 1), Point(2, 2)));
  segs2.push_back(X_monotone_subcurve(Point(2, 2), Point(3, 1)));
  segs2.push_back(X_monotone_subcurve(Point(3, 1), Point(4, 0)));
  X_monotone_curve pc2 = ctr(segs2.begin(), segs2.end());

  insert_non_intersecting_curve(arr, pc1);
  insert_non_intersecting_curve(arr, pc2);

  print_arrangement_size(arr);          // print the arrangement size

  return 0;
}
