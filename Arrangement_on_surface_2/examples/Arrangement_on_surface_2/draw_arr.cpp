#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/draw_arrangement_2.h>

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using Traits = CGAL::Arr_segment_traits_2<Kernel>;
using Point = Traits::Point_2;
using X_monotone_curve = Traits::X_monotone_curve_2;
using Arrangement_2 = CGAL::Arrangement_2<Traits>;

int main() {
  Traits traits;
  Arrangement_2 arr(&traits);
  auto ctr_xcv = traits.construct_x_monotone_curve_2_object();

  CGAL::insert(arr, ctr_xcv(Point(-2,-2), Point(2,-2)));
  CGAL::insert(arr, ctr_xcv(Point(2,-2), Point(2,2)));
  CGAL::insert(arr, ctr_xcv(Point(2,2), Point(-2,2)));
  CGAL::insert(arr, ctr_xcv(Point(-2,2), Point(-2,-2)));

  CGAL::insert(arr, ctr_xcv(Point(-1,-1), Point(1,-1)));
  CGAL::insert(arr, ctr_xcv(Point(1,-1), Point(1,1)));
  CGAL::insert(arr, ctr_xcv(Point(1,1), Point(-1,1)));
  CGAL::insert(arr, ctr_xcv(Point(-1,1), Point(-1,-1)));

  CGAL::insert(arr, ctr_xcv(Point(-2,-2), Point(-2,-4)));
  CGAL::insert(arr, ctr_xcv(Point(2,-2), Point(4,-2)));

  CGAL::insert(arr, ctr_xcv(Point(0,0), Point(0,-3)));

  CGAL::draw(arr);

  return EXIT_SUCCESS;
}
