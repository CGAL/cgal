// Testing the do_intersect function

#include <CGAL/Quotient.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>

#include <list>
#include <iostream>

using Number_type = CGAL::Quotient<int>;
using Kernel = CGAL::Simple_cartesian<Number_type>;
using Traits_2 = CGAL::Arr_segment_traits_2<Kernel>;
using Point_2 = Traits_2::Point_2;
using Segment_2 = Traits_2::X_monotone_curve_2;
using Arrangement_2 = CGAL::Arrangement_2<Traits_2>;
using Halfedge_handle = Arrangement_2::Halfedge_handle;

int main () {
  Arrangement_2 arr;
  using Tt = Arrangement_2::Topology_traits;
  Tt::Default_point_location_strategy def_pl(arr);

  Segment_2 segs[] = {
    Segment_2(Point_2(-2, -2), Point_2(-1, -1)),
    Segment_2(Point_2(-1, 1), Point_2(0, 1)),
    Segment_2(Point_2(-1, 0), Point_2(0, 0))
  };

  bool expected_intersect[] = {false, true, true};

  insert(arr, Segment_2(Point_2(0, 0), Point_2(2, 0)));
  insert(arr, Segment_2(Point_2(2, 0), Point_2(2, 2)));
  insert(arr, Segment_2(Point_2(2, 2), Point_2(0, 2)));
  insert(arr, Segment_2(Point_2(0, 2), Point_2(0, 0)));

  size_t k = 0;
  for (const auto& seg : segs) {
    bool do_inter_0 = do_intersect(arr, seg);
    bool do_inter_1 = do_intersect(arr, seg, def_pl, std::true_type());
    bool do_inter_2 = do_intersect(arr, seg, def_pl, std::false_type());

    std::cout << "Segment: " << segs[k] << std::endl;
    std::cout << "  Expected: " << expected_intersect[k] << std::endl;
    std::cout << "  Actual auto: " << do_inter_0 << std::endl;
    std::cout << "  Actual x-monotone curve: " << do_inter_1 << std::endl;
    std::cout << "  Actual curve: " << do_inter_2 << std::endl;

    if (expected_intersect[k] != do_inter_0) return -1;
    if (expected_intersect[k] != do_inter_1) return -1;
    if (expected_intersect[k] != do_inter_2) return -1;
    ++k;
  }

  return 0;
}
