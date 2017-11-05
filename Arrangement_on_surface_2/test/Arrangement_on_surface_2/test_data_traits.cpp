#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arr_curve_data_traits_2.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel    Kernel;
typedef Kernel::Point_2                                      Point_2;
typedef Kernel::Segment_2                                    Segment_2;
typedef CGAL::Arr_segment_traits_2<Kernel>                   Segment_traits_2;
typedef CGAL::Arr_curve_data_traits_2<Segment_traits_2, int> Traits_2;
typedef Traits_2::X_monotone_curve_2                         X_monotone_curve_2;

int main()
{
  Point_2 p1 = Point_2(0.0, 0.0);
  Point_2 p2 = Point_2(1.0, 0.0);
  Point_2 p3 = Point_2(2.0, 0.0);
  Segment_2 s1 = Segment_2(p1, p2);
  Segment_2 s2 = Segment_2(p2, p3);

  Traits_2 traits;
  X_monotone_curve_2 c1 = X_monotone_curve_2(s1, 1);
  X_monotone_curve_2 fc1 = traits.construct_opposite_2_object()(c1);

  if (c1.data() !=  fc1.data()) return 1;

  X_monotone_curve_2 c2 = X_monotone_curve_2(s2, 1);
  X_monotone_curve_2 mc;
  if (!traits.are_mergeable_2_object()(c1, c2)) return 1;
  traits.merge_2_object()(c1, c2, mc);
  if (c1.data() !=  mc.data()) return 1;

  X_monotone_curve_2 sc1, sc2;
  traits.split_2_object()(mc, p2, sc1, sc2);
  if (mc.data() !=  sc1.data()) return 1;
  if (sc1.data() !=  sc2.data()) return 1;

  return 0;
}
