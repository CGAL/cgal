// Testing the compute zone function

#include <list>
#include <iostream>

#include <CGAL/Quotient.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>

typedef CGAL::Quotient<int>                           Number_type;
typedef CGAL::Simple_cartesian<Number_type>           Kernel;
typedef CGAL::Arr_segment_traits_2<Kernel>            Traits_2;
typedef Traits_2::Point_2                             Point_2;
typedef Traits_2::X_monotone_curve_2                  Segment_2;
typedef CGAL::Arrangement_2<Traits_2>                 Arrangement_2;
typedef Arrangement_2::Halfedge_handle                Halfedge_handle;

#define N_SEGMENTS 3

int main ()
{
  Arrangement_2   arr;
  Segment_2       segs[N_SEGMENTS];
  std::size_t     zone_expected_comp[N_SEGMENTS];
  int             k;

  segs[0] = Segment_2 (Point_2 (-2, -2), Point_2 (-1, -1));
  segs[1] = Segment_2 (Point_2 (-1, 1), Point_2 (1, 1));
  segs[2] = Segment_2 (Point_2 (0, 0), Point_2 (3, 0));

  zone_expected_comp[0] = 1;
  zone_expected_comp[1] = 3;
  zone_expected_comp[2] = 4;

  insert(arr, Segment_2(Point_2(0, 0), Point_2(2, 0)));
  insert(arr, Segment_2(Point_2(2, 0), Point_2(2, 2)));
  insert(arr, Segment_2(Point_2(2, 2), Point_2(0, 2)));
  insert(arr, Segment_2(Point_2(0, 2), Point_2(0, 0)));

  for (k = 0; k < N_SEGMENTS; k++)
  {
    std::list<CGAL::Object> zone_elems;
    zone(arr, segs[k], std::back_inserter(zone_elems));
    std::size_t zone_actual_comp = zone_elems.size();

    std::cout << "Segment: " << segs[k];
    std::cout << "        Expected: " << zone_expected_comp[k];
    std::cout << "        Actual: " << zone_actual_comp << std::endl;

    if (zone_expected_comp[k] != zone_actual_comp)
      return (1);
  }

  return (0);
}
