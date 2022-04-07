// Testing the do_intersect function

#include <CGAL/Quotient.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>

#include <list>
#include <iostream>

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
  bool            expected_intersect[N_SEGMENTS];
  int             k;

  segs[0] = Segment_2 (Point_2 (-2, -2), Point_2 (-1, -1));
  segs[1] = Segment_2 (Point_2 (-1, 1), Point_2 (0, 1));
  segs[2] = Segment_2 (Point_2 (-1, 0), Point_2 (0, 0));

  expected_intersect[0] = false;
  expected_intersect[1] = true;
  expected_intersect[2] = true;

  insert(arr, Segment_2(Point_2(0, 0), Point_2(2, 0)));
  insert(arr, Segment_2(Point_2(2, 0), Point_2(2, 2)));
  insert(arr, Segment_2(Point_2(2, 2), Point_2(0, 2)));
  insert(arr, Segment_2(Point_2(0, 2), Point_2(0, 0)));

  for (k = 0; k < N_SEGMENTS; k++)
  {
    bool do_inter = do_intersect(arr, segs[k]);

    std::cout << "Segment: " << segs[k];
    std::cout << "        Expected: " << expected_intersect[k];
    std::cout << "        Actual: " << do_inter << std::endl;

    if (expected_intersect[k] != do_inter)
      return (1);
  }

  return (0);
}
