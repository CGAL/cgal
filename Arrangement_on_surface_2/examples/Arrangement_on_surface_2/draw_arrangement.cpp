#include <CGAL/Cartesian.h>
#include <CGAL/Quotient.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_naive_point_location.h>
#include <CGAL/draw_arrangement_2.h>

typedef CGAL::Quotient<CGAL::MP_Float>                Number_type;
typedef CGAL::Cartesian<Number_type>                  Kernel;
typedef CGAL::Arr_segment_traits_2<Kernel>            Traits_2;
typedef Traits_2::Point_2                             Point_2;
typedef Traits_2::X_monotone_curve_2                  Segment_2;
typedef CGAL::Arrangement_2<Traits_2>                 Arrangement_2;

int main ()
{
  // Construct the arrangement of five intersecting segments.
  Arrangement_2 arr;
  Segment_2     S1 [5];

  S1[0]=Segment_2(Point_2 (1, 2.5), Point_2 (4, 5));
  S1[1]=Segment_2(Point_2 (1, 2.5), Point_2 (6, 2.5));
  S1[2]=Segment_2(Point_2 (1, 2.5), Point_2 (4, 0));  
  S1[3]=Segment_2(Point_2 (4, 5), Point_2 (6, 2.5));
  S1[4]=Segment_2(Point_2 (4, 0), Point_2 (6, 2.5));

  insert_non_intersecting_curves(arr, S1, S1 + 5);

  // Perform an incremental insertion of a single overlapping segment.
  insert(arr, Segment_2 (Point_2 (0, 2.5), Point_2 (4, 2.5)));

  // Aggregately insert an additional set of five segments.
  Segment_2     S2 [5];
  S2[0]=Segment_2(Point_2 (0, 4), Point_2 (6, 5));
  S2[1]=Segment_2(Point_2 (0, 3), Point_2 (6, 4));
  S2[2]=Segment_2(Point_2 (0, 2), Point_2 (6, 1));
  S2[3]=Segment_2(Point_2 (0, 1), Point_2 (6, 0));
  S2[4]=Segment_2(Point_2 (6, 1), Point_2 (6, 4));

  insert(arr, S2, S2 + 5);

  // Draw the arrangement.
  CGAL::draw(arr);
  
  return EXIT_SUCCESS;
}
