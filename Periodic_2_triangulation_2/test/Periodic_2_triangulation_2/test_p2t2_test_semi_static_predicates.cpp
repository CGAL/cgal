// Author(s)     : Nico Kruithof  <Nico@nghk.nl>

#include "./types.h"
#include <map>

const int N = 10;

void test_orientation()
{
  Gt traits;
  traits.set_domain(Gt::Iso_rectangle_2(0, 0, 1, 1));

  Gt::Offset_2 o0(0, 0);
  Gt::Offset_2 o1(0, 1);

  /// Near degenerate points, which cause the predicate to fail if not filtered
  Point p0(0.5 + (0.4999 / N) * 2, 0.5 + (0.4999 / N) * -5);
  Point p1(0.5 + (0.4999 / N) * 4, 0.5 + (0.4999 / N) * -4);
  Point p2(0.5 + (0.4999 / N) * 6, 0.5 + (0.4999 / N) * -3);

  CGAL_assertion(traits.orientation_2_object()(p0, p1, p2) == 1);
  CGAL_assertion(traits.orientation_2_object()(p2, p0, p1) == 1);
  CGAL_assertion(traits.orientation_2_object()(p1, p2, p0) == 1);
  CGAL_assertion(traits.orientation_2_object()(p0, p2, p1) == -1);
  CGAL_assertion(traits.orientation_2_object()(p1, p0, p2) == -1);
  CGAL_assertion(traits.orientation_2_object()(p2, p1, p0) == -1);

  CGAL_assertion(traits.orientation_2_object()(p0, p1, p2) ==
                 traits.orientation_2_object()(p0, p1, p2, o0, o0, o0));
  CGAL_assertion(traits.orientation_2_object()(p2, p0, p1) ==
                 traits.orientation_2_object()(p2, p0, p1, o0, o0, o0));
  CGAL_assertion(traits.orientation_2_object()(p1, p2, p0) ==
                 traits.orientation_2_object()(p1, p2, p0, o0, o0, o0));
  CGAL_assertion(traits.orientation_2_object()(p0, p2, p1) ==
                 traits.orientation_2_object()(p0, p2, p1, o0, o0, o0));
  CGAL_assertion(traits.orientation_2_object()(p1, p0, p2) ==
                 traits.orientation_2_object()(p1, p0, p2, o0, o0, o0));
  CGAL_assertion(traits.orientation_2_object()(p2, p1, p0) ==
                 traits.orientation_2_object()(p2, p1, p0, o0, o0, o0));

  CGAL_assertion(traits.orientation_2_object()(p0, p1, p2) ==
                 traits.orientation_2_object()(p0, p1, p2, o1, o1, o1));
  CGAL_assertion(traits.orientation_2_object()(p2, p0, p1) ==
                 traits.orientation_2_object()(p2, p0, p1, o1, o1, o1));
  CGAL_assertion(traits.orientation_2_object()(p1, p2, p0) ==
                 traits.orientation_2_object()(p1, p2, p0, o1, o1, o1));
  CGAL_assertion(traits.orientation_2_object()(p0, p2, p1) ==
                 traits.orientation_2_object()(p0, p2, p1, o1, o1, o1));
  CGAL_assertion(traits.orientation_2_object()(p1, p0, p2) ==
                 traits.orientation_2_object()(p1, p0, p2, o1, o1, o1));
  CGAL_assertion(traits.orientation_2_object()(p2, p1, p0) ==
                 traits.orientation_2_object()(p2, p1, p0, o1, o1, o1));
}

void test_in_circle()
{
  Gt traits;
  traits.set_domain(Gt::Iso_rectangle_2(0, 0, 1, 1));

  Gt::Offset_2 o0(0, 0);
  Gt::Offset_2 o1(0, 1);

  /// Near degenerate points, which cause the predicate to fail if not filtered
  /// On the circle with center (0.4999, 0.4999) and radius 5
  Point p0( 5 - 0.4999,   -0.4999);
  Point p1( 3 - 0.4999,  4 - 0.4999);
  Point p2(-4 - 0.4999,  3 - 0.4999);
  Point p3( 4 - 0.4999, -3 - 0.4999);

  CGAL_assertion(traits.side_of_oriented_circle_2_object()(p0, p1, p2, p3) == 1);
  CGAL_assertion(traits.side_of_oriented_circle_2_object()(p2, p0, p1, p3) == 1);
  CGAL_assertion(traits.side_of_oriented_circle_2_object()(p1, p2, p0, p3) == 1);
  CGAL_assertion(traits.side_of_oriented_circle_2_object()(p0, p2, p1, p3) == -1);
  CGAL_assertion(traits.side_of_oriented_circle_2_object()(p1, p0, p2, p3) == -1);
  CGAL_assertion(traits.side_of_oriented_circle_2_object()(p2, p1, p0, p3) == -1);

  CGAL_assertion(traits.side_of_oriented_circle_2_object()(p0, p1, p2, p3) ==
                 traits.side_of_oriented_circle_2_object()(p0, p1, p2, p3, o0, o0, o0, o0));
  CGAL_assertion(traits.side_of_oriented_circle_2_object()(p2, p0, p1, p3) ==
                 traits.side_of_oriented_circle_2_object()(p2, p0, p1, p3, o0, o0, o0, o0));
  CGAL_assertion(traits.side_of_oriented_circle_2_object()(p1, p2, p0, p3) ==
                 traits.side_of_oriented_circle_2_object()(p1, p2, p0, p3, o0, o0, o0, o0));
  CGAL_assertion(traits.side_of_oriented_circle_2_object()(p0, p2, p1, p3) ==
                 traits.side_of_oriented_circle_2_object()(p0, p2, p1, p3, o0, o0, o0, o0));
  CGAL_assertion(traits.side_of_oriented_circle_2_object()(p1, p0, p2, p3) ==
                 traits.side_of_oriented_circle_2_object()(p1, p0, p2, p3, o0, o0, o0, o0));
  CGAL_assertion(traits.side_of_oriented_circle_2_object()(p2, p1, p0, p3) ==
                 traits.side_of_oriented_circle_2_object()(p2, p1, p0, p3, o0, o0, o0, o0));

  CGAL_assertion(traits.side_of_oriented_circle_2_object()(p0, p1, p2, p3) ==
                 traits.side_of_oriented_circle_2_object()(p0, p1, p2, p3, o1, o1, o1, o1));
  CGAL_assertion(traits.side_of_oriented_circle_2_object()(p2, p0, p1, p3) ==
                 traits.side_of_oriented_circle_2_object()(p2, p0, p1, p3, o1, o1, o1, o1));
  CGAL_assertion(traits.side_of_oriented_circle_2_object()(p1, p2, p0, p3) ==
                 traits.side_of_oriented_circle_2_object()(p1, p2, p0, p3, o1, o1, o1, o1));
  CGAL_assertion(traits.side_of_oriented_circle_2_object()(p0, p2, p1, p3) ==
                 traits.side_of_oriented_circle_2_object()(p0, p2, p1, p3, o1, o1, o1, o1));
  CGAL_assertion(traits.side_of_oriented_circle_2_object()(p1, p0, p2, p3) ==
                 traits.side_of_oriented_circle_2_object()(p1, p0, p2, p3, o1, o1, o1, o1));
  CGAL_assertion(traits.side_of_oriented_circle_2_object()(p2, p1, p0, p3) ==
                 traits.side_of_oriented_circle_2_object()(p2, p1, p0, p3, o1, o1, o1, o1));
}

int main()
{
  test_orientation();
  test_in_circle();

  return 0;
}
