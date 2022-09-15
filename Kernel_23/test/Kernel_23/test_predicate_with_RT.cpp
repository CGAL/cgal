#define CGAL_NO_MPZF_DIVISION_OPERATOR 1

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Mpzf.h>

template <typename R>
void test(const R& rep) {
  using Point_3 = typename R::Point_3;
  using Segment_3 = typename R::Segment_3;
  using Line_3 = typename R::Line_3;

  auto construct_point = rep.construct_point_3_object();
  Point_3 p2 = construct_point(CGAL::ORIGIN);
  Point_3 p3 = construct_point(1,1,1);
  Point_3 p4 = construct_point(1,1,2);
  Point_3 p5 = construct_point(1,2,3);
  Point_3 p6 = construct_point(4,2,1);

  auto construct_segment = rep.construct_segment_3_object();
  Segment_3 s2 = construct_segment(p2,p3), s1 = s2;

  auto construct_line = rep.construct_line_3_object();
  Line_3 l2 = construct_line(p5,p6);

  auto compare_distance = rep.compare_distance_3_object();
  // compare_distance(p2, p2, p2);
  compare_distance(p2, s2, p2);
  // compare_distance(p2, l2, p2);
}

int main()
{
  test(CGAL::Simple_cartesian<CGAL::Mpzf>());
  return 0;
}
