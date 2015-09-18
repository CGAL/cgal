#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/intersections.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel K;
typedef K::Point_2 Point_2;
typedef K::Segment_2 Segment_2;
typedef K::Line_2 Line_2;
typedef K::Intersect_2 Intersect_2;

int main()
{
  Segment_2 seg(Point_2(0,0), Point_2(1,1));
  Line_2 lin(1,0,0);

  if (result) {
    if (const Segment_2* s = boost::get<Segment_2>(&*result)) {
      // handle segment
    } else {
      const Point_2* p = boost::get<Point_2 >(&*result);
      // handle point
    }
  }
  return 0;
}
