#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/intersections.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel K;
typedef K::Point_2 Point_2;
typedef K::Segment_2 Segment_2;
typedef K::Line_2 Line_2;
typedef K::Intersect_2 Intersect_2;

int main()
{
  Segment_2 seg(Point_2(0,0), Point_2(2,2));
  Line_2 lin(1,-1,0);

  const auto result = intersection(seg, lin);
  if (result) {
    if (const Segment_2* s = std::get_if<Segment_2>(&*result)) {
      std::cout << *s << std::endl;
    } else {
      const Point_2* p = std::get_if<Point_2 >(&*result);
      std::cout << *p << std::endl;
    }
  }
  return 0;
}
