#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/intersections.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel K;
typedef K::Point_2 Point_2;
typedef K::Segment_2 Segment_2;
typedef K::Line_2 Line_2;
typedef K::Intersect_2 Intersect_2;

struct Intersection_visitor {
  void operator()(const Point_2& p) const
  {
    std::cout << p << std::endl;
  }

  void operator()(const Segment_2& s) const
  {
    std::cout << s << std::endl;
  }
};

int main()
{
  Segment_2 seg(Point_2(0,0), Point_2(1,1));
  Line_2 lin(1,-1,0);

  const auto result = intersection(seg, lin);
  if (result) {
    std::visit(Intersection_visitor(), *result);
  } else {
    // no intersection
  }

  return 0;
}
