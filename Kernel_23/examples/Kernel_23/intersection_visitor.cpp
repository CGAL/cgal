#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/intersections.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel K;
typedef K::Point_2 Point_2;
typedef K::Segment_2 Segment_2;
typedef K::Line_2 Line_2;
typedef K::Intersect_2 Intersect_2;

struct Intersection_visitor {
  typedef void result_type;

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

  // with C++11 support
  // auto result = intersection(seg, lin);
  // without C++11
  CGAL::cpp11::result_of<Intersect_2(Segment_2, Line_2)>::type
    result = intersection(seg, lin);
  if (result) {
    boost::apply_visitor(Intersection_visitor(), *result);
  } else {
    // no intersection
  }

  return 0;
}
