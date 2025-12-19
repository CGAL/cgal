
#include <cassert>

template <class R>
bool
_test_depth(const R& )
{
  typedef typename R::Point_3 Point_3;
  typedef typename R::Segment_3 Segment_3;

  Point_3 p(CGAL::ORIGIN), q(1, 1, 1) , r(1, 0, 0);
  Segment_3 s0(p, q), s1(p, r);

  Point_3 m = CGAL::midpoint(p,q);

  auto result = CGAL::intersection(s0, s1);
  const Point_3* ip = std::get_if<Point_3>(&*result);

  assert(CGAL::depth(p) == 0);
  assert(CGAL::depth(q) == 0);
  assert(CGAL::depth(m) == 1);
  assert(CGAL::depth(s0) == 1);
  assert(CGAL::depth(s1) == 1);
  assert(CGAL::depth(*ip) == 3);
  assert(CGAL::depth(m.x()) == 2);
  return true;
}
