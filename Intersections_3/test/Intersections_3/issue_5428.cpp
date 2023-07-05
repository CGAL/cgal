#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel EPICK;

typedef EPICK::Point_3 Point;

int main(int, char**)
{
  std::vector<Point> pts;
  pts.emplace_back(-10,-10,-10);
  pts.emplace_back(10,10,10);

  EPICK::Plane_3 pl(Point(0, 0, 0), Point(1, 0, 0), Point(0, 1, 0));
  CGAL::Bbox_3 cub = CGAL::bbox_3(pts.begin(), pts.end());

  auto result = intersection(cub, pl);

  const std::vector<Point>* res = std::get_if <std::vector<Point> >(&*result);
  for(const Point& p : *res)
    std::cout << p << std::endl;

  assert(res->size() == 4);

  return 0;
}
