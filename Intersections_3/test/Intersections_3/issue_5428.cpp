#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel EPICK;

typedef EPICK::Point_3 Point;

int main()
{
  EPICK::Plane_3 pl(Point(0, 0, 0), Point(1, 0, 0), Point(0, 1, 0));

  std::vector<Point> pts;
  pts.push_back(Point(-10,-10,-10));
  pts.push_back(Point(10,10,10));
  CGAL::Bbox_3 cub = CGAL::bbox_3(pts.begin(), pts.end());

  auto result = intersection(cub,pl);

  const std::vector<Point>* res = boost::get <std::vector<Point> >(&*result);
  for (Point p : *res) {
    std::cout << p << std::endl;
  }

  return 0;
}
