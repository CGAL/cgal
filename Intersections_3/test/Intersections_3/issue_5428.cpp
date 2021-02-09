#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel EPECK;
typedef CGAL::Exact_predicates_inexact_constructions_kernel EPICK;

typedef CGAL::Cartesian_converter<EPICK,EPECK> IK_to_EK;

typedef EPICK::Point_3 Point;

int main()
{
  IK_to_EK to_exact;

  //EPICK::Plane_3 pl(Point(0,0,10),Point(1,0,10),Point(0,1,10));

  // EPICK::Plane_3 pl(0.265189, 0.902464, 0.33946, -2.47551);

  EPICK::Plane_3 pl(Point(1, 1, 1), Point(1, 2, 1), Point(1, 2, 2));

  EPECK::Plane_3 epl = to_exact(pl);

  std::vector<Point> pts;
  pts.push_back(Point(1,1,1));
  pts.push_back(Point(2,2,2));
  CGAL::Bbox_3 cub = CGAL::bbox_3(pts.begin(), pts.end());

  EPECK::Iso_cuboid_3 ecub = to_exact(cub);

  auto result = intersection(cub,pl);

  //auto result2 = intersection(ecub,epl);

  const std::vector<Point>* res = boost::get <std::vector<Point> >(&*result);
  for (Point p : *res) {
	  std::cout << p << std::endl;
  }

  return 0;
}
