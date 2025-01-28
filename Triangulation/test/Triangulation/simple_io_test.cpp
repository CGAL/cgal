#include <CGAL/Epick_d.h>
#include <CGAL/Delaunay_triangulation.h>
#include <sstream>

int main()
{
  typedef CGAL::Delaunay_triangulation<CGAL::Epick_d<CGAL::Dimension_tag<2>>> T;
  T dt1(2), dt2(2);

  std::vector<T::Point> points;
  points.emplace_back(1,0);
  points.emplace_back(0,1);
  points.emplace_back(2,2);
  dt1.insert(points.begin(), points.end());

  std::stringstream f;
  f << dt1 << std::endl;
  std::cout << f.str();
  f >> dt2;
  assert(dt2.is_valid(true));
}
