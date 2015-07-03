#include <iostream>
#include <fstream>
#include <algorithm>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Advancing_front_surface_reconstruction.h>
#include <CGAL/tuple.h>

typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_3  Point_3;

typedef CGAL::cpp11::array<std::size_t,3> Facet;

namespace std {
std::ostream&
operator<<(std::ostream& os, const Facet& f)
{
  os << "3 " << f[0] << " " << f[1] << " " << f[2];
  return os;
}

}

void fct(const char* fname)
{
  std::ifstream in(fname);
  std::vector<Point_3> points;
  std::vector<Facet> facets;

  std::copy(std::istream_iterator<Point_3>(in),
            std::istream_iterator<Point_3>(),
            std::back_inserter(points));

  CGAL::advancing_front_surface_reconstruction(points.begin(),
                                               points.end(),
                                               std::back_inserter(facets));

  std::cout << "OFF\n" << points.size() << " " << facets.size() << " 0\n";
  std::copy(points.begin(),
            points.end(),
            std::ostream_iterator<Point_3>(std::cout, "\n"));
  std::copy(facets.begin(),
            facets.end(),
            std::ostream_iterator<Facet>(std::cout, "\n"));
}

int main()
{
  fct("data/point.xyz");
  fct("data/segment.xyz");
  fct("data/triangle.xyz");
  fct("data/planar.xyz");
  return 0;
}
