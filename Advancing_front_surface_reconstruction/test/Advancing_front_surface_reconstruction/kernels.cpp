#include <iostream>
#include <fstream>
#include <algorithm>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Advancing_front_surface_reconstruction.h>

typedef CGAL::Simple_cartesian<double> K;


typedef std::array<std::size_t,3> Facet;

namespace std {
std::ostream&
operator<<(std::ostream& os, const Facet& f)
{
  os << "3 " << f[0] << " " << f[1] << " " << f[2];
  return os;
}

}

template <typename K>
void fct(const char* fname)
{
  typedef typename K::Point_3  Point_3;
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
  {
    typedef CGAL::Simple_cartesian<float> K;
    fct<K>("data/planar.xyz");
  }
  {
    typedef CGAL::Simple_cartesian<double> K;
    fct<K>("data/planar.xyz");
  }
  {
    typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
    fct<K>("data/planar.xyz");
  }
  {
    typedef CGAL::Exact_predicates_exact_constructions_kernel K;
    fct<K>("data/planar.xyz");
  }
  return 0;
}
