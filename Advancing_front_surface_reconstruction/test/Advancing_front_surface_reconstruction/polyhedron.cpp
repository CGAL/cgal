#include <iostream>
#include <fstream>
#include <algorithm>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Advancing_front_surface_reconstruction.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_3  Point_3;

typedef CGAL::Polyhedron_3<K> Polyhedron;

typedef CGAL::cpp11::array<std::size_t,3> Facet;

namespace std {
std::ostream&
operator<<(std::ostream& os, const Facet& f)
{
  os << "3 " << f[0] << " " << f[1] << " " << f[2];
  return os;
}

}


int main()
{
  Polyhedron polyhedron;
  std::ifstream in("data/planar.xyz");
  std::vector<Point_3> points;
  std::vector<Facet> facets;

  std::copy(std::istream_iterator<Point_3>(in),
            std::istream_iterator<Point_3>(),
            std::back_inserter(points));

  CGAL::advancing_front_surface_reconstruction(points.begin(),
                                               points.end(),
                                               polyhedron);

  std::cout << polyhedron << std::endl;

  return 0;
}
