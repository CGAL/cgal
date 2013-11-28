

#include <iostream>
#include <algorithm>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Advancing_front_surface_reconstruction.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3  Point_3;
typedef std::vector<Point_3> V;

typedef CGAL::Triple<std::size_t,std::size_t,std::size_t> Facet;
typedef std::vector<Facet> F;

namespace std {
std::ostream& 
operator<<(std::ostream& os, const Facet f)
{
  os << "3 " << f.first << " " << f.second << " " << f.third << std::endl;
  return os;
}

}

int main()
{
  V points;
  F facets;

  std::copy(std::istream_iterator<Point_3>(std::cin), std::istream_iterator<Point_3>(), std::back_inserter(points));
  CGAL::advancing_front_surface_reconstruction(points.begin(), points.end(), std::back_inserter(facets));
  std::cout << "OFF\n" << points.size() << " " << facets.size() << " 0\n";
  std::copy(points.begin(), points.end(), std::ostream_iterator<Point_3>(std::cout, "\n"));
  std::copy(facets.begin(), facets.end(), std::ostream_iterator<Facet>(std::cout));
  return 0;
}
