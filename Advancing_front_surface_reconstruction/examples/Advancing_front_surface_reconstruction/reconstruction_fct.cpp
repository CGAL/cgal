#include <iostream>
#include <fstream>
#include <algorithm>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Advancing_front_surface_reconstruction.h>
#include <CGAL/tuple.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3  Point_3;

typedef CGAL::cpp11::tuple<std::size_t,std::size_t,std::size_t> Facet;

namespace std {
  std::ostream& 
  operator<<(std::ostream& os, const Facet& f)
  {
    os << "3 " << get<0>(f) << " " << get<1>(f) << " " << get<2>(f);
    return os;
  }
}

struct Perimeter {

  double bound;

  Perimeter(double bound)
    : bound(bound)
  {}

  bool operator()(const Point_3& p, const Point_3& q, const Point_3& r) const
  {
    double d  = sqrt(squared_distance(p,q));
    if(d>bound) return true;
    d += sqrt(squared_distance(p,r)) ;
    if(d>bound) return true;
    d+= sqrt(squared_distance(q,r));
    return d>bound;
  }
};

int main(int argc, char* argv[])
{
  std::ifstream in((argc>1)?argv[1]:"data/half.xyz"); 
  std::vector<Point_3> points;
  std::vector<Facet> facets;
  
  std::copy(std::istream_iterator<Point_3>(in), 
            std::istream_iterator<Point_3>(), 
            std::back_inserter(points));
  
  Perimeter perimeter(0.5);
  CGAL::advancing_front_surface_reconstruction(points.begin(),
                                               points.end(),
                                               std::back_inserter(facets),
                                               perimeter);

  std::cout << "OFF\n" << points.size() << " " << facets.size() << " 0\n";
  std::copy(points.begin(),
            points.end(),
            std::ostream_iterator<Point_3>(std::cout, "\n"));
  std::copy(facets.begin(),
            facets.end(),
            std::ostream_iterator<Facet>(std::cout, "\n"));
  
  return 0;
}
