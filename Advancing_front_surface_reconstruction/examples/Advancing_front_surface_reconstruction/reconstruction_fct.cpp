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

struct Perimeter {

  double bound;

  Perimeter(double bound)
    : bound(bound)
  {}

  // The point type that will be injected here will be
  // CGAL::Exact_predicates_inexact_constructions_kernel::Point_3
  template <typename Point>
  bool operator()(const Point& p, const Point& q, const Point& r) const
  {
    // bound == 0 is better than bound < infinity
    // as it avoids the distance computations
    if(bound == 0){
      return false;
    }
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

  Perimeter perimeter(0);
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
