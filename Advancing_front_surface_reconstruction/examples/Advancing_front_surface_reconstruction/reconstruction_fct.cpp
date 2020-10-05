#include <iostream>
#include <fstream>
#include <algorithm>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Advancing_front_surface_reconstruction.h>
#include <CGAL/tuple.h>
#include <boost/lexical_cast.hpp>

typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_3  Point_3;

typedef std::array<std::size_t,3> Facet;

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

  template <typename AdvancingFront, typename Cell_handle>
  double operator() (const AdvancingFront& adv, Cell_handle& c,
                     const int& index) const
  {
    // bound == 0 is better than bound < infinity
    // as it avoids the distance computations
    if(bound == 0){
      return adv.smallest_radius_delaunay_sphere (c, index);
    }

    // If perimeter > bound, return infinity so that facet is not used
    double d  = 0;
    d = sqrt(squared_distance(c->vertex((index+1)%4)->point(),
                              c->vertex((index+2)%4)->point()));
    if(d>bound) return adv.infinity();
    d += sqrt(squared_distance(c->vertex((index+2)%4)->point(),
                               c->vertex((index+3)%4)->point()));
    if(d>bound) return adv.infinity();
    d += sqrt(squared_distance(c->vertex((index+1)%4)->point(),
                               c->vertex((index+3)%4)->point()));
    if(d>bound) return adv.infinity();

    // Otherwise, return usual priority value: smallest radius of
    // delaunay sphere
    return adv.smallest_radius_delaunay_sphere (c, index);
  }
};

int main(int argc, char* argv[])
{
  std::ifstream in((argc>1)?argv[1]:"data/half.xyz");
  double per = (argc>2)?boost::lexical_cast<double>(argv[2]):0;
  double radius_ratio_bound = (argc>3)?boost::lexical_cast<double>(argv[3]):5.0;

  std::vector<Point_3> points;
  std::vector<Facet> facets;

  std::copy(std::istream_iterator<Point_3>(in),
            std::istream_iterator<Point_3>(),
            std::back_inserter(points));

  Perimeter perimeter(per);
  CGAL::advancing_front_surface_reconstruction(points.begin(),
                                               points.end(),
                                               std::back_inserter(facets),
                                               perimeter,
                                               radius_ratio_bound);

  std::cout << "OFF\n" << points.size() << " " << facets.size() << " 0\n";
  std::copy(points.begin(),
            points.end(),
            std::ostream_iterator<Point_3>(std::cout, "\n"));
  std::copy(facets.begin(),
            facets.end(),
            std::ostream_iterator<Facet>(std::cout, "\n"));

  return 0;
}
