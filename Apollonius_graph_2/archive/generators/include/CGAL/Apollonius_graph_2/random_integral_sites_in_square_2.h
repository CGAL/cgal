#ifndef CGAL_APOLLONIUS_GRAPH_2_RANDOM_INTEGRAL_SITES_IN_SQUARE_2_H
#define CGAL_APOLLONIUS_GRAPH_2_RANDOM_INTEGRAL_SITES_IN_SQUARE_2_H 1

#include <CGAL/basic.h>
#include <CGAL/random_integer.h>

namespace CGAL {

// creates a random site with x, y in [-M,M]x[-M,M]
// where M = 2^b - 1 and r in [0,R], R = 2^B - 1

template<class Site, class R>
class Random_integral_sites_in_square_2
{
public:
  typedef Site  Site_2;
  typedef R     Random;

public:
  Random_integral_sites_in_square_2(unsigned int b)
    : b_(b), B_(b), r_(0) {
    CGAL_precondition( b >= 0 && b <= 52 );
  }

  Random_integral_sites_in_square_2(unsigned int b, int seed)
    : b_(b), B_(b), r_(seed) {
    CGAL_precondition( b >= 0 && b <= 52 );
  }

  Random_integral_sites_in_square_2(unsigned int b, unsigned int B, int seed)
    : b_(b), B_(B), r_(seed) {
    CGAL_precondition( b >= 0 && b <= 52 );
    CGAL_precondition( B >= 0 && B <= 52 );
  }

  Site_2 operator*()
  {
    double x = random_integer(r_, b_, true);
    double y = random_integer(r_, b_, true);
    double w = random_integer(r_, B_, false);
    CGAL_assertion( w >= 0 );

    typename Site_2::Point_2 p(x, y);
    return Site_2(p, w);
  }

private:
  unsigned int b_, B_;
  Random r_;
};

} //namespace CGAL

#endif // CGAL_APOLLONIUS_GRAPH_2_RANDOM_INTEGRAL_SITES_IN_SQUARE_2_H
