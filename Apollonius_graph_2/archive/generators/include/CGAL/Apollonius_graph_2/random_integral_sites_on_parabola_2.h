#ifndef CGAL_APOLLONIUS_GRAPH_2_RANDOM_INTEGRAL_SITES_ON_PARABOLA_2_H
#define CGAL_APOLLONIUS_GRAPH_2_RANDOM_INTEGRAL_SITES_ON_PARABOLA_2_H 1

#include <CGAL/basic.h>
#include <CGAL/random_integer.h>

namespace CGAL {

// creates a random site of the form {(t, t^2), t^2}, where t is in
// the range [-M, M], M = 2^b - 1

template<class Site, class R>
class Random_integral_sites_on_parabola_2
{
public:
  typedef Site  Site_2;
  typedef R     Random;

public:
  Random_integral_sites_on_parabola_2(unsigned int b)
    : b_(b), p_(0), r_(0) {
    CGAL_precondition( b >= 0 && b <= 26 );
  }

  Random_integral_sites_on_parabola_2(unsigned int b, int seed)
    : b_(b), p_(0), r_(seed) {
    CGAL_precondition( b >= 0 && b <= 26 );
  }

  Random_integral_sites_on_parabola_2(unsigned int b, unsigned int p,
				      int seed)
    : b_(b), p_(p), r_(seed) {
    CGAL_precondition( b >= 0 && b <= 26 );
    CGAL_precondition( p >= 0 && p <= 26 );
  }

  Site_2 operator*()
  {
    double t = random_integer(r_, b_, true);
    double t2 = t * t;
    double perb = random_integer(r_, p_, true);

    typename Site_2::Point_2 p(t, t2);
    return Site_2(p, t2 + perb);
  }

private:
  unsigned int b_, p_;
  Random r_;
};

} //namespace CGAL

#endif // CGAL_APOLLONIUS_GRAPH_2_RANDOM_INTEGRAL_SITES_ON_PARABOLA_2_H
