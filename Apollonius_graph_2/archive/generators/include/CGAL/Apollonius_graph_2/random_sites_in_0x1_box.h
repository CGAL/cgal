#ifndef CGAL_APOLLONIUS_GRAPH_2_RANDOM_SITES_IN_0X1_BOX_H
#define CGAL_APOLLONIUS_GRAPH_2_RANDOM_SITES_IN_0X1_BOX_H 1

#include <CGAL/basic.h>
#include <CGAL/Random.h>

namespace CGAL {

template<class Site>
class Random_sites_in_0x1_box
{
public:
  typedef Site    Site_2;
  typedef Site_2  result_type;

private:
  typedef typename Site_2::Point_2  Point_2;

public:
  Random_sites_in_0x1_box(double rmax = 0.125, int seed = 0)
    : rmax_(rmax), random_(seed)  {}

  Site_2 operator*() {
    double x = random_.get_double(0, 1);
    double y = random_.get_double(0, 1);
    double w = random_.get_double(0, rmax_);
    return Site_2(Point_2(x,y),w);
  }

private:
  double rmax_;
  Random random_;
};

} //namespace CGAL

#endif // CGAL_APOLLONIUS_GRAPH_2_RANDOM_SITES_IN_0X1_BOX_H
