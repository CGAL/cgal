#include <CGAL/Polytope_distance_d.h>
#include <CGAL/Polytope_distance_d_traits_2.h>
#include <CGAL/Simple_cartesian.h>
#include <iostream>
#ifdef CGAL_USE_GMP
#include <CGAL/Gmpzf.h>
typedef CGAL::Gmpzf ET;
#else
#include <CGAL/MP_Float.h>
typedef CGAL::MP_Float ET;
#endif

typedef CGAL::Simple_cartesian<double>    K;
typedef K::Point_2                        Point;
typedef CGAL::Polytope_distance_d_traits_2<K, ET, double>
                                          Traits;
typedef CGAL::Polytope_distance_d<Traits> Polytope_distance;

int main(void)
{
  Point P[3] = {Point(0.6075439453125,  0.826324462890625),
                Point(0.853607177734375, 0.414520263671875),
                Point(0.376190185546875, 0.5203857421875)};
  Point Q[3] = {Point(0.54795532226562504, 0.31597900390625),
                Point(0.46067504882812499, 0.879180908203125),
                Point(1.234539794921875, 0.92095947265625)};

  Polytope_distance pd(P, P+3, Q, Q+3);

  std::cout << "Squared distance: " <<
    CGAL::to_double(pd.squared_distance_numerator()) /
    CGAL::to_double(pd.squared_distance_denominator()) << std::endl;

  std::cout << "Number of support points: " <<
    pd.number_of_support_points() << std::endl;
  return 0;
}
