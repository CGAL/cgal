#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polytope_distance_d.h>
#include <CGAL/Polytope_distance_d_traits_2.h>


typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_2 Point;
typedef CGAL::Polytope_distance_d_traits_2<K> Traits;

typedef CGAL::Polytope_distance_d<Traits> Dist;


int main()
{
  Point P[3] = { Point(0,0), Point(1,1), Point(1,0) };
  Point Q[3] = { Point(10,10), Point(11,11), Point(11,10) };
  Dist dist(P,P+3, Q, Q+3);
  return 0;

} 
