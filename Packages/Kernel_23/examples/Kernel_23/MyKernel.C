#include <CGAL/basic.h>
#include <CGAl/Filtered_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/squared_distance_2.h>

#include "MyKernel.h"


typedef MyKernel<double>      MK;
typedef CGAL::Filtered_kernel<MK> K;
typedef CGAL::Delaunay_triangulation_2<K> Delaunay_triangulation_2;


typedef K::Point_2         Point;
typedef K::Segment_2       Segment;
typedef K::Ray_2           Ray;
typedef K::Line_2          Line;
typedef K::Triangle_2      Triangle;
typedef K::Iso_rectangle_2 Iso_rectangle;


int main()
{
  Point a(0,0), b(1,0), c(1,1), d(0,1);

  Delaunay_triangulation_2 dt;
  dt.insert(a);

  K::Orientation_2 orientation;
  orientation(a,b,c);
  Point p(1,2), q;
  p.c() = 7812;
  q.c() = 13;
  std::cout << p << std::endl;

  K::Compute_squared_distance_2 squared_distance;

  std::cout << "squared_distance(a, b) == " << squared_distance(a, b) << std::endl;

  Segment s1(p,q), s2(a, c);
  K::Intersect_2 intersection;

  CGAL::Object o = intersection(s1, s2);

  p.cartesian_begin();
  Line l1(a,b), l2(p, q);
  intersection(l1, l2);


  intersection(s1, l1);
  
  Ray r1(d,b), r2(d,c);
  intersection(r1, r2);
  
  intersection(r1, l1);

  squared_distance(r1, r2);
  squared_distance(r1, l2);
  squared_distance(r1, s2);

  Triangle t1(a,b,c), t2(a,c,d);
  intersection(t1, t2);

  intersection(t1, s1);

  intersection(t1, l1);
  intersection(t1, r1);
  
  Iso_rectangle i1(a,c), i2(d,p);
  intersection(i1, i2);
  intersection(i1, s1);
  intersection(i1, r1);
  intersection(i1, l1);
  
  t1.orientation();

  std::cout << s1.source() << std::endl;
  return 0;
}
