#include <CGAL/basic.h>
#include <CGAL/Filtered_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/squared_distance_2.h>
#include <cassert>
#include "MyKernel.h"
#include "MyPointC2_iostream.h"

typedef MyKernel<double>                   MK;
typedef CGAL::Filtered_kernel_adaptor<MK>  K;
typedef CGAL::Delaunay_triangulation_2<K>  Delaunay_triangulation_2;


typedef K::Point_2         Point;
typedef K::Segment_2       Segment;
typedef K::Ray_2           Ray;
typedef K::Line_2          Line;
typedef K::Triangle_2      Triangle;
typedef K::Iso_rectangle_2 Iso_rectangle;

const int RED= 1;
const int BLACK=2;

int main()
{
  Point a(0,0), b(1,0), c(1,1), d(0,1);
  a.color()=RED;
  b.color()=BLACK;
  d.color()=RED;

  Delaunay_triangulation_2 dt;
  dt.insert(a);

  K::Orientation_2 orientation;
  orientation(a,b,c);

  Point p(1,2), q;
  p.color() = RED;
  q.color() = BLACK;
  std::cout << p << std::endl;

  K::Compute_squared_distance_2 squared_distance;

  std::cout << "squared_distance(a, b) == "
            << squared_distance(a, b) << std::endl;

  Segment s1(p,q), s2(a, c);

  K::Construct_midpoint_2 construct_midpoint_2;

  Point mp = construct_midpoint_2(p,q);

  std::cout << "midpoint(" << p << " , " << q << ") == " << mp << std::endl;
  assert(s1.source().color() == RED);


  K::Intersect_2 intersection;

  CGAL::cpp11::result_of<K::Intersect_2(Segment, Segment)>::type
    intersect = intersection(s1, s2);

  K::Construct_cartesian_const_iterator_2 construct_it;
  K::Cartesian_const_iterator_2  cit = construct_it(a);
  assert(*cit == a.x());

  cit = construct_it(a,0);

  cit--;
  assert(*cit == a.y());

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
  intersection(t1, l1);

  intersection(t1, s1);

  intersection(t1, r1);

  Iso_rectangle i1(a,c), i2(d,p);
  intersection(i1, i2);
  intersection(i1, s1);

  intersection(i1, r1);
  intersection(i1, l1);

  t1.orientation();

  std::cout << s1.source() << std::endl;

  std::cout << t1.bbox() << std::endl;

  std::cout << "done" << std::endl;
  return 0;
}
