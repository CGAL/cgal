#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Projection_traits_xy_3.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Epik;
typedef CGAL::Projection_traits_xy_3<Epik> K;

typedef K::Orientation_2 Orientation_2;
typedef K::Intersect_2 Intersect_2;

typedef K::Point_2 Point_2;
typedef K::Segment_2 Segment_2;
typedef CGAL::Object Object;

int main()
{
  Point_2 p(0,0,0), q(1,1,1), r(1,0,0), s(0,1,1), t(0,1,0), u(0,1,-1);
  Point_2 v(0.5, 0, 0), w(0.5,1,0);

  Segment_2 pq(p,q), rs(r,s), rt(r,t), ru(r,u), vw(v,w);

  Point_2 pqrs, pqrt, pqru, pqvw;

  Object o = Intersect_2()(pq,rs);
  assert(assign(pqrs,o));
  assert(pqrs == Point_2(0.5, 0.5, 0.5));

  o = Intersect_2()(pq,rt);
  assert(assign(pqrt,o));
  assert(pqrt == Point_2(0.5, 0.5, 0.25));

 o = Intersect_2()(pq,ru);
 assert(assign(pqru,o));
 assert(pqru == Point_2(0.5, 0.5, 0));

 o = Intersect_2()(pq,vw);
 assert(assign(pqvw,o));
 assert(pqvw == Point_2(0.5, 0.5, 0.25));

 std::cerr << "done" << std::endl;
 return 0;
}
