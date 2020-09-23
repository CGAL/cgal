#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>


typedef CGAL::Exact_predicates_exact_constructions_kernel Epeck;
typedef CGAL::Exact_predicates_inexact_constructions_kernel Epick;

template<typename K>
void test_P_Cub()
{
  typedef CGAL::Point_3<K> P;
  typedef CGAL::Iso_cuboid_3<K> Cub;
  P p(0.0,1.0,0.5);
  Cub cub(P(-1.0,-1.0,-1.0), P(1.0,1.0,1.0));
  assert(CGAL::do_intersect(p, cub));
  CGAL::Object o = CGAL::intersection(p,cub);
  P res;
  assert(assign(res, o));
  assert(res == p);
}
template<typename K>
void test_P_L()
{
  typedef CGAL::Point_3<K> P;
  typedef CGAL::Line_3<K> Line;
  P p(0.99,0.99,0.99);
  Line line(P(-1.0,-1.0,-1.0), P(1.0,1.0,1.0));
  assert(CGAL::do_intersect(p, line));
  CGAL::Object o = CGAL::intersection(p,line);
  P res;
  assert(assign(res, o));
  assert(res == p);
}
template<typename K>
void test_P_R()
{
  typedef CGAL::Point_3<K> P;
  typedef CGAL::Ray_3<K> Ray;
  P p(0.99,0.99,0.99);
  Ray ray(P(-1.0,-1.0,-1.0), P(1.0,1.0,1.0));
  assert(CGAL::do_intersect(p, ray));
  CGAL::Object o = CGAL::intersection(p,ray);
  P res;
  assert(assign(res, o));
  assert(res == p);
}
template<typename K>
void test_P_S()
{
  typedef CGAL::Point_3<K> P;
  typedef CGAL::Segment_3<K> S;
  P p(0.99,0.99,0.99);
  S s(P(-1.0,-1.0,-1.0), P(1.0,1.0,1.0));
  assert(CGAL::do_intersect(p, s));
  CGAL::Object o = CGAL::intersection(p,s);
  P res;
  assert(assign(res, o));
  assert(res == p);
}
template<typename K>
void test_P_P()
{
  typedef CGAL::Point_3<K> P;
  P p(0.99,0.99,0.99);
  assert(CGAL::do_intersect(p, p));
  CGAL::Object o = CGAL::intersection(p,p);
  P res;
  assert(assign(res, o));
  assert(res == p);
}
template<typename K>
void test_P_Pl()
{
  typedef CGAL::Point_3<K> P;
  typedef CGAL::Plane_3<K> Pl;
  P p(0.99,0.99,0.99);
  Pl pl(P(-1.0,-1.0,-1.0), P(1.0,1.0,1.0), P(0.0,0.0,0.0));
  assert(CGAL::do_intersect(p, pl));
  CGAL::Object o = CGAL::intersection(p,pl);
  P res;
  assert(assign(res, o));
  assert(res == p);
}

template<typename K>
void test_P_Tet()
{
  typedef CGAL::Point_3<K> P;
  typedef CGAL::Tetrahedron_3<K> T;
  P p(0,0,0), q(1,0,0), r(1,1,0), s(0,0,1);
  T t(p,q,r,s);
  P q0(0.1, 0.1, 0.1), q1(10,10,10);

  assert(CGAL::do_intersect(p,t));
  assert(CGAL::do_intersect(t,p));

  assert(CGAL::do_intersect(q0,t));
  assert(CGAL::do_intersect(t,q0));

  assert(! CGAL::do_intersect(q1,t));
  assert(! CGAL::do_intersect(t,q1));
}

int main()
{
  test_P_Cub<Epick>();
  test_P_Cub<Epeck>();

  test_P_L<Epick>();
  test_P_L<Epeck>();

  test_P_R<Epick>();
  test_P_R<Epeck>();
  test_P_S<Epick>();
  test_P_S<Epeck>();

  test_P_P<Epick>();
  test_P_P<Epeck>();

  test_P_Pl<Epick>();
  test_P_Pl<Epeck>();

  test_P_Tet<Epick>();
  test_P_Tet<Epeck>();
}
#include <iostream>



