// Tests of `do_intersect()` between `Ball_3` (a solid ball) and other 3D objects.
//
// Every predicate is checked against an independent criterion: an object intersects
// the ball iff its squared distance to the center of the ball is not larger than the
// squared radius. It is also checked that the answer is `true` whenever the boundary
// sphere intersects the object, and that objects strictly inside the ball (which do not
// intersect the sphere) do intersect the ball.

#include <CGAL/intersection_3.h>
#include <CGAL/squared_distance_3.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Homogeneous.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_rational.h>
#include <CGAL/Exact_integer.h>
#include <CGAL/Random.h>

#include <CGAL/Cartesian_converter.h>

#include <iostream>
#include <cassert>
#include <type_traits>

// The expected answers are computed with an exact kernel `EK`
template <typename K, typename EK = K>
struct Ball_3_do_intersect_tester
{
  typedef typename K::FT              FT;
  typedef typename K::Point_3         P;
  typedef typename K::Vector_3        V;
  typedef typename K::Ball_3          B;
  typedef typename K::Sphere_3        Sph;
  typedef typename K::Line_3          L;
  typedef typename K::Ray_3           R;
  typedef typename K::Segment_3       S;
  typedef typename K::Plane_3         Pl;
  typedef typename K::Triangle_3      Tr;
  typedef typename K::Tetrahedron_3   Tet;
  typedef typename K::Iso_cuboid_3    Cub;

  template <typename T>
  static decltype(auto) X(const T& t)
  {
    if constexpr(std::is_same_v<K, EK>)
      return t;
    else
      return CGAL::Cartesian_converter<K, EK>()(t);
  }

  CGAL::Random& rnd;
  const int N = 500;
  const int M = 10; // coordinates in [-M, M]

  explicit Ball_3_do_intersect_tester(CGAL::Random& r) : rnd(r) { }

  P p(int x, int y, int z) const { return P(x, y, z); }
  P random_point() const { return P(rnd.get_int(-M, M+1), rnd.get_int(-M, M+1), rnd.get_int(-M, M+1)); }
  FT random_sq_radius() const { const int r = rnd.get_int(0, M); return FT(r*r + rnd.get_int(0, 3)); }

  template <typename O1, typename O2>
  void check(const O1& o1, const O2& o2, bool expected) const
  {
    assert(CGAL::do_intersect(o1, o2) == expected);
    assert(CGAL::do_intersect(o2, o1) == expected);
  }

  // Squared distance from a point to an iso-cuboid (0 if inside)
  static FT sq_dist(const P& c, const Cub& ic)
  {
    FT d(0);
    for(int i=0; i<3; ++i)
    {
      if(c[i] < (ic.min)()[i]) d += CGAL::square((ic.min)()[i] - c[i]);
      else if(c[i] > (ic.max)()[i]) d += CGAL::square(c[i] - (ic.max)()[i]);
    }
    return d;
  }

  void B_P()
  {
    std::cout << "Ball - Point\n";
    const B b(p(0,0,0), 4);
    check(b, p(0,0,0), true);   // center
    check(b, p(1,1,0), true);   // strictly inside: not on the sphere
    check(b, p(2,0,0), true);   // on the boundary
    check(b, p(2,1,0), false);
    assert(!CGAL::do_intersect(Sph(p(0,0,0), 4), p(1,1,0)));

    for(int i=0; i<N; ++i)
    {
      const B ball(random_point(), random_sq_radius());
      const P q = random_point();
      check(ball, q, CGAL::squared_distance(X(ball.center()), X(q)) <= X(ball.squared_radius()));
    }
  }

  void B_L_R_S_Pl()
  {
    std::cout << "Ball - Line/Ray/Segment/Plane\n";
    const B b(p(0,0,0), 4);
    check(b, S(p(-1,0,0), p(1,0,0)), true); // strictly inside the ball
    assert(!CGAL::do_intersect(Sph(p(0,0,0), 4), S(p(-1,0,0), p(1,0,0))));
    check(b, S(p(2,0,0), p(5,0,0)), true);  // touches the boundary
    check(b, S(p(3,0,0), p(5,0,0)), false);
    check(b, S(p(-5,1,0), p(5,1,0)), true); // crosses
    check(b, R(p(0,0,0), p(1,0,0)), true);  // starts at the center
    check(b, R(p(3,0,0), p(4,0,0)), false);
    check(b, L(p(0,2,0), p(1,2,0)), true);  // tangent
    check(b, L(p(0,3,0), p(1,3,0)), false);
    check(b, Pl(1, 0, 0, -2), true);        // tangent
    check(b, Pl(1, 0, 0, -3), false);

    for(int i=0; i<N; ++i)
    {
      const B ball(random_point(), random_sq_radius());
      const P a = random_point(), c = random_point(), d = random_point();
      if(a == c)
        continue;
      const FT sr = ball.squared_radius();
      check(ball, S(a, c), CGAL::squared_distance(X(ball.center()), X(S(a, c))) <= X(sr));
      check(ball, R(a, c), CGAL::squared_distance(X(ball.center()), X(R(a, c))) <= X(sr));
      check(ball, L(a, c), CGAL::squared_distance(X(ball.center()), X(L(a, c))) <= X(sr));
      if(!CGAL::collinear(a, c, d))
        check(ball, Pl(a, c, d), CGAL::squared_distance(X(ball.center()), X(Pl(a, c, d))) <= X(sr));
    }
  }

  void B_Tr_Tet()
  {
    std::cout << "Ball - Triangle/Tetrahedron\n";
    const B b(p(0,0,0), 100);
    const Tr small_tr(p(0,0,0), p(1,0,0), p(0,1,0));
    const Tet small_tet(p(0,0,0), p(1,0,0), p(0,1,0), p(0,0,1));
    check(b, small_tr, true); // strictly inside the ball
    check(b, small_tet, true);
    assert(!CGAL::do_intersect(Sph(p(0,0,0), 100), small_tr));
    check(B(p(1,1,1), 1), Tet(p(-5,-5,-5), p(20,0,0), p(0,20,0), p(0,0,20)), true); // ball inside the tetrahedron
    check(B(p(0,0,0), 1), Tet(p(-50,-50,-50), p(200,0,0), p(0,200,0), p(0,0,200)), true);
    check(B(p(30,30,30), 1), Tet(p(0,0,0), p(1,0,0), p(0,1,0), p(0,0,1)), false);

    for(int i=0; i<N; ++i)
    {
      const B ball(random_point(), random_sq_radius());
      const P a = random_point(), c = random_point(), d = random_point(), e = random_point();
      const FT sr = ball.squared_radius();
      if(!CGAL::collinear(a, c, d))
        check(ball, Tr(a, c, d), CGAL::squared_distance(X(ball.center()), X(Tr(a, c, d))) <= X(sr));
      if(CGAL::orientation(a, c, d, e) != CGAL::COPLANAR && sr != FT(0))
        check(ball, Tet(a, c, d, e), CGAL::squared_distance(X(ball.center()), X(Tet(a, c, d, e))) <= X(sr));
    }
  }

  void B_Cub_Bbox()
  {
    std::cout << "Ball - Iso_cuboid/Bbox\n";
    const B b(p(0,0,0), 100);
    const Cub small_cub(p(0,0,0), p(1,1,1));
    check(b, small_cub, true); // strictly inside the ball
    check(b, small_cub.bbox(), true);
    assert(!CGAL::do_intersect(Sph(p(0,0,0), 100), small_cub));
    check(B(p(0,0,0), 1), Cub(p(-5,-5,-5), p(5,5,5)), true); // ball inside the cuboid
    check(B(p(5,5,11), 1), Cub(p(0,0,0), p(10,10,10)), true); // touches a face
    check(B(p(11,11,11), 2), Cub(p(0,0,0), p(10,10,10)), false);

    for(int i=0; i<N; ++i)
    {
      const B ball(random_point(), random_sq_radius());
      const Cub ic(random_point(), random_point());
      const bool expected = sq_dist(ball.center(), ic) <= ball.squared_radius();
      check(ball, ic, expected);
      check(ball, ic.bbox(), expected);
    }
  }

  void B_B()
  {
    std::cout << "Ball - Ball\n";
    check(B(p(0,0,0), 1), B(p(0,0,0), 4), true);   // nested, concentric
    check(B(p(0,0,0), 1), B(p(1,0,0), 16), true);  // nested
    assert(!CGAL::do_intersect(Sph(p(0,0,0), 1), Sph(p(1,0,0), 16)));
    check(B(p(0,0,0), 1), B(p(3,0,0), 4), true);   // tangent
    check(B(p(0,0,0), 1), B(p(4,0,0), 4), false);
    check(B(p(0,0,0), 0), B(p(0,0,0), 0), true);   // two identical points

    for(int i=0; i<N; ++i)
    {
      // integer radii, so that the exact criterion does not need a square root
      const int r1 = rnd.get_int(0, M), r2 = rnd.get_int(0, M);
      const B b1(random_point(), FT(r1*r1)), b2(random_point(), FT(r2*r2));
      check(b1, b2, CGAL::squared_distance(b1.center(), b2.center()) <= FT((r1+r2)*(r1+r2)));
    }
  }

  // A ball intersects everything its boundary sphere intersects
  void B_vs_Sph()
  {
    std::cout << "Ball contains Sphere\n";
    for(int i=0; i<N; ++i)
    {
      const P c = random_point();
      const FT sr = random_sq_radius();
      const B ball(c, sr);
      const Sph sph(c, sr);
      const P a = random_point(), d = random_point(), e = random_point();
      if(a == d || CGAL::collinear(a, d, e))
        continue;
      if(CGAL::do_intersect(sph, S(a, d)))  assert(CGAL::do_intersect(ball, S(a, d)));
      if(CGAL::do_intersect(sph, Tr(a, d, e))) assert(CGAL::do_intersect(ball, Tr(a, d, e)));
      if(CGAL::do_intersect(sph, Cub(a, d))) assert(CGAL::do_intersect(ball, Cub(a, d)));
      if(CGAL::do_intersect(sph, a)) assert(CGAL::do_intersect(ball, a));
    }
  }

  void run()
  {
    B_P();
    B_L_R_S_Pl();
    B_Tr_Tet();
    B_Cub_Bbox();
    B_B();
    B_vs_Sph();
  }
};

int main(int, char**)
{
  CGAL::Random r;
  std::cout << "random seed = " << r.get_seed() << std::endl;

  std::cout << " |||||||| Test Simple_cartesian<Exact_rational> ||||||||" << std::endl;
  Ball_3_do_intersect_tester< CGAL::Simple_cartesian<CGAL::Exact_rational> >(r).run();

  std::cout << " |||||||| Test CGAL::Homogeneous<Exact_integer> ||||||||" << std::endl;
  Ball_3_do_intersect_tester< CGAL::Homogeneous<CGAL::Exact_integer> >(r).run();

  std::cout << " |||||||| Test EPICK ||||||||" << std::endl;
  Ball_3_do_intersect_tester< CGAL::Epick, CGAL::Epeck >(r).run();

  std::cout << " |||||||| Test EPECK ||||||||" << std::endl;
  Ball_3_do_intersect_tester< CGAL::Epeck >(r).run();

  std::cout << "OK!" << std::endl;
  return EXIT_SUCCESS;
}
