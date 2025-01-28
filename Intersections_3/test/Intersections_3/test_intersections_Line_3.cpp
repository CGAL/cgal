// 3D intersection tests.

// We want to check that no division is performed for interface macro Do_intersect_3_RT
#define CGAL_NO_MPZF_DIVISION_OPERATOR

// Line_3_X with X < Line (lexicographically) are tested in other files

#include <CGAL/Intersections_3/Line_3_Line_3.h>
#include <CGAL/Intersections_3/Line_3_Plane_3.h>
#include <CGAL/Intersections_3/Line_3_Point_3.h>
#include <CGAL/Intersections_3/Line_3_Ray_3.h>
#include <CGAL/Intersections_3/Line_3_Segment_3.h>
#include <CGAL/Intersections_3/Line_3_Sphere_3.h>
#include <CGAL/Intersections_3/Line_3_Tetrahedron_3.h>
#include <CGAL/Intersections_3/Line_3_Triangle_3.h>

// just to cross verify
#include <CGAL/Intersections_3/Plane_3_Segment_3.h>
#include <CGAL/Intersections_3/Plane_3_Ray_3.h>
#include <CGAL/Intersections_3/Segment_3_Segment_3.h>
#include <CGAL/Intersections_3/Segment_3_Tetrahedron_3.h>
#include <CGAL/Intersections_3/Segment_3_Triangle_3.h>

#include <CGAL/Vector_3.h>

#include "intersection_test_helper.h"

#include <CGAL/MP_Float.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Homogeneous.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <iostream>
#include <cassert>

template <typename K>
struct Line_3_intersection_tester
  : public Intersection_3_tester<K>
{
  typedef Intersection_3_tester<K>  Base;

  typedef typename K::FT            FT;

  typedef CGAL::Iso_cuboid_3<K>     Cub;
  typedef CGAL::Line_3<K>           L;
  typedef CGAL::Point_3<K>          P;
  typedef CGAL::Plane_3<K>          Pl;
  typedef CGAL::Ray_3<K>            R;
  typedef CGAL::Segment_3<K>        S;
  typedef CGAL::Sphere_3<K>         Sph;
  typedef CGAL::Tetrahedron_3<K>    Tet;
  typedef CGAL::Triangle_3<K>       Tr;

  typedef CGAL::Vector_3<K>         V;

  using Base::p;
  using Base::pl;
  using Base::random_point;

  using Base::check_do_intersect;
  using Base::check_do_not_intersect;
  using Base::check_intersection;
  using Base::check_no_intersection;

private:
  int N = 1000;

public:
  Line_3_intersection_tester(CGAL::Random& r,
                             const bool has_exact_p = false,
                             const bool has_exact_c = false)
    : Base(r, has_exact_p, has_exact_c)
  { }

public:
  void L_L()
  {
    std::cout << "Line - Line\n";

    for(int vx=0; vx<4; ++vx)
    {
      for(int vy=1; vy<4; ++vy)
      {
        for(int vz=0; vz<4; ++vz)
        {
          if(vx == 0 && vy == 0 && vz == 0)
            continue;

          const int a = vx;
          const int b = vy;
          const int c = vz;

          const L l =   L(p(0,0,0), p( a,  b, c));
          const L l_1 = L(p(0,0,0), p(-a, -b,-c));
          const L l_2 = L(p(0,0,0), p(-a, -b, 7));
          const L l_3 = L(p(1,0,0), p( a,  b, c));
          const L l_4 = L(p(1,0,0), p( a+1,b, c));

          const CGAL::Object obj1 = CGAL::intersection(l, l_1);
          const CGAL::Object obj2 = CGAL::intersection(l, l_2);
          const CGAL::Object obj3 = CGAL::intersection(l, l_3);
          const CGAL::Object obj4 = CGAL::intersection(l, l_4);

          check_do_intersect(l, l_1);
          check_do_intersect(l, l_2);
          check_do_intersect(l, l_3);
          check_do_not_intersect(l, l_4);

          Base::template check_intersection<L>(l, l_1);
          L interl1;
          assert(assign(interl1, obj1));

          check_intersection(l, l_2, P(0,0,0));
          P interp2;
          assert(assign(interp2, obj2));
          assert(interp2 == P(0,0,0));

          check_intersection(l, l_3, P(a,b,c));
          P interp3;
          assert(assign(interp3, obj3));
          assert(interp3 == P(a,b,c));

          check_no_intersection(l, l_4);
          assert(obj4.is_empty());
        }
      }
    }
  }

  void L_Pl()
  {
    std::cout << "Line - Plane\n";

    // No intersection
    check_no_intersection         (pl(0, 0, 1,-2), L(p(1,1,1), p(2,3,1)));

    // Point intersection
    check_intersection            (L(p(1,1,1), p(2,3,4)), pl(1, 1, 1, 0),
                                   P(0.5, 0, -0.5));
    check_intersection            (L(p(1,1,1), p(2,3,-1)), pl(1, 0, 1, 3),
                                   P(6,11,-9));
    Base::template check_intersection<P>(L(p(1,1,1), p(2,3,4)), pl(1, 2, 4, 7));
    Base::template check_intersection<P>(L(p(-1,-1,-1), p(2,3,4)), pl(1, 2, 4, 7),
                                         pl(1, 2, 4, 7), S(p(-1,-1,-1), p(2,3,4)));
    Base::template check_intersection<P>(L(p(2,3,4), p(-1,-1,-1)), pl(1, 2, 4, 7),
                                         pl(1, 2, 4, 7), S(p(-1,-1,-1), p(2,3,4)));

    if(this->has_exact_c) // because Plane and Line inner coeffs are actually constructions...
    {
      for(int i=0; i<N; ++i)
      {
        P pl0 = random_point(), pl1 = random_point(), pl2 = random_point();
        P l0 = random_point(), l1 = random_point();

        Pl pl(pl0, pl1, pl2);
        P pl3 = pl0 + FT(this->r.get_double()) * V(pl1 - pl0) + FT(this->r.get_double()) * V(pl1 - pl0);
        if(pl.has_on(l1))
          Base::check_intersection(L(pl3, l1), pl, L(pl3, l1)); // both points on the plane
        else
          Base::check_intersection(L(pl3, l1), pl, pl3); // single point on the plane

        if(pl.oriented_side(l0) != pl.oriented_side(l1)) // l0 xor l1 on pl is fine
        {
          Base::template check_intersection<P>(L(l0, l1), pl, pl, S(l0, l1));
          Base::template check_intersection<P>(L(l0, l1), pl, pl, R(l0, l1));
        }
      }
    }

    // Line intersection
    Base::template check_intersection<L>(L(p(1,1,1), p(2,3,1)), pl(0, 0, 1,-1));

    if(this->has_exact_c)
    {
      for(int i=0; i<N; ++i)
      {
        P pl0 = random_point(), pl1 = random_point(), pl2 = random_point();
        if(CGAL::collinear(pl0, pl1, pl2))
          continue;

        Pl pl(pl0, pl1, pl2);
        L l(pl0 + (pl2 - pl1), pl0 + (pl1 - pl2));
        check_intersection(l, pl, l);
      }
    }
  }

  void L_P()
  {
    std::cout << "Line - Point\n";

    // No intersection
    check_no_intersection(L(p(1,1,1), p(2,3,4)), p(0,0,0));
    check_no_intersection(L(p(1,1,1), p(2,3,4)), P(3,5,7+1e-8));
    check_no_intersection(L(p(2,3,4), p(1,1,1)), P(3,5,7+1e-8));

    // Extremity
    check_intersection   (L(p(1,1,1), p(2,3,4)), p(1,1,1), p(1,1,1));
    check_intersection   (L(p(1,1,1), p(2,3,4)), p(2,3,4), p(2,3,4));

    // Point on the line
    check_intersection   (L(p(1,1,1), p(2,3,4)), p(3,5,7), p(3,5,7));
    check_intersection   (L(p(1,1,1), p(2,3,4)), P(1.5,2,2.5), P(1.5,2,2.5));

    if(this->has_exact_c)
    {
      for(int i=0; i<N; ++i)
      {
        P l0 = random_point(), l1 = random_point();
        if(l0 == l1)
          continue;

        L l(l0, l1);
        P pt = l.point(double(i)/N);
        check_intersection(l, pt, pt);

        P pt2 = random_point();
        if(CGAL::collinear(l.point(0), l.point(1), pt2))
          check_intersection(l, pt2, pt2);
        else
          check_no_intersection(l, pt2);
      }
    }
  }

  void L_R()
  {
    std::cout << "Line - Ray\n";

    // No intersection
    check_no_intersection(L(p(0,0,0),p(1,0,0)), R(p(3,0,1),p(6,0,1)));
    check_no_intersection(L(p(0,0,0),p(1,0,0)), R(p(0,2,0),p(0,4,0)));
    check_no_intersection(L(p(0,0,0),p(1,0,0)), R(p(6,2,0),p(5,4,0)));
    check_no_intersection(L(p(0,0,0),p(0,1,0)), R(p(1,-1,0),p(1,0,0)));
    check_no_intersection(L(p(0,-10,0),p(0,-9,0)), R(p(1,-1,0),p(2,0,0)));
    check_no_intersection(L(p(0,-10,0),p(0,0,0)), R(p(1,-1,0),p(2,0,0)));
    check_no_intersection(L(p(0,0,0),p(0,1,0)), R(p(1,-1,0),p(2,0,0)));

    // Point intersection
    check_intersection   (L(p(0,0,0),p(1,0,0)), R(p(3,0,0),p(6,4,0)),
                          P(3,0,0));
    check_intersection   (L(p(0,0,0),p(1,0,0)), R(p(5,-2,0),p(5,-1,0)),
                          P(5,0,0));
    check_intersection   (L(p(0,0,0),p(1,0,0)), R(p(0,-2,0),p(0,-1,0)),
                          P(0,0,0));

    // Ray Intersection
    check_intersection   (L(p(0,0,0),p(1,0,0)), R(p(3,0,0),p(6,0,0)),
                          R(P(3,0,0),P(6,0,0)));
  }

  void L_S()
  {
    std::cout << "Line - Segment\n";

    // No intersection
    check_no_intersection(L(p(0,1,0),p(0,-1,0)), S(p(4,2,-1),p(8,4,3)));
    check_no_intersection(L(p(3,1,9),p(7,-1,2)), S(p(4,2,-1),p(8,4,3)));

    // Point
    check_intersection(L(p(0,1,0),p(0,-1,0)), S(p(0,1,0),p(1,-5,-1)),
                       P(0,1,0));
    check_intersection(L(p(0,1,0),p(0,-1,0)), S(p(0,0,0),p(12,0,0)),
                       P(0,0,0));
    check_intersection(L(p(0,1,0),p(0,-1,0)), S(p(-4,3,0),p(4,-3,0)),
                       P(0,0,0));

    // Segment
    check_intersection(L(p(0,0,0),p(1,3,7)), S(p(3,9,21),p(6,18,42)),
                       S(p(3,9,21),p(6,18,42)));
    check_intersection(L(p(0,0,0),p(1,0,0)), S(p(0,0,0),p(-9,0,0)),
                       S(P(0,0,0),P(-9,0,0)));

    for(int i=0; i<N; ++i)
    {
      P l0 = random_point(), l1 = random_point();
      P m = CGAL::midpoint(l0, l1);
      P s0 = m + V(CGAL::ORIGIN, P(1e-5, 1e-5, 1e-5)), s1 = m + (random_point() - s0);

      if(l0 == l1 || s0 == s1)
        continue;

      if(CGAL::do_intersect(S(l0, l1), S(s0, s1)))
        check_intersection(L(l0, l1), S(s0, s1), S(l0, l1), S(s0, s1));
    }
  }

  void L_Sph()
  {
    std::cout << "Line - Sphere\n";

    // No intersection
    check_do_not_intersect(L(p(10,4,9), p(9,-4,8)), Sph(p(1,0,0), p(0,1,0), p(0,0,1), p(-1,0,0)));
    check_do_not_intersect(L(p(-1,10,-3), p(-1,-4,5)), Sph(p(1,0,0), p(0,1,0), p(0,0,1), p(-1,0,0)));
    check_do_not_intersect(L(p(-1,-2,-3), p(5,1,8)), Sph(p(1,0,0), p(0,1,0), p(0,0,1), p(-1,0,0)));

    // Point
    check_do_intersect(L(p(-1,8,-10), p(-1,-4,5)), Sph(p(1,0,0), p(0,1,0), p(0,0,1), p(-1,0,0)));
    check_do_intersect(L(p(-8,4,9), p(19,-8,-18)), Sph(p(10,4,9), p(0,1,0), p(0,0,1), p(-1,0,0)));

    // Two points
    check_do_intersect(L(p(-7,8,0), p(3,-4,1)), Sph(p(1,0,0), p(0,1,0), p(0,0,1), p(-1,0,0)));
    check_do_intersect(L(p(-1,-2,-3), p(2,1,3)), Sph(p(1,0,0), p(0,1,0), p(0,0,1), p(-1,0,0)));
  }

  void L_Tet()
  {
    std::cout << "Line - Tetrahedron\n";

    Tet tet(p(0,0,0), p(0,1,0), p(1,0,0), p(0,0,1));
    check_no_intersection(L(p(5,0,0), p(5,1,0)), tet);

    // on edge
    check_intersection(L(p(0,2,0), p(0,3,0)), tet,
                       S(p(0,0,0), p(0,1,0)));

    // through vertex and out from face
    check_intersection(L(p(0,1,0), P(0.25,0,0.25)), tet,
                       S(P(0,1,0), P(0.25,0,0.25)));

    // through a vertex only
    check_intersection(L(p(1,1,0), p(-1,1,0)), tet,
                       P(0,1,0));

    // through 2 faces
    check_intersection(L(P(0.25,0.25,0), P(0.25,0,0.25)), tet,
                       S(P(0.25,0,0.25), P(0.25,0.25,0)));

    // through one edge
    check_intersection(L(P(0.5,-0.5,0.5), P(0.5,0.5,0.5)), tet,
                       P(0.5,0,0.5));

    // in a single face through 2 edges
    check_intersection(L(P(0,0.5,0), P(0.5,0.5,0)), tet,
                             S(P(0,0.5,0), P(0.5, 0.5, 0)));

    // in a single face through 1 vertex and 1 edge
    check_intersection(L(p(-1,1,0), p(1,0,0)), tet,
                             S(P(0,0.5,0), P(1, 0, 0)));

    for(int i=0; i<N; ++i)
    {
      P tet0 = random_point(), tet1 = random_point(), tet2 = random_point(), tet3 = random_point();

      Tet tet(tet0, tet1, tet2, tet3);
      if(tet.is_degenerate())
        continue;

      if(tet.orientation() == CGAL::NEGATIVE)
        tet = Tet(tet1, tet0, tet2, tet3);

      P l0 = tet[0] - CGAL::cross_product(V(tet[0], tet[1]), V(tet[0], tet[2]));
      P l1 = tet[3] + CGAL::cross_product(V(tet[3], tet[1]), V(tet[3], tet[2]));

      assert(tet.has_on_unbounded_side(l0) && tet.has_on_unbounded_side(l1));

      const S s {l0, l1};
      if(s.is_degenerate())
        continue;

      if(CGAL::do_intersect(s, tet))
        check_intersection(L(l0, l1), tet, tet, s);
    }
  }

  void L_Tr()
  {
    std::cout << "Line - Triangle\n";

    // No intersection
    Tr tr(p(0,0,0), p(1,0,1), p(1,1,0));

    check_no_intersection(L(p(5,0,0), p(5,1,0)), tr);

    // on edge
    check_intersection(L(p(-1,0,-1), p(3,0,3)), tr,
                       S(p(0,0,0), p(1,0,1)));

    // through a vertex only
    check_intersection(L(p(-2,-1,0), p(4,2,0)), tr,
                       P(0,0,0));

    // through vertex and out from face
    Base::template check_intersection<S>(L(p(-3,-3,0), p(-4,-4,0)), Tr(p(0,0,0), p(1,0,0), p(0,1,0)));

    // through one edge
    Base::template check_intersection<P>(L(p(0,1,0), p(1,0,0)), tr);

    // within the face
    Base::template check_intersection<S>(L(P(0.25,0.25,0), P(0.3,0.5,0)), Tr(p(0,0,0), p(1,0,0), p(0,1,0)));

    for(int i=0; i<N; ++i)
    {
      P tr0 = random_point(), tr1 = random_point(), tr2 = random_point();

      P l0 = tr0 + CGAL::cross_product(V(tr0, tr1), V(tr0, tr2));
      P l1 = tr2 + CGAL::cross_product(V(tr2, tr1), V(tr2, tr0));

      S s{l0, l1};
      Tr tr{tr0, tr1, tr2};

      if(s.is_degenerate() || tr.is_degenerate())
        continue;

      if(CGAL::do_intersect(s, tr))
        check_intersection(L(l0, l1), tr, tr, s);
    }
  }

  void run()
  {
    std::cout << "3D Line Intersection tests\n";

    L_L();
    L_Pl();
    L_P();
    L_R();
    L_S();
    L_Sph();
    L_Tet();
    L_Tr();
  }
};

int main(int, char**)
{
  std::cout.precision(17);
  std::cerr.precision(17);

  CGAL::Random r;
  std::cout << "random seed = " << r.get_seed() << std::endl;

  std::cout << " |||||||| Test Simple_cartesian<double> ||||||||" << std::endl;
  Line_3_intersection_tester< CGAL::Simple_cartesian<double> >(r).run();

  std::cout << " |||||||| Test CGAL::Homogeneous<CGAL::MP_Float> ||||||||" << std::endl;
  Line_3_intersection_tester< CGAL::Homogeneous<CGAL::MP_Float> >(r).run();

  std::cout << " |||||||| Test EPICK ||||||||" << std::endl;
  Line_3_intersection_tester< CGAL::Epick >(r, true /*exact predicates*/).run();

  std::cout << " |||||||| Test EPECK ||||||||" << std::endl;
  Line_3_intersection_tester< CGAL::Epeck >(r, true /*exact predicates*/, true /*exact constructions*/).run();

  std::cout << " |||||||| Test CGAL::Homogeneous<CGAL::Epeck_ft> ||||||||" << std::endl;
  Line_3_intersection_tester< CGAL::Homogeneous<CGAL::Epeck_ft> >(r, true /*exact predicates*/, true /*exact constructions*/).run();

  std::cout << "OK!" << std::endl;
}
