// 3D intersection tests.

// We want to check that no division is performed for interface macro Do_intersect_3_RT
#define CGAL_NO_MPZF_DIVISION_OPERATOR

// Plane_3_X with X < Plane (lexicographically) are tested in other files

#include <CGAL/Intersections_3/Plane_3_Plane_3.h>
#include <CGAL/Intersections_3/Plane_3_Point_3.h>
#include <CGAL/Intersections_3/Plane_3_Ray_3.h>
#include <CGAL/Intersections_3/Plane_3_Segment_3.h>
#include <CGAL/Intersections_3/Plane_3_Sphere_3.h>
#include <CGAL/Intersections_3/Plane_3_Tetrahedron_3.h>
#include <CGAL/Intersections_3/Plane_3_Triangle_3.h>

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
struct Plane_3_intersection_tester
  : public Intersection_3_tester<K>
{
  typedef Intersection_3_tester<K>  Base;

  typedef typename K::FT            FT;

  typedef CGAL::Circle_3<K>         C;
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
  Plane_3_intersection_tester(CGAL::Random& r,
                              const bool has_exact_p = false,
                              const bool has_exact_c = false)
    : Base(r, has_exact_p, has_exact_c)
  { }

public:
  void Pl_Pl()
  {
    std::cout << "Plane - Plane\n";

    // No intersection
    check_no_intersection(pl(2,1,3,4), pl(6,3,9,3));

    // Self
    check_intersection(pl(2,1,3,4), pl(2,1,3,4), pl(2,1,3,4));
    Base::template check_intersection<Pl>(Pl(p(0,1,8), p(4,3,7), p(2,9,1)), Pl(p(4,3,7), p(0,1,8), p(2,9,1)));

    // Line (directed...)
    check_intersection(pl(0,0,1,0), pl(0,1,0,0), L(P(0,0,0), P(-1,0,0)), false /*do_opposite*/);
    check_intersection(pl(0,1,0,0), pl(0,0,1,0), L(P(-85,0,0), P(0,0,0)), false /*do_opposite*/);
    check_intersection(pl(2,3,7,5), pl(9,7,1,3), L(P(2,-3,0), P(-3908,5182,-1105)), false /*do_opposite*/);
    check_intersection(pl(9,7,1,3), pl(2,3,7,5), L(P(-3908,5182,-1105), P(2,-3,0)), false /*do_opposite*/);

    // generic
    for(int i=0; i<N; ++i)
    {
      P pl0 = random_point(), pl1 = random_point(), pl2 = random_point();
      if(CGAL::collinear(pl0, pl1, pl2))
        continue;

      Pl pln(pl0, pl1, pl2);
      P q = random_point();
      Pl pln2(q,pl2,pl1);
      if(pln.has_on(q))
      {
        check_intersection(pln, pln2, pln, false /*do_opposite*/);
      }
      else
      {
        auto res = CGAL::intersection(pln, pln2);
        assert(res);
        L l;
        assert(assign(l, res));

        if(this->has_exact_c)
        {
          assert(l.has_on(pl1));
          assert(l.has_on(pl2));
        }
      }
    }

    // Plane
    Base::template check_intersection<Pl>(pl(0,0,1,1), pl(0,0,3,3));
    Base::template check_intersection<Pl>(pl(2,1,3,4), pl(6,3,9,12));
  }

  void Pl_P()
  {
    std::cout << "Plane - Point\n";

    // no intersection
    check_no_intersection(Pl(p(5,3,2), p(1,1,1), p(0,0,0)), P(99,98,99));

    check_intersection(Pl(p(5,3,2), p(1,1,1), p(0,0,0)), P(99,99,99),
                       P(99,99,99));

    // Homogeneous projections are broken
    if(this->has_exact_c && std::is_same<typename K::Kernel_tag, CGAL::Cartesian_tag>::value)
    {
      for(int i=0; i<N; ++i)
      {
        const P pl0 = random_point(), pl1 = random_point(), pl2 = random_point();
        if(CGAL::collinear(pl0, pl1, pl2))
          continue;

        const Pl pln(pl0, pl1, pl2);
        const P& plnpt = pln.point();
        check_intersection(pln, plnpt, plnpt);

        const P q = random_point();
        const P pq = pln.projection(q);
        check_intersection(pln, pq, pq);
        if(q != pq)
          check_no_intersection(pln, q);
      }
    }
  }

  void Pl_R()
  {
    std::cout << "Plane - Ray\n";

    // No intersection
    check_no_intersection  (pl(1,1,1,0), R(p(1,1,1), p(2,3,4)));
    check_no_intersection  (pl(0,0,1,-2), R(p(1,1,1), p(2,3,1)));
    check_no_intersection  (pl(1,2,4,7), R(p(1,1,1), p(2,3,4)));

    // Point
    check_intersection     (pl(1,0,1,3), R(p(1,1,1), p(2,3,-1)),
                            p(6,11,-9));
    check_intersection     (pl(0,0,1,0), R(p(1,1,-1), p(-1,-1,1)),
                            P(0,0,0));
    check_intersection     (pl(0,0,1,0), R(p(7,1,0), p(83,1,-4)),
                            p(7,1,0));
    check_intersection     (pl(0,0,1,0), R(p(12,6,-4), p(7,25,0)),
                            P(7,25,0));

    // Ray
    Base::template check_intersection<R>  (pl(0,0,1,-1), R(p(1,1,1), p(2,3,1)));
    check_intersection(pl(0,0,1,-1), R(p(1,1,1), p(2,3,1)), R(p(1,1,1), p(2,3,1)));
  }

  void Pl_S()
  {
    std::cout << "Plane - Segment\n";

    // No intersection
    check_no_intersection  (pl(1,1,1,0), S(p(1,1,1), p(2,3,4)));
    check_no_intersection  (pl(0,0,1,-2), S(p(1,1,1), p(2,3,1)));
    check_no_intersection  (pl(1,0,1,3), S(p(1,1,1), p(2,3,-1)));
    check_no_intersection  (pl(1,2,4,7), S(p(1,1,1), p(2,3,4)));

    // Point
    check_intersection     (pl(0,0,1,0), S(p(1,1,-1), p(-1,-1,1)),
                            P(0,0,0));
    check_intersection     (pl(0,0,1,0), S(p(7,1,0), p(83,1,-4)),
                            p(7,1,0));
    check_intersection     (pl(0,0,1,0), S(p(12,6,-4), p(7,25,0)),
                            p(7,25,0));

    // Segment
    Base::template check_intersection<S>(pl(0,0,1,-1), S(p(1,1,1), p(2,3,1)));
  }

  // see also circle_other.cpp
  void Pl_Sph()
  {
    std::cout << "Plane - Sphere\n";

    // No intersection
    check_no_intersection(pl(1,1,1,0), Sph(p(4,3,1), 1));

    // Point
    check_intersection  (pl(0,10,0,0), Sph(p(4,-2,1), 4), p(4,0,1));
    check_intersection  (pl(1,1,0,0), Sph(p(4,4,-10), 32), p(0,0,-10));

    // generic
    check_intersection  (pl(0,0,1,-1), Sph(p(0,0,2), 4), C(pl(0,0,1,-1), Sph(p(0,0,2), 4)));

    if(this->has_exact_c)
    {
      for(int i=0; i<N; ++i)
      {
        const P pl0 = random_point(), pl1 = random_point(), pl2 = random_point();
        if(CGAL::collinear(pl0, pl1, pl2))
          continue;

        const Pl pln(pl0, pl1, pl2);

        const P sph4 = random_point();
        if(pln.has_on(sph4))
          continue;

        Sph sph(pl0, pl1, pl2, sph4);
        check_intersection(pln, sph, C(pln, sph));

        const P sph5 = random_point();
        if(pln.has_on(sph5))
          continue;

        Tet tet{pl0, pl1, sph5, sph4};
        if(tet.is_degenerate())
          continue;

        sph = Sph(pl0, pl1, sph5, sph4);
        check_intersection(pln, sph, C(pln, sph));

        const P sph6 = random_point();
        if(pln.has_on(sph6))
          continue;

        // The sphere constructor will assert on degenerate positions
        tet = Tet{pl0, sph6, sph5, sph4};
        if(tet.is_degenerate())
          continue;

        sph = Sph(pl0, sph6, sph5, sph4);
        check_intersection(pln, sph, C(pln, sph));
      }
    }
  }

  void Pl_Tet()
  {
    std::cout << "Plane - Tetrahedron\n";

    Tet tet(p(0,0,0), p(0,1,0), p(1,0,0), p(0,0,1));

    // No intersection
    check_no_intersection(tet, Pl(p(2,0,0),
                                  p(2,1,0),
                                  p(2,0,1)));

    // Point
    check_intersection(tet, Pl(p(-1,1,12), p(0,1,-50), P(0.5,1,-0.5)),
                       p(0,1,0));

    // Segment
    check_intersection(tet, Pl(p(0,1,0), p(1,0,0), P(0.5,0,-0.5)),
                       S(p(0,1,0), p(1,0,0)));

    // Triangle
    check_intersection(tet, Pl(p(0,2,0), p(0,0,0), p(0,0,1)),
                       Tr(p(0,0,0), p(0,1,0), p(0,0,1)));

    if(this->has_exact_c)
    {
      check_intersection(tet, Pl(P(0,0.5,0), P(1,0.5,-5), P(0.5,0.5,0.5)),
                         Tr(P(0,0.5,0), P(0.5,0.5,0),  P(0,0.5,0.5)));

      Pl pln(P(0,0.9,0), P(0.9,0,0), P(0.9,0.01,0.06));

      // Don't have the right values to test further.
      auto res = CGAL::intersection(tet, pln);
      const std::vector<P>* poly = std::get_if<std::vector<P> >(&*res);
      assert(poly != nullptr);
      assert(poly->size() == 4);
    }

    for(int i=0; i<N; ++i)
    {
      const P pl0 = random_point(), pl1 = random_point(), pl2 = random_point();
      if(CGAL::collinear(pl0, pl1, pl2))
        continue;

      Pl pln(pl0, pl1, pl2);

      const P tet0 = random_point(), tet1 = random_point(), tet2 = random_point(), tet3 = random_point();
      Tet t(tet0, tet1, tet2, tet3);
      if(t.is_degenerate())
        continue;

      std::set<CGAL::Orientation> os;
      os.insert(pln.oriented_side(tet0));
      os.insert(pln.oriented_side(tet1));
      os.insert(pln.oriented_side(tet2));
      os.insert(pln.oriented_side(tet3));
      if(os.size() == 1 && (*os.begin() != CGAL::ZERO))
        check_no_intersection(pln, t);
      else
        check_do_intersect(pln, t);
    }
  }

  void Pl_Tr()
  {
    std::cout << "Plane - Triangle\n";

    // No intersection

    // Point
    check_intersection(Pl(p(0,0,0), p(12,0,0), p(0,11,0)), Tr(p(0,0,0), p(1,0,1), p(0,1,1)),
                       p(0,0,0));

    // Segment
    check_intersection(Pl(p(0,0,0), p(12,0,0), p(0,11,0)), Tr(p(0,0,0), p(1,0,0), p(0,1,1)),
                       S(p(0,0,0), p(1,0,0)));
    check_intersection(Pl(p(0,0,0), p(12,0,0), p(0,11,0)), Tr(p(1,0,-1), p(-1,0,-1), p(0,0,1)),
                       S(P(0.5,0,0), P(-0.5,0,0)));

    // Triangle
    check_intersection(Pl(p(0,0,0), p(12,0,0), p(0,11,0)), Tr(p(0,0,0), p(1,0,0), p(0,1,0)),
                       Tr(p(0,0,0), p(1,0,0), p(0,1,0)));

    for(int i=0; i<N; ++i)
    {
      const P pl0 = random_point(), pl1 = random_point(), pl2 = random_point();
      if(CGAL::collinear(pl0, pl1, pl2))
        continue;

      Pl pln(pl0, pl1, pl2);
      Tr tr(pl1, pl0, pl2);
      check_intersection(pln, tr, tr);

      const P pl3 = random_point(), pl4 = random_point();
      if(CGAL::collinear(pl0, pl1, pl3))
        continue;

      tr = Tr(pl0, pl1, pl3);
      if(pln.has_on(pl3))
      {
        check_intersection(pln, tr, tr);
     }
      else
      {
        check_intersection(pln, tr, S(pl0, pl1));

        if(CGAL::collinear(pl0, pl3, pl4))
          continue;

        if(pln.oriented_side(pl3) == pln.oriented_side(pl4))
          check_intersection(pln, Tr(pl0, pl3, pl4), pl0);
        else
          Base::template check_intersection<S>(pln, Tr(pl0, pl3, pl4));
      }
    }
  }

  void Pl_Pl_Pl()
  {
    std::cout << "Plane - Plane - Plane\n";

    Pl pl1(1,0,0,0);
    Pl pl2(0,1,0,0);
    Pl pl3(0,0,1,0);
    Pl pl4(1,0,0,1); // pl4 is parallel to pl1.

    // No intersection
    check_no_intersection(pl1, pl2, pl4);
    check_no_intersection(pl1, pl4, pl2);
    check_no_intersection(pl4, pl2, pl1);

    check_no_intersection(Pl(p(0,0,0),p(1,1,1),p(1,1,0)),
                          Pl(p(0,0,0),p(1,-1,1),p(1,-1,-2)),
                          Pl(1,0,0,-5));

    // Intersection in a line.
    Pl pl5(1,1,0,0); // pl1, pl2, pl5 intersect in the line l.
    L l;
    assert(CGAL::do_intersect(pl1, pl2, pl5));
    CGAL::Object o4 = CGAL::intersection(pl1, pl2, pl5);
    assert(assign(l, o4));
    assert(l == L(P(0,0,0), P(0,0,1)));

    // Intersection in a plane.
    assert(CGAL::do_intersect(pl1, pl1, pl1));
    CGAL::Object o5 = CGAL::intersection(pl1, pl1, pl1);
    Pl pln;
    assert(assign(pln, o5));
    assert(pln == pl1);

    // Generic intersection.
    assert(CGAL::do_intersect(pl1, pl2, pl3));
    CGAL::Object o = CGAL::intersection(pl1, pl2, pl3);
    P pt;
    assert(assign(pt, o));
    assert(pt == P(0,0,0));
  }

  void run()
  {
    std::cout << "3D Plane Intersection tests\n";

    Pl_Pl();
    Pl_P();
    Pl_R();
    Pl_S();
    Pl_Sph();
    Pl_Tet();
    Pl_Tr();

    Pl_Pl_Pl();
  }
};

int main(int, char**)
{
  std::cout.precision(17);
  std::cerr.precision(17);

  CGAL::Set_ieee_double_precision pfr;

  CGAL::Random r;
  std::cout << "random seed = " << r.get_seed() << std::endl;

  std::cout << " |||||||| Test Simple_cartesian<double> ||||||||" << std::endl;
  Plane_3_intersection_tester< CGAL::Simple_cartesian<double> >(r).run();

  // Homogeneous is broken for projection and Pln-Sphere
//  std::cout << " |||||||| Test CGAL::Homogeneous<CGAL::MP_Float> ||||||||" << std::endl;
//  Plane_3_intersection_tester< CGAL::Homogeneous<CGAL::MP_Float> >(r).run();

  std::cout << " |||||||| Test EPICK ||||||||" << std::endl;
  Plane_3_intersection_tester< CGAL::Epick >(r, true /*exact predicates*/).run();

  std::cout << " |||||||| Test EPECK ||||||||" << std::endl;
  Plane_3_intersection_tester< CGAL::Epeck >(r, true /*exact predicates*/, true /*exact constructions*/).run();

  // See above
//  std::cout << " |||||||| Test CGAL::Homogeneous<CGAL::Epeck_ft> ||||||||" << std::endl;
//  Plane_3_intersection_tester< CGAL::Homogeneous<CGAL::Epeck_ft> >(r, true /*exact predicates*/, true /*exact constructions*/).run();

  std::cout << "OK!" << std::endl;
}
