// 3D intersection tests.

// We want to check that no division is performed for interface macro Do_intersect_3_RT
#define CGAL_NO_MPZF_DIVISION_OPERATOR

// Sphere_3_X with X < Sphere (lexicographically) are tested in other files

#include <CGAL/Intersections_3/Segment_3_Segment_3.h>
#include <CGAL/Intersections_3/Segment_3_Sphere_3.h>
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
struct Segment_3_intersection_tester
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
  Segment_3_intersection_tester(CGAL::Random& r,
                                const bool has_exact_p = false,
                                const bool has_exact_c = false)
    : Base(r, has_exact_p, has_exact_c)
  { }

public:
  void S_S()
  {
    // see segment_segment.cpp
  }

  void S_Sph()
  {
    std::cout << "Segment - Sphere\n";

    // No intersection
    check_do_not_intersect(S(p(4,1,9), p(4,6,-1)), Sph(p(9,2,4), 1));
    check_do_not_intersect(S(p(-2,3,4), p(1,2,-1)), Sph(p(9,2,4), 10000));

    check_do_intersect(S(p(0,1,0), p(2,1,0)), Sph(p(1,0,0), 1));
    check_do_intersect(S(p(-2,3,4), p(1,2,-1)), Sph(p(9,2,4), 100));
    check_do_intersect(S(p(-2,3,4), p(3,2,-5)), Sph(p(9,2,4), 100));

    for(int i=0; i<N; ++i)
    {
      P s0 = random_point(), s1 = random_point();
      P c = random_point();
      if(s0 == s1 || c == s0 || c == s1)
        continue;

      check_do_intersect(S(s0, s1), Sph(c, CGAL::squared_distance(c, s1)));

      FT sqr = this->r.get_int(this->m + 1, this->M);
      Sph sph(c, sqr);
      if(sph.oriented_side(s0) != sph.oriented_side(s1))
        check_do_intersect(sph, S(s0, s1));
    }
  }

  void S_Tet()
  {
    std::cout << "Segment - Tetrahedron\n";

    Tet tet(p(0,0,0), p(0,1,0), p(1,0,0), p(0,0,1));

    // No intersection
    check_no_intersection(S(p(5,0,0), p(5,1,0)), tet);

    // Point
    check_intersection(S(p(0,1,7), p(-7,-1,0)), Tet(p(0,1,7), p(9,4,-2), p(3,-7,2), p(7,-8,1)),
                       p(0,1,7));
    check_intersection(S(p(4,3,2), p(-7,-1,0)), Tet(p(0,1,6), p(8,5,-2), p(3,-7,2), p(7,-8,1)),
                       p(4,3,2));
    check_intersection(S(p(3,2,2), p(-7,-1,3)), Tet(p(4,1,6), p(2,-1,-2), p(3,6,2), p(7,-8,1)),
                       p(3,2,2));

    // Segment
    check_intersection(S(p(0,1,7), p(-7,-1,0)), Tet(p(0,1,7), p(-7,-1,0), p(3,-7,2), p(7,-8,1)),
                       S(p(0,1,7), p(-7,-1,0)));
    check_intersection(S(p(0,1,7), p(-8,-1,1)), Tet(p(0,1,7), p(-4,0,4), p(3,-7,2), p(7,-8,1)),
                       S(p(0,1,7), p(-4,0,4)));
    check_intersection(S(p(0,1,7), p(3,-7,2)), Tet(p(0,1,7), p(-4,0,4), p(3,-7,2), p(7,-8,1)),
                       S(p(0,1,7), p(3,-7,2)));

    FT fourth = FT(1)/FT(4);
    FT fifth = FT(1)/FT(5);
    FT tenth = FT(1)/FT(10);

    if(this->has_exact_c)
    {
      check_intersection(tet, S(p(0,2,0), p(0,-2,0)),
                         S(p(0,1,0), p(0,0,0)));
      check_intersection(tet, S(p(0,1,0), P(fourth,0,fourth)),
                         S(p(0,1,0), P(fourth,0,fourth)));
      check_intersection(tet, S(p(2,1,0), p(-2,1,0)), p(0,1,0));

      check_intersection(tet, S(P(tenth,tenth,tenth), P(fifth,fifth,fifth)),
                         S(P(tenth,tenth,tenth), P(fifth,fifth,fifth)));
      typename K::FT third = FT(1)/FT(3);
      check_intersection(tet, S(P(fourth,fourth,fourth), p(2,2,2)),
                         S(P(fourth,fourth,fourth), P(third, third, third)));
    }
    else
    {
      CGAL::intersection(tet, S(p(0,2,0), p(0,-2,0)));
      CGAL::intersection(tet, S(p(0,1,0), P(fourth,0,fourth)));
      CGAL::intersection(tet, S(p(2,1,0), p(-2,1,0)));
      CGAL::intersection(tet, S(P(tenth,tenth,tenth), P(fifth,fifth,fifth)));
      CGAL::intersection(tet, S(P(fourth,fourth,fourth), p(2,2,2)));
    }
  }

  void S_Tr()
  {
    std::cout << "Segment - Triangle\n";

    // No intersection
    check_no_intersection(S(p(1,4,7), p(2,3,9)), Tr(p(-1,9,8), p(-5,3,4), p(-0,8,9)));

    // Point
    check_intersection(S(p(4,2,7), p(1,3,7)), Tr(p(1,3,7), p(3,9,-1), p(-7,6,3)),
                       p(1,3,7));
    check_intersection(S(p(1,3,9), p(5,1,7)), Tr(p(3,2,8), p(-1,9,7), p(-2,2,8)),
                       p(3,2,8));

    // Segment
    check_intersection(S(p(3,2,7), p(1,-2,-9)), Tr(p(3,2,7), p(1,-2,-9), p(4,3,9)),
                       S(p(3,2,7), p(1,-2,-9)));
    check_intersection(S(p(3,2,7), p(1,-2,-9)), Tr(p(5,6,23), p(1,-2,-9), p(4,3,9)),
                       S(p(3,2,7), p(1,-2,-9)));
    check_intersection(S(p(3,2,7), p(1,-2,-9)), Tr(p(5,6,23), p(-1,-6,-25), p(4,3,9)),
                       S(p(3,2,7), p(1,-2,-9)));
    check_intersection(S(p(1,3,9), p(5,1,7)), Tr(p(4,0,6), p(2,4,10), p(-2,2,8)),
                       S(p(1,3,9), p(3,2,8)));

    Base::template check_intersection<S>(S(p(1,3,7), p(5,1,7)), Tr(p(4,0,7), p(2,4,7), p(-2,2,7)));
    Base::template check_intersection<S>(S(p(1,3,7), p(5,1,7)), Tr(p(4,0,7), p(2,4,7), p(1,2,7)));
  }

  void run()
  {
    std::cout << "3D Segment Intersection tests\n";

    S_S();
    S_Sph();
    S_Tet();
    S_Tr();
  }
};

int main(int, char**)
{
  std::cout.precision(17);
  std::cerr.precision(17);

  CGAL::Random r;
  std::cout << "random seed = " << r.get_seed() << std::endl;

  std::cout << " |||||||| Test Simple_cartesian<double> ||||||||" << std::endl;
  Segment_3_intersection_tester< CGAL::Simple_cartesian<double> >(r).run();

  std::cout << " |||||||| Test CGAL::Homogeneous<CGAL::MP_Float> ||||||||" << std::endl;
  Segment_3_intersection_tester< CGAL::Homogeneous<CGAL::MP_Float> >(r).run();

  std::cout << " |||||||| Test EPICK ||||||||" << std::endl;
  Segment_3_intersection_tester< CGAL::Epick >(r, true /*exact predicates*/).run();

  std::cout << " |||||||| Test EPECK ||||||||" << std::endl;
  Segment_3_intersection_tester< CGAL::Epeck >(r, true /*exact predicates*/, true /*exact constructions*/).run();

  std::cout << " |||||||| Test CGAL::Homogeneous<CGAL::Epeck_ft> ||||||||" << std::endl;
  Segment_3_intersection_tester< CGAL::Homogeneous<CGAL::Epeck_ft> >(r, true /*exact predicates*/, true /*exact constructions*/).run();

  std::cout << "OK!" << std::endl;
}
