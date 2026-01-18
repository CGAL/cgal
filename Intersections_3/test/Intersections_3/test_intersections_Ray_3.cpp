// 3D intersection tests.

// We want to check that no division is performed for interface macro Do_intersect_3_RT
#define CGAL_NO_MPZF_DIVISION_OPERATOR

// Sphere_3_X with X < Sphere (lexicographically) are tested in other files

#include <CGAL/Intersections_3/Ray_3_Ray_3.h>
#include <CGAL/Intersections_3/Ray_3_Segment_3.h>
#include <CGAL/Intersections_3/Ray_3_Sphere_3.h>
#include <CGAL/Intersections_3/Ray_3_Tetrahedron_3.h>
#include <CGAL/Intersections_3/Ray_3_Triangle_3.h>
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
struct Ray_3_intersection_tester
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
  Ray_3_intersection_tester(CGAL::Random& r,
                            const bool has_exact_p = false,
                            const bool has_exact_c = false)
    : Base(r, has_exact_p, has_exact_c)
  { }

public:
  void R_R()
  {
    std::cout << "Ray - Ray\n";

    // No intersection
    check_no_intersection  (R(p(0,0,0), p(1,0,0)), R(p(-2,0,0), p(-3,0,0)));
    check_no_intersection  (R(p(0,0,0), p(1,0,0)), R(p(1,-1,0), p(1,-2,0)));
    check_no_intersection  (R(p(0,0,0), p(1,0,0)), R(p(0,-1,0), p(0,-2,0)));
    check_no_intersection  (R(p(0,0,0), p(1,0,0)), R(p(0,1,0), p(0,2,0)));
    check_no_intersection  (R(p(0,0,0), p(1,0,0)), R(p(-1,0,0), p(-1,-1,0)));

    check_no_intersection  (R(p(1,-1,0), p(2,0,0)), R(p(2,-1,0), p(3,0,0)));
    check_no_intersection  (R(p(1,-1,0), p(2,0,0)), R(p(2,-1,0), p(3,-1,0)));
    check_no_intersection  (R(p(1,-1,0), p(2,0,0)), R(p(2,-1,0), p(3,-2,0)));
    check_no_intersection  (R(p(0,0,0), p(0,1,0)), R(p(0,-1,0), p(1,-1,0)));
    check_no_intersection  (R(p(0,0,0), p(0,1,0)), R(p(-1,-3,0),p(2,0,0)));
    check_no_intersection  (R(p(0,0,0), p(0,1,0)), R(p(-2,-4,0),p(-1,-3,0)));
    check_no_intersection  (R(p(0,0,0), p(0,1,0)), R(p(1,-1,0), p(2,0,0)));

    // Point
    check_intersection     (R(p(0,0,0), p(1,0,0)), R(p(0,0,0), p(-1,0,0)),
                            p(0,0,0));
    check_intersection     (R(p(0,0,0), p(1,0,0)), R(p(0,0,0), p(0,1,0)),
                            p(0,0,0));
    check_intersection     (R(p(0,0,0), p(1,0,0)), R(p(5,0,0), p(5,1,0)),
                            p(5,0,0));
    check_intersection     (R(p(0,0,0), p(1,0,0)), R(p(1,0,0), p(1,1,0)),
                            p(1,0,0));
    check_intersection     (R(p(0,0,0), p(1,0,0)), R(p(1,-1,0), p(1,1,0)),
                            p(1,0,0));
    check_intersection     (R(p(0,0,0), p(1,0,0)), R(p(1,-2,0), p(1,-1,0)),
                            p(1,0,0));

    check_intersection     (R(p(0,0,0), p(1,0,0)), R(p(1,-2,0), p(1,-1,0)),
                            p(1,0,0));


    // Segment
    check_intersection     (R(p(0,0,0), p(1,0,0)), R(p(2,0,0), p(-3,0,0)),
                            S(p(0,0,0), p(2,0,0)), false);
    check_intersection     (R(p(2,0,0), p(-3,0,0)), R(p(0,0,0), p(1,0,0)),
                            S(p(2,0,0), p(0,0,0)), false);

    // Ray
    check_intersection     (R(p(0,0,0), p(1,0,0)), R(p(-1,0,0), p(3,0,0)),
                            R(p(0,0,0), p(1,0,0)));
    check_intersection     (R(p(0,0,0), p(1,0,0)), R(p(2,0,0), p(3,0,0)),
                            R(p(2,0,0), p(6,0,0)));
  }

  void R_S()
  {
    std::cout << "Ray - Segment\n";

    // No intersection
    check_no_intersection  (S(p(-3,0,0), p(3,0,0)), R(p(0,-1,0), p(0,-2,0)));
    check_no_intersection  (S(p(-3,0,0), p(3,0,0)), R(p(-3,-1,0), p(-3,-2,0)));
    check_no_intersection  (S(p(-3,0,0), p(3,0,0)), R(p(3,-1,0), p(3,-2,0)));
    check_no_intersection  (S(p(-3,0,0), p(3,0,0)), R(p(4,0,0), p(5,0,0)));
    check_no_intersection  (S(p(-3,0,0), p(3,0,0)), R(p(-5,0,0), p(-8,0,0)));

    // Point
    check_intersection     (S(p(-3,0,0), p(3,0,0)), R(p(0,-2,0), p(0,-1,0)),
                            p(0,0,0));
    check_intersection     (S(p(-3,0,0), p(3,0,0)), R(p(-3,-2,0), p(-3,-1,0)),
                            p(-3,0,0));
    check_intersection     (S(p(-3,0,0), p(3,0,0)), R(p(3,-2,0), p(3,-1,0)),
                            p(3,0,0));
    check_intersection     (S(p(-3,0,0), p(3,0,0)), R(p(0,0,0), p(0,-1,0)),
                            p(0,0,0));
    check_intersection     (S(p(-3,0,0), p(3,0,0)), R(p(-3,0,0), p(-3,-1,0)),
                            p(-3,0,0));
    check_intersection     (S(p(-3,0,0), p(3,0,0)), R(p(3,0,0), p(3,-1,0)),
                            p(3,0,0));
    check_intersection     (S(p(-3,0,0), p(3,0,0)), R(p(0,-2,0), p(0,0,0)),
                            p(0,0,0));
    check_intersection     (S(p(-3,0,0), p(3,0,0)), R(p(-3,-2,0), p(-3,0,0)),
                            p(-3,0,0));
    check_intersection     (S(p(-3,0,0), p(3,0,0)), R(p(3,-2,0), p(3,0,0)),
                            p(3,0,0));
    check_intersection     (S(p(-3,0,0), p(3,0,0)), R(p(-3,0,0), p(-4,0,0)),
                            p(-3,0,0));
    check_intersection     (S(p(-3,0,0), p(3,0,0)), R(p(3,0,0), p(4,0,0)),
                            p(3,0,0));

    // Segment
    check_intersection     (S(p(-3,0,0), p(3,0,0)), R(p(-4,0,0), p(-2,0,0)),
                            S(p(-3,0,0), p(3,0,0)));
    check_intersection     (S(p(-3,0,0), p(3,0,0)), R(p(-4,0,0), p(7,0,0)),
                            S(p(-3,0,0), p(3,0,0)));
    check_intersection     (S(p(-3,0,0), p(3,0,0)), R(p(-3,0,0), p(-2,0,0)),
                            S(p(-3,0,0), p(3,0,0)));
    check_intersection     (S(p(-3,0,0), p(3,0,0)), R(p(-3,0,0), p(7,0,0)),
                            S(p(-3,0,0), p(3,0,0)));
    check_intersection     (S(p(-3,0,0), p(3,0,0)), R(p(0,0,0), p(-2,0,0)),
                            S(p(0,0,0), p(-3,0,0)));
    check_intersection     (S(p(-3,0,0), p(3,0,0)), R(p(0,0,0), p(2,0,0)),
                            S(p(0,0,0), p(3,0,0)));
  }

  void R_Sph()
  {
    std::cout << "Ray - Sphere\n";

    // No intersection
    check_do_not_intersect(R(p(0,1,7), p(1,2,1)), Sph(p(0,8,4), 4));

    for(int i=0; i<N; ++i)
    {
      P c = random_point(), q = random_point();
      while(c==q)
        q = random_point();
      Sph sph(c, CGAL::squared_distance(c, q));

      // single point
      if(this->has_exact_c)
      {
        Pl pln(q, V(c, q)); // plane tangent to the sphere at q
        V v1 = pln.base1();
        check_do_intersect(R(q, q + v1), sph);
      }

      // maybe two points
      P r = random_point();
      if(q == r)
        continue;

      check_do_intersect(R(q, r), sph);
    }
  }

  void R_Tet()
  {
    std::cout << "Ray - Tetrahedron\n";

    Tet tet(p(0,0,0), p(0,1,0), p(1,0,0), p(0,0,1));

    check_no_intersection (R(p(5,0,0), p(5,1,0)), tet);

    check_intersection (R(p(0,2,0), p(0,-2,0)), tet,
                        S(p(0,1,0), p(0,0,0)));
    check_intersection (R(p(0,1,0), P(0.25,0,0.25)), tet,
                        S(p(0,1,0), P(0.25,0,0.25)));
    check_intersection (R(p(2,1,0), p(-2,1,0)), tet,
                        p(0,1,0));

    typename K::FT third = FT(1)/FT(3);
    typename K::FT fifth = FT(1)/FT(5);
    typename K::FT tenth = FT(1)/FT(10);

    if(this->has_exact_c)
    {
      check_intersection (tet, R(P(tenth,tenth,tenth), P(fifth,fifth,fifth)),
                          S(P(tenth,tenth,tenth), P(third, third, third)));
      check_intersection (tet, R(P(0.25,0.25,0.25), p(2,2,2)),
                          S(P(0.25,0.25,0.25), P(third, third, third)));
    }
    else
    {
      CGAL::intersection(tet, R(p(0,2,0), p(0,-2,0)));
      CGAL::intersection(tet, R(p(0,1,0), P(0.25,0,0.25)));
      CGAL::intersection(tet, R(p(2,1,0), p(-2,1,0)));
      CGAL::intersection(tet, R(P(tenth,tenth,tenth), P(fifth,fifth,fifth)));
      CGAL::intersection(tet, R(P(0.25,0.25,0.25), p(2,2,2)));
    }
  }

  void R_Tr()
  {
    std::cout << "Ray - Triangle\n";

    check_no_intersection(R(p(0,0,0), p(1,1,1)), Tr(p(-4, 2, 2), p(3, 7, -1), p(0, 9, -4)));
    check_no_intersection(R(p(0,0,0), p(1,1,1)), Tr(p(-1, -1, -1), p(3, 7, -1), p(0, 9, -4)));

    // Point
    check_intersection(R(p(0,0,0), p(1,1,1)), Tr(p(0, 0, 0), p(3, 7, -1), p(0, 9, -4)),
                       p(0, 0, 0));
    check_intersection(R(p(0,0,0), p(1,1,1)), Tr(p(2, 2, 2), p(3, 7, -1), p(0, 9, -4)),
                       p(2, 2, 2));
    check_intersection(R(p(0,0,0), p(4,-2,6)), Tr(p(2, -1, 3), p(3, 7, -1), p(0, 9, -4)),
                       p(2, -1, 3));
    check_intersection(R(p(0,0,0), p(2,-1,3)), Tr(p(4, -2, 6), p(3, 7, -2), p(0, 8, -4)),
                       p(4, -2, 6));

    Base::template check_intersection<P>(R(p(0,0,0), p(1,1,1)), Tr(p(2, 4, 2), p(5, -7, -1), p(0, -9, 4)));

    // Segment
    check_intersection(R(p(0,0,1), p(1,1,1)), Tr(p(3, 3, 1), p(7, 7, 1), p(0, -9, -4)),
                       S(p(3, 3, 1), p(7,7,1)));
    check_intersection(R(p(0,0,1), p(1,1,1)), Tr(p(7, 7, 1), p(3, 3, 1), p(0, -9, -4)),
                       S(p(3, 3, 1), p(7,7,1)));
    check_intersection(R(p(3,3,1), p(7,7,1)), Tr(p(4, 4, 1), p(2, 2, 1), p(0, -9, -4)),
                       S(p(3, 3, 1), p(4,4,1)));
  }

  void run()
  {
    std::cout << "3D Ray Intersection tests\n";

    R_R();
    R_S();
    R_Sph();
    R_Tet();
    R_Tr();
  }
};

int main(int, char**)
{
  std::cout.precision(17);
  std::cerr.precision(17);

  CGAL::Random r;
  std::cout << "random seed = " << r.get_seed() << std::endl;

  std::cout << " |||||||| Test Simple_cartesian<double> ||||||||" << std::endl;
  Ray_3_intersection_tester< CGAL::Simple_cartesian<double> >(r).run();

  std::cout << " |||||||| Test CGAL::Homogeneous<CGAL::MP_Float> ||||||||" << std::endl;
  Ray_3_intersection_tester< CGAL::Homogeneous<CGAL::MP_Float> >(r).run();

  std::cout << " |||||||| Test EPICK ||||||||" << std::endl;
  Ray_3_intersection_tester< CGAL::Epick >(r, true /*exact predicates*/).run();

  std::cout << " |||||||| Test EPECK ||||||||" << std::endl;
  Ray_3_intersection_tester< CGAL::Epeck >(r, true /*exact predicates*/, true /*exact constructions*/).run();

  std::cout << " |||||||| Test CGAL::Homogeneous<CGAL::Epeck_ft> ||||||||" << std::endl;
  Ray_3_intersection_tester< CGAL::Homogeneous<CGAL::Epeck_ft> >(r, true /*exact predicates*/, true /*exact constructions*/).run();

  std::cout << "OK!" << std::endl;
}
