// 3D intersection tests.

// We want to check that no division is performed for interface macro Do_intersect_3_RT
#define CGAL_NO_MPZF_DIVISION_OPERATOR

// Point_3_X with X < Point (lexicographically) are tested in other files

#include <CGAL/Intersections_3/Point_3_Point_3.h>
#include <CGAL/Intersections_3/Point_3_Ray_3.h>
#include <CGAL/Intersections_3/Point_3_Segment_3.h>
#include <CGAL/Intersections_3/Point_3_Sphere_3.h>
#include <CGAL/Intersections_3/Point_3_Tetrahedron_3.h>
#include <CGAL/Intersections_3/Point_3_Triangle_3.h>

#include <CGAL/Vector_3.h>

#include "intersection_test_helper.h"

#include <CGAL/MP_Float.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Homogeneous.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/representation_tags.h>

#include <iostream>
#include <cassert>
#include <type_traits>

template <typename K>
struct Point_3_intersection_tester
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
  Point_3_intersection_tester(CGAL::Random& r,
                              const bool has_exact_p = false,
                              const bool has_exact_c = false)
    : Base(r, has_exact_p, has_exact_c)
  { }

public:
  void P_P()
  {
    std::cout << "Point - Point\n";

    P p(0.99,0.99,0.99);
    check_intersection(p,p,p);

    P q = random_point();
    check_no_intersection(p,q);
  }

  void P_R()
  {
    std::cout << "Point - Ray\n";

    P pt(1,1,1);
    R ray(p(-1,-1,-1), p(3,3,3));
    check_intersection(ray, pt, pt);

    if(this->has_exact_c && std::is_same<typename K::Kernel_tag, CGAL::Cartesian_tag>::value)
    {
      for(int i=0; i<N; ++i)
      {
        P rp0 = random_point(), rp1 = random_point();
        if(rp0 != rp1)
        {
          R r(rp0, rp1);
          L l(rp0, rp1);
          P q = random_point();
          P pq = l.projection(q);

          assert(CGAL::collinear(rp0, rp1, pq));
          if(pq == rp0 || !CGAL::collinear_are_ordered_along_line(pq, rp0, rp1))
            check_intersection(r, pq, pq);
          else
            check_no_intersection(r, pq);
        }
      }
    }
  }

  void P_S()
  {
    std::cout << "Point - Segment\n";

    P pt(1,1,1);
    S s(p(-1,-1,-1), p(3,3,3));
    check_intersection(s, pt, pt);

    // @fixme Homogeneous' projection is broken
    if(this->has_exact_c && std::is_same<typename K::Kernel_tag, CGAL::Cartesian_tag>::value)
    {
      for(int i=0; i<N; ++i)
      {
        P rp0 = random_point(), rp1 = random_point();
        if(rp0 != rp1)
        {
          s = S(rp0, rp1);
          L l(rp0, rp1);
          P q = random_point(), pq = l.projection(q);
          assert(CGAL::collinear(rp0, rp1, pq));
          if(pq == rp0 || pq == rp1 || CGAL::collinear_are_strictly_ordered_along_line(rp0, pq, rp1))
            check_intersection(s, pq, pq);
          else
            check_no_intersection(s, pq);
        }
      }
    }
  }

  void P_Sph()
  {
    std::cout << "Point - Segment\n";

    P pt(1,1,1);
    Sph sph(p(0,0,0), 3);

    check_no_intersection(p(0,0,0), sph);

    check_intersection(pt, sph, pt);

    for(int i=0; i<N; ++i)
    {
      P c = random_point(), q = random_point();
      if(c != q)
      {
        sph = Sph(c, CGAL::squared_distance(c, q));
        check_intersection(sph, q, q);
      }
    }
  }

  void P_Tet()
  {
    std::cout << "Point - Tetrahedron\n";

    Tet t(p(0,0,0), p(10,0,0),p(10,10,0),p(0,0,10));
    P q0 = p(1, 1, 1), q1 = p(11,11,11);

    check_no_intersection(t, q1);

    for(int i=0;i<4;++i)
      check_intersection(t[i], t, t[i]);

    check_intersection(q0, t, q0);

    for(int i=0; i<N; ++i)
    {
      P tet0 = random_point(), tet1 = random_point(), tet2 = random_point(), tet3 = random_point();
      Tet tet(tet0, tet1, tet2, tet3);
      if(tet.is_degenerate())
        continue;

      if(tet.orientation() == CGAL::NEGATIVE)
        tet = Tet(tet1, tet0, tet2, tet3);

      P q = random_point();

      // not using bounded_side intentionally, otherwise it's a tautology
      if(Pl(tet[0], tet[1], tet[2]).oriented_side(q) != CGAL::NEGATIVE &&
         Pl(tet[0], tet[3], tet[1]).oriented_side(q) != CGAL::NEGATIVE &&
         Pl(tet[2], tet[1], tet[3]).oriented_side(q) != CGAL::NEGATIVE &&
         Pl(tet[0], tet[2], tet[3]).oriented_side(q) != CGAL::NEGATIVE)
        check_intersection(tet, q, q);
      else
        check_no_intersection(tet, q);
    }
  }

  void P_Tr()
  {
    std::cout << "Point - Triangle\n";

    check_no_intersection(Tr(p(0,0,0), p(1,0,0), p(0,1,0)), p(0,0,1));

    check_intersection(Tr(p(0,0,0), p(2,0,0), p(0,2,0)), p(1,0,0), p(1,0,0));
    check_intersection(Tr(p(0,0,0), p(2,0,0), p(0,2,0)), p(1,1,0), p(1,1,0));
    check_intersection(Tr(p(0,0,0), p(4,0,0), p(0,4,0)), p(1,1,0), p(1,1,0));

    if(this->has_exact_c && std::is_same<typename K::Kernel_tag, CGAL::Cartesian_tag>::value)
    {
      for(int i=0; i<N; ++i)
      {
        P tr0 = random_point(), tr1 = random_point(), tr2 = random_point();
        if(CGAL::collinear(tr0, tr1, tr2))
          continue;

        Tr tr(tr0, tr1, tr2);
        Pl pln(tr0, tr1, tr2);
        P q = random_point(), pq = pln.projection(q);
        if(tr.has_on(pq))
          check_intersection(tr, pq, pq);
        else
          check_no_intersection(tr, pq);
      }
    }
  }

  void run()
  {
    std::cout << "3D Point Intersection tests\n";

    P_P();
    P_R();
    P_S();
    P_Sph();
    P_Tet();
    P_Tr();
  }
};

int main(int, char**)
{
  std::cout.precision(17);
  std::cerr.precision(17);

  CGAL::Random r;
  std::cout << "random seed = " << r.get_seed() << std::endl;

  std::cout << " |||||||| Test Simple_cartesian<double> ||||||||" << std::endl;
  Point_3_intersection_tester< CGAL::Simple_cartesian<double> >(r).run();

  std::cout << " |||||||| Test CGAL::Homogeneous<CGAL::MP_Float> ||||||||" << std::endl;
  Point_3_intersection_tester< CGAL::Homogeneous<CGAL::MP_Float> >(r).run();

  std::cout << " |||||||| Test EPICK ||||||||" << std::endl;
  Point_3_intersection_tester< CGAL::Epick >(r, true /*exact predicates*/).run();

  std::cout << " |||||||| Test EPECK ||||||||" << std::endl;
  Point_3_intersection_tester< CGAL::Epeck >(r, true /*exact predicates*/, true /*exact constructions*/).run();

  std::cout << " |||||||| Test CGAL::Homogeneous<CGAL::Epeck_ft> ||||||||" << std::endl;
  Point_3_intersection_tester< CGAL::Homogeneous<CGAL::Epeck_ft> >(r, true /*exact predicates*/, true /*exact constructions*/).run();

  std::cout << "OK!" << std::endl;
}
