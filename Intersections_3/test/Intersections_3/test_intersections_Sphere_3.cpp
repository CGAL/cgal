// 3D intersection tests.

// We want to check that no division is performed for interface macro Do_intersect_3_RT
#define CGAL_NO_MPZF_DIVISION_OPERATOR

// Sphere_3_X with X < Sphere (lexicographically) are tested in other files

#include <CGAL/Intersections_3/Sphere_3_Sphere_3.h>
#include <CGAL/Intersections_3/Sphere_3_Tetrahedron_3.h>
#include <CGAL/Intersections_3/Sphere_3_Triangle_3.h>

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
struct Sphere_3_intersection_tester
  : public Intersection_3_tester<K>
{
  typedef Intersection_3_tester<K>  Base;

  typedef typename K::FT            FT;

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
  Sphere_3_intersection_tester(CGAL::Random& r,
                               const bool has_exact_p = false,
                               const bool has_exact_c = false)
    : Base(r, has_exact_p, has_exact_c)
  { }

public:
  void Sph_Sph()
  {
    std::cout << "Sphere - Sphere\n";

    // No intersection
    check_no_intersection(Sph(p(0,0,0), 4), Sph(p(6,3,8), 9));

    for(int i=0; i<N; ++i)
    {
      P c0 = random_point(), c1 = random_point();

      const FT r0 = this->r.get_int(this->m + 1, this->M), r1 = this->r.get_int(this->m + 1, this->M);
      if(r0 == r1)
        check_intersection(Sph(c0, r0), Sph(c0, r1), Sph(c0, r0));
      else
        check_no_intersection(Sph(c0, r0), Sph(c0, r1));

      const FT sqd = CGAL::squared_distance(c0, c1);
      if(sqd <= 8)
        continue;

      const int r = int(CGAL::to_double((sqd / 4)) - 1);
      check_no_intersection(Sph(c0, r), Sph(c1, r));
    }
  }

  void Sph_Tet()
  {
    std::cout << "Sphere - Tetrahedron\n";

    // No intersection
    check_do_not_intersect(Sph(p(0,0,0), 3), Tet(p(6,3,8), p(7,2,8), p(1,9,7), p(2,3,4)));
    check_do_not_intersect(Sph(p(0,0,0), 10000), Tet(p(6,3,8), p(7,2,8), p(1,9,7), p(2,3,4)));

    // Point
    check_do_intersect(Sph(p(0,0,0), 4), Tet(p(2,0,0), p(9,3,7), p(4,9,7), p(8,-3,4))); // vertex
    check_do_intersect(Sph(p(0,0,0), 4), Tet(p(2,2,0), p(2,-2,0), p(4,9,7), p(8,-3,4))); // edge
    check_do_intersect(Sph(p(0,0,0), 4), Tet(p(2,2,3), p(2,-2,-3), p(4,2,-3), p(8,-3,4))); // face

    for(int i=0; i<N; ++i)
    {
      P tet0 = random_point(), tet1 = random_point(), tet2 = random_point(), tet3 = random_point();
      Tet tet(tet0, tet1, tet2, tet3);
      if(tet.is_degenerate())
        continue;

      P c = random_point();
      FT sqm, sqM;
      sqm = (std::min)((std::min)(CGAL::squared_distance(c, tet[0]), CGAL::squared_distance(c, tet[1])),
                       (std::min)(CGAL::squared_distance(c, tet[2]), CGAL::squared_distance(c, tet[3])));
      sqM = (std::max)((std::max)(CGAL::squared_distance(c, tet[0]), CGAL::squared_distance(c, tet[1])),
                       (std::max)(CGAL::squared_distance(c, tet[2]), CGAL::squared_distance(c, tet[3])));
      FT r = (sqm + sqM) / 2;
      if(r == 0)
        continue;

      check_do_intersect(Sph(c, r), tet);
    }
  }

  void Sph_Tr()
  {
    std::cout << "Sphere - Triangle\n";

    // No intersection
    check_do_not_intersect(Sph(p(0,0,0), 3), Tr(p(6,3,8), p(7,2,8), p(1,9,7)));
    check_do_not_intersect(Sph(p(0,0,0), 10000), Tr(p(6,3,8), p(7,2,8), p(1,9,7)));

    check_do_intersect(Sph(p(1, 2, 3), 100), Tr(p(-15, -18, 3), p(-15, 13, 1), p(15, 0, 2)));
    check_do_intersect(Sph(p(1, 2, 3), 100), Tr(p(-150, -180, 3), p(-150, 130, 1), p(150, 0, 2)));

    // Generic
    for(int i=0; i<N; ++i)
    {
      const P tr0 = random_point(), tr1 = random_point(), tr2 = random_point();
      const Tr tr(tr0, tr1, tr2);
      if(tr.is_degenerate())
        continue;

      const P c = random_point();

      if(this->has_exact_c)
        check_do_intersect(Sph(c, CGAL::squared_distance(c, tr)), tr);

      const FT r = this->r.get_int(this->m + 1, this->M);
      Sph sph(c, r);

      if(sph.oriented_side(tr0) != sph.oriented_side(tr1) ||
         sph.oriented_side(tr0) != sph.oriented_side(tr2) ||
         sph.oriented_side(tr2) != sph.oriented_side(tr1))
        check_do_intersect(sph, tr);
    }
  }

  void run()
  {
    std::cout << "3D Sphere Intersection tests\n";

    Sph_Sph();
    Sph_Tet();
    Sph_Tr();
  }
};

int main(int, char**)
{
  std::cout.precision(17);
  std::cerr.precision(17);

  CGAL::Random r;
  std::cout << "random seed = " << r.get_seed() << std::endl;

  std::cout << " |||||||| Test Simple_cartesian<double> ||||||||" << std::endl;
  Sphere_3_intersection_tester< CGAL::Simple_cartesian<double> >(r).run();

  std::cout << " |||||||| Test CGAL::Homogeneous<CGAL::MP_Float> ||||||||" << std::endl;
  Sphere_3_intersection_tester< CGAL::Homogeneous<CGAL::MP_Float> >(r).run();

  std::cout << " |||||||| Test EPICK ||||||||" << std::endl;
  Sphere_3_intersection_tester< CGAL::Epick >(r, true /*exact predicates*/).run();

  std::cout << " |||||||| Test EPECK ||||||||" << std::endl;
  Sphere_3_intersection_tester< CGAL::Epeck >(r, true /*exact predicates*/, true /*exact constructions*/).run();

  std::cout << " |||||||| Test CGAL::Homogeneous<CGAL::Epeck_ft> ||||||||" << std::endl;
  Sphere_3_intersection_tester< CGAL::Homogeneous<CGAL::Epeck_ft> >(r, true /*exact predicates*/, true /*exact constructions*/).run();

  std::cout << "OK!" << std::endl;
}
