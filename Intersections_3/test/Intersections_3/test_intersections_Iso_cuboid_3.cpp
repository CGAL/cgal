// 3D intersection tests.

// We want to check that no division is performed for interface macro Do_intersect_3_RT
#define CGAL_NO_MPZF_DIVISION_OPERATOR

// Iso_cuboid_3_X with X < Iso_cuboid (lexicographically) are tested in other files

#include <CGAL/Intersections_3/Iso_cuboid_3_Iso_cuboid_3.h>
#include <CGAL/Intersections_3/Iso_cuboid_3_Line_3.h>
#include <CGAL/Intersections_3/Iso_cuboid_3_Plane_3.h>
#include <CGAL/Intersections_3/Iso_cuboid_3_Point_3.h>
#include <CGAL/Intersections_3/Iso_cuboid_3_Ray_3.h>
#include <CGAL/Intersections_3/Iso_cuboid_3_Segment_3.h>
#include <CGAL/Intersections_3/Iso_cuboid_3_Sphere_3.h>
#include <CGAL/Intersections_3/Iso_cuboid_3_Tetrahedron_3.h>
#include <CGAL/Intersections_3/Iso_cuboid_3_Triangle_3.h>

#include <CGAL/squared_distance_3.h>

#include <CGAL/Vector_3.h>

#include "intersection_test_helper.h"

#include <CGAL/MP_Float.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Homogeneous.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/number_utils.h>

#include <cassert>
#include <iostream>

template <typename K>
struct Iso_cuboid_3_intersection_tester
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

  typedef std::vector<P>            Pol;

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
  Iso_cuboid_3_intersection_tester(CGAL::Random& r,
                                   const bool has_exact_p = false,
                                   const bool has_exact_c = false)
    : Base(r, has_exact_p, has_exact_c)
  { }

public:
  void Cub_Cub()
  {
    std::cout << "Iso_cuboid - Iso_cuboid\n";

    // No intersection
    check_no_intersection(Cub(p(0,1,2), p(0,1,2)), Cub(p(4,8,6), p(5,9,7))); // degenerate cuboids
    check_no_intersection(Cub(p(4,8,6), p(5,9,7)), Cub(p(0,4,1), p(0,4,1)));
    check_no_intersection(Cub(p(4,8,6), p(4,8,6)), Cub(p(0,4,1), p(0,4,1)));

    check_no_intersection(Cub(p(-7, 6, 1), p(7092, 71, 58)), Cub(p(-758, 98725, 43), p(17, 9025473, 47)));
    check_no_intersection(Cub(p(-73, 6, 1), p(-70, 71, 58)), Cub(p(8, -98725, 43), p(17, 9025473, 47)));
    check_no_intersection(Cub(p(-7, 6, 1), p(7092, 71, 58)), Cub(p(-758, -98725, -47), p(17, 9025473, -43)));

    // Sharing a point
    check_intersection   (Cub(p(0,1,2), p(7,6,4)), Cub(p(7,1,4), p(7,1,4)), // degenerate cuboids
                          Cub(p(7,1,4), p(7,1,4)));
    check_intersection   (Cub(p(0,1,2), p(7,8,9)), Cub(p(5,1,2), p(5,1,2)), // degenerate cuboids
                          Cub(p(5,1,2), p(5,1,2)));
    check_intersection   (Cub(p(0,1,2), p(0,1,2)), Cub(p(0,1,2), p(0,1,2)), // degenerate cuboids
                          Cub(p(0,1,2), p(0,1,2)));

    check_intersection   (Cub(p(1,3,2), p(4,8,6)), Cub(p(4,8,6), p(5,9,7)),
                          Cub(p(4,8,6), p(4,8,6)));
    check_intersection   (Cub(p(1,3,2), p(4,8,6)), Cub(p(-3,-6,-4), p(1,3,2)),
                          Cub(p(1,3,2), p(1,3,2)));
    check_intersection   (Cub(p(1,3,2), p(4,8,6)), Cub(p(4,0,0), p(6,3,2)),
                          Cub(p(4,3,2), p(4,3,2)));
    check_intersection   (Cub(p(1,3,2), p(4,8,5)), Cub(p(4,8,1), p(9,10,2)),
                          Cub(p(4,8,2), p(4,8,2)));

    // Sharing an edge
    check_intersection   (Cub(p(1,3,2), p(4,8,5)), Cub(p(4,3,0), p(7,8,2)),
                          Cub(p(4,3,2), p(4,8,2)));
    check_intersection   (Cub(p(1,3,2), p(4,8,5)), Cub(p(4,4,0), p(7,6,2)),
                          Cub(p(4,4,2), p(4,6,2)));

    // Sharing a face
    check_intersection   (Cub(p(1,3,2), p(4,8,5)), Cub(p(2,4,0), p(3,6,2)), // face within a face
                          Cub(p(2,4,2), p(3,6,2)));
    check_intersection   (Cub(p(1,3,2), p(4,8,5)), Cub(p(2,4,0), p(4,8,2)), // face within a face
                          Cub(p(2,4,2), p(4,8,2)));
    check_intersection   (Cub(p(1,3,2), p(4,8,5)), Cub(p(1,8,2), p(4,8,5)), // full face
                          Cub(p(1,8,2), p(4,8,5)));
    check_intersection   (Cub(p(1,3,2), p(4,8,5)), Cub(p(1,8,2), p(10,12,8)), // full face
                          Cub(p(1,8,2), p(4,8,5)));

    // Proper intersection
    check_intersection   (Cub(p(-7, 6, 1), p(7092, 71, 58)), Cub(p(758, 50, 43), p(758, 50, 43)), // degenerate cuboids
                          Cub(p(758, 50, 43), p(758, 50, 43)));
    check_intersection   (Cub(p(-7, 6, 1), p(7092, 71, 58)), Cub(p(-758, -98725, 43), p(17, 9025473, 47)),
                          Cub(p(-7, 6, 43), p(17, 71, 47)));
  }

  void Cub_L()
  {
    std::cout << "Iso_cuboid - Line\n";

    // No intersection
    check_no_intersection(L(p(0, 0, 0), p(1, 0, 3)), Cub(p(2, 1, 1), p(3, 5, 8)));
    check_no_intersection(L(p(4, 0, 0), p(4, 1, 3)), Cub(p(2, 1, 1), p(3, 5, 8)));

    // through vertex
    check_intersection   (Cub(p(-7, -8, -9), p(-1, 2, -4)), L(p(-8,1,-9), p(-6, 3, -9)),
                          p(-7,2,-9));

    // Segment
    check_intersection   (Cub(p(1, 1, 1), p(3, 5, 8)), L(p(0, 0, 3), p(1, 2, 3)),
                          S(P(1, 2, 3), P(2.5, 5, 3)));
    check_intersection   (Cub(p(1, 1, 1), p(3, 5, 8)), L(p(1, 0, 0), p(1, 2, 3)),
                          S(P(1, 1, 1.5), P(1, 5, 7.5)));
    check_intersection   (Cub(p(1, 1, 1), p(3, 5, 8)), L(p(0, 0, 0), p(1, 2, 3)),
                          S(P(1, 2, 3), P(2.5, 5, 7.5)));

    for(int i=0; i<N; ++i)
    {
      int mx = this->r.get_int(this->m, this->M);
      int my = this->r.get_int(this->m, this->M);
      int mz = this->r.get_int(this->m, this->M);
      int Mx = this->r.get_int(this->m, this->M);
      int My = this->r.get_int(this->m, this->M);
      int Mz = this->r.get_int(this->m, this->M);

      if(abs(mx - Mx) < 2) Mx = mx + 2;
      if(abs(my - My) < 2) My = my + 2;
      if(abs(mz - Mz) < 2) Mz = mz + 2;
      if(mx > Mx) std::swap(mx, Mx);
      if(my > My) std::swap(my, My);
      if(mz > Mz) std::swap(mz, Mz);

      P mp = p(mx, my, mz), Mp = p(Mx, My, Mz);
      Cub cub(mp, Mp);

      P q1 = p(this->r.get_int(mx - this->M, mx - 1), this->r.get_int(my, My), this->r.get_int(mz, Mz)),
        q2 = p(this->r.get_int(Mx + 1, Mx + this->M), this->r.get_int(my, My), this->r.get_int(mz, Mz));
      assert(!cub.has_on_bounded_side(q1) && !cub.has_on_bounded_side(q2));

      L l(q1, q2);
      Base::template check_intersection<S>(cub, l);
      Base::template check_intersection<S>(cub, l, S(q1, q2), cub);

      q1 = p(this->r.get_int(mx, Mx), this->r.get_int(my - this->M, my), this->r.get_int(mz, Mz)),
      q2 = p(this->r.get_int(mx+1, Mx), this->r.get_int(my+1, My), this->r.get_int(mz+1, Mz));
      assert(cub.has_on_unbounded_side(q1) && cub.has_on_bounded_side(q2));

      l = L(q1, q2);
      Base::template check_intersection<S>(cub, l);
    }
  }

  void Cub_Pl()
  {
    std::cout << "Iso_cuboid - Plane\n";

    Cub cub(p(1,1,1), p(2,2,2));

    // no intersection
    check_no_intersection(cub, Pl(p(3,0,0), p(3,1,0), p(3,0,1)));
    check_no_intersection(Cub(p(1,1,1), p(1,1,1)), Pl(p(3,0,0), p(3,1,0), p(3,0,1)));

    // vertex
    check_intersection(cub, Pl(0.5, -0.5, -0.5, 0),
                       p(2,1,1));
    check_intersection(Cub(p(1, 2, 3), p(5,9,5)), Pl(p(4,-2,3), p(1,2,3), p(3,5,-8)),
                       p(1,2,3));

    // edge
    check_intersection(cub, Pl(p(1,1,1), p(1,2,1), P(1.5,0,0)),
                       S(p(1,1,1), p(1,2,1)));

    if(this->has_exact_c)
    {
      // face
      auto res = CGAL::intersection(cub, Pl(p(1,1,1), p(1,2,1), p(1,2,2)));

      const std::vector<P>* poly = std::get_if<std::vector<P> >(&*res);
      assert(poly != nullptr);
      assert(poly->size() == 4);
      for(const P& pt : *poly)
      {
        assert(pt.x() == 1);
      }
      res = CGAL::intersection(cub, Pl(p(1,1,1), p(1,2,1), p(2,2,2)));

      poly = std::get_if<std::vector<P> >(&*res);
      assert(poly != nullptr);
      assert(poly->size() == 4);
      for(const P& pt : *poly) {
        assert(cub.has_on_boundary(pt));
      }

      // other edge
      Pl pln(p(1,1,1), p(1,2,1), P(1.5, 1, 2));
      res = CGAL::intersection(cub, pln);
      poly = std::get_if<std::vector<P> >(&*res);
      assert(poly != nullptr);
      assert(poly->size() == 4);
      for(const P& pt : *poly) {
        assert(pln.has_on(pt) && cub.has_on_boundary(pt));
      }

      // triangle
      check_intersection(cub, Pl(P(2, 1.66, 2),
                                 P(1.66,2,2),
                                 P(2,2,1.66)),
                         Tr(P(1.66,2,2),
                            P(2, 2, 1.66),
                            P(2,1.66,2)));

      // random
      pln = Pl(0.265189, 0.902464, 0.33946, -2.47551);
      res = CGAL::intersection(cub, pln);
      poly = std::get_if<std::vector<P> >(&*res);
      assert(poly != nullptr);
      assert(poly->size() == 5);
      for(const P& pt : *poly) {
        assert(pln.has_on(pt) && cub.has_on_boundary(pt));
      }
    }
    else
    {
      // face
      auto res = CGAL::intersection(cub, Pl(p(1,1,1), p(1,2,1), p(1,2,2)));
      CGAL_USE(res);

      res = CGAL::intersection(cub, Pl(p(1,1,1), p(1,2,1), p(2,2,2)));
      CGAL_USE(res);

      CGAL::intersection(cub, Pl(0.5, -0.5, -0.5, 0));
      CGAL::intersection(cub, Pl(P(2, 1.66, 2),
                                 P(1.66,2,2),
                                 P(2,2,1.66)));

      Pl pln(0.265189, 0.902464, 0.33946, -2.47551);
      res = CGAL::intersection(cub, pln);
      CGAL_USE(res);
    }
  }

  void Cub_P()
  {
    std::cout << "Iso_cuboid - Point\n";

    // Outside
    check_no_intersection(Cub(p(0,1,2), p(0,1,2)), p(4,8,6)); // degenerate
    check_no_intersection(Cub(p(0,1,2), p(3,7,5)), p(4,8,6));

    // On vertex
    check_intersection(Cub(p(0,4,2), p(0,4,2)), p(0,4,2), p(0,4,2)); // degenerate
    check_intersection(Cub(p(0,1,2), p(8,7,2)), p(8,7,2), p(8,7,2));
    check_intersection(Cub(p(0,1,2), p(8,4,8)), p(8,1,8), p(8,1,8));
    check_intersection(Cub(p(0,1,2), p(0,4,2)), p(0,4,2), p(0,4,2));

    // On edge
    check_intersection(Cub(p(6,3,-2), p(6,3,2)), p(6,3,2), p(6,3,2)); // degenerate
    check_intersection(Cub(p(6,3,-2), p(6,3,2)), p(6,3,0), p(6,3,0)); // degenerate
    check_intersection(Cub(p(0,1,2), p(0,4,5)), p(0,3,2), p(0,3,2));
    check_intersection(Cub(p(0,1,2), p(0,4,5)), p(0,1,3), p(0,1,3));

    // On face
    check_intersection(Cub(p(0,1,2), p(0,4,5)), p(0,3,3), p(0,3,3)); // degenerate
    check_intersection(Cub(p(0,1,2), p(3,4,5)), p(2,4,4), p(2,4,4));

    // Within
    check_intersection(Cub(p(0,1,2), p(10,11,12)), p(4,8,6), p(4,8,6));
  }

  void Cub_R()
  {
    std::cout << "Iso_cuboid - Ray\n";

    check_no_intersection  (Cub(p(2,1,1), p(3,5,8)), R(p(0,0,0), p(1,0,3)));
    check_no_intersection  (Cub(p(2,1,1), p(3,5,8)), R(p(4,0,0), p(4,1,3)));

    // point
    check_intersection     (Cub(p(-7, -8, -9), p(-1, 2, -4)), R(p(-7, 2, -4), p(-9, -3, 0)),
                            p(-7,2,-4));

    // segment
    check_intersection     (Cub(p(1, 4, 8), p(2, 9, 10)), R(p(-1, -6, 4), p(0, -1, 6)),
                            S(p(1, 4, 8), P(2,9,10)));

    check_intersection     (Cub(p(-7, -8, -9), p(-1, 2, -4)), R(p(-3, 1, -5), p(-2, -5, -7)),
                            S(p(-3, 1, -5), P(-1.5, -8, -8)));
    check_intersection     (Cub(p(1, 1, 1), p(3, 5, 8)), R(p(0, 0, 3), p(1, 2, 3)),
                            S(p(1, 2, 3), P(2.5, 5, 3)));
    check_intersection     (Cub(p(1, 1, 1), p(3, 5, 8)), R(p(1, 0, 0), p(1, 2, 3)),
                            S(P(1, 1,1.5), P(1, 5,7.5)));
    check_intersection     (Cub(p(1, 1, 1), p(3, 5, 8)), R(p(0, 0, 0), p(1, 2, 3)),
                            S(p(1, 2, 3), P(2.5, 5, 7.5)));

    Base::template check_intersection<S>(Cub(p(2, 1, 1), p(3, 5, 8)), R(p(0, 2, 0), p(1, 2, 3)));
  }

  void Cub_S()
  {
    std::cout << "Segment - Iso_cuboid\n";

    // no intersection
    check_no_intersection  (Cub(p(2, 1, 1), p(3, 5, 8)), S(p(0, 2, 0), p(1, 2, 3)));
    check_no_intersection  (Cub(p(2, 1, 1), p(3, 5, 8)), S(p(0, 0, 0), p(1, 0, 3)));
    check_no_intersection  (Cub(p(2, 1, 1), p(3, 5, 8)), S(p(4, 0, 0), p(4, 1, 3)));

    // point
    check_intersection     (Cub(p(1, 1, 1), p(3, 5, 8)), S(p(0, 0, 3), p(1, 2, 3)),
                            p(1, 2, 3));
    check_intersection     (Cub(p(1, 1, 1), p(3, 5, 8)), S(p(0, 0, 0), p(1, 2, 3)),
                            p(1, 2, 3));

    // segment
    check_intersection     (Cub(p(-7, -8, -9), p(-1, 2, -4)), S(p(-3, 1, -5), p(-2, -5, -7)),
                            S(p(-3, 1, -5), p(-2, -5, -7)));
    check_intersection     (Cub(p(1, 1, 1), p(3, 5, 8)), S(p(1, 0, 0), p(1, 2, 3)),
                            S(P(1, 1, 1.5), p(1, 2, 3)));

    for(int i=0; i<N; ++i)
    {
      int mx = this->r.get_int(this->m, this->M);
      int my = this->r.get_int(this->m, this->M);
      int mz = this->r.get_int(this->m, this->M);
      int Mx = this->r.get_int(this->m, this->M);
      int My = this->r.get_int(this->m, this->M);
      int Mz = this->r.get_int(this->m, this->M);

      if(mx == Mx) ++Mx;
      if(my == My) ++My;
      if(mz == Mz) ++Mz;
      if(mx > Mx) std::swap(mx, Mx);
      if(my > My) std::swap(my, My);
      if(mz > Mz) std::swap(mz, Mz);

      P mp = p(mx, my, mz), Mp = p(Mx, My, Mz);
      Cub cub(mp, Mp);

      P s0 = p(this->r.get_int(mx, Mx), this->r.get_int(my - this->M, my), this->r.get_int(mz, Mz)),
        s1 = p(this->r.get_int(mx, Mx), this->r.get_int(My, My + this->M), this->r.get_int(mz, Mz));
      assert(!cub.has_on_bounded_side(s0) && !cub.has_on_bounded_side(s1));

      if(CGAL::do_intersect(cub, S(s0, s1)))
      {
        check_intersection(cub, S(s0, s1), L(s0, s1), cub);
        check_intersection(cub, S(s0, s1), R(s0, s1), cub);
      }

      s1 = p(this->r.get_int(mx, Mx), this->r.get_int(my, My), this->r.get_int(mz, Mz));
      assert(!cub.has_on_unbounded_side(s1));

      check_do_intersect(cub, S(s0, s1));
      check_intersection(cub, S(s1, s0), R(s1, s0), cub);
    }
  }

  void Cub_Sph()
  {
    std::cout << "Iso_cuboid - Sphere\n";

    // no intersection
    check_do_not_intersect(Cub(p(0,0,0), p(0,0,0)), Sph(p(1,0,0), 0));

    // degenerate
    check_do_intersect(Cub(p(0,0,0), p(0,0,0)), Sph(p(0,0,0), 0));
    check_do_intersect(Cub(p(0,0,0), p(0,0,0)), Sph(p(0,0,0), 1));
    check_do_intersect(Cub(p(0,0,0), p(1,1,1)), Sph(p(0,0,0), 0));

    // tangent
    check_do_intersect(Cub(p(0,0,0), p(0,0,0)), Sph(p(1,0,0), 1)); // at vertex
    check_do_intersect(Cub(p(0,0,0), p(2,2,2)), Sph(p(1,-1,0), 1)); // at edge
    check_do_intersect(Cub(p(0,0,0), p(2,2,2)), Sph(p(1,1,-1), 1)); // at face

    for(int i=0; i<N; ++i)
    {
      int mx = this->r.get_int(this->m, this->M);
      int my = this->r.get_int(this->m, this->M);
      int mz = this->r.get_int(this->m, this->M);
      int Mx = this->r.get_int(this->m, this->M);
      int My = this->r.get_int(this->m, this->M);
      int Mz = this->r.get_int(this->m, this->M);

      if(mx == Mx) ++Mx;
      if(my == My) ++My;
      if(mz == Mz) ++Mz;
      if(mx > Mx) std::swap(mx, Mx);
      if(my > My) std::swap(my, My);
      if(mz > Mz) std::swap(mz, Mz);

      P mp = p(mx, my, mz), Mp = p(Mx, My, Mz);
      Cub cub(mp, Mp);

      P c = p(this->r.get_int(this->m, this->M), this->r.get_int(this->m, this->M), this->r.get_int(this->m, this->M));

      // smallest distance to the iso cuboid, manually...
      FT sq_r = CGAL::squared_distance(c, Tr(p(mx, my, mz), p(Mx, my, mz), p(mx, my, Mz)));
      sq_r = (std::min)(sq_r, CGAL::squared_distance(c, Tr(p(Mx, my, mz), p(Mx, my, Mz), p(mx, my, Mz))));
      sq_r = (std::min)(sq_r, CGAL::squared_distance(c, Tr(p(Mx, my, mz), p(Mx, my, Mz), p(Mx, My, Mz))));
      sq_r = (std::min)(sq_r, CGAL::squared_distance(c, Tr(p(Mx, my, mz), p(Mx, My, mz), p(Mx, My, Mz))));
      sq_r = (std::min)(sq_r, CGAL::squared_distance(c, Tr(p(mx, my, mz), p(Mx, my, mz), p(Mx, My, mz))));
      sq_r = (std::min)(sq_r, CGAL::squared_distance(c, Tr(p(mx, my, mz), p(Mx, My, mz), p(mx, My, mz))));
      sq_r = (std::min)(sq_r, CGAL::squared_distance(c, Tr(p(mx, my, mz), p(mx, My, mz), p(mx, My, Mz))));
      sq_r = (std::min)(sq_r, CGAL::squared_distance(c, Tr(p(mx, my, mz), p(mx, my, Mz), p(mx, My, Mz))));
      sq_r = (std::min)(sq_r, CGAL::squared_distance(c, Tr(p(mx, my, Mz), p(Mx, my, Mz), p(mx, My, Mz))));
      sq_r = (std::min)(sq_r, CGAL::squared_distance(c, Tr(p(mx, My, Mz), p(Mx, My, Mz), p(Mx, my, Mz))));
      sq_r = (std::min)(sq_r, CGAL::squared_distance(c, Tr(p(mx, My, mz), p(mx, My, Mz), p(Mx, My, Mz))));
      sq_r = (std::min)(sq_r, CGAL::squared_distance(c, Tr(p(mx, My, mz), p(Mx, My, mz), p(Mx, My, Mz))));

      check_do_intersect(cub, Sph(c, sq_r));
    }

    // generic inside
    for(int i=0; i<N; ++i)
    {
      int mx = this->r.get_int(this->m, this->M);
      int my = this->r.get_int(this->m, this->M);
      int mz = this->r.get_int(this->m, this->M);
      int Mx = this->r.get_int(this->m, this->M);
      int My = this->r.get_int(this->m, this->M);
      int Mz = this->r.get_int(this->m, this->M);

      if(mx == Mx) ++Mx;
      if(my == My) ++My;
      if(mz == Mz) ++Mz;
      if(mx > Mx) std::swap(mx, Mx);
      if(my > My) std::swap(my, My);
      if(mz > Mz) std::swap(mz, Mz);

      P mp = p(mx, my, mz), Mp = p(Mx, My, Mz);
      Cub cub(mp, Mp);

      P c = p(this->r.get_int(mx, Mx), this->r.get_int(my, My), this->r.get_int(mz, Mz));
      int ms = (std::min)((std::min)(Mx - mx, My - my), Mz - mz) + 1;

      check_do_intersect(cub, Sph(c, CGAL::square(ms)));
    }
  }

  void Cub_Tet()
  {
    std::cout << "Iso_cuboid - Sphere\n";

    // no intersection
    check_do_not_intersect(Cub(p(0,0,0), p(1,1,1)), Tet(p(2,4,1), p(1, 2, 0), p(2,3,1), p(4,5,6)));

    // point
    check_do_intersect(Cub(p(0,0,0), p(1,1,1)), Tet(p(2,4,1), p(1,1,0), p(2,3,1), p(4,5,6))
                       /*p(1,1,1)*/); // shared vertex
    check_do_intersect(Cub(p(0,0,0), p(2,1,1)), Tet(p(2,4,1), p(1,1,0), p(2,3,1), p(4,5,6))
                       /*p(1,1,1)*/); // on edge
    check_do_intersect(Cub(p(0,0,0), p(2,2,2)), Tet(p(2,4,1), p(1,2,1), p(2,3,1), p(4,5,6))
                       /*p(1,1,2)*/); // on face

    // edge
    check_do_intersect(Cub(p(0,0,0), p(1,1,1)), Tet(p(0,1,0), p(1,1,0), p(2,3,1), p(4,5,6))
                       /*S(p(0,1,0), p(1,1,0))*/);
    check_do_intersect(Cub(p(0,0,0), p(2,2,2)), Tet(p(0,2,0), p(1,2,0), p(2,3,1), p(4,5,6))
                       /*S(p(0,1,0), p(1,1,0))*/);
    check_do_intersect(Cub(p(0,0,0), p(2,2,2)), Tet(p(0,2,0), p(2,2,2), p(2,3,1), p(4,5,6))
                       /*S(p(0,2,0), p(2,2,0))*/);
    check_do_intersect(Cub(p(0,0,0), p(2,2,2)), Tet(p(-1,2,-1), p(3,2,3), p(2,3,1), p(4,5,6))
                       /*S(p(0,2,0), p(2,2,0))*/);

    // face
    check_do_intersect(Cub(p(0,0,0), p(1,1,1)), Tet(p(0,1,0), p(1,1,0), p(1,1,1), p(4,5,6))
                       /*Tr(p(0,1,0), p(1,1,0), p(1,1,1))*/);
    check_do_intersect(Cub(p(0,0,0), p(1,1,1)), Tet(p(0,1,0), p(1,1,0), p(4,5,6), p(1,1,1))
                       /*Tr(p(0,1,0), p(1,1,0), p(1,1,1))*/);

    // generic inside
    for(int i=0; i<N; ++i)
    {
      int mx = this->r.get_int(this->m, this->M);
      int my = this->r.get_int(this->m, this->M);
      int mz = this->r.get_int(this->m, this->M);
      int Mx = this->r.get_int(this->m, this->M);
      int My = this->r.get_int(this->m, this->M);
      int Mz = this->r.get_int(this->m, this->M);

      if(mx == Mx) ++Mx;
      if(my == My) ++My;
      if(mz == Mz) ++Mz;
      if(mx > Mx) std::swap(mx, Mx);
      if(my > My) std::swap(my, My);
      if(mz > Mz) std::swap(mz, Mz);

      P mp = p(mx, my, mz), Mp = p(Mx, My, Mz);
      Cub cub(mp, Mp);

      P c = p(this->r.get_int(mx, Mx), this->r.get_int(my, My), this->r.get_int(mz, Mz));
      P q1 = p(this->r.get_int(this->m, this->M), this->r.get_int(this->m, this->M), this->r.get_int(this->m, this->M)),
        q2 = p(this->r.get_int(this->m, this->M), this->r.get_int(this->m, this->M), this->r.get_int(this->m, this->M)),
        q3 = p(this->r.get_int(this->m, this->M), this->r.get_int(this->m, this->M), this->r.get_int(this->m, this->M));

      Tet tet(c, q1, q2, q3);
      if(!tet.is_degenerate())
        check_do_intersect(cub, tet);

      int id = this->r.get_int(0, 8);
      Tet teti(q1, cub[id], q2, q3);
      if(!teti.is_degenerate())
        check_do_intersect(cub, teti);

      int jd = this->r.get_int(0, 8);
      Tet tetj(q1, q2, cub[id], cub[jd]);
      if(!tetj.is_degenerate())
        check_do_intersect(cub, tetj);
    }
  }

  void Cub_Tr()
  {
    typedef typename CGAL::Intersection_traits<K, Tr, Cub>::result_type Res;

    std::cout << "Iso_cuboid - Triangle\n";

    // no intersection
    Cub cub(p(1,1,1), p(2,2,2));
    check_no_intersection(cub, Tr(P(1.1,2,0), p(2,3,1), p(4,5,6)));

    // tr adj to a cuboid vertex
    check_intersection(cub, Tr(P(1, 0.5, 0.5), p(3, 2, 1), p(3, 1, 2)), p(2,1,1));

    // tr adj to a point on a cuboid edge
    check_intersection(cub, Tr(P(1, 0.5, 0.5), p(3, 2, 1), p(3, 1, 2)), p(2,1,1));

    // tr adj to a point on a cuboid face
    check_intersection(cub, Tr(P(1, 1.5, 1.5), p(0, 0, 0), p(-4, 3, 1)), P(1, 1.5, 1.5));

    // tr adj to an edge
    check_intersection(cub, Tr(P(2, 1.5, 2), p(5, 6, 7), p(4, 7, 6)), P(2, 1.5, 2));

    // tr sharing an edge
    check_intersection(cub, Tr(P(2, 1.5, 2), P(2, 2.5, 2), p(4, 7, 6)),
                       S(P(2, 1.5, 2), p(2, 2, 2)));

    // tr sharing part of an edge
    check_intersection(cub, Tr(P(2, 1.5, 2), p(5, 6, 7), p(4, 7, 6)), P(2, 1.5, 2));

    // tr in a face
    check_intersection(cub, Tr(P(1, 1.1, 1), P(1, 1.5, 1), P(1, 1, 1.1)),
                       Tr(P(1, 1.1, 1), P(1, 1.5, 1), P(1, 1, 1.1)));

    // face in a tr
    Tr tr(p(-3, -3, 1), p(3, -3, 1), P(1.5, 6, 1));
    Res res = CGAL::intersection(cub, tr);
    Pol* poly = std::get_if<std::vector<P> >(&*res);
    assert(poly != nullptr);
    assert(poly->size() == 4);
    if(this->has_exact_c)
    {
      for(const P& pt : *poly)
        assert(tr.has_on(pt) && cub.has_on_boundary(pt));
    }

    // tr inside
    check_intersection(cub, Tr(P(1.1,1.1,1.1), P(1.8,1.8,1.8), P(1.5,1.8,1.1)),
                       Tr(P(1.1,1.1,1.1), P(1.8,1.8,1.8), P(1.5,1.8,1.1)));

    for(int i=0; i<N; ++i)
    {
      int mx = this->r.get_int(this->m, this->M);
      int my = this->r.get_int(this->m, this->M);
      int mz = this->r.get_int(this->m, this->M);
      int Mx = this->r.get_int(this->m, this->M);
      int My = this->r.get_int(this->m, this->M);
      int Mz = this->r.get_int(this->m, this->M);

      if(mx == Mx) ++Mx;
      if(my == My) ++My;
      if(mz == Mz) ++Mz;
      if(mx > Mx) std::swap(mx, Mx);
      if(my > My) std::swap(my, My);
      if(mz > Mz) std::swap(mz, Mz);

      P mp = p(mx, my, mz), Mp = p(Mx, My, Mz);
      Cub rcub(mp, Mp);

      P q1 = p(this->r.get_int(mx, Mx), this->r.get_int(my, My), this->r.get_int(mz, Mz)),
        q2 = p(this->r.get_int(mx, Mx), this->r.get_int(my, My), this->r.get_int(mz, Mz)),
        q3 = p(this->r.get_int(mx, Mx), this->r.get_int(my, My), this->r.get_int(mz, Mz));
      Tr rtr(q1, q2, q3);

      check_intersection(rcub, rtr, rtr);
    }

    // tr through
    tr = Tr(p(2, 4, 2), P(1, 3.5, -0.5), p(1, -1, 1));
    res = CGAL::intersection(cub, tr);
    poly = std::get_if<std::vector<P> >(&*res);
    assert(poly != nullptr);
    assert(poly->size() == 4);
    if(this->has_exact_c)
    {
      for(const P& pt : *poly) {
        assert(tr.has_on(pt) && cub.has_on_boundary(pt));
      }
    }

    // cutting in half along diagonal (intersection == triangle)
    check_intersection(cub, Tr(p(1, 1, 1), p(2, 2, 2), p(2, 2, 1)),
                       Tr(p(1, 1, 1), p(2, 2, 2), p(2, 2, 1)));

    // cutting in half along diagonal (intersection included in triangle)
    tr = Tr(p(1, 1, 10), p(10, 10, 1), p(1, 1, 1));
    res = CGAL::intersection(cub, tr);
    poly = std::get_if<std::vector<P> >(&*res);
    assert(poly != nullptr);
    assert(poly->size() == 4);
    if(this->has_exact_c)
    {
      for(const P& pt : *poly)
        assert(tr.has_on(pt) && cub.has_on_boundary(pt));
    }

    // 6 points intersection
    tr = Tr(P(18.66, -5.4, -11.33), P(-2.41, -7.33, 19.75), P(-10.29, 20.15, -10.33));
    res = CGAL::intersection(cub, tr);
    poly = std::get_if<std::vector<P> >(&*res);
    assert(poly != nullptr);
    assert(poly->size() == 6);
    if(this->has_exact_c)
    {
      for(const P& pt : *poly)
        assert(tr.has_on(pt) && cub.has_on_boundary(pt));
    }

    // triangle clipping a cuboid corner
    tr = Tr(P(1.02, 1.33, 0.62), P(1.95, 2.54, 0.95), P(0.79, 2.36, 1.92));
    res = CGAL::intersection(cub, tr);
    Tr* tr_res = std::get_if<Tr>(&*res);
    assert(tr_res != nullptr);
    if(this->has_exact_c)
    {
      assert(cub.has_on_boundary((*tr_res)[0]));
      assert(cub.has_on_boundary((*tr_res)[1]));
      assert(cub.has_on_boundary((*tr_res)[2]));
    }
  }

  void run()
  {
    std::cout << "3D Line Intersection tests\n";

    Cub_Cub();
    Cub_L();
    Cub_Pl();
    Cub_P();
    Cub_R();
    Cub_S();
    Cub_Sph();
    Cub_Tet();
    Cub_Tr();
  }
};

int main(int, char**)
{
  std::cout.precision(17);
  std::cerr.precision(17);

  CGAL::Random r;
  std::cout << "random seed = " << r.get_seed() << std::endl;

  std::cout << " |||||||| Test Simple_cartesian<double> ||||||||" << std::endl;
  Iso_cuboid_3_intersection_tester< CGAL::Simple_cartesian<double> >(r).run();

  std::cout << " |||||||| Test CGAL::Homogeneous<CGAL::MP_Float> ||||||||" << std::endl;
  Iso_cuboid_3_intersection_tester< CGAL::Homogeneous<CGAL::MP_Float> >(r).run();

  std::cout << " |||||||| Test EPICK ||||||||" << std::endl;
  Iso_cuboid_3_intersection_tester< CGAL::Epick >(r, true /*exact predicates*/).run();

  std::cout << " |||||||| Test EPECK ||||||||" << std::endl;
  Iso_cuboid_3_intersection_tester< CGAL::Epeck >(r, true /*exact predicates*/, true /*exact constructions*/).run();

  std::cout << " |||||||| Test CGAL::Homogeneous<CGAL::Epeck_ft> ||||||||" << std::endl;
  Iso_cuboid_3_intersection_tester< CGAL::Homogeneous<CGAL::Epeck_ft> >(r, true /*exact predicates*/, true /*exact constructions*/).run();

  std::cout << "OK!" << std::endl;
}
