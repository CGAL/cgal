// 3D intersection tests.

// We want to check that no division is performed for interface macro Do_intersect_3_RT
#define CGAL_NO_MPZF_DIVISION_OPERATOR

// Tetrahedron_3_X with X < Tetrahedron (lexicographically) are tested in other files

#include <CGAL/Intersections_3/Tetrahedron_3_Tetrahedron_3.h>
#include <CGAL/Intersections_3/Tetrahedron_3_Triangle_3.h>

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
struct Tetrahedron_3_intersection_tester
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

  typedef std::vector<P>            Poly;

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
  Tetrahedron_3_intersection_tester(CGAL::Random& r,
                                    const bool has_exact_p = false,
                                    const bool has_exact_c = false)
    : Base(r, has_exact_p, has_exact_c)
  { }

public:
  void Tet_Tet()
  {
    std::cout << "Tetrahedron - Tetrahedron\n";

    for(int i=0; i<N; ++i)
    {
      P tet00 = random_point(), tet01 = random_point(), tet02 = random_point(), tet03 = random_point();
      Tet tet0(tet00, tet01, tet02, tet03);
      if(tet0.is_degenerate())
        continue;

      P tet10 = random_point(), tet11 = random_point(), tet12 = random_point(), tet13 = random_point();

      // ---
      Tet tet1(tet00, tet11, tet12, tet13);
      if(tet1.is_degenerate())
        continue;

      check_do_intersect(tet0, tet1);

      // ---
      Tet tet2(tet00, tet11, tet02, tet13);
      if(tet2.is_degenerate())
        continue;

      check_do_intersect(tet0, tet2);

      // ---
      Tet tet3(tet00, tet03, tet02, tet13);
      if(tet3.is_degenerate())
        continue;

      check_do_intersect(tet0, tet3);

      // ---
      Tet tet4(tet10, tet13, tet12, tet13);
      if(tet4.is_degenerate())
        continue;

      std::set<CGAL::Bounded_side> os;
      os.insert(tet0.bounded_side(tet10));
      os.insert(tet0.bounded_side(tet11));
      os.insert(tet0.bounded_side(tet12));
      os.insert(tet0.bounded_side(tet13));

      if(os.size() != 1 || (*os.begin() != CGAL::ON_UNBOUNDED_SIDE))
        check_do_intersect(tet0, tet4);
    }
  }

  void Tet_Tr()
  {
    Tet tet(p(0,0,0), p(0,1,0), p(1,0,0), p(0,0,1));

    // no intersection
    check_no_intersection(tet, Tr(p(2,0,0), p(2,1,0), p(2,0,1)));

    // tr share a vertex, triangle on face
    check_intersection(tet, Tr(p(0,1,0), P(0.1,0.1,0), P(0.5,0.1,0)),
                       Tr(p(0,1,0), P(0.1,0.1,0), P(0.5,0.1,0)));

    // tr vertex adjacent to a vertex
    check_intersection(tet, Tr(p(0,1,0), p(3,1,0), p(0,3,0)),
                       p(0,1,0));

    // tr edge adjacent to a vertex
    check_intersection(tet, Tr(p(-1,1,0), p(3,1,0), p(0,3,0)),
                       p(0,1,0));

    // tr adjacent to a vertex, outside
    check_intersection(tet, Tr(p(-1,1,12), p(0,1,-50), P(0.5,1,-0.5)),
                       p(0,1,0));

    // only one edge intersecting
    check_intersection(tet, Tr(P(-2, 0.5, 0.5), P(-2,0.75,1), P(1.5, 0.5, 0.5)),
                       P(0,0.5,0.5));

    // edge shared, 3rd point outside
    check_intersection(tet, Tr(p(0,1,0), p(1,0,0), P(0.5,0,-100)),
                       S(p(0,1,0), p(1,0,0)));
    check_intersection(tet, Tr(P(0.75,0.25,0), p(10,10,10), P(0.25,0.75,0)),
                       S(P(0.75,0.25,0), P(0.25,0.75,0)));

    // shared edge, 3rd point inside
    check_intersection(tet, Tr(p(0,1,0), p(1,0,0), P(0.25,0.25,0.25)),
                       Tr(p(0,1,0), p(1,0,0), P(0.25,0.25,0.25)));

    // tr edge containing a tetrahedron edge
    check_intersection(tet, Tr(p(0,-1,0), p(0,3,0), p(1,1,-1)), S(p(0,0,0), p(0,1,0)));

    // tr edge intersecting a tet face
    check_intersection(tet, Tr(P(-0.5,0,0.5), P(2,0,0.5), P(0.5,-2,0)),
                       S(P(0,0,0.5), P(0.5,0,0.5)));

    // tr face containing a tetrahedron edge
    check_intersection(Tet(p(0,1,0), p(0,0,1), p(0,1,1), p(1, 0, 1)),
                       Tr(p(-2,0,-2), P(0,0,5), p(10,0,0)),
                       S(p(0,0,1), p(1,0,1)));

    // tr inside a face
    check_intersection(tet, Tr(P(0,0.9,0), P(0,0.1,0), P(0,0,0.9)),
                       Tr(P(0,0.9,0), P(0,0.1,0), P(0,0,0.9)));

    // face inside tr
    check_intersection(tet, Tr(p(0,2,0), p(0,-2,0), p(0,0,2)),
                       Tr(p(0,1,0), p(0,0,0), p(0,0,1)));

    // on the same plane, polygonal intersection
    Base::template check_intersection<Poly>(tet, Tr(p(0,-2,0), p(1,1,0), p(1,2,0)));

    // traversing triangles
    Base::template check_intersection<Tr>(tet, Tr(P(-2, 0.5, 0.25), P(-2, 0.75, 0.6), P(1.5, 0.5, 0.25)));

    // tr share an edge, third point on face
    Base::template check_intersection<Tr>(tet, Tr(p(0,1,0), P(0.1,0.1,0), P(0.9,0.1,0)));

    // sharing a vertex + edge/face intersection
    Base::template check_intersection<Poly>(tet, Tr(p(0,0,0), P(0.5,-1,0), P(0.5,1,0)));

    // vertex on edge & triangle inside, triangle intersection
    Base::template check_intersection<Tr>(tet, Tr(P(0.5,0,0), P(0.5,0.25,0.25), P(0.5,-0.5,0)));

    // vertex on edge & triangle inside, double segment
    Base::template check_intersection<Tr>(tet, Tr(P(0,0,0.5), P(-1,1.5,0.25), P(1.5,-1,0.25)));

    // vertex on edge & triangle inside, double segment non-incident
    Base::template check_intersection<Poly>(tet, Tr(P(0.25,0,0.25), P(-1,0.5,0.25), P(1.5,0.5,0.25)));

    // vertex on face, triangle outside & point intersection
    Base::check_intersection(tet, Tr(P(-1,1,0.25), P(-1,0,0.25), P(0,0.25,0.25)),
                             P(0, 0.25, 0.25));
    Base::check_intersection(tet, Tr(P(-1,0,0.25), P(-1,1,0.25), P(0,0.25,0.25)),
                             P(0, 0.25, 0.25));
    Base::check_intersection(tet, Tr(P(0,0.25,0.25), P(-1,1,0.25), P(-1,0,0.25)),
                             P(0, 0.25, 0.25));

    // vertex on face, triangle outside & segment intersection
    Base::check_intersection(tet, Tr(P(0.5,0,-0.25), P(0.5,0,0.25), P(0.5,-0.5,0)),
                             S(P(0.5, 0, 0.25), P(0.5, 0, 0)));

    // vertex on face & inside, triangle intersection
    Base::template check_intersection<Tr>(tet, Tr(P(0.5,0,0), P(0.5,0.25,0.25), P(0.5,-0.5,0)));

    // vertex on face, triangle intersection
    Base::template check_intersection<Tr>(tet, Tr(P(0, 0.25, 0.25), P(0.5, -0.25, 0.25), P(0.5, -0.25, 0.5)));

    // vertex inside, segment-triangle intersection
    Base::template check_intersection<Tr>(tet, Tr(P(0.2, 0.15, 0.3), P(-2, 0.6, 0.15), P(-2, 0.12, 0.15)));

    // vertex inside, double segment intersection
    Base::template check_intersection<Poly>(tet, Tr(P(0.2, 0.15, 0.3), P(-1, 0.15, -2), P(-1, 0.15, 2)));

    // two vertices inside
    Base::template check_intersection<Poly>(tet, Tr(P(0.1, 0.5, 0.1), P(0.3,0.1,0.1), P(4,0.3,0.9)));

    if(this->has_exact_c)
    {
      // sharing a point on edge
      check_intersection(tet, Tr(P(-2, 1, 0.5), P(-2, 0.75, 1.6), P(2, 0, 0.5)),
                         P(0, 0.5, 0.5));

      Tr tr(p(-2,2,0), p(2,2,0), P(0.25,0.25,0));
      auto res = CGAL::intersection(tet, tr);
      const Poly* poly = std::get_if<Poly>(&*res);
      assert(poly != nullptr);
      assert(poly->size() == 4);
      for(const P& pt : *poly) {
        assert(tet.has_on_boundary(pt) && tr.has_on(pt));
      }

      tr = Tr(P(0.45, 0.20, 0.1), P(0.1, 0.20, 0.5), P(-0.5, 0.25, -0.5));
      res = CGAL::intersection(tet, tr);
      const Poly* inter = std::get_if<Poly>(&*res);
      assert(inter != nullptr);
      assert(inter->size() == 5);
      for(const P& pt : *inter) {
        assert((tet.has_on_bounded_side(pt) || tet.has_on_boundary(pt)) && tr.has_on(pt));
      }
    }

    for(int i=0; i<N; ++i)
    {
      P tet0 = random_point(), tet1 = random_point(), tet2 = random_point(), tet3 = random_point();
      Tet tet(tet0, tet1, tet2, tet3);
      if(tet.is_degenerate())
        continue;

      P tr0 = random_point(), tr1 = random_point(), tr2 = random_point();

      // ---
      Tr tr(tet0, tr1, tr2);
      if(tr.is_degenerate())
        continue;

      check_do_intersect(tet, tr);

      // ---
      tr = Tr(tet0, tr1, tet2);
      if(tr.is_degenerate())
        continue;

      check_do_intersect(tet, tr);

      // ---
      tr = Tr(tr0, tr1, tr2);
      if(tr.is_degenerate())
        continue;

      std::set<CGAL::Bounded_side> os;
      os.insert(tet.bounded_side(tr0));
      os.insert(tet.bounded_side(tr1));
      os.insert(tet.bounded_side(tr2));

      // conservative
      if(os.size() != 1 || (*os.begin() != CGAL::ON_UNBOUNDED_SIDE))
      {
        check_do_intersect(tet, tr);

        if(!this->has_exact_c)
          continue;

        auto res = CGAL::intersection(tet, tr);

        std::vector<P> points;
        if(const P* pt = std::get_if<P>(&*res))
        {
          points.push_back(*pt);
        }
        else if(const S* s = std::get_if<S>(&*res))
        {
          points.push_back(s->source());
          points.push_back(s->target());
        }
        else if(const Tr* itr = std::get_if<Tr>(&*res))
        {
          points.push_back(itr->operator[](0));
          points.push_back(itr->operator[](1));
          points.push_back(itr->operator[](2));
        }
        else if(const Poly* poly = std::get_if<Poly>(&*res))
        {
          points = *poly;

//          assert(CGAL::is_simple_2(poly->begin(), poly->end(),
//                                   CGAL::Projection_traits_3(tr.supporting_plane().orthogonal_vector())));
        }

        for(const P& pt : points) {
          assert(!tet.has_on_unbounded_side(pt) && tr.has_on(pt));
        }
      }
    }
  }

  void issue_6777()
  {
    Tr tri(P(0.191630, -0.331630, -0.370000), P(-0.124185, -0.385815, -0.185000), P(-0.0700000, -0.0700000, 0.00000));
    Tet tet(P(0, -1, 0), P(-1, 0, 0), P(0, 0, 0), P(0, 0, -1));
    auto res = intersection(tri, tet);
    assert(res != std::nullopt);
    const std::vector<P> *vps = std::get_if<std::vector<P>>(&*res);
    assert(vps!=nullptr);
  }

  void run()
  {
    std::cout << "3D Tetrahedron Intersection tests\n";

    Tet_Tet();
    Tet_Tr();
    issue_6777();
  }
};

int main(int, char**)
{
  std::cout.precision(17);
  std::cerr.precision(17);

  CGAL::Random r;
  std::cout << "random seed = " << r.get_seed() << std::endl;

  std::cout << " |||||||| Test Simple_cartesian<double> ||||||||" << std::endl;
  Tetrahedron_3_intersection_tester< CGAL::Simple_cartesian<double> >(r).run();

  std::cout << " |||||||| Test CGAL::Homogeneous<CGAL::MP_Float> ||||||||" << std::endl;
  Tetrahedron_3_intersection_tester< CGAL::Homogeneous<CGAL::MP_Float> >(r).run();

  std::cout << " |||||||| Test EPICK ||||||||" << std::endl;
  Tetrahedron_3_intersection_tester< CGAL::Epick >(r, true /*exact predicates*/).run();

  std::cout << " |||||||| Test EPECK ||||||||" << std::endl;
  Tetrahedron_3_intersection_tester< CGAL::Epeck >(r, true /*exact predicates*/, true /*exact constructions*/).run();

  std::cout << " |||||||| Test CGAL::Homogeneous<CGAL::Epeck_ft> ||||||||" << std::endl;
  Tetrahedron_3_intersection_tester< CGAL::Homogeneous<CGAL::Epeck_ft> >(r, true /*exact predicates*/, true /*exact constructions*/).run();

  std::cout << "OK!" << std::endl;
}
