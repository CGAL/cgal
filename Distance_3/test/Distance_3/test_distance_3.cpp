#include <CGAL/Simple_cartesian.h>
#include <CGAL/Simple_homogeneous.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_integer.h>
#include <CGAL/Homogeneous.h>

#include <CGAL/squared_distance_3.h>

#include <CGAL/Random.h>
#include <CGAL/Timer.h>

// #define CGAL_USE_GTE_AS_SANITY_CHECK
#ifdef CGAL_USE_GTE_AS_SANITY_CHECK
#include <Mathematics/DistTriangle3Triangle3.h>
#include <Mathematics/DistSegmentSegment.h>
#endif

#include <cassert>
#include <iostream>

struct randomint
{
  randomint() ;
  int get() const { return sequence[cur]; }
  int next() {
    cur = (cur + 1) % 11;
    return get();
  }

private:
  int sequence[11];
  int cur;
};

inline randomint::randomint()
{
  cur = 0;
  sequence[0] = 19;
  sequence[1] = 5;
  sequence[2] = 17;
  sequence[3] = 13;
  sequence[4] = 29;
  sequence[5] = 2;
  sequence[6] = 23;
  sequence[7] = 31;
  sequence[8] = 3;
  sequence[9] = 37;
  sequence[10] = 11;
}

randomint ri;

template <typename K>
struct Test
{
  typedef typename K::RT              RT;
  typedef typename K::FT              FT;
  typedef typename K::Point_3         P;
  typedef typename K::Segment_3       S;
  typedef typename K::Vector_3        V;
  typedef typename K::Ray_3           R;
  typedef typename K::Line_3          L;
  typedef typename K::Triangle_3      T;
  typedef typename K::Plane_3         Pl;
  typedef typename K::Tetrahedron_3   Tet;
  typedef typename K::Iso_cuboid_3    Cub;

private:
  CGAL::Random& r;
  const double epsilon = 1e-14;
  int N = 10;
  double m = 0, M = 1;

public:
   Test(CGAL::Random& r, const double epsilon) : r(r), epsilon(epsilon) { }

private:
  inline RT to_nt(int d) const { return RT(d); }

  P p(int x, int y, int z)
  {
    int w = ri.next();
    return P(to_nt(x*w), to_nt(y*w), to_nt(z*w), to_nt(w));
  }

  P random_point() const
  {
    return P(FT(r.get_double(m, M)), FT(r.get_double(m, M)), FT(r.get_double(m, M)));
  }

  Pl pl(int a, int b, int c, int d)
  {
    int w = ri.next();
    return Pl(to_nt(a*w), to_nt(b*w), to_nt(c*w), to_nt(d*w));
  }

private:
  template <typename Type>
  bool are_equal(const Type& t1, const Type& t2, const bool verbose = true)
  {
    const FT diff = CGAL::abs(t1 - t2);
    if(diff > std::numeric_limits<FT>::epsilon() &&
       diff > epsilon * (CGAL::abs(t1) + CGAL::abs(t2)))
    {
      if(verbose)
      {
        std::cerr << "Approximate comparison failed (t1|t2): got " << t1 << " but expected " << t2 << std::endl;
        std::cerr << "Diff: " << CGAL::abs(t1 - t2) << " vs eps: " <<  epsilon * (CGAL::abs(t1) + CGAL::abs(t2)) << std::endl;
      }
      return false;
    }

    return true;
  }

  template <typename O1, typename O2>
  void check_ss_distance(const O1& o1, const O2& o2)
  {
    FT res = CGAL::squared_distance(o1, o2);
    FT asd = compute_squared_distance_interval_between_segments(o1.source(), o1.target(),
                                                                o2.source(), o2.target(), K());

    assert(res == asd);

    std::cout << "input: " << o1.source() << " " << o1.target() << " " << o2.source() << " " << o2.target() << std::endl;
    std::cout << "result (old) = " << res << std::endl;
    std::cout << "result (new) = " << asd << std::endl;
  }

  void do_intersect_check(const P&, const P&) { }

  template <typename O2>
  void do_intersect_check(const P& p, const O2& o2)
  {
    if(!o2.is_degenerate() && CGAL::do_intersect(p, o2))
    {
      assert(are_equal(CGAL::squared_distance(p, o2), FT(0)));
      assert(are_equal(CGAL::squared_distance(o2, p), FT(0)));
    }
  }

  template <typename O1, typename O2>
  void do_intersect_check(const O1& o1, const O2& o2)
  {
    if(!o1.is_degenerate() && !o2.is_degenerate() && CGAL::do_intersect(o1, o2))
    {
      assert(are_equal(CGAL::squared_distance(o1, o2), FT(0)));
      assert(are_equal(CGAL::squared_distance(o2, o1), FT(0)));
    }
  }

  template <typename O1, typename O2>
  void check_squared_distance(const O1& o1, const O2& o2, const FT& expected_result)
  {
    const FT res_o1o2 = CGAL::squared_distance(o1, o2);
    const FT res_o2o1 = CGAL::squared_distance(o2, o1);

    assert(are_equal(res_o1o2, expected_result));
    assert(are_equal(res_o2o1, expected_result));

    do_intersect_check(o1, o2);
  }

  template <typename O1, typename O2>
  void check_squared_distance_with_bound(const O1& o1, const O2& o2, const FT& ubound)
  {
    const FT res_o1o2 = CGAL::squared_distance(o1, o2);
    const FT res_o2o1 = CGAL::squared_distance(o2, o1);

    do_intersect_check(o1, o2);

    assert(res_o1o2 <= ubound);
    assert(res_o2o1 <= ubound);
  }

private:
  void P_P()
  {
    std::cout << "Point - Point" << std::endl;
    check_squared_distance(p(0, 0, 0), p(0, 0, 0), 0);
    check_squared_distance(p(3, 5, 7), p(0, 0, 0), 83);
  }

  void P_S()
  {
    std::cout << "Point - Segment" << std::endl;
    check_squared_distance(p(0, 1, 2), S{p(-3, 0, 0), p( 2, 0, 0)}, 5);
    check_squared_distance(p(0, 1, 2), S{p( 3, 0, 0), p( 2, 0, 0)}, 9);
    check_squared_distance(p(0, 1, 2), S{p( 2, 0, 0), p( 3, 0, 0)}, 9);
    check_squared_distance(p(6, 1, 2), S{p( 2, 0, 0), p( 3, 0, 0)}, 14);
  }

  void P_T()
  {
    std::cout << "Point - Triangle" << std::endl;
    check_squared_distance(p(0, 1, 2), T{p(0, 0, 0), p(2, 0, 0), p(0, 2, 0)}, 4);

    T t(p(0,0,0), p(3,0,0), p(3,3,0));
    check_squared_distance (p(-1, -1, 0), t, 2);
    check_squared_distance (p(-1,  0, 0), t, 1);
    check_squared_distance (p(0, 0, 0), t, 0);
    check_squared_distance (p(1, 0, 0), t, 0);
    check_squared_distance (p(4, 0, 0), t, 1);
    check_squared_distance (p(1, -1, 0), t, 1);
    check_squared_distance (p(1, 1, 1), t, 1);
    check_squared_distance (p(1, 0, 1), t, 1);
    check_squared_distance (p(0, 0, 1), t, 1);

    // Degenerate
    check_squared_distance (p(1, 2, 3), T(p(4,3,2), p(4,3,2), p(4,3,2)), squared_distance(p(1, 2, 3), p(4,3,2)));
    check_squared_distance (p(1, 2, 3), T(p(4,3,2), p(10,12,3), p(4,3,2)), squared_distance(p(1, 2, 3), p(4,3,2)));
    check_squared_distance (p(0, 0, 0), T(p(4,3,2), p(4,-3,-2), p(4,3,2)), squared_distance(p(0, 0, 0), p(4,0,0)));

    // On the triangle
    check_squared_distance (p(7, 1, -5), T(p(2,9,8), p(-4,-3,-5), p(7, 1, -5)), 0); // vertex
    check_squared_distance (p(7, 1, -5), T(p(14,2,-10), p(-7,-1,5), p(8, 3, -1)), 0); // edge
    check_squared_distance (p(1, 4, -3), T(p(0,-8,-3), p(-5,14,-3), p(10, 1, -3)), 0); // face

    // General
    check_squared_distance (p(-15, 1, 0), T(p(-10, 1, 0), p(0,0,0), p(10,0,0)), 25);
    check_squared_distance (p(-5, 0, 0), T(p(-10, 1, 0), p(0,0,0), p(10,0,0)), squared_distance(p(-5, 0, 0), S(p(-10, 1, 0), p(0,0,0))));
    check_squared_distance (p(0, -3, 0), T(p(-10, 1, 0), p(0,0,0), p(10,0,0)), 9);
    check_squared_distance (p(3, -3, 0), T(p(-10, 1, 0), p(0,0,0), p(10,0,0)), squared_distance(p(3, -3, 0), S(p(0,0,0), p(10,0,0))));
    check_squared_distance (p(16, 1, 1), T(p(-10, 1, 0), p(0,0,0), p(10,0,0)), 38);
    check_squared_distance (p(5, 5, 2), T(p(-10, 1, 0), p(0,0,0), p(10,0,0)), squared_distance(p(5, 5, 2), S(p(10,0,0), p(-10,1,0))));

    for(int i=0; i<N; ++i)
    {
      P p0 = random_point();
      P p1 = random_point();
      P p2 = random_point();
      P q = random_point();

      check_squared_distance_with_bound(q, T(p0, p1, p2), squared_distance(q, S(p0, p1)));
      check_squared_distance_with_bound(q, T(p0, p1, p2), squared_distance(q, S(p1, p2)));
      check_squared_distance_with_bound(q, T(p0, p1, p2), squared_distance(q, S(p2, p0)));
    }
  }

  void P_Tet()
  {
    std::cout << "Point - Tetrahedron\n";
    check_squared_distance (p(0, 0, 0), Tet(p(0, 0, 0), p( 1, 0, 0), p( 0, 1, 0), p( 0, 0, 1)), 0);
    check_squared_distance (p(0, 0, 2), Tet(p(0, 0, 0), p( 1, 0, 0), p( 0, 1, 0), p( 0, 0, 1)), 1);
    check_squared_distance (p(0, 0, -1), Tet(p(0, 0, 0), p( 1, 0, 0), p( 0, 1, 0), p( 0, 0, 1)), 1);
    check_squared_distance (p(5, 0, 0), Tet(p(0, 0, 0), p( 1, 0, 0), p( 0, 1, 0), p( 4, 0, 1)), 2);
  }

  void S_S()
  {
    std::cout << "Segment - Segment" << std::endl;

    // coplanar segments (hardcoded)
    FT z(std::sqrt(2.));
    P p0{-1, 0, z};
    P p1{ 1, 0, z};

    // translations of (0, -1, z) -- (0, 1, z) to have all variations of x&y (<0, [0;1]; >1) in the code
    for(int j=-2;j<4; j+=2)
    {
      for(int k=-3;k<3; k+=2)
      {
        P p2{j,   k, z};
        P p3{j, k+2, z};

        // to explicit the expected distances
        if(j == -2 && k == -3)
          check_squared_distance(S{p0, p1}, S{p2, p3}, CGAL::squared_distance(p3, p0));
        else if(j == -2 && k == -1)
          check_squared_distance(S{p0, p1}, S{p2, p3}, 1);
        else if(j == -2 && k == 1)
          check_squared_distance(S{p0, p1}, S{p2, p3}, CGAL::squared_distance(p2, p0));
        else if(j == 0 && k == -3)
          check_squared_distance(S{p0, p1}, S{p2, p3}, 1);
        else if(j == 0 && k == -1)
          check_squared_distance(S{p0, p1}, S{p2, p3}, 0);
        else if(j == 0 && k == 1)
          check_squared_distance(S{p0, p1}, S{p2, p3}, 1);
        else if(j == 2 && k == -3)
          check_squared_distance(S{p0, p1}, S{p2, p3}, CGAL::squared_distance(p3, p1));
        else if(j == 2 && k == -1)
          check_squared_distance(S{p0, p1}, S{p2, p3}, 1);
        else if(j == 2 && k == 1)
          check_squared_distance(S{p0, p1}, S{p2, p3}, CGAL::squared_distance(p2, p1));
      }
    }

    for(int i=0; i<N; ++i)
    {
      P p0 = random_point();
      P p1 = random_point();
      P p2 = random_point();
      P p3 = random_point();
      p0 = CGAL::midpoint(p0, p1);
      p1 = p0 + FT(0.1) * V{p1 - p0};
      p2 = p2 + V{p2 - CGAL::ORIGIN} / CGAL::approximate_sqrt(CGAL::square(p2.x()) + CGAL::square(p2.y()) + CGAL::square(p2.z()) + 3);
      p3 = p3 + V{p3 - CGAL::ORIGIN} * FT(std::cos(1.3));

      // degenerate inputs
      check_squared_distance(S{p0, p0}, S{p0, p0}, 0); // both degen
      check_squared_distance(S{p3, p3}, S{p3, p3}, 0); // both degen
      check_squared_distance(S{p0, p0}, S{p0, p1}, 0); // left degen + common extremity (left)
      check_squared_distance(S{p0, p0}, S{p1, p0}, 0); // left degen + common extremity (right)
      check_squared_distance(S{p0, p0}, S{p2, p3}, CGAL::squared_distance(p0, S(p2, p3))); // left degen

      // common extremities
      check_squared_distance(S{p2, p3}, S{p2, p3}, 0); // equal segments
      check_squared_distance(S{p3, p2}, S{p2, p3}, 0); // equal segments but opposite dirs
      check_squared_distance(S{p2, p3}, S{p2, p1}, 0); // common generic (p2 common)
      check_squared_distance(S{p2, p3}, S{p1, p2}, 0); // common generic (p2 common)
      check_squared_distance(S{p2, p3}, S{p3, p1}, 0); // common generic (p3 common)
      check_squared_distance(S{p2, p3}, S{p1, p3}, 0); // common generic (p3 common)

      // collinear segments
      const double lambda_4 = r.get_double(0, 1);
      const P p4 = p2 + FT(lambda_4) * V{p3 - p2};
      const double lambda_5 = r.get_double(0, 1);
      const P p5 = p2 + FT(lambda_5) * V{p3 - p2};

      // [p2;p3) fully contains [p4;p5]
      check_squared_distance(S{p2, p3}, S{p4, p5}, 0);
      check_squared_distance(S{p2, p3}, S{p5, p4}, 0);
      check_squared_distance(S{p3, p2}, S{p4, p5}, 0);
      check_squared_distance(S{p3, p2}, S{p5, p4}, 0);

      const double lambda_6 = r.get_double(0, 1);
      const P p6 = p3 + FT(lambda_6) * V{p3 - p2};
      // [p2;p3] overlaps [p5;p6]
      check_squared_distance(S{p2, p3}, S{p6, p5}, 0);
      check_squared_distance(S{p2, p3}, S{p5, p6}, 0);
      check_squared_distance(S{p3, p2}, S{p6, p5}, 0);
      check_squared_distance(S{p3, p2}, S{p5, p6}, 0);

      const double lambda_7 = r.get_double(1, 2);
      const P p7 = p3 + FT(lambda_7) * V{p3 - p2};

      // [p2;p3] disjoint && left of [p6;p7]
      check_squared_distance(S{p2, p3}, S{p6, p7}, CGAL::squared_distance(p3, p6));
      check_squared_distance(S{p2, p3}, S{p7, p6}, CGAL::squared_distance(p3, p6));
      check_squared_distance(S{p3, p2}, S{p6, p7}, CGAL::squared_distance(p3, p6));
      check_squared_distance(S{p3, p2}, S{p7, p6}, CGAL::squared_distance(p3, p6));

      // Generic collinear
      const double lambda_8 = r.get_double(-M, M);
      const P p8 = p2 + FT(lambda_8) * V{p3 - p2};
      const double lambda_9 = r.get_double(-M, M);
      const P p9 = p2 + FT(lambda_9) * V{p3 - p2};

      S s89(p8, p9);
      S s32(p3, p2);
      FT result;
      if(!s89.is_degenerate() && !s32.is_degenerate()) // for do_intersect...
      {
        if(CGAL::do_intersect(s89, s32))
          result = 0;
        else
          result = (std::min)(CGAL::squared_distance(p2, p8),
                     (std::min)(CGAL::squared_distance(p2, p9),
                       (std::min)(CGAL::squared_distance(p3, p8),
                                  CGAL::squared_distance(p3, p9))));

#ifdef CGAL_USE_GTE_AS_SANITY_CHECK
        gte::DCPQuery<FT, gte::Segment<3, FT>, gte::Segment<3, FT> > GTE_SS_checker;
        gte::Segment<3, FT> gte_s1{{p8.x(), p8.y(), p8.z()}, {p9.x(), p9.y(), p9.z()}};
        gte::Segment<3, FT> gte_s2{{p3.x(), p3.y(), p3.z()}, {p2.x(), p2.y(), p2.z()}};
        auto gte_res = GTE_SS_checker(gte_s1, gte_s2);
        std::cout << "-------" << std::endl;
        std::cout << "old: " << CGAL::internal::squared_distance_old(s89, s32, K()) << std::endl;
        std::cout << "dist (GTE) : " << gte_res.sqrDistance << std::endl;
#endif

        // Because do_intersect() with constructions is likely to return 'false' even for overlaps
        assert(are_equal(CGAL::squared_distance(s89, s32), result, false /*verbose*/) ||
               are_equal(CGAL::squared_distance(s32, s89), FT(0)));
      }

      // completely generic
      S s1{p0, p1}, s2{p2, p3};
      do_intersect_check(s1, s2);

#ifdef CGAL_USE_GTE_AS_SANITY_CHECK
      gte::DCPQuery<FT, gte::Segment<3, FT>, gte::Segment<3, FT> > GTE_SS_checker;
      gte::Segment<3, FT> gte_s1{{p0.x(), p0.y(), p0.z()}, {p1.x(), p1.y(), p1.z()}};
      gte::Segment<3, FT> gte_s2{{p2.x(), p2.y(), p2.z()}, {p3.x(), p3.y(), p3.z()}};
      auto gte_res = GTE_SS_checker(gte_s1, gte_s2);

      std::cout << "dist (CGAL) : " << CGAL::squared_distance(s1, s2) << std::endl;
      std::cout << "dist (GTE) : " << gte_res.sqrDistance << std::endl;
      assert(are_equal(CGAL::squared_distance(s1, s2), gte_res.sqrDistance));
#endif
    }

    // a few brute force checks: sample a segments and use squared_distance(P3, S3)
    for(int i=0; i<10; ++i)
    {
      P p0 = random_point();
      P p1 = random_point();
      P p2 = random_point();
      P p3 = random_point();

      S s01{p0, p1};
      S s23{p2, p3};

      FT upper_bound = CGAL::squared_distance(p0, p2);

      V p01 = V{p0, p1} / FT(N);
      for(int l=0; l<N; ++l)
      {
        P tp = p0 + FT(l) * p01;
        FT sqd = CGAL::squared_distance(tp, s23);
        if(sqd < upper_bound)
          upper_bound = sqd;
      }

      // bit ugly, but if constructions are not exact, building `tp` introduces some error
      if(epsilon != 0)
        upper_bound *= (1 + 1e-10);

      check_squared_distance_with_bound(s01, s23, upper_bound);
    }
  }

  void P_R()
  {
    // Note : the value is not verified by hand
    std::cout << "Point - Ray" << std::endl;
    check_squared_distance_with_bound(p( -8, -7,  0), R{p(23, -27, 2), p( -17, 16, 2)}, 86.368512613);
  }

  void R_R()
  {
    // Note : the values are not verified by hand
    std::cout << "Ray - Ray" << std::endl;
    check_squared_distance_with_bound(R{p( 0, 0, 30), p(  0, 30, 30)}, R{p(100, -100, 0), p( 200,  1, 0)}, 20899.504975002);
    check_squared_distance(R{p( 1, 0,  0), p(  0,  0,  0)}, R{p(  1,    3, 3), p(   0,  0, 3)}, 9);
    check_squared_distance(R{p( 0, 0,  0), p(  1,  0,  0)}, R{p(  0,    0, 2), p(  -1,  0, 2)}, 4);
  }

  void S_R()
  {
    // Note : the values are not verified by hand
    std::cout << "Segment - Ray" << std::endl;
    check_squared_distance_with_bound(S{p( 0, 0, 30), p(  0, 30, 30)}, R{p(100, -100, 0), p( 200,  1, 0)}, 20899.504975002);
  }

  void R_L()
  {
    // Note : the values are not verified by hand
    std::cout << "Ray - Line" << std::endl;
    check_squared_distance_with_bound(R{p( 0, 0, 30), p(  0, 30, 30)}, L{p(100, -100, 0), p( 200,  1, 0)}, 20899.504975002);
    check_squared_distance(R{p(10, 0,  0), p( 20,  0,  0)}, L{p(  0,    0, 3), p(   0,  3, 3)}, 109);
    check_squared_distance(R{p( 1, 0,  0), p(  0,  0,  0)}, L{p(  1,    3, 3), p(   0,  0, 3)}, 9);
    check_squared_distance(R{p( 0, 0,  0), p(  1,  0,  0)}, L{p(  0,    0, 2), p(  -1,  0, 2)}, 4);
  }

  void P_L()
  {
    std::cout << "Point - Line" << std::endl;
    check_squared_distance(p(  0,  1,  2), L{p(  2,    0, 0), p(   3,  0, 0)}, 5);
    check_squared_distance(p(  0,  0,  2), L{p(  0,    0, 0), p(   1,  2, 0)}, 4);
  }

  void S_L()
  {
    // Note : the values are not verified by hand
    std::cout << "Segment - Line" << std::endl;
    check_squared_distance(S{p( 1, 0,  0), p(  0,  0,  0)}, L{p(  1,    3, 3), p(   0,  0, 3)}, 9);
    check_squared_distance(S{p(-90, 0,  0), p(-10,  0,  0)}, L{p(  0,    0, 3), p(   0,  3, 3)}, 109);
    check_squared_distance(S{p(  0, 0,  0), p(  1,  0,  0)}, L{p(  0,    0, 2), p(  -1,  0, 2)}, 4);
  }

  void L_L()
  {
    // Note : the values are not verified by hand
    std::cout << "Line - Line" << std::endl;
    check_squared_distance(L{p(-10, 0,  0), p(-90,  0,  0)}, L{p(  0,    0, 3), p(   0,  3, 3)}, 9);
    check_squared_distance(L{p(  1, 0,  0), p(  0,  0,  0)}, L{p(  1,    3, 3), p(   0,  0, 3)}, 9);
    check_squared_distance(L{p(  0, 0,  0), p(  1,  0,  0)}, L{p(  0,    0, 2), p(  -1,  0, 2)}, 4);
  }

  void P_Pl()
  {
    std::cout << "Point - Plane" << std::endl;
    check_squared_distance(p(2, 5,  3), Pl(0, 1, 0, 0), 25);
  }

  void S_Pl()
  {
    std::cout << "Segment - Plane" << std::endl;
    check_squared_distance(S{p(2, -3,  3), p( 3,-7, 4)}, pl(0, 1, 0, 0), 9);
  }

  void R_Pl()
  {
    std::cout << "Ray - Plane" << std::endl;
    check_squared_distance(R{p(2, -4,  3), p( 3,-4, 4)}, Pl(0, 1, 0, 0), 16);
    check_squared_distance(R{p(2, -4,  3), p( 3, 4, 4)}, Pl(0, 1, 0, 0), 0);
    check_squared_distance(R{p(2, -4,  3), p( 3,-8, 4)}, Pl(0, 1, 0, 0), 16);
  }

  void L_Pl()
  {
    std::cout << "Line - Plane" << std::endl;
    check_squared_distance(L{p(2, -4,  3), p( 3,-4, 4)}, Pl(0, 1, 0, 0), 16);
    check_squared_distance(L{p(2, -4,  3), p( 3, 4, 4)}, Pl(0, 1, 0, 0), 0);
    check_squared_distance(L{p(2, -4,  3), p( 3,-8, 4)}, Pl(0, 1, 0, 0), 0);
  }

  void Pl_Pl()
  {
    std::cout << "Plane - Plane" << std::endl;
    Pl p1(0, 1, 0, 0);
    typename K::Vector_3 v = -p1.orthogonal_vector();
    v /= CGAL::approximate_sqrt(v.squared_length());
    Pl p2 = Pl(0,-1,0,6);
    check_squared_distance(p1,p2, 36);
    check_squared_distance(Pl(-2, 1, 1, 0), Pl(2, 1, 3, 0), 0);
  }

  void T_T()
  {
    std::cout << "Triangle - Triangle" << std::endl;

    // min between vertices (hardcoded)
    check_squared_distance(T{p(0,0,0), p(1,0,0), p(0,1,0)}, T{p(0,0,2), p(-1,0,2), p(0,-1,2)}, 4);
    check_squared_distance(T{p(0,0,0), p(1,0,0), p(0,1,0)}, T{p(-1,0,2), p(0,0,2), p(0,-1,2)}, 4);

    check_squared_distance(T{p(1,2,3), P{FT(4.2),FT(5.3),-6}, p(7,-8,9)},
                           T{P{FT(10.1), 12, -10}, p(15, 14, -12), p(19, 45, -20)},
                           CGAL::squared_distance(P{FT(4.2),FT(5.3),-6}, P{FT(10.1), 12, -10}));

    // min vertex-edge (hardcoded)
    check_squared_distance(T{p(0,0,0), p(1,0,0), p(0,1,0)}, T{p(1,1,0), p(2,1,0), p(1,2,0)}, 0.5);
    check_squared_distance(T{p(0,0,0), p(2,0,0), p(0,2,0)}, T{p(0,-1,1), p(2,0,1), p(2,-1,1)}, 1);

    for(int i=0; i<N; ++i)
    {
      P p0 = random_point();
      P p1 = random_point();
      P p2 = random_point();
      P p3 = random_point();
      P p4 = random_point();
      P p5 = random_point();

      // these are still exact with EPECK
      p0 = CGAL::midpoint(p0, p1);
      p1 = p0 + FT(0.1) * V{p1 - p0};
      p2 = p2 + V{p2 - p0} / FT(CGAL_PI);

      // this is still exact with EPECK_with_sqrt
      p4 = p4 + V{p4 - CGAL::ORIGIN} / CGAL::approximate_sqrt(CGAL::square(p4.x()) + CGAL::square(p4.y()) + CGAL::square(p4.z()) + 3);

      p5 = p5 + V{p5 - CGAL::ORIGIN} * FT(std::cos(1.3));

      // degenerate inputs
      check_squared_distance(T{p3, p3, p3}, T{p3, p3, p3}, 0); // both degen
      check_squared_distance(T{p0, p0, p0}, T{p3, p3, p3}, CGAL::squared_distance(p0, p3)); // both degen

      check_squared_distance(T{p0, p0, p0}, T{p0, p0, p3}, 0); // single degen and common edge
      check_squared_distance(T{p0, p0, p0}, T{p3, p0, p0}, 0);
      check_squared_distance(T{p0, p0, p0}, T{p0, p3, p0}, 0);

      check_squared_distance(T{p0, p0, p0}, T{p0, p3, p4}, 0); // single degen and common vertex
      check_squared_distance(T{p0, p0, p0}, T{p3, p0, p4}, 0);
      check_squared_distance(T{p0, p0, p0}, T{p3, p4, p0}, 0);

      // degen into point & degen into segment
      check_squared_distance(T{p1, p1, p1}, T{p4, p3, p3}, CGAL::squared_distance(p1, S{p3, p4}));
      check_squared_distance(T{p5, p5, p5}, T{p3, p3, p4}, CGAL::squared_distance(p5, S{p3, p4}));

      // both degen into segment
      check_squared_distance(T{p0, p1, p0}, T{p3, p3, p4}, CGAL::squared_distance(S{p0, p1}, S{p3, p4}));
      check_squared_distance(T{p5, p5, p4}, T{p4, p3, p3}, CGAL::squared_distance(S{p5, p4}, S{p3, p4}));

      // common vertex
      check_squared_distance(T{p0, p1, p2}, T{p0, p3, p4}, 0);
      check_squared_distance(T{p0, p1, p2}, T{p3, p0, p4}, 0);
      check_squared_distance(T{p0, p1, p2}, T{p4, p3, p0}, 0);
      check_squared_distance(T{p0, p1, p2}, T{p1, p3, p4}, 0);
      check_squared_distance(T{p0, p1, p2}, T{p3, p1, p4}, 0);
      check_squared_distance(T{p0, p1, p2}, T{p4, p3, p1}, 0);
      check_squared_distance(T{p0, p1, p2}, T{p2, p3, p4}, 0);
      check_squared_distance(T{p0, p1, p2}, T{p3, p2, p4}, 0);
      check_squared_distance(T{p0, p1, p2}, T{p4, p3, p2}, 0);

      // common edge
      check_squared_distance(T{p0, p1, p2}, T{p0, p1, p4}, 0);
      check_squared_distance(T{p0, p1, p2}, T{p1, p0, p4}, 0);
      check_squared_distance(T{p0, p1, p2}, T{p4, p0, p1}, 0);
      check_squared_distance(T{p0, p1, p2}, T{p4, p1, p0}, 0);
      check_squared_distance(T{p0, p1, p2}, T{p0, p4, p1}, 0);
      check_squared_distance(T{p0, p1, p2}, T{p1, p4, p2}, 0);

      check_squared_distance(T{p0, p1, p2}, T{p2, p1, p4}, 0);
      check_squared_distance(T{p0, p1, p2}, T{p1, p2, p4}, 0);
      check_squared_distance(T{p0, p1, p2}, T{p4, p2, p1}, 0);
      check_squared_distance(T{p0, p1, p2}, T{p4, p1, p2}, 0);
      check_squared_distance(T{p0, p1, p2}, T{p2, p4, p1}, 0);
      check_squared_distance(T{p0, p1, p2}, T{p1, p4, p2}, 0);

      check_squared_distance(T{p0, p1, p2}, T{p0, p2, p4}, 0);
      check_squared_distance(T{p0, p1, p2}, T{p2, p0, p4}, 0);
      check_squared_distance(T{p0, p1, p2}, T{p4, p0, p2}, 0);
      check_squared_distance(T{p0, p1, p2}, T{p4, p2, p0}, 0);
      check_squared_distance(T{p0, p1, p2}, T{p0, p4, p2}, 0);
      check_squared_distance(T{p0, p1, p2}, T{p2, p4, p0}, 0);

      // same triangle
      check_squared_distance(T{p0, p1, p2}, T{p0, p1, p2}, 0);
      check_squared_distance(T{p0, p1, p2}, T{p1, p2, p0}, 0);
      check_squared_distance(T{p0, p1, p2}, T{p2, p0, p1}, 0);
      check_squared_distance(T{p0, p1, p2}, T{p2, p1, p0}, 0);
      check_squared_distance(T{p0, p1, p2}, T{p0, p2, p1}, 0);
      check_squared_distance(T{p0, p1, p2}, T{p1, p0, p2}, 0);

      // vertex on triangle
      double lambda = r.get_double(0, 1);
      double mu = r.get_double(0, 1 - lambda);
      const P bp = CGAL::barycenter(p0, lambda, p1, mu, p2, 1 - lambda - mu);
      check_squared_distance(T{p0, p1, p2}, T{bp, p3, p4}, 0);
      check_squared_distance(T{p0, p1, p2}, T{p3, bp, p4}, 0);
      check_squared_distance(T{p0, p1, p2}, T{p3, p4, bp}, 0);

      // edge on triangle
      lambda = r.get_double(0, 1);
      mu = r.get_double(0, 1 - lambda);
      P bp2 = CGAL::barycenter(p0, lambda, p1, mu, p2, 1 - lambda - mu);
      check_squared_distance(T{p0, p1, p2}, T{bp, bp2, p4}, 0);
      check_squared_distance(T{p0, p1, p2}, T{bp2, bp, p4}, 0);
      check_squared_distance(T{p0, p1, p2}, T{bp, p4, bp2}, 0);
      check_squared_distance(T{p0, p1, p2}, T{bp2, p4, bp}, 0);
      check_squared_distance(T{p0, p1, p2}, T{p4, bp, bp2}, 0);
      check_squared_distance(T{p0, p1, p2}, T{p4, bp2, bp}, 0);

      // part of the edge crossing the triangle
      lambda = r.get_double(-1, 1);
      mu = r.get_double(-1, 1);
      bp2 = CGAL::barycenter(p0, lambda, p1, mu, p2, 1 - lambda - mu);
      check_squared_distance(T{p0, p1, p2}, T{bp, bp2, p4}, 0);
      check_squared_distance(T{p0, p1, p2}, T{bp2, bp, p4}, 0);
      check_squared_distance(T{p0, p1, p2}, T{bp, p4, bp2}, 0);
      check_squared_distance(T{p0, p1, p2}, T{bp2, p4, bp}, 0);
      check_squared_distance(T{p0, p1, p2}, T{p4, bp, bp2}, 0);
      check_squared_distance(T{p0, p1, p2}, T{p4, bp2, bp}, 0);

      // generic triangles
      T tr1{p0, p1, p2}, tr2{p3, p4, p5};
      do_intersect_check(tr1, tr2);

#ifdef CGAL_USE_GTE_AS_SANITY_CHECK
      gte::DCPQuery<FT, gte::Triangle3<FT>, gte::Triangle3<FT> > GTE_TT_checker;
      gte::Triangle3<FT> gte_tr1{{p0.x(), p0.y(), p0.z()}, {p1.x(), p1.y(), p1.z()}, {p2.x(), p2.y(), p2.z()}};
      gte::Triangle3<FT> gte_tr2{{p3.x(), p3.y(), p3.z()}, {p4.x(), p4.y(), p4.z()}, {p5.x(), p5.y(), p5.z()}};
      auto gte_res = GTE_TT_checker(gte_tr1, gte_tr2);

      std::cout << "dist (CGAL) : " << CGAL::squared_distance(tr1, tr2) << std::endl;
      std::cout << "dist (GTE) : " << gte_res.sqrDistance << std::endl;
      std::cout << "diff CGAL GTE : " << (gte_res.sqrDistance - CGAL::squared_distance(tr1, tr2)) << std::endl;

      // don't assert on purpose, GTE has slightly (10^-30 different results, even with an exact NT)
      are_equal(CGAL::squared_distance(tr1, tr2), gte_res.sqrDistance);
#endif
    }
  }

public:
  void run()
  {
    std::cout << "Kernel: " << typeid(K).name() << std::endl;

    P_P();
    P_S();
    P_R();
    P_L();
    P_T();
    P_Pl();
    P_Tet();

    S_S();
    S_R();
    S_L();
    S_Pl();

    R_R();
    R_L();
    R_Pl();

    L_L();
    L_Pl();

    T_T();

    Pl_Pl();
  }
};

int main()
{
  std::cout.precision(17);
  std::cerr.precision(17);

  std::cout << "3D Distance tests" << std::endl;

  CGAL::Random r;
  std::cout << "random seed = " << r.get_seed() << std::endl;

  // @todo Some tests are too difficult for these kernels
//  Test<CGAL::Simple_cartesian<double> >(r, 1e-5).run();
//  Test<CGAL::Simple_homogeneous<double> >(r, 1e-5).run();
//  Test<CGAL::Simple_cartesian<CGAL::Interval_nt<true> > >(r, 1e-5).run();

  Test<CGAL::Homogeneous<CGAL::Exact_integer> >(r, 0).run();

  const double epick_eps = 10 * std::numeric_limits<double>::epsilon();
  Test<CGAL::Exact_predicates_inexact_constructions_kernel>(r, epick_eps).run();

  Test<CGAL::Exact_predicates_exact_constructions_kernel>(r, 0).run();

  std::cout << "Done!" << std::endl;

  return EXIT_SUCCESS;
}
