// 2D intersection tests.

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Homogeneous.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/FPU.h>
#include <CGAL/Number_types/internal/Exact_type_selector.h>
#include <CGAL/intersection_2.h>

#include <vector>
#include <iostream>
#include <cassert>

const double epsilon = 0.001;

struct randomint
{
  randomint() ;
  int        get() const { return sequence[cur]; }
  int next() { cur = (cur+1)%11; return get();}
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

inline double to_nt(int d)
{
    return double(d);
}

template <typename K >
struct Test
{
  typedef typename K::Point_2               P;
  typedef typename K::Line_2                L;
  typedef typename K::Segment_2             S;
  typedef typename K::Ray_2                 R;
  typedef typename K::Triangle_2            T;
  typedef typename K::Iso_rectangle_2       Rec;
  typedef typename K::Circle_2              C;

  typedef CGAL::Bbox_2                      B;
  typedef std::vector<P>                    Pol;

  Test(const K& k = K()) : kernel(k) { }

  template <typename Type >
  bool approx_equal_nt(const Type &t1, const Type &t2)
  {
    if (t1 == t2)
      return true;

    if (CGAL::abs(t1 - t2) / (CGAL::max)(CGAL::abs(t1), CGAL::abs(t2)) < epsilon)
      return true;

    std::cout << " Approximate comparison failed between : " << t1 << "  and  " << t2 << "" << std::endl;
    return false;
  }

  template <typename Type >
  bool approx_equal(const Type&t1, const Type&t2)
  {
    return t1 == t2;
    // we need approx equal to check approx kernels, but maybe we should only test with exact kernels
    // (approx kernels were useful before, when the text output was checked by diff ?)
    // idea : test containment with intervals ?  or use some "epsilon double"?
    // I need to convert the text output to exact rationals in the source...
    // Well, for now the current scheme works.
  }

  bool approx_equal(const P & p, const P & q)
  {
    return approx_equal_nt(p.x(), q.x()) && approx_equal_nt(p.y(), q.y());
  }

  bool approx_equal(const Pol & p, const Pol & q)
  {
    if (p.size() != q.size())
      return false;
    if (CGAL::Intersections::internal::Is_cw<K,
        typename CGAL::Algebraic_structure_traits<typename K::FT>::Is_exact>()(p))
      return false;
    for(typename Pol::const_iterator itp = p.begin(), itq = q.begin(); itp != p.end(); ++itp, ++itq)
      if (!approx_equal(*itp, *itq))
        return false;

    return true;
  }

  bool approx_equal(const T& p, const T& q)
  {
    if(p.size() != q.size())
      return false;

    std::vector<P> vec { p[0], p[1], p[2] };

    if (CGAL::Intersections::internal::Is_cw<K,
        typename CGAL::Algebraic_structure_traits<typename K::FT>::Is_exact>(vec))
      return false;

    return approx_equal(p[0], q[0]) && approx_equal(p[1], q[1]) && approx_equal(p[2], q[2]);
  }

  template <typename O1, typename O2>
  void check_no_do_intersect(const O1& o1, const O2& o2)
  {
    assert(!CGAL::do_intersect(o1, o2));
    assert(!CGAL::do_intersect(o2, o1));

    typename K::Do_intersect_2 do_2 = kernel.do_intersect_2_object();
    assert(!do_2(o1, o2));
    assert(!do_2(o2, o1));
  }

  template <typename O1, typename O2>
  void check_no_intersection(const O1& o1, const O2& o2)
  {
    check_no_do_intersect(o1, o2);

    assert(!CGAL::intersection(o2, o1));

    typename K::Intersect_2 i_2 = kernel.intersect_2_object();
    assert(!i_2(o1, o2));
    assert(!i_2(o2, o1));
  }

  template <typename O1, typename O2 >
  void check_do_intersect(const O1& o1, const O2& o2)
  {
    assert(CGAL::do_intersect(o1, o2));
    assert(CGAL::do_intersect(o2, o1));
  }

  template <typename Res, typename O1, typename O2 >
  void check_intersection(const O1& o1, const O2& o2)
  {
    check_do_intersect(o1, o2);

    Res tmp;
    assert(CGAL::assign(tmp, CGAL::intersection(o1, o2)));
    assert(CGAL::assign(tmp, CGAL::intersection(o2, o1)));
  }

  template <typename Res, typename O1, typename O2 >
  void check_intersection(const O1& o1, const O2& o2, const Res& result, bool do_opposite=true)
  {
    check_do_intersect(o1, o2);

    Res tmp;

    assert(CGAL::assign(tmp, CGAL::intersection(o1, o2)));
    assert(approx_equal(tmp, result));
    if (do_opposite)
    {
      assert(CGAL::assign(tmp, CGAL::intersection(o2, o1)));
      assert(approx_equal(tmp, result));
    }
  }

  template <typename O >
  void check_intersection(const O& o)
  {
    return check_intersection(o, o, o);
  }

  P p(int x, int y)
  {
    int w = ri.next();
    return P(to_nt(x*w), to_nt(y*w), to_nt(w));
  }

  void B_C()
  {
    std::cout << "Bbox - Circle" << std::endl;

    // no intersection
    check_no_do_intersect  (p(0, 0).bbox() + p( 2,  3).bbox(), C(p( 8, 9), 3));
    check_no_do_intersect  (p(8, 9).bbox() + p( 9, 10).bbox(), C(p( 9, 9), 9)); // circle containing the bbox

    // point intersection
    check_do_intersect     (p(0, 0).bbox() + p( 2,  3).bbox(), C(p(-1, 0), 1));
    check_do_intersect     (p(0, 0).bbox() + p( 2,  3).bbox(), C(p( 1, 1), 1));

    // generic intersection
    check_do_intersect     (p(0, 0).bbox() + p( 2,  3).bbox(), C(p(-1, 0), 2));
    check_do_intersect     (p(0, 0).bbox() + p( 2,  0).bbox(), C(p( 1, 0), 1));
    check_do_intersect     (p(0, 0).bbox() + p(10, 10).bbox(), C(p( 3, 2), 3));
  }

  void B_L()
  {
    std::cout << "Bbox - Line" << std::endl;

    // no intersection
    check_no_do_intersect  (p(0,0).bbox()                , L(p( 1,1), p( 0,1)));
    check_no_do_intersect  (p(0,0).bbox() + p(4,5).bbox(), L(p(-1,1), p(-1,5)));

    // point intersection
    check_do_intersect     (p( 0, 0).bbox()                , L(p(-1,1), p( 1,-1)));
    check_do_intersect     (p(-1,-1).bbox() + p(4,2).bbox(), L(p( 1,5), p(-3,-1)));

    // segment intersection
    check_do_intersect     (p(0,0).bbox() + p(4,5).bbox(), L(p(-1, 5), p(7,5)));
    check_do_intersect     (p(0,0).bbox() + p(4,5).bbox(), L(p(-1, 5), p(2,5)));
    check_do_intersect     (p(0,0).bbox() + p(4,5).bbox(), L(p( 4,-3), p(4,1)));
  }

  void B_P()
  {
    std::cout << "Bbox - Point" << std::endl;

    // no intersection
    check_no_intersection  (p(1,2).bbox() + p(4,6).bbox(), p(8,9)); // point is outside the bbox

    // point intersection
    check_intersection     (p(1,2).bbox()                , p(1,2), p(1,2)); // degenerate bbox (0d)
    check_intersection     (p(1,2).bbox() + p(4,2).bbox(), p(2,2), p(2,2)); // degenerate bbox (1d)
    check_intersection     (p(1,2).bbox() + p(4,6).bbox(), p(1,6), p(1,6)); // point is a bbox corner
    check_intersection     (p(1,2).bbox() + p(4,6).bbox(), p(3,6), p(3,6)); // point is on a bbox edge
    check_intersection     (p(1,2).bbox() + p(4,6).bbox(), p(3,3), p(3,3)); // point is within the bbox
  }

  void B_R()
  {
    std::cout << "Bbox - Ray" << std::endl;

    // no intersection
    check_no_do_intersect  (p(0,0).bbox()                , R(p( 1,1), p( 0,1)));
    check_no_do_intersect  (p(0,0).bbox() + p(4,5).bbox(), R(p(-1,1), p(-1,5)));

    // point intersection
    check_do_intersect     (p( 0, 0).bbox()                , R(p(-1,1), p( 1,-1)));
    check_do_intersect     (p(-1,-1).bbox() + p(4,2).bbox(), R(p( 1,5), p(-3,-1)));

    // segment intersection
    check_do_intersect     (p(0,0).bbox() + p(4,5).bbox(), R(p(-1, 5), p(7,5)));
    check_do_intersect     (p(0,0).bbox() + p(4,5).bbox(), R(p(-1, 5), p(2,5)));
    check_do_intersect     (p(0,0).bbox() + p(4,5).bbox(), R(p( 4,-3), p(4,1)));
  }

  void C_C()
  {
    std::cout << "Circle - Circle" << std::endl;

    // no intersection
    check_no_do_intersect  (C(p(13, 14),  5), C(p( 3, 2), 1));
    check_no_do_intersect  (C(p( 2,  3),  9), C(p( 3, 4), 1)); // one contains the other

    // point intersection
    check_do_intersect     (C(p(-1, -4),  0), C(p(-1, -4), 0));
    check_do_intersect     (C(p( 3,  4), 25), C(p( 6,  8), 0));
    check_do_intersect     (C(p( 3,  4),  1), C(p( 3,  2), 1));

    // more than point intersection
    check_do_intersect     (C(p( 2,  5),  9), C(p(-3, 2), 9));
    check_do_intersect     (C(p( 1,  2),  1), C(p( 1, 2), 1)); // same cicle twice
  }

  void C_Rec()
  {
    std::cout << "Circle - Iso_rectangle" << std::endl;

    // no intersection
    check_no_do_intersect  (C(p( 2,  1),  0), Rec(p( 3,  2), p( 4,  6)));
    check_no_do_intersect  (C(p(13, 14),  5), Rec(p( 3,  2), p( 4,  6)));
    check_no_do_intersect  (C(p(13, 14),  9), Rec(p(12, 13), p(14, 15))); // circle contains the box

    // point intersection
    check_do_intersect     (C(p(-3, -2),  0), Rec(p(-3, -2), p(14,  2))); // vertex of the rectangle
    check_do_intersect     (C(p( 0,  0),  4), Rec(p( 0,  2), p(12,  9))); // vertex of the rectangle
    check_do_intersect     (C(p(-4,  3),  0), Rec(p(-6,  2), p( 7,  3))); // point on an edge
    check_do_intersect     (C(p( 0,  0),  4), Rec(p(-2,  2), p(13,  4))); // point on an edge

    // more than point intersection
    check_do_intersect     (C(p( 5,  5),  4), Rec(p( 3,  3), p( 7,  7)));
    check_do_intersect     (C(p(10, 10),  2), Rec(p( 0,  0), p(10, 10)));
    check_do_intersect     (C(p(13, 14),  3), Rec(p( 1,  1), p(30, 32))); // rectangle contains the circle
    check_do_intersect     (C(p( 0,  0), 25), Rec(p(-3, -4), p( 3,  4))); // inscribed circle of the box

    // intersection with all vertices outside of the circle
    check_do_intersect     (C(p( 0,  0),  9), Rec(p(-10, 2), p(10, 13)));
  }

  void C_L()
  {
    std::cout << "Circle - Line" << std::endl;

    // no intersection
    check_no_do_intersect  (C(p( 2, 8), 6), L(p(-3, -2), p( 2,  4)));
    check_no_do_intersect  (C(p( 2, 8), 6), L(p(-3, 22), p( 2, 14)));

    // point intersection
    check_do_intersect     (C(p( 3, 4), 0), L(p(-3,  8), p( 6,  2)));
    check_do_intersect     (C(p( 4, 3), 4), L(p( 6, -7), p( 6, -2)));
    check_do_intersect     (C(p( 4, 3), 4), L(p( 6, -7), p( 6, -9)));

    // two points intersection
    check_do_intersect     (C(p(-3, 1), 7), L(p(-1, -3), p(-6,  7)));
  }

  void C_P()
  {
    std::cout << "Circle - Point" << std::endl;

    // no intersection
    check_no_intersection  (C(p(13, 14),  5)                , p(13, 14));
    check_no_intersection  (C(p( 2,  9),  4)                , p(11, 12));

    // point intersection
    check_intersection     (C(p( 3,  4), 16)                , p( 7, 4), p(7, 4));
    check_intersection     (C(p( 0,  5), p(5,  0), p(-5, 0)), p( 0, -5), p(0, -5));
//    check_intersection     (C(p( 3,  4), p(2,  6), p( 1, 5)), p( 1, 5), p(1, 5));
    check_intersection<P>  (C(p( 0,  0), 25)                , p( 5, 0));
  }

  void C_S()
  {
    std::cout << "Circle - Segment" << std::endl;

    check_no_do_intersect  (C(p( 2, 8), 6), S(p(-3, -2), p( 2,  -4)));
    check_no_do_intersect  (C(p( 2, 8), 6), S(p(-3, 22), p( 2, 34)));
    check_no_do_intersect  (C(p( 4, 16), 18), S(p(5, 7), p(6, 8)));

    check_do_intersect     (C(p( 3, 4), 0), S(p(-3,  8), p( 6,  2)));
    check_do_intersect     (C(p( 4, 3), 4), S(p( 6, -7), p( 5, 2)));

    check_do_intersect     (C(p(-3, 1), 7), S(p(-1, -3), p(-6,  7)));

    check_do_intersect     (C(p(1, 1), 4), S(p(-1, 1), p(-6,  -8)));
    check_do_intersect     (C(p(1, 1), 4), S(p(-1, 4), p(-1,  -4)));
  }

  void C_R()
  {
    std::cout << "Circle - Ray" << std::endl;

    check_no_do_intersect  (C(p( 2, 8), 6), R(p(-3, -2), p( 2,  -4)));
    check_no_do_intersect  (C(p( 2, 8), 6), R(p(-3, 22), p( 2, 34)));

    check_do_intersect     (C(p( 3, 4), 0), R(p(-3,  8), p( 6,  2)));
    check_do_intersect     (C(p( 4, 3), 4), R(p( 6, -7), p( 5, 2)));

    check_do_intersect     (C(p(-3, 1), 7), R(p(-1, -3), p(-6,  7)));
  }

  void C_T()
  {
    std::cout << "Circle - Triangle" << std::endl;
    check_no_do_intersect  (C(p( 2, 8), 6), T(p( 6,  0), p( 8,  0), p(8, 2)));
    check_no_do_intersect  (C(p( 0, 0), 9), T(p( 1,  1), p( 1, -1), p(0, 0)));

    check_do_intersect     (C(p( 3, 4), 0), T(p(-3,  8), p( 6, 2), p(4,6)));
    check_do_intersect     (C(p( 4, 3), 4), T(p( 6, -7), p( 5, 2), p(2,-3)));

    check_do_intersect     (C(p( 4, 3), 4), T(p( 6, -7), p( 4, 6), p(2,-3)));
  }

  void L_L()
  {
    std::cout << "Line - Line" << std::endl;

    // no intersection
    check_no_intersection  (L(p(0, 0), p(10,10)), L(p(8,7), p(1, 0)));

    // point intersection
    check_intersection     (L(p(0, 0), p(10, 0)), L(p( 1, 7), p(  1,  -2)), P(1,0));
    check_intersection     (L(p(0,-1), p(10, 0)), L(p( 2, 1), p(  8,  -6)), P(3.42105,-0.657895));
    check_intersection<P>  (L(p(0, 0), p( 0, 4)), L(p(-1,-3), p(-10, -10)));

    // line intersection
    check_intersection     (L(p( 1, 1), p( 5, 5)));
    check_intersection<L>  (L(p( 1, 1), p( 5, 5)), L(p(3, 3), p( 7,  7))); // ps0 < ps1 < pt0 < pt1
    check_intersection<L>  (L(p( 5, 5), p( 1, 1)), L(p(3, 3), p( 7,  7))); // L0 & L1 have opposite directions
    check_intersection<L>  (L(p( 0, 0), p(10, 0)), L(p(0, 0), p(10,  0))); // ps0 == ps1 < pt0 == pt1
    check_intersection<L>  (L(p(10, 0), p( 0, 0)), L(p(0, 0), p(10,  0))); // ps1 == pt0 < ps0 == pt1
    check_intersection<L>  (L(p( 0, 0), p(10, 0)), L(p(1, 0), p( 8,  0))); // ps0 < ps1 < pt1 < pt0
    check_intersection<L>  (L(p( 0, 0), p(10, 0)), L(p(8, 0), p( 1,  0))); // L0 & L1 have opposite directions
  }

  void L_P()
  {
    std::cout << "Line - Point" << std::endl;

    // no intersection
    check_no_intersection  (L(p(-3,  0), p( 7, 10)), p(9, 11));

    // point intersection
    check_intersection     (L(p(-3,  0), p(-3, 10)), p(-3, 2), p(-3, 2));
    check_intersection     (L(p(-3, 10), p(-3,  7)), p(-3, 4), p(-3, 4));
    check_intersection     (L(p(-3,  4), p( 6,  4)), p( 7, 4), p( 7, 4));
    check_intersection     (L(p(-3,  4), p(-6,  4)), p( 9, 4), p( 9, 4));
    check_intersection     (L(p( 3, -1), p( 6,  1)), p(12, 5), p(12, 5));
  }

  void S_S()
  {
    std::cout << "Segment - Segment" << std::endl;

    // no intersection
    check_no_intersection  (S(p(29,  16), p( 28,   9)), S(p( 30,  12), p( 29,   6)));
    check_no_intersection  (S(p( 0,   0), p(103567,9826)), S(p(10000,3782), p(76250,83736)));

    // point intersection
    check_intersection     (S(p( 0,   0), p( 10,   0)), S(p(  1,   6), p(  1,  -3)), P(1, 0));
    check_intersection     (S(p( 0,  -1), p( 10,   0)), S(p(  2,   1), p(  8,  -6)), P(3.42105, -0.657895));
    check_intersection     (S(p( 0,   0), p( 10,   0)), S(p(  0,   0), p(  4,  -7)), P(0, 0)); // meeting at an extremity

    // segment intersection
    check_intersection     (S(p( 0,   0), p( 10,   0)));
    check_intersection     (S(p( 0,   0), p( 10,   0)), S(p(  1,   0), p(  8,   0)), S(P(  1,   0), P(  8,   0)));
    check_intersection     (S(p(68, 106), p(192, 106)), S(p(150, 106), p(255, 106)), S(P(150, 106), P(192, 106)));
    check_intersection     (S(p( 1,  10), p(  1,   2)), S(p(  1,   7), p(  1,   3)), S(P(  1,   3), P(  1,   7)));

    // exact point intersection
    check_intersection     (S(p( 3,   0), p( 3,   10)), S(p(  3,   3), p(  5,   10)), p(  3,   3));
    check_intersection     (S(p( 3,   0), p( 3,   10)), S(p(  5,   10), p(  3,   3)), p(  3,   3));
    check_intersection     (S(p( 3,   0), p( 3,   10)), S(p(  3,   3), p(  -5,   10)), p(  3,   3));
    check_intersection     (S(p( 3,   0), p( 3,   10)), S(p(  -5,   10), p(  3,   3)), p(  3,   3));
    check_intersection     (S(p( 0,   0), p( 44,   44)), S(p(  44,   44), p(  55,   55)), p(  44,   44));
    check_intersection     (S(p( 0,   0), p( 44,   44)), S(p(  55,   55), p(  44,   44)), p(  44,   44));
    check_intersection     (S(p( 44,   44), p( 0,   0)), S(p(  44,   44), p(  55,   55)), p(  44,   44));
    check_intersection     (S(p( 44,   44), p( 0,   0)), S(p(  55,   55), p(  44,   44)), p(  44,   44));
    check_intersection     (S(p( 0,   0), p( -44,   -44)), S(p(  -44,   -44), p(  -55,   -55)), p(  -44,   -44));
    check_intersection     (S(p( 0,   0), p( -44,   -44)), S(p(  -55,   -55), p(  -44,   -44)), p(  -44,   -44));
    check_intersection     (S(p( -44,   -44), p( 0,   0)), S(p(  -44,   -44), p(  -55,   -55)), p(  -44,   -44));
    check_intersection     (S(p( -44,   -44), p( 0,   0)), S(p(  -55,   -55), p(  -44,   -44)), p(  -44,   -44));

    // more segment intersection (containment)
    check_intersection     (S(p(0,0), p(4,4)), S(p(1,1), p(2,2)), S(p(1,1), p(2,2)));
    check_intersection     (S(p(0,0), p(4,4)), S(p(2,2), p(1,1)), S(p(1,1), p(2,2)));
    check_intersection     (S(p(4,4), p(0,0)), S(p(1,1), p(2,2)), S(p(1,1), p(2,2)));
    check_intersection     (S(p(4,4), p(0,0)), S(p(2,2), p(1,1)), S(p(1,1), p(2,2)));

    // more segment intersection (overlap)
    check_intersection     (S(p( 0,0), p( 2,2)), S(p(1,1), p(4,4)), S(p(1,1), p(2,2)));
    check_intersection     (S(p( 0,0), p( 2,2)), S(p(4,4), p(1,1)), S(p(1,1), p(2,2)));
    check_intersection     (S(p( 2,2), p( 0,0)), S(p(1,1), p(4,4)), S(p(1,1), p(2,2)));
    check_intersection     (S(p( 2,2), p( 0,0)), S(p(4,4), p(1,1)), S(p(1,1), p(2,2)));
    check_intersection     (S(p( 0,0), p( -2,-2)), S(p(-1,-1), p(-4,-4)), S(p(-2,-2),p(-1,-1)));
    check_intersection     (S(p( 0,0), p( -2,-2)), S(p(-4,-4), p(-1,-1)), S(p(-2,-2),p(-1,-1)));
    check_intersection     (S(p( -2,-2), p( 0,0)), S(p(-1,-1), p(-4,-4)), S(p(-2,-2),p(-1,-1)));
    check_intersection     (S(p( -2,-2), p( 0,0)), S(p(-4,-4), p(-1,-1)), S(p(-2,-2),p(-1,-1)));

    // more segment intersection (one common point)
    check_intersection     (S(p( 0,0), p( 2,2)), S(p(1,1), p(2,2)), S(p(1,1), p(2,2)));
    check_intersection     (S(p( 0,0), p( 2,2)), S(p(2,2), p(1,1)), S(p(1,1), p(2,2)));
    check_intersection     (S(p( 2,2), p( 0,0)), S(p(1,1), p(2,2)), S(p(1,1), p(2,2)));
    check_intersection     (S(p( 2,2), p( 0,0)), S(p(2,2), p(1,1)), S(p(1,1), p(2,2)));
    check_intersection     (S(p( 0,0), p( -2,-2)), S(p(-1,-1), p(-2,-2)), S(p(-2,-2),p(-1,-1)));
    check_intersection     (S(p( 0,0), p( -2,-2)), S(p(-2,-2), p(-1,-1)), S(p(-2,-2),p(-1,-1)));
    check_intersection     (S(p( -2,-2), p( 0,0)), S(p(-1,-1), p(-2,-2)), S(p(-2,-2),p(-1,-1)));
    check_intersection     (S(p( -2,-2), p( 0,0)), S(p(-2,-2), p(-1,-1)), S(p(-2,-2),p(-1,-1)));

    // more segment intersection (two identical points)
    check_intersection     (S(p( 1,1), p( 2,2)), S(p(1,1), p(2,2)), S(p(1,1), p(2,2)));
    check_intersection     (S(p( 1,1), p( 2,2)), S(p(2,2), p(1,1)), S(p(1,1), p(2,2)));
    check_intersection     (S(p( 2,2), p( 1,1)), S(p(1,1), p(2,2)), S(p(1,1), p(2,2)));
    check_intersection     (S(p( 2,2), p( 1,1)), S(p(2,2), p(1,1)), S(p(1,1), p(2,2)));
    check_intersection     (S(p( -1,-1), p( -2,-2)), S(p(-1,-1), p(-2,-2)), S(p(-2,-2),p(-1,-1)));
    check_intersection     (S(p( -1,-1), p( -2,-2)), S(p(-2,-2), p(-1,-1)), S(p(-2,-2),p(-1,-1)));
    check_intersection     (S(p( -2,-2), p( -1,-1)), S(p(-1,-1), p(-2,-2)), S(p(-2,-2),p(-1,-1)));
    check_intersection     (S(p( -2,-2), p( -1,-1)), S(p(-2,-2), p(-1,-1)), S(p(-2,-2),p(-1,-1)));

    //check determinism of crossing point
    P p1a(-122.37323046264295, 37.7435274415764);
    P p1b(-122.3711959178425,  37.74348027376899);
    P p2a(-122.37130249711004, 37.74203327688176);
    P p2b(-122.3722247426892,  37.74401427059434);
    std::set<double> ds;
    auto test = [&ds](S s1, S s2)
    {
      P i = std::get<P>(*CGAL::intersection(s1,s2));
      ds.insert(CGAL::to_double(i.x())); ds.insert(CGAL::to_double(i.y()));
      assert(ds.size()==2);
    };
    test(S(p1a,p1b), S(p2a,p2b));
    test(S(p1a,p1b), S(p2b,p2a));
    test(S(p1b,p1a), S(p2b,p2a));
    test(S(p1b,p1a), S(p2a,p2b));
    test(S(p2a,p2b), S(p1a,p1b));
    test(S(p2b,p2a), S(p1a,p1b));
    test(S(p2b,p2a), S(p1b,p1a));
    test(S(p2a,p2b), S(p1b,p1a));
  }

  void R_R()
  {
    std::cout << "Ray - Ray" << std::endl;

    // no intersection
    check_no_intersection  (R(p( 3,   4), p(  5,   7)), R(p(  2,   0), p(  2,   4)));

    // point intersection
    check_intersection     (R(p( 2,  -1), p(  2,   1)), R(p(  1,   3), p(  2,   3)), P(2, 3));
    check_intersection     (R(p( 0,  -1), p( 10,   0)), R(p(  2,   1), p(  8,  -6)), P(3.42105, -0.657895));
    check_intersection     (R(p( 0,   0), p( 10,   0)), R(p(  1,   6), p(  1,  -3)), P(1, 0));
    check_intersection     (R(p( 0,   0), p( 10,   0)), R(p(  15,  0), p( 18,  30)), P(15, 0)); // R1's source on R0's path
    check_intersection     (R(p( 0,   0), p( 10,   0)), R(p(  0,   0), p( -3,   4)), P(0, 0)); // same source but different directions

    // segment intersection (same supporting line, different directions)
    check_intersection<S>  (R(p( 2,   4), p(  6,   1)), R(p(  14, -5), p(10, -2)));
    check_intersection<S>  (R(p( 2,   4), p( 10,  -2)), R(p(   6,  1), p(-2,  7)));

    // ray intersection
    check_intersection     (R(p( 0,   0), p( 10,   0)));
    check_intersection<R>  (R(p( 0,   0), p( 10,   0)), R(p(  -1,  0), p(0,   0))); // R0 'runs' into R1's source
  }

  void R_S()
  {
    std::cout << "Ray - Segment" << std::endl;

    // no intersection
    check_no_intersection  (S(p( 2,  -1), p(  2,   1)), R(p(  1,   3), p(  2,   3)));

    // point intersection
    check_intersection     (S(p( 0,  -1), p( 10,   0)), R(p(  2,   1), p(  8,  -6)), P(3.42105, -0.657895));
    check_intersection     (S(p( 0,   0), p( 10,   0)), R(p(  1,   6), p(  1,  -3)), P(1, 0));
    check_intersection     (S(p( 0,   0), p( 10,   0)), R(p(  1,   6), p(  1,  -3)), P(1, 0));
    check_intersection     (S(p( 0,   0), p( 10,   0)), R(p(  0,   0), p(-10,   4)), P(0, 0)); // start of ray is exactly on the segment
    check_intersection     (S(p( 0,   0), p( 10,   0)), R(p(  4,   0), p(-10,   4)), P(4, 0)); // start of ray is a segment extremity

    // segment intersection
    check_intersection     (S(p(  0,   0), p(  1,  0)), R(p(  0,   0), p( 10,   0)), S(P(0, 0), P(1,0)));
    check_intersection<S>  (S(p(  3,   7), p(  2,  5)), R(p(  1,   3), p(  4,   9)));
  }

  void L_R()
  {
    std::cout << "Line - Ray" << std::endl;

    // no intersection
    check_no_intersection  (L(p( 2,  -1), p(  2,   1)), R(p(  1,   -3), p(  1,   3)));

    // point intersection
    check_intersection     (L(p( 2,  -1), p(  2,   1)), R(p(  1,   3), p(  2,   3)), P(2, 3));
    check_intersection     (L(p( 0,  -1), p( 10,   0)), R(p(  2,   1), p(  8,  -6)), P(3.42105, -0.657895));
    check_intersection     (L(p( 0,   0), p( 10,   0)), R(p(  1,   6), p(  1,  -3)), P(1, 0));
    check_intersection     (L(p( 0,   0), p(  4,   5)), R(p(  1,   8), p(  7,   2)), P(4, 5)); // ray starts on the line

    // ray intersection
    check_intersection<R>  (L(p( 2,  -1), p(  2,   1)), R(p(  2,   -3), p(  2,   3)));
    check_intersection<R>  (L(p( 2,  -1), p(  2,   1)), R(p(  2,    3), p(  2,  -3))); // opposite direction
  }

  void L_S()
  {
    std::cout << "Line - Segment" << std::endl;

    // no intersection
    check_no_intersection  (S(p( 2,  -1), p(  2,   1)), L(p(  1,   3), p(  2,   3)));

    // point intersection
    check_intersection     (S(p( 0,  -1), p( 10,   0)), L(p(  2,   1), p(  8,  -6)), P(3.42105, -0.657895));
    check_intersection     (S(p( 0,   0), p( 10,   0)), L(p(  1,   6), p(  1,  -3)), P(1, 0));
    check_intersection     (S(p( 0,   0), p( 10,   0)), L(p(  1,   6), p(  2,  12)), P(0, 0));

    // segment intersection
    check_intersection<S>  (S(p(-3,   5), p( 12,   1)), L(p( 12,   1), p( 27,  -3)));
    check_intersection<S>  (S(p(-3,   5), p( 12,   1)), L(p(-18,   9), p( 27,  -3)));
    check_intersection<S>  (S(p(-3,   5), p( 12,   1)), L(p( 27,  -3), p(-18,   9)));
  }

  void T_T()
  {
    std::cout << "Triangle - Triangle" << std::endl;

    // no intersection
    check_no_intersection  (T(p( -10,-10), p(  0,  10), p( 20, -5)), T(p(   90, -10), p(100,  10), p(120, -5)));

    // point intersection
    check_intersection     (T(p( -10,  0), p( 10,   0), p(0,  3)), T(p( -12,   3), p( 12,   3), p(1, 5)), P(0, 3)); // intersection on an edge
    check_intersection     (T(p( -25, -4), p( 13, -18), p(0,  3)), T(p( -12,   5), p( 12,   5), p(0, 3)), P(0, 3)); // intersection at a vertex

    // segment intersection
    check_intersection<S>  (T(p( -10,  0), p( 10,   0), p(0, -3)), T(p(  -8,   0), p( 12,   0), p(1, 5)));
    check_intersection<S>  (T(p( -10,  0), p( 10,   0), p(0, -3)), T(p(  -8,   0), p(  8,   0), p(1, 5)));
    check_intersection<S>  (T(p( -10,  0), p( 10,   0), p(0, -3)), T(p( -10,   0), p( 10,   0), p(1, 5)));
    check_intersection<S>  (T(p( -10,  0), p( 10,   0), p(0, -3)), T(p(  10,   0), p(-10,   0), p(1, 5)));
    check_intersection<S>  (T(p( -10,  0), p( 10,   0), p(0, -3)), T(p( -10,   0), p( 10,   0), p(1, 5)));
    check_intersection<S>  (T(p( -10,  0), p( 10,   0), p(0, -3)), T(p( -12,   0), p( 12,   0), p(1, 5)));

    // triangle intersection
    check_intersection<T>  (T(p(   0, 10), p(-10, -10), p( 20,  -5)), T(p( -10, -10), p(  0,  10), p( 20, -5)));
    check_intersection<T>  (T(p( -12,  1), p(  5,   3), p( -7, -15)), T(p(  29,  -2), p(  0, -13), p(1, 21)));

    // polygon intersection
    Pol pol0;
    pol0.push_back(P( 8, 4 ));
    pol0.push_back(P( 0, 10 ));
    pol0.push_back(P( -5.11111, -0.222222 ));
    pol0.push_back(P(-6, -4));
    check_intersection     (T(p(   0, 10), p(-10, -10), p( 20, -5)), T(p(   2,  30), p( -6,  -4), p(15, 8)), pol0, false);
    check_intersection     (T(p( -10,-10), p(  0,  10), p( 20, -5)), T(p(   2,  30), p( -6,  -4), p(15, 8)), pol0, false);

    Pol pol2;
    pol2.push_back(P( 10.2222, 2.33333 ));
    pol2.push_back(P( 1.96923, 8.52308 ));
    pol2.push_back(P( -0.680851, 8.6383 ));
    pol2.push_back(P( -5.94872, -1.89744 ));
    pol2.push_back(P( -3.96178, -8.99363 ));
    pol2.push_back(P( 3.5, -7.75 ));
    check_intersection     (T(p( -10,-10), p(  0,  10), p( 20, -5)), T(p(   -9,   9), p( 14,   8), p(-2, -16)), pol2, false);
  }

  void L_T()
  {
    std::cout << "Line - Triangle" << std::endl;

    // no intersection
    check_no_intersection  (L(p(-10,   0), p( 10,   0)), T(p(  -12,   3), p(  12,   3), p(1,  5)));

    // point intersection
    check_intersection<P>  (L(p( -8,  30), p(  8,  30)), T(p(    2,  30), p(  14,   2), p(-7, -2)));
    check_intersection<P>  (L(p( -8,  30), p( -7,  30)), T(p(    2,  30), p(  14,   2), p(-7, -2)));

    // segment intersection
    check_intersection<S>  (L(p( -1,  -1), p(  0,  -1)), T(p(  -10,   0), p(  10,   0), p(0, -4)));
    check_intersection<S>  (L(p( -1,  -1), p(  0,  -1)), T(p(  -10,   0), p(  15,   2), p(0, -4)));
    check_intersection<S>  (L(p(-10,   0), p( 10,   0)), T(p(   -8,   0), p(   8,   0), p(1,  5)));
    check_intersection<S>  (L(p(-10,   0), p( 10,   0)), T(p(  -12,   0), p(  12,   0), p(1,  5)));
    check_intersection<S>  (L(p(  0,  10), p(-10, -10)), T(p(    2,  30), p(  -6,  -4), p(15, 8)));
    check_intersection<S>  (L(p(-12,   1), p(  5,   3)), T(p(   29,  -2), p(   0, -13), p( 1, 21)));
    check_intersection<S>  (L(p(-10, -10), p(  0,  10)), T(p(    2,  30), p(  -6,  -4), p(15,  8)));
    check_intersection<S>  (L(p(-10, -10), p(  0,  10)), T(p(   -9,   9), p(  14,   8), p(-2,-16)));
  }

  void R_T()
  {
    std::cout << "Ray - Triangle" << std::endl;

    // no intersection
    check_no_intersection  (R(p(-10,   0), p( 10,   0)), T(p(  -12,   3), p(  12,   3), p(1,  5)));

    // point intersection
    check_intersection<P>  (R(p(-2, -16), p(  4,  -20)), T(p(   -9,   9), p(  14,   8), p(-2, -16)));
    check_intersection<P>  (R(p(-8, -1),  p(  -8, -12)), T(p(  -12,   2), p(  10,   3), p(-4,  -4)));
    check_intersection<P>  (R(p(-8, 30),  p(   4,  30)), T(p(    2,  30), p(  14,   2), p(-7,  -2)));
    check_intersection<P>  (R(p(-8, 30),  p(  -7,  30)), T(p(    2,  30), p(  14,   2), p(-7,  -2)));

    // segment intersection
    check_intersection<S>  (R(p( -1,  -1), p(  0,  -1)), T(p(  -10,   0), p(  10,   0), p(0, -4)));
    check_intersection<S>  (R(p( -1,  -1), p(  0,  -1)), T(p(  -10,   0), p(  15,   2), p(0, -4)));
    check_intersection<S>  (R(p(-10,   0), p( 10,   0)), T(p(   -8,   0), p(   8,   0), p(1,  5)));
    check_intersection<S>  (R(p(-10,   0), p( 10,   0)), T(p(  -12,   0), p(  12,   0), p(1,  5)));
    check_intersection<S>  (R(p(  0,  10), p(-10, -10)), T(p(    2,  30), p(  -6,  -4), p(15, 8)));
    check_intersection<S>  (R(p(-12,   1), p(  5,   3)), T(p(   29,  -2), p(   0, -13), p( 1, 21)));
    check_intersection<S>  (R(p(-10, -10), p(  0,  10)), T(p(    2,  30), p(  -6,  -4), p(15,  8)));
    check_intersection<S>  (R(p(-10, -10), p(  0,  10)), T(p(   -9,   9), p(  14,   8), p(-2,-16)));
  }

  void S_T()
  {
    std::cout << "Segment - Triangle" << std::endl;

    // no intersection
    check_no_intersection  (S(p(-10, -10), p(  0,  10)), T(p(   90, -10), p( 100,  10), p(120, -5)));
    check_no_intersection  (S(p(-10,   0), p( 10,   0)), T(p(  -12,   3), p(  12,   3), p(1,  5)));

    // point intersection
    check_intersection<P>  (S(p(-2, -16), p(  4,  -20)), T(p(   -9,   9), p(  14,   8), p(-2, -16)));
    check_intersection<P>  (S(p(-8,  -1), p(  -8, -12)), T(p(  -12,   2), p(  10,   3), p(-4,  -4)));
    check_intersection<P>  (S(p(-8, -12), p(  -8,  -1)), T(p(  -12,   2), p(  10,   3), p(-4,  -4)));
    check_intersection<P>  (S(p(-8,  30), p(   1,  30)), T(p(   -2,  30), p(  14,   2), p(-7,  -2)));

    // segment intersection
    check_intersection<S>  (S(p( -1,  -1), p(  0,  -1)), T(p(  -10,   0), p(  10,   0), p(0, -4)));
    check_intersection<S>  (S(p( -1,  -1), p(  0,  -1)), T(p(  -10,   0), p(  15,   2), p(0, -4)));
    check_intersection<S>  (S(p(-10,   0), p( 10,   0)), T(p(   -8,   0), p(   8,   0), p(1,  5)));
    check_intersection<S>  (S(p(-10,   0), p( 10,   0)), T(p(  -12,   0), p(  12,   0), p(1,  5)));
    check_intersection<S>  (S(p(  0,  10), p(-10, -10)), T(p(    2,  30), p(  -6,  -4), p(15, 8)));
    check_intersection<S>  (S(p(-12,   1), p(  5,   3)), T(p(   29,  -2), p(   0, -13), p( 1, 21)));
    check_intersection<S>  (S(p(-10, -10), p(  0,  10)), T(p(    2,  30), p(  -6,  -4), p(15,  8)));
    check_intersection<S>  (S(p(-10, -10), p(  0,  10)), T(p(   -9,   9), p(  14,   8), p(-2,-16)));
  }

  void P_P()
  {
    std::cout << "Point - Point" << std::endl;

    check_no_intersection<P>(p(8, 4), p(-4, 8));
    check_intersection<P>   (p(8, 4), p( 8, 4));
  }

  void P_R()
  {
    std::cout << "Point - Ray" << std::endl;

    // no intersection
    check_no_intersection   (R(p(-1, -1), p(0, 12)), p( 9, -1));

    // point intersection
    check_intersection      (R(p(-3,  0), p(7, 10)), p(-3,  0), p(-3,  0)); // origin of the ray
    check_intersection      (R(p(-3,  0), p(7, 10)), p( 9, 12), p( 9, 12));
  }

  void P_S()
  {
    std::cout << "Point - Segment" << std::endl;

    // no intersection
    check_no_intersection   (S(p(-1, -1), p( 0, 12)), p( 9, -1));
    check_no_intersection   (S(p(-2, -3), p( 0,  4)), p( 4, 18));

    // point intersection
    check_intersection      (S(p(-2, -3), p( 4, 18)), p( 0,  4), p( 0,  4));
    check_intersection      (S(p(-3,  0), p( 7, 10)), p(-3,  0), p(-3,  0));
    check_intersection      (S(p( 3,  2), p(-2, 11)), p(-2, 11), p(-2, 11));
  }

  void P_T()
  {
    std::cout << "Point - Triangle" << std::endl;

    // no intersection
    check_no_intersection  (p(  8,   6), T(p(    4,   0), p(  12,   4), p(-4,  8)));

    // point intersection
    check_intersection<P>  (p(  8,   4), T(p(    4,   0), p(  12,   4), p(-4,  8)));
    check_intersection<P>  (p(  8,   5), T(p(    4,   0), p(  12,   4), p(-4,  8)));
    check_intersection<P>  (p(  4,   0), T(p(    4,   0), p(  12,   4), p(-4,  8)));
    check_intersection<P>  (p(  12,  4), T(p(    4,   0), p(  12,   4), p(-4,  8)));
    check_intersection<P>  (p(  -4,  8), T(p(    4,   0), p(  12,   4), p(-4,  8)));
  }

  void Rec_L()
  {
    std::cout << "Iso_rectangle - Line" << std::endl;

    // no intersection
    check_no_intersection  (L(p( 18,  6), p( 16,  4)), Rec(p( 2,  0), p(6,  3)));

    // point intersection
    check_intersection     (L(p( -1,  0), p( 4,   5)), Rec(p( 0,  0), p(1,  1)), P(0, 1));
    check_intersection<P>  (L(p( -5, 10), p(-1, -12)), Rec(p(-3, -1), p(2,  14)));

    // segment intersection
    check_intersection<S>  (L(p( 18,  6), p( 0,   0)), Rec(p( 2,  0), p(6,  3)));
    check_intersection<S>  (L(p( 18,  6), p( 0,   0)), Rec(p( 2,  0), p(6,  3)));
    check_intersection<S>  (L(p( 2,  14), p( 2, -14)), Rec(p( 2,  0), p(6,  3)));
    check_intersection<S>  (L(p( 6,   1), p( 6,   2)), Rec(p( 2,  0), p(6,  3)));
    check_intersection<S>  (L(p(-1,   3), p(-2,   3)), Rec(p( 2,  0), p(6,  3)));
  }

  void Rec_P()
  {
    std::cout << "Iso_rectangle - Point" << std::endl;

    // no intersection
    check_no_intersection  (Rec(p(-2, -6), p( 6, 3)), p( 7,  2));
    check_no_intersection  (Rec(p(-2, -6), p( 6, 3)), p(-2, -7));

    // point intersection
    check_intersection     (Rec(p(-1,  4), p(-1, 4)), p(-1, 4), p(-1, 4)); // degenerate rectangle (0d)
    check_intersection     (Rec(p(-2,  4), p(-2, 7)), p(-2, 6), p(-2, 6)); // degenerate rectangle (1d)
    check_intersection     (Rec(p(-2,  4), p(-2, 7)), p(-2, 7), p(-2, 7)); // degenerate rectangle (1d)
    check_intersection     (Rec(p(-3,  0), p( 4, 2)), p(-3, 2), p(-3, 2)); // on vertex
    check_intersection     (Rec(p( 7,  8), p( 9, 9)), p( 8, 9), p( 8, 9)); // on edge
    check_intersection     (Rec(p(-2,  0), p( 6, 7)), p( 1, 1), p( 1, 1)); // within
  }

  void Rec_R()
  {
    std::cout << "Iso_rectangle - Ray" << std::endl;

    // no intersection
    check_no_intersection  (R(p( 18,  6), p( 16,  4)), Rec(p( 2,  0), p(6,  3)));

    // point intersection
    check_intersection     (R(p( -1,  0), p( 4,   5)), Rec(p( 0,  0), p(1,  1)), P(0, 1));
    check_intersection     (R(p(  0,  2), p(-4,  14)), Rec(p( 0,  0), p(5,  2)), P(0, 2));
    check_intersection<P>  (R(p( -5, 10), p(-1, -12)), Rec(p(-3, -1), p(2,  14)));

    // segment intersection
    check_intersection<S>  (R(p( 18,  6), p( 0,   0)), Rec(p( 2,  0), p(6,  3)));
    check_intersection<S>  (R(p( 2,  14), p( 2, -14)), Rec(p( 2,  0), p(6,  3)));
    check_intersection<S>  (R(p( 2,   1), p( 2,   4)), Rec(p( 2,  0), p(6,  3)));
    check_intersection<S>  (R(p(-2,   3), p(-1,   3)), Rec(p( 2,  0), p(6,  3)));
  }

  void Rec_S()
  {
    std::cout << "Iso_rectangle - Segment" << std::endl;

    // no intersection
    check_no_intersection  (S(p( 18,  6), p( 16,  4)), Rec(p( 2,  0), p(6,  3)));

    // point intersection
    check_intersection     (S(p( -1,  0), p( 4,   5)), Rec(p( 0,  0), p(1,  1)), P(0, 1));
    check_intersection     (S(p(  0,  2), p(-4,  14)), Rec(p( 0,  0), p(5,  2)), P(0, 2));
    check_intersection<P>  (S(p( -5, 10), p(-1, -12)), Rec(p(-3, -1), p(2,  14)));

    // segment intersection
    check_intersection<S>  (S(p( 18,  6), p( 0,   0)), Rec(p( 2,  0), p(6,  3)));
    check_intersection<S>  (S(p( 2,  14), p( 2, -14)), Rec(p( 2,  0), p(6,  3)));
    check_intersection<S>  (S(p( 6,   1), p( 6,   2)), Rec(p( 2,  0), p(6,  3)));
    check_intersection<S>  (S(p(-2,   3), p( 6,   3)), Rec(p( 2,  0), p(6,  3)));
  }

  void Rec_Rec()
  {
    std::cout << "Iso_rectangle - Iso_rectangle" << std::endl;

    // no intersection
    check_no_intersection  (Rec(p( -4, -12), p(12, 23)), Rec(p( -4,  24), p(  5, 26)));
    check_no_intersection  (Rec(p( -4, -12), p(12, 23)), Rec(p( 13, -24), p( 14, 15)));

    // point intersection
    check_intersection<Rec> (Rec(p( 10,  12), p(30, 40)), Rec(p(  30,   40), p( 31, 42))/*, p(30, 40)*/);
    check_intersection<Rec> (Rec(p( 10,  12), p(30, 40)), Rec(p(  30,  -13), p( 33, 12))/*, p(30, 12)*/);

    // segment intersection
    check_intersection<Rec>  (Rec(p( 3,  5), p(4, 6)), Rec(p( 2, 4), p( 6, 5)));
    check_intersection<Rec>  (Rec(p( 3,  5), p(4, 6)), Rec(p( 1, 1), p( 3, 8)));
    check_intersection<Rec>  (Rec(p( 3,  5), p(9, 9)), Rec(p( 1, 4), p( 3, 8)));

    // Iso rectangle intersection
    check_intersection     (Rec(p( 10,  12), p(30, 40)));
    check_intersection     (Rec(p( 10,  12), p(30, 40)), Rec(p(  25,   40), p( 26,  103)), Rec(P(25, 40), P(26, 40)));
  }

  void Rec_T()
  {
    std::cout << "Iso_rectangle - Triangle" << std::endl;

    // no intersection
    check_no_intersection  (Rec(p( 10,  12), p(30, 40)), T(p(    4,   0), p(  12,   4), p(-4,  8)));
    check_no_intersection  (Rec(p(  2,   0), p( 3,  1)), T(p(    0,   3), p(   3,   3), p( 0,  0)));

    // point intersection
    check_intersection<P>  (Rec(p( 0,  0), p(1, 1)), T(p(  -1,  0), p(  0,  0), p( 0,  -1))); // intersection at a vertex
    check_intersection<P>  (Rec(p( 0,  0), p(1, 1)), T(p(   0,  0), p( -1,  0), p( 0,  -1))); // inversed triangle
    check_intersection<P>  (Rec(p( 0,  0), p(1, 1)), T(p(  -1,  0), p(  2,  3), p(-4,   6))); // intersection on an edge of the triangle
    check_intersection     (Rec(p( 0,  0), p(3, 3)), T(p( -10,  0), p(  0,  2), p(-1,   4)), p(0, 2)); // intersection on an edge of the iso rectangle

    // segment intersection
    check_intersection<S>  (Rec(p( 0,  0), p(3, 3)), T(p( -10,  0), p(  0,  0), p(  0, 3)));
    check_intersection<S>  (Rec(p(-2,  2), p(3, 3)), T(p( -15, 12), p( -1,  3), p(  2, 3)));
    check_intersection<S>  (Rec(p(-2,  2), p(3, 3)), T(p( -15, 12), p( -4,  3), p(  2, 3)));
    check_intersection<S>  (Rec(p(-2,  2), p(3, 3)), T(p( -15, 12), p( -4,  3), p( 15, 3)));

    // triangle intersection
    check_intersection<T>  (Rec(p( 0,  0), p(1, 1)), T(p(  -1,  0), p( -1,  2), p(  2,  2)));
    check_intersection<T>  (Rec(p( 0,  0), p(1, 1)), T(p(  -1,  0), p(  2,  2), p( -1,  2)));
    check_intersection<T>  (Rec(p( 0,  0), p(2, 2)), T(p(   0,  0), p(  1,  0), p(  0,  1)));
    check_intersection<T>  (Rec(p( 0,  0), p(3, 3)), T(p(   1,  1), p(  2,  1), p(  1,  2)));

    // polygon intersection
    check_intersection<Pol>(Rec(p(   0,   0), p(  1,   1)), T(p( -1, -2), p( -1,   2), p( 5,   2)));
    check_intersection<Pol>(Rec(p( 100, 100), p(200, 200)), T(p(150, 50), p(250, 170), p(50, 170)));
  }

  void run()
  {
    std::cout << "2D Intersection tests with Kernel: " << typeid(K).name() << std::endl;

    B_C();
    B_L();
    B_P();
    B_R();

    C_C();
    C_Rec();
    C_L();
    C_P();
    C_S();
    C_R();
    C_T();

    Rec_Rec();
    Rec_L();
    Rec_P();
    Rec_R();
    Rec_S();
    Rec_T();

    L_L();
    L_P();
    L_R();
    L_S();
    L_T();

    P_P();
    P_R();
    P_S();
    P_T();

    R_R();
    R_S();
    R_T();

    S_S();
    S_T();

    T_T();
  }

private:
  K kernel;
};

int main()
{
  CGAL::Set_ieee_double_precision pfr;

  Test< CGAL::Simple_cartesian<CGAL::internal::Exact_field_selector<double>::Type > >().run();
  Test< CGAL::Cartesian<double>   >().run();
  Test< CGAL::Homogeneous<CGAL::internal::Exact_field_selector<double>::Type > >().run();
  Test< CGAL::Exact_predicates_inexact_constructions_kernel >().run();
  Test< CGAL::Exact_predicates_exact_constructions_kernel >().run();
}
