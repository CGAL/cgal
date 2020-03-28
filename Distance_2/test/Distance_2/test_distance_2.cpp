// 2D distance tests.

#ifdef NDEBUG
#undef NDEBUG //this testsuite requires NDEBUG to be not defined
#endif


#include <CGAL/Simple_cartesian.h>
#include <CGAL/Simple_homogeneous.h>

#include <vector>
#include <iostream>
#include <cassert>

const double epsilon = 0.001;

struct randomint {
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

template < typename K >
struct Test {

  typedef typename K::FT              FT;
  typedef CGAL::Point_2< K >          P;
  typedef CGAL::Line_2< K >           L;
  typedef CGAL::Segment_2< K >        S;
  typedef CGAL::Ray_2< K >            R;
  typedef CGAL::Triangle_2< K >       T;
  typedef CGAL::Iso_rectangle_2< K >  Rec;


  template < typename Type >
  bool approx_equal_nt(const Type &t1, const Type &t2)
  {
        if (t1 == t2)
                return true;
        if (CGAL::abs(t1 - t2) / (CGAL::max)(CGAL::abs(t1), CGAL::abs(t2)) < epsilon)
                return true;
        std::cout << " Approximate comparison failed between : " << t1 << "  and  " << t2 << "\n";
        return false;
  }

  template < typename O1, typename O2 >
  void check_squared_distance(const O1& o1, const O2& o2, const FT& result)
  {
        assert(approx_equal_nt(CGAL::squared_distance(o1, o2), result));
        assert(approx_equal_nt(CGAL::squared_distance(o2, o1), result));
  }


  P p(int x, int y)
  {
    int w = ri.next();
    return P(to_nt(x*w), to_nt(y*w), to_nt(w));
  }

  void P_P()
  {
    std::cout << "Point - Point\n";
    check_squared_distance (p(0, 0), p(0, 0), 0);
    check_squared_distance (p(1, 1), p(0, 0), 2);
  }

  void P_S()
  {
    std::cout << "Point - Segment\n";
    check_squared_distance (p(1, 1), S(p(4, 0), p(-3, -1)), 2);
  }

  void S_S()
  {
    std::cout << "Segment - Segment\n";
    check_squared_distance (S(p(2,     0), p(0,   2)), S(p(1, 1),      p(4, 0)), 0);
    check_squared_distance (S(p(10,    0), p(0,  10)), S(p(6, 6),      p(20, 0)), 2);
    check_squared_distance (S(p(-10, -13), p(0,  10)), S(p(10, 5),     p(70, -30)), 124.642);
    check_squared_distance (S(p( 0,    0), p(6,  -2)), S(p(-1, 4),     p(8, 6)), 16.9882);
    check_squared_distance (S(p(-8,   -7), p( 11, 6)), S(p(  23, -27), p(-17, 16)), 0);
    check_squared_distance (S(p( 5,    0), p(  0, 0)), S(p(   1,   1), p(  2,  1)), 1);
    check_squared_distance (S(p( 0,    0), p(  5, 0)), S(p(   6,   1), p(  8,  1)), 2);
    check_squared_distance (S(p( 0,    0), p(  0,-3)), S(p(   1,   4), p(  1,  7)), 17);
    check_squared_distance (S(p( 0,    0), p(  0, 0)), S(p(   8,   1), p(  6,  1)), 37);
    check_squared_distance (S(p( 0,    0), p(  5, 0)), S(p(   8,   1), p(  6,  1)),  2);
    check_squared_distance (S(p( 0,    0), p(  5, 0)), S(p(   8,   0), p( 12,  0)),  9);
    check_squared_distance (S(p( 0,    0), p( 50, 0)), S(p(  80,  10), p(120, 11)), 1000);
    check_squared_distance (S(p( 0,    0), p(  1, 0)), S(p(   2,   1), p(  2, -1)),  1);
    check_squared_distance (S(p( 4,    0), p( -3,-1)), S(p(   1,   1), p(  2, 11)),  2);
    check_squared_distance (S(p( 3,    4), p(  7, 7)), S(p(   7,   0), p(  6,  5)),  1);
    check_squared_distance (S(p(-1,    1), p(  3, 4)), S(p(   7,   0), p(  6,  5)),  9.84615);
    check_squared_distance (S(p( 0,    0), p(  5, 0)), S(p(   1,   1), p(  6,  1)),  1);
    check_squared_distance (S(p( 0,    0), p(  5, 0)), S(p(   1,   1), p(  2,  1)),  1);
    check_squared_distance (S(p( 5,    0), p(  8, 0)), S(p(   1,   1), p(  2,  1)),  10);
    check_squared_distance (S(p( 5,    0), p(  0, 0)), S(p(   1,   1), p(  2,  1)),  1);
  }

  void P_R()
  {
    std::cout << "Point - Ray\n";
    check_squared_distance (p(0, 0), R(p(30,10), p(  25,   9)), 15.3846);
    check_squared_distance (p(1, 1), R(p( 2, 1), p(  10,   1)), 1);
    check_squared_distance (p(1, 1), R(p( 2, 1), p( -10,   1)), 0);
    check_squared_distance (p(1, 1), R(p(-2, 1), p( -10,   1)), 9);
    check_squared_distance (p(1, 1), R(p(-2, 1), p(  10,   1)), 0);
    check_squared_distance (p(0, 0), R(p( 2, 1), p(  10,   1)), 5);
    check_squared_distance (p(0, 0), R(p( 2, 1), p( -10,   1)), 1);
    check_squared_distance (p(0, 0), R(p(-2, 1), p( -10,   1)), 5);
    check_squared_distance (p(0, 0), R(p(-2, 1), p(  10,   1)), 1);
    check_squared_distance (p(0, 0), R(p( 3, 1), p(  10,   3)), 10);
  }

  void R_R()
  {
    std::cout << "Ray - Ray\n";
    check_squared_distance (R(p(0, 0), p(10, 0)), R(p( 9, 9), p( 20, 20)), 81);
    check_squared_distance (R(p(10, 0), p(20, 0)), R(p( 9, 9), p( 20, 20)), 82);
    check_squared_distance (R(p(0, 0), p(10, 0)), R(p( 11, 11), p( 20, 20)), 121);
    check_squared_distance (R(p(11, 11), p(20, 20)), R(p( 0, 0), p( 10, 0)), 121);
    check_squared_distance (R(p(2, 0), p(0, 2)), R(p( 1, 1), p( 4, 0)), 0);
    check_squared_distance (R(p(10, 0), p(0, 10)), R(p( 6, 6), p( 20, 0)), 2);
    check_squared_distance (R(p(-10, -13), p(0, 10)), R(p( 10, 5), p( 70, -30)), 124.642);
    check_squared_distance (R(p( 0,    0), p(30, -10)), R(p( -5, 20), p( 40, 30)), 424.706);
    check_squared_distance (R(p( 0,    0), p( 1,   0)), R(p( -1,  1), p(  1,  1)), 1);
    check_squared_distance (R(p( 3,    4), p( 7,   7)), R(p(  7,  0), p(  6,  5)), 0);
    check_squared_distance (R(p(-1,    1), p( 3,   4)), R(p(  7,  0), p(  6,  5)), 0);
    check_squared_distance (R(p( 0,    0), p( 5,   0)), R(p(  1,  1), p(  6,  1)), 1);
    check_squared_distance (R(p( 0,    0), p( 5,   0)), R(p(  1,  1), p(  2,  1)), 1);
    check_squared_distance (R(p( 5,    0), p( 8,   0)), R(p(  1,  1), p(  2,  1)), 1);
    check_squared_distance (R(p( 5,    0), p( 0,   0)), R(p(  1,  1), p(  2,  1)), 1);
    check_squared_distance (R(p( 0,    0), p( 5,   0)), R(p(  6,  1), p(  8,  1)), 1);
    check_squared_distance (R(p( 0,    0), p( 0,  -3)), R(p(  1,  4), p(  1,  7)), 17);
    check_squared_distance (R(p( 0,    0), p( 5,   0)), R(p(  8,  1), p(  6,  1)), 1);
    check_squared_distance (R(p( 0,    0), p( 1,   0)), R(p( -1,  1), p( -2,  1)), 2);
    check_squared_distance (R(p( 0,    0), p( 1,   0)), R(p(  1,  1), p( -2,  1)), 1);
    check_squared_distance (R(p( 0,    0), p( 1,   0)), R(p( -1, -1), p(  1,  1)), 0);
    check_squared_distance (R(p(-1,   -1), p( 1,   1)), R(p(  0,  0), p(  1,  0)), 0);
    check_squared_distance (R(p(-1,   -1), p( 1,   1)), R(p( -2,  0), p(  0,  0)), 0);
    check_squared_distance (R(p(-8,   -7), p(11,   6)), R(p( 23,-27), p(-17, 16)), 0);
    check_squared_distance (R(p( 0,    0), p( 1,   0)), R(p(  2,  1), p(  2, -1)), 0);
    check_squared_distance (R(p( 4,    0), p(-3,  -1)), R(p(  1,  1), p(  2, 11)), 2);
  }

  void R_S()
  {
    std::cout << "Ray - Segment\n";
    check_squared_distance (R(p(2, 0), p( 0, 2)), S(p( 1, 1), p( 4, 0)), 0);
    check_squared_distance (R(p(10, 0), p( 0, 10)), S(p( 6, 6), p( 20, 0)), 2);
    check_squared_distance (R(p(-10, -13), p( 0, 10)), S(p( 10, 5), p( 70, -30)), 124.642);
    check_squared_distance (R(p(  0,   0), p( 6, -2)), S(p( -1, 4), p( 8, 6)), 16.9882);
    check_squared_distance (R(p( -8,  -7), p(11,  6)), S(p( 23, -27), p( -17, 16)), 0);
    check_squared_distance (R(p(  5,   0), p( 0,  0)), S(p(  1,   1), p(   2,  1)), 1);
    check_squared_distance (R(p(  0,   0), p( 5,  0)), S(p(  6,   1), p(   8,  1)), 1);
    check_squared_distance (R(p(  0,   0), p( 0, -3)), S(p(  1,   4), p(   1,  7)), 17);
    check_squared_distance (R(p(  0,   0), p( 5,  0)), S(p(  8,   1), p(   6,  1)), 1);
    check_squared_distance (R(p(  8,   0), p(12,  0)), S(p(  0,   0), p(   5,  0)), 9);
    check_squared_distance (R(p(  0,   0), p( 1,  0)), S(p(  2,   1), p(   2, -1)), 0);
    check_squared_distance (R(p(  4,   0), p(-3, -1)), S(p(  1,   1), p(   2, 11)), 2);
    check_squared_distance (R(p(  3,   4), p( 7,  7)), S(p(  7,   0), p(   6,  5)), 1);
    check_squared_distance (R(p( -1,   1), p( 3,  4)), S(p(  7,   0), p(   6,  5)), 1);
    check_squared_distance (R(p(  0,   0), p( 5,  0)), S(p(  1,   1), p(   6,  1)), 1);
    check_squared_distance (R(p(  0,   0), p( 5,  0)), S(p(  1,   1), p(   2,  1)), 1);
    check_squared_distance (R(p(  5,   0), p( 8,  0)), S(p(  1,   1), p(   2,  1)), 10);
    check_squared_distance (R(p(  5,   0), p( 0,  0)), S(p(  1,   1), p(   2,  1)), 1);
  }

  void R_L()
  {
    std::cout << "Ray - Line\n";
    check_squared_distance (R(p( 10,  0), p(  11,  10)), L(p(  0,  0), p(  0,  10)), 100);
    check_squared_distance (R(p( 10,  0), p(  11, -10)), L(p(  0,  0), p(  0,  10)), 100);
    check_squared_distance (R(p( 10,  0), p(  11,  10)), L(p(  0, -100), p(  0, -90)), 100);
    check_squared_distance (R(p( 10,  0), p(  11, -10)), L(p(  0, -100), p(  0, -90)), 100);
    check_squared_distance (R(p( 10,  0), p(  11, -10)), L(p(  0, 90), p(  0, 100)), 100);
    check_squared_distance (R(p( 10,  0), p(  11,  10)), L(p(  0,100), p(  0,  90)), 100);
    check_squared_distance (R(p( 10,  0), p(   9,  10)), L(p(  0,  0), p(  0,  10)), 0);
    check_squared_distance (R(p( 10,  0), p(   9, -10)), L(p(  0,  0), p(  0,  10)), 0);
    check_squared_distance (R(p(  1,  0), p(   1,   1)), L(p(  0,  0), p(  0,   1)), 1);
    check_squared_distance (R(p(  1,  0), p(   1,  -1)), L(p(  0,  0), p(  0,   1)), 1);
    check_squared_distance (R(p(  1,  0), p(   1,   1)), L(p(  0, -10), p(  0, -9)), 1);
    check_squared_distance (R(p(  1,  0), p(   1,  -1)), L(p(  0, -10), p(  0, -9)), 1);
    check_squared_distance (R(p(  1,  0), p(   1,   1)), L(p(  0,   9), p(  0, 10)), 1);
    check_squared_distance (R(p(  1,  0), p(   1,  -1)), L(p(  0,   9), p(  0, 10)), 1);
    check_squared_distance (R(p(  1,  0), p(   1,   1)), L(p(  0,  10), p(  0,  9)), 1);
  }

  void P_L()
  {
    std::cout << "Point - Line\n";
    check_squared_distance (p( 1,  1), L(p(  4,  0), p(  -3,  -1)), 2);
  }

  void L_S()
  {
    std::cout << "Line - Segment\n";
    check_squared_distance (L(p( 0,  0), p( 1,  0)), S(p(  2,  2), p(  3,  3)), 4);
    check_squared_distance (L(p( 0,  0), p( 1,  0)), S(p(  2,  2), p(  3,  1)), 1);
    check_squared_distance (L(p( 0,  0), p( 1,  0)), S(p(  2,  2), p(  3, -1)), 0);
    check_squared_distance (L(p( 0,  0), p( 1,  0)), S(p(  2,  2), p(  3,  2)), 4);
    check_squared_distance (L(p( 0,  0), p( 1,  0)), S(p(  2,  2), p( -3,  2)), 4);
  }

  void L_L()
  {
    std::cout << "Line - Line\n";
    check_squared_distance (L(p( 0,  0), p( 1,  0)), L(p(  2,  2), p(  3,  3)), 0);
    check_squared_distance (L(p( 0,  0), p( 1,  0)), L(p(  2,  2), p(  3,  1)), 0);
    check_squared_distance (L(p( 0,  0), p( 1,  0)), L(p(  2,  2), p(  3,  2)), 4);
    check_squared_distance (L(p( 0,  0), p( 1,  0)), L(p(  2,  2), p( -3,  2)), 4);
  }

  void P_T()
  {
    std::cout << "Point - Triangle\n";
    check_squared_distance (p( 10,  0), T(p(  5,  2), p( 18, 57), p(  2,56)), 29);
    check_squared_distance (p(  5,101), T(p(  0,  0), p(  5,  1), p( 10, 0)), 10000);
  }

  void L_T()
  {
    std::cout << "Line - Triangle\n";
    check_squared_distance (L(p( 0,  0), p( 1,  0)), T(p(  5,  2), p( 18, 57), p(  2,56)), 4);
  }

  void R_T()
  {
    std::cout << "Ray - Triangle\n";
    check_squared_distance (R(p( 1,  3), p( 0, 11)), T(p(  5,  2), p(  8, 57), p(  2,26)), 14.7846);
  }

  void S_T()
  {
    std::cout << "Segment - Triangle\n";
    check_squared_distance (S(p( 60,  0), p( 6,  0)), T(p(  5,  2), p( 18, 57), p(  2,56)), 5);
  }

  void T_T()
  {
    std::cout << "Triangle - Triangle\n";
    check_squared_distance (T(p(  7,  1), p(100, 15), p(103,  7)), T(p(  5,  2), p( 18, 57), p(  2,56)), 5);
  }

  void run()
  {
    std::cout << "2D Distance tests\n";
    P_P();
    P_S();
    S_S();
    P_R();
    R_R();
    R_S();
    R_L();
    P_L();
    L_S();
    L_L();
    P_T();
    L_T();
    R_T();
    S_T();
    T_T();
  }

};

int main()
{
        Test< CGAL::Simple_cartesian<double>   >().run();
        Test< CGAL::Simple_homogeneous<double> >().run();
        // TODO : test more kernels.
}
