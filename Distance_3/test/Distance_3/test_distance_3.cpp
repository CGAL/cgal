// 3D distance tests.

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
  typedef CGAL::Point_3< K >          P;
  typedef CGAL::Line_3< K >           L;
  typedef CGAL::Segment_3< K >        S;
  typedef CGAL::Ray_3< K >            R;
  typedef CGAL::Triangle_3< K >       T;
  typedef CGAL::Plane_3< K >          Pl;
  typedef CGAL::Iso_cuboid_3< K >     Cub;
  typedef CGAL::Tetrahedron_3< K >    Tet;


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


  P p(int x, int y, int z)
  {
    int w = ri.next();
    return P(to_nt(x*w), to_nt(y*w), to_nt(z*w), to_nt(w));
  }

  Pl pl(int a, int b, int c, int d)
  {
    int w = ri.next();
    return Pl(to_nt(a*w), to_nt(b*w), to_nt(c*w), to_nt(d*w));
  }


  void P_P()
  {
    std::cout << "Point - Point\n";
    check_squared_distance (p(0, 0, 0), p(0, 0, 0), 0);
    check_squared_distance (p(1, 1, 1), p(0, 0, 0), 3);
  }

  void P_S()
  {
    // Note : the values are not verified by hand
    std::cout << "Point - Segment\n";
    check_squared_distance (p(0, 1, 2), S(p(-3, 0, 0), p( 2, 0, 0)), 5);
    check_squared_distance (p(0, 1, 2), S(p( 3, 0, 0), p( 2, 0, 0)), 9);
    check_squared_distance (p(0, 1, 2), S(p( 2, 0, 0), p( 3, 0, 0)), 9);
    check_squared_distance (p(6, 1, 2), S(p( 2, 0, 0), p( 3, 0, 0)), 14);
  }

  void P_T()
  {
    std::cout << "Point - Triangle\n";
    check_squared_distance (p(0, 1, 2), T(p(0, 0, 0), p( 2, 0, 0), p( 0, 2, 0)), 4);

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
    std::cout << "Segment - Segment\n";
    check_squared_distance (S(p( -8, -7,  0), p( 11, 6, 0)), S(p(23, -27, 2), p( -17, 16, 2)), 4);
    check_squared_distance (S(p(  0,  0,  0), p(  5, 0, 0)), S(p( 1,   1, 2), p(   6,  1, 2)), 5);
    check_squared_distance (S(p(  0,  0,  0), p(  5, 0, 0)), S(p( 1,   1, 2), p(   2,  1, 2)), 5);
    check_squared_distance (S(p(  5,  0,  0), p(  8, 0, 0)), S(p( 1,   1, 2), p(   2,  1, 2)), 14);
    check_squared_distance (S(p(  5,  0,  0), p(  0, 0, 0)), S(p( 1,   1, 2), p(   2,  1, 2)), 5);
    check_squared_distance (S(p(  0,  0,  0), p(  5, 0, 0)), S(p( 6,   1, 2), p(   8,  1, 2)), 6);
    check_squared_distance (S(p(  0,  0,  0), p(  0,-3, 0)), S(p( 1,   4, 2), p(   1,  7, 2)), 21);
    check_squared_distance (S(p(  0,  0,  0), p(  5, 0, 0)), S(p( 8,   1, 2), p(   6,  1, 2)), 6);
    check_squared_distance (S(p(  0,  0,  0), p(  0, 0, 0)), S(p( 8,   1, 2), p(   6,  1, 2)), 41);
    check_squared_distance (S(p(  0,  0,  0), p(  1, 0, 0)), S(p( 2,   1, 2), p(   2, -1, 2)), 5);
    check_squared_distance (S(p(  2,  0,  0), p(  0, 2, 0)), S(p( 1,   1, 4), p(   4,  0, 4)), 16);
    check_squared_distance (S(p( 10,  0,  0), p(  0,10, 0)), S(p( 6,   6,20), p(  20,  0,20)), 402);
    check_squared_distance (S(p(-10,-13,  0), p(  0,10, 0)), S(p(10,   5,20), p(  70,-30,20)), 524.642);
    check_squared_distance (S(p(  0,  0,  0), p(30,-10, 0)), S(p(-5,  20,20), p(  40, 30,20)), 824.706);
    check_squared_distance (S(p(  4,  0,  0), p(-3, -1, 0)), S(p( 1,   1, 2), p(   2, 11, 2)), 6);
    check_squared_distance (S(p(  3,  4,  0), p( 7,  7, 0)), S(p( 7,   0, 2), p(   6,  5, 2)), 5);
    check_squared_distance (S(p( -1,  1,  0), p( 3,  4, 0)), S(p( 7,   0, 2), p(   6,  5, 2)), 13.8462);
  }

  void P_R()
  {
    // Note : the value is not verified by hand
    std::cout << "Point - Ray\n";
    check_squared_distance (p( -8, -7,  0), R(p(23, -27, 2), p( -17, 16, 2)), 86.3685);
  }

  void R_R()
  {
    // Note : the values are not verified by hand
    std::cout << "Ray - Ray\n";
    check_squared_distance (R(p( 0, 0, 30), p(  0, 30, 30)), R(p(100, -100, 0), p( 200,  1, 0)), 20899.5);
    check_squared_distance (R(p( 1, 0,  0), p(  0,  0,  0)), R(p(  1,    3, 3), p(   0,  0, 3)), 9);
    check_squared_distance (R(p( 0, 0,  0), p(  1,  0,  0)), R(p(  0,    0, 2), p(  -1,  0, 2)), 4);
  }

  void S_R()
  {
    // Note : the values are not verified by hand
    std::cout << "Segment - Ray\n";
    check_squared_distance (S(p( 0, 0, 30), p(  0, 30, 30)), R(p(100, -100, 0), p( 200,  1, 0)), 20899.5);
  }

  void R_L()
  {
    // Note : the values are not verified by hand
    std::cout << "Ray - Line\n";
    check_squared_distance (R(p(10, 0,  0), p( 20,  0,  0)), L(p(  0,    0, 3), p(   0,  3, 3)), 109);
    check_squared_distance (R(p( 0, 0, 30), p(  0, 30, 30)), L(p(100, -100, 0), p( 200,  1, 0)), 20899.5);
    check_squared_distance (R(p( 1, 0,  0), p(  0,  0,  0)), L(p(  1,    3, 3), p(   0,  0, 3)), 9);
    check_squared_distance (R(p( 0, 0,  0), p(  1,  0,  0)), L(p(  0,    0, 2), p(  -1,  0, 2)), 4);
  }

  void P_L()
  {
    std::cout << "Point - Line\n";
    check_squared_distance (p(  0,  1,  2), L(p(  2,    0, 0), p(   3,  0, 0)), 5);
    check_squared_distance (p(  0,  0,  2), L(p(  0,    0, 0), p(   1,  2, 0)), 4);
  }

  void S_L()
  {
    // Note : the values are not verified by hand
    std::cout << "Segment - Line\n";
    check_squared_distance (S(p( 1, 0,  0), p(  0,  0,  0)), L(p(  1,    3, 3), p(   0,  0, 3)), 9);
    check_squared_distance (S(p(-90, 0,  0), p(-10,  0,  0)), L(p(  0,    0, 3), p(   0,  3, 3)), 109);
    check_squared_distance (S(p(  0, 0,  0), p(  1,  0,  0)), L(p(  0,    0, 2), p(  -1,  0, 2)), 4);
  }

  void L_L()
  {
    // Note : the values are not verified by hand
    std::cout << "Line - Line\n";
    check_squared_distance (L(p(-10, 0,  0), p(-90,  0,  0)), L(p(  0,    0, 3), p(   0,  3, 3)), 9);
    check_squared_distance (L(p(  1, 0,  0), p(  0,  0,  0)), L(p(  1,    3, 3), p(   0,  0, 3)), 9);
    check_squared_distance (L(p(  0, 0,  0), p(  1,  0,  0)), L(p(  0,    0, 2), p(  -1,  0, 2)), 4);
  }

  void P_Pl()
  {
    std::cout << "Point - Plane\n";
    check_squared_distance (p(2, 5,  3), Pl(0, 1, 0, 0), 25);
  }

  void S_Pl()
  {
    std::cout << "Segment - Plane\n";
    check_squared_distance (S(p(2, -3,  3), p( 3,-7, 4)), pl(0, 1, 0, 0), 9);
  }

  void R_Pl()
  {
    std::cout << "Ray - Plane\n";
    check_squared_distance (R(p(2, -4,  3), p( 3,-4, 4)), Pl(0, 1, 0, 0), 16);
    check_squared_distance (R(p(2, -4,  3), p( 3, 4, 4)), Pl(0, 1, 0, 0), 0);
    check_squared_distance (R(p(2, -4,  3), p( 3,-8, 4)), Pl(0, 1, 0, 0), 16);
  }

  void L_Pl()
  {
    std::cout << "Line - Plane\n";
    check_squared_distance (L(p(2, -4,  3), p( 3,-4, 4)), Pl(0, 1, 0, 0), 16);
    check_squared_distance (L(p(2, -4,  3), p( 3, 4, 4)), Pl(0, 1, 0, 0), 0);
    check_squared_distance (L(p(2, -4,  3), p( 3,-8, 4)), Pl(0, 1, 0, 0), 0);
  }

  void Pl_Pl()
  {
    std::cout << "Plane - Plane\n";
    Pl p1(0, 1, 0, 0);
    typename K::Vector_3 v = -p1.orthogonal_vector();
    v /= CGAL::sqrt(v.squared_length());
    Pl p2 = Pl(0,-1,0,6);
    check_squared_distance (p1,p2, 36);
    check_squared_distance (Pl(-2, 1, 1, 0), Pl(2, 1, 3, 0), 0);
  }

  void run()
  {
    std::cout << "3D Distance tests\n";
    P_P();
    P_S();
    P_T();
    P_Tet();
    S_S();
    P_R();
    R_R();
    S_R();
    R_L();
    P_L();
    S_L();
    L_L();
    P_Pl();
    S_Pl();
    R_Pl();
    L_Pl();
    Pl_Pl();
  }

};

int main()
{
        Test< CGAL::Simple_cartesian<double>   >().run();
        Test< CGAL::Simple_homogeneous<double> >().run();
        // TODO : test more kernels.
}
