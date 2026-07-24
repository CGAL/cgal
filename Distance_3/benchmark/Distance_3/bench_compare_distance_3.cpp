#include <CGAL/Simple_cartesian.h>
#include <CGAL/Simple_homogeneous.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_integer.h>
#include <CGAL/Homogeneous.h>

#include <CGAL/squared_distance_3.h>

#include <CGAL/Cartesian_converter.h>

#include <CGAL/Random.h>
#include <CGAL/Timer.h>
#include <CGAL/Real_timer.h>

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
  typedef typename K::Comparison_result Comparison_result;
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
  double m = 0, M = 1;

public:
   Test(CGAL::Random& r) : r(r) { }

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
  void P_P(int N, FT d2)
  {
    std::cout << "Point - Point" << std::endl;
    CGAL::Real_timer t;
    t.start();
    for(int i=0; i<N; ++i)
    {
      P p0 = random_point();
      P p1 = random_point();
      K().compare_squared_distance_3_object()(p0, p1, d2);
    }
    t.stop();
    std::cout << t.time() << " sec" << std::endl;
  }

  void P_S(int N, FT d2)
  {
    std::cout << "Point - Segment" << std::endl;
    CGAL::Real_timer t;
    t.start();
    for(int i=0; i<N; ++i)
    {
      P p0 = random_point();
      P p1 = random_point();
      P p2 = random_point();
      K().compare_squared_distance_3_object()(p0, S(p1,p2), d2);
    }
    t.stop();
    std::cout << t.time() << " sec" << std::endl;
  }

  void P_T(int N, FT d2)
  {
    std::cout << "Point - Triangle" << std::endl;
    CGAL::Real_timer t;
    t.start();
    for(int i=0; i<N; ++i)
    {
      P p0 = random_point();
      P p1 = random_point();
      P p2 = random_point();
      P q = random_point();

      K().compare_squared_distance_3_object()(q, T(p0, p1, p2), d2);
    }
    t.stop();
    std::cout << t.time() << " sec" << std::endl;
  }

  void P_Tet(int N, FT d2)
  {
    std::cout << "Point - Tetrahedron\n";
    CGAL::Real_timer t;
    t.start();
    for(int i=0; i<N; ++i)
    {
      P p0 = random_point();
      P p1 = random_point();
      P p2 = random_point();
      P p3 = random_point();
      P q = random_point();

      K().compare_squared_distance_3_object()(q, Tet(p0, p1, p2, p3), d2);
    }
    t.stop();
    std::cout << t.time() << " sec" << std::endl;
  }

  void S_S(int N, FT d2)
  {
    std::cout << "Segment - Segment" << std::endl;
    CGAL::Real_timer t;
    t.start();
    for(int i=0; i<N; ++i)
    {
      P p0 = random_point();
      P p1 = random_point();
      P q0 = random_point();
      P q1 = random_point();

      K().compare_squared_distance_3_object()(S(p0, p1), S(q0, q1), d2);
    }
    t.stop();
    std::cout << t.time() << " sec" << std::endl;

  }

  void S_L(int N, FT d2)
  {
    std::cout << "Segment - Line" << std::endl;
    CGAL::Real_timer t;
    t.start();
    for(int i=0; i<N; ++i)
    {
      P p0 = random_point();
      P p1 = random_point();
      P q0 = random_point();
      P q1 = random_point();

      K().compare_squared_distance_3_object()(S(p0, p1), L(q0, q1), d2);
    }
    t.stop();
    std::cout << t.time() << " sec" << std::endl;

  }

  void T_T(int N, FT d2)
  {
    std::cout << "Triangle - Triangle with distance " << d2 << std::endl;
    CGAL::Real_timer t;
    t.start();
    for(int i=0; i<N; ++i)
    {
      P p0 = random_point();
      P p1 = random_point();
      P p2 = random_point();
      P p3 = random_point();
      P p4 = random_point();
      P p5 = random_point();

      // // these are still exact with EPECK
      // p0 = CGAL::midpoint(p0, p1);
      // p1 = p0 + FT(FT(0.1)) * V{p1 - p0};
      // p2 = p2 + V{p2 - p0} / FT(CGAL_PI);

      // // this is still exact with EPECK_with_sqrt
      // p4 = p4 + V{p4 - CGAL::ORIGIN} / CGAL::approximate_sqrt(CGAL::square(p4.x()) + CGAL::square(p4.y()) + CGAL::square(p4.z()) + 3);

      // p5 = p5 + V{p5 - CGAL::ORIGIN} * FT(std::cos(1.3));
      // generic triangles
      T tr1{p0, p1, p2}, tr2{p3, p4, p5};
      K().compare_squared_distance_3_object()(tr1, tr2, d2);
    }
    t.stop();
    std::cout << t.time() << " sec" << std::endl;
  }

public:
  void run()
  {
    std::cout << "Kernel: " << typeid(K).name() << std::endl;
    auto test_3=[]<typename F>(F& f){

    };
    P_P(10000000, FT(0.1));
    P_S(1000000, FT(0.1));

    std::cout << std::endl;
    P_T(500000, FT(10));
    P_T(500000, FT(0.1));
    P_T(500000, FT(0.001));
    std::cout << std::endl;

    P_Tet(200000, FT(0.1));

    std::cout << std::endl;
    S_S(500000, FT(10));
    S_S(500000, FT(0.1));
    S_S(500000, FT(0.001));
    std::cout << std::endl;
    S_L(500000, FT(0.1));

    std::cout << std::endl;
    T_T(500000, FT(10));
    T_T(500000, FT(0.1));
    T_T(500000, FT(0.001));
    std::cout << std::endl;
    std::cout << std::endl;
  }
};

int main(int argc, char** argv)
{
  std::cout.precision(17);
  std::cerr.precision(17);

  std::cout << "3D Distance tests" << std::endl;

  CGAL::Random rp;
  CGAL::Random r(argc==1?rp.get_seed():std::stoi(argv[1]));
  std::cout << "random seed = " << r.get_seed() << std::endl;

 Test<CGAL::Simple_cartesian<double> >(r).run();
 Test<CGAL::Simple_homogeneous<double> >(r).run();
//  Test<CGAL::Simple_cartesian<CGAL::Interval_nt<true> > >(r).run();

  // Test<CGAL::Homogeneous<CGAL::Exact_integer> >(r).run();

  Test<CGAL::Exact_predicates_inexact_constructions_kernel>(r).run();

  Test<CGAL::Exact_predicates_exact_constructions_kernel>(r).run();

  std::cout << "Done!" << std::endl;

  return EXIT_SUCCESS;
}
