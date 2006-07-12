// Benchmark for kernel constructions.
// Sylvain Pion, September 2000.

// Select the kernel you want to benchmark :

#include <CGAL/Cartesian.h>
// #include <CGAL/Simple_Cartesian_v2.h>
// #include <CGAL/Simple_cartesian.h>
// #include <CGAL/Homogeneous.h>

#include <iostream>
#include <CGAL/Timer.h>

// NB : with doubles, I've got reliability problems of the benchmarks,
//      probably because of alignment...
//      float is a little bit better it seems, but still...
// typedef int FT;
typedef double FT;

// The second line is the direct Point, not the wrapper provided by the _2
// classes.  They are indeed more efficient, as shown by the benchmark, and
// this is what is actually provided by the advanced kernel !
// So I guess it's worth making this kernel actually work !
// And what about making it the default, on compiler that support it ?
//
// More benchs :
// - triangulation.
// - orientation.
#ifdef CGAL_CARTESIAN_H
  typedef CGAL::Cartesian<FT>::Point_2 Point2;
  // typedef CGAL::PointC2<CGAL::Cartesian<FT> > Point2;
#elif defined CGAL_SIMPLE_CARTESIAN_H
  // typedef CGAL::Simple_cartesian<FT>::Point_2 Point2;
  typedef CGAL::PointS2<FT> Point2;
#elif defined CGAL_SIMPLE_CARTESIAN_V2_H
  typedef CGAL::Simple_Cartesian_v2<FT>::Point_2 Point2;
  // typedef CGAL::PointC2<CGAL::Simple_Cartesian_v2<FT> > Point2;
#elif defined CGAL_HOMOGENEOUS_H
  // typedef CGAL::Homogeneous<FT>::Point_2 Point2;
  typedef CGAL::PointH2<CGAL::Quotient<FT>, FT> Point2;
#endif

namespace CGAL {

#ifdef CGAL_CARTESIAN_H
inline
Point2
midpoint1(const Point2 &p, const Point2 &q)
{
  FT x,y;
  midpointC2(p.x(),p.y(),q.x(),q.y(),x,y);
  return Point2(x,y);
}

inline
Point2
midpoint2(const Point2 &p, const Point2 &q)
{
  // return Point2((p.x_ref() + q.x_ref())/FT(2),(p.y_ref() + q.y_ref())/FT(2));
  return Point2((p.x() + q.x())/FT(2), (p.y() + q.y())/FT(2));
}

#if 0
inline
Point2
midpoint3(const Point2 &p, const Point2 &q)
{
  Point2 M;
  midpointC2(p.x(),p.y(),q.x(),q.y(),M.x_ref_non_const(),M.y_ref_non_const());
  return M;
}
#endif // 0
#endif

} // namespace CGAL

Point2 A(5.7, 7.5);  // How is it possible when FT=int ?????
Point2 B(3.4, 4.5);  // I don't even get a warning !

int main()
{
  int i, loops=10000000;
  CGAL::Timer t;
  double dt;

#if 0

// Those timings were made on my laptop, which is now not mine anymore,
// so I need to make them again to be able to make useful comparisons...
// Maybe automazing the process would be useful to test on different
// platforms...

  //                   mp   / mp1  / mp2  / mp3

  // Cartesian       : 3.44 / 2.71 / 2.67 / 3.5
  // PointC2         : 2.27 / 2.26 / 2.17 / 2.78
  // Advanced kernel : 2.25 / 2.26 / 2.17 / 2.78
  // SimpleCartesian : 1.23 (1.21) (= without the wrapper classes)
  // Homogeneous     : 4.46 (3.47)

  dt = t.time(); t.start();
  for (i=0; i<loops; i++)
    Point2 C = CGAL::midpoint(A,B);

#else

  // Cartesian       : 4.13 / 3.68 / 3.63 / 4.65
  // PointC2         : 3.29 / 3.29 / 3.16 / 3.5
  // Advanced kernel : 3.29 / 3.29 / 3.16 / 3.51
  // SimpleCartesian : 1.32 (1.21)
  // Homogeneous     : 5.23 (4.22)

  Point2 C;
  dt = t.time(); t.start();
  for (i=0; i<loops; i++)
    C = CGAL::midpoint(A,B);

#endif

  t.stop();
  std::cout << "time = " << t.time()-dt << std::endl;

  return 0;
}

