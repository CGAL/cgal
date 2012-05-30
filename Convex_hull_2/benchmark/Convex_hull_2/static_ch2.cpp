#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/ch_jarvis.h>
#include <CGAL/ch_eddy.h>
#include <CGAL/Timer.h>

#include <CGAL/point_generators_2.h>
#include <CGAL/algorithm.h>
#include <boost/math/special_functions/next.hpp>

#include <iostream>
#include <vector>


#define bench(METHOD) \
{\
  unsigned run=0;\
  CGAL::Timer time;\
  do{\
    result.clear();\
    time.start();\
    METHOD( points.begin(), points.end(), std::back_inserter(result) );\
    time.stop();\
  }while(++run<repeat+1);\
  std::cout << result.size() << " points on the convex hull using "<< #METHOD << "; Done in "<< time.time() << "s\n";\
}

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point_2;
typedef std::vector<Point_2> Points;
typedef CGAL::Creator_uniform_2<double,Point_2>  Creator;

int main(int argc, char** argv)
{
  unsigned nbpts=100000;
  unsigned repeat=0;
  unsigned seed=0;
  
  if (argc>1)  nbpts=atoi(argv[1]);
  if (argc>2)  repeat=atoi(argv[2]);
  if (argc>3)  seed=atoi(argv[3]);
  
  Points points, result;
  
  CGAL::Random r(seed);
  CGAL::Random_points_in_disc_2<Point_2,Creator> g( 150.0,r);
  CGAL::cpp0x::copy_n( g, nbpts, std::back_inserter(points));
  
  std::cout << "seed is " << seed << "; using " << nbpts << " pts; on " << repeat+1 <<  " run(s).\n";   
  
  bench(CGAL::convex_hull_2)
  //bench(CGAL::ch_akl_toussaint)
  //bench(CGAL::ch_bykat)
  //bench(CGAL::ch_eddy)
  //bench(CGAL::ch_graham_andrew)
  //bench(CGAL::ch_jarvis)
}
