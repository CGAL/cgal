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


#define bench(METHOD,CONTAINER) \
{\
  std::size_t previous=0;\
  unsigned run=0;\
  CGAL::Timer time;\
  do{\
    result.clear();\
    time.start();\
    METHOD( CONTAINER.begin(), CONTAINER.end(), std::back_inserter(result) );\
    time.stop();\
    if( previous!=0 && previous!=result.size()) std::cerr << "error got different result" << std::endl;\
    previous=result.size();\
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
  std::copy_n( g, nbpts, std::back_inserter(points));

  //the following code is for testing when there is only two extreme points, affine hull is 2D
  /*
  CGAL::Bbox_2 bbox=points.begin()->bbox();
  for (Points::iterator it=points.begin();it!=points.end();++it)
    bbox=bbox+it->bbox();

  points.push_back( Point_2(bbox.xmin()-1,bbox.ymin()-1) );
  points.push_back( Point_2(bbox.xmax()+1,bbox.ymax()+1) );
  */

  //the following code is for testing when there is only three extreme points
  /*
  CGAL::Bbox_2 bbox=points.begin()->bbox();
  for (Points::iterator it=points.begin();it!=points.end();++it)
    bbox=bbox+it->bbox();

  points.push_back( Point_2(bbox.xmin()-1,bbox.ymin()-1) );
  points.push_back( Point_2(bbox.xmax()+1,bbox.ymax()+1) );
  points.push_back( Point_2(bbox.xmax(),bbox.ymax()+2) );
  */

  //the following code is for testing when there is only two extreme points, affine hull is 1D
  /*
  points.clear();
  for (unsigned i=0;i<nbpts;++i)
    points.push_back(Point_2(i,i));
  */





  std::cout << "seed is " << seed << "; using " << nbpts << " pts; on " << repeat+1 <<  " run(s).\n";

  std::cout << "Using vector" << std::endl;
  bench(CGAL::convex_hull_2,points)
  //bench(CGAL::ch_akl_toussaint,points)
  //bench(CGAL::ch_bykat,points)
  //bench(CGAL::ch_eddy,points)
  //bench(CGAL::ch_graham_andrew,points)
  //bench(CGAL::ch_jarvis,points)

  {
    std::list<Point_2> pt_list;
    std::copy(points.begin(),points.end(),std::back_inserter(pt_list));
    std::cout << "Using list" << std::endl;
    bench(CGAL::convex_hull_2,pt_list)
  }
}
