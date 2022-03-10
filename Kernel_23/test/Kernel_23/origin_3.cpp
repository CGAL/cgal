//#define CGAL_PROFILE
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Timer.h>
#include <vector>
#include <iostream>
typedef CGAL::Exact_predicates_exact_constructions_kernel K;
typedef K::Point_3 Point_3;
typedef CGAL::Timer Timer;

int main(int argc, char* argv[] )
{
  std::ifstream ifs((argc>1)? argv[1]:CGAL::data_file_path("points/cube.xyz"));

  std::vector<Point_3> points;
  Point_3 p;

  while(ifs >> p){
    points.push_back(p);
  }

  const int N = points.size()-3;

  const K::Orientation_3 orientation = K().orientation_3_object();

  int positive = 0;
  Timer t;
  t.start();
  {

    for(int k = 0; k < 100; ++k)
    for(int i = 0; i < N; ++i){
       Point_3 o(CGAL::ORIGIN);
      if(orientation(o, points[i], points[i+1], points[i+2]) == CGAL::POSITIVE){
        ++positive;
      }
    }
  }
  t.stop();

  std::cout << t.time() << " sec." << std::endl;

  t.reset();
  t.start();
  {
      for (int k = 0; k < 100; ++k)
    for(int i = 0; i < N; ++i){
      if(orientation(CGAL::ORIGIN, points[i], points[i+1], points[i+2]) == CGAL::POSITIVE){
        --positive;
     }
    }
  }
  t.stop();

  assert(positive == 0);
  std::cout << t.time() << " sec." << std::endl;
  return 0;
}
