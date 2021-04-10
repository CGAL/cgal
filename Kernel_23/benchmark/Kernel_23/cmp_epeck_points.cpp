#define TEST_IDENTICAL

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Timer.h>

#include <vector>
#include <algorithm>

typedef CGAL::Exact_predicates_exact_constructions_kernel K;
typedef K::Point_3 Point_3;
typedef K::Vector_3 Vector_3;


int main(){

  int N = 10000000;

  std::vector<Point_3> points(N);

  int count = 0;
  for(int i = 0; i < N; i+=10){
    points[i] = Point_3((N-i) * CGAL_PI , 0, 0) + Vector_3(i * (CGAL_PI/2.0), 0, 0);
    if(points[i].approx().x().inf() != points[i].approx().x().sup()){
      ++count;
    }
    for(int j = 1; j < 10; j++){
      points[i+j] =  points[i];
    }
  }
  std::cout << count << " points approx().x().inf() != approx().x().sup()" << std::endl;
  CGAL::Timer t;
  t.start();
  std::sort(points.begin(), points.end());
  std::cout << t.time() << " sec." << std::endl;
  return 0;
}
