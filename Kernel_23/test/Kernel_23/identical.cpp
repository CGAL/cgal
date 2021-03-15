#define CGAL_IDENTICAL

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Timer.h>

#include <vector>
#include <algorithm>

typedef CGAL::Exact_predicates_exact_constructions_kernel K;
typedef K::Point_3 Point_3;


int main(){

  int N = 10000000;

  std::vector<Point_3> points(N);

  for(int i = 0; i < N; i+=10){
    points[i] = Point_3((N-i) * CGAL_PI , 0, 0);
    for(int j = 1; j < 10; j++){
      points[i+j] =  points[i];
    }
  }
  CGAL::Timer t;
  t.start();
  std::sort(points.begin(), points.end());
  std::cout << t.time() << " sec." << std::endl;
  return 0;
}
