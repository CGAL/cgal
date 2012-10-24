#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#define BENCH_CLASS_LOCAL
//#define BENCH_STACK
//#define BENCH_VECTOR
//#define BENCH_SMALL_VECTOR
#define BENCH_ARRAY

#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Timer.h>


#include <iostream>
#include <fstream>
#include <string>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel  K;
typedef CGAL::Delaunay_triangulation_2<K> Delaunay;
typedef K::Point_2                                     Point;



int main(int argc, char **argv)
{
  int n;
  std::cin >> n;
  std::vector<Point> points;
  points.reserve(n);
  std::copy(std::istream_iterator<Point>(std::cin), std::istream_iterator<Point>(), std::back_inserter(points));
  
  CGAL::Timer t;
  t.start();
  Delaunay delaunay;
  delaunay.insert(points.begin(), points.end());
  t.stop();
  std::cerr << t.time() << " seconds (";
#ifdef CGAL_TRIANGULATION_2_USE_OLD_PROPAGATING_FLIP 
  std::cerr << "Old implementation)\n";
#elif not defined(BENCH_CLASS_LOCAL)
  std::cerr << "Thread-local vector)\n";
#elif defined (BENCH_STACK)
  std::cerr << "Local stack)\n";
#elif defined (BENCH_VECTOR)
  std::cerr << "Local vector)\n";
#elif defined (BENCH_SMALL_VECTOR)
  std::cerr << "Local \"home-made\" small vector)\n";
#elif defined (BENCH_ARRAY)
  std::cerr << "Local fixed size array)\n";
#else
#  error
#endif
            
 return 0;
}
