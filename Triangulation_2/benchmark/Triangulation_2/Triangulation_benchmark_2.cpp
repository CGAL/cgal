#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

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
  std::cout << t.time() << " seconds\n";
            
 return 0;
}
