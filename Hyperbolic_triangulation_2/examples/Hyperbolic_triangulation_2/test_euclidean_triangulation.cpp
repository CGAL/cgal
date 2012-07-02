#include <fstream>

// CGAL headers
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Delaunay_triangulation_2.h>

#include <CGAL/Timer.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef K::Point_2 Point_2;

typedef CGAL::Delaunay_triangulation_2<K> Dt;

int main()
{
  CGAL::Timer timer;
  
  std::vector<Point_2> pts;
  
  std::ifstream f("points.cin");
  
  Point_2 p;
  while(f >> p) {
    pts.push_back(p);
  }
  f.close();
  
  timer.start();
  
  Dt dt = Dt();
  
  dt.insert(pts.begin(), pts.end());
  
  timer.stop();
  
  std::cout << "R^2" << std::endl;
  std::cout << "Number of points: " << dt.number_of_vertices() << std::endl;
  std::cout << "Time: " << timer.time() << std::endl;
  
  timer.reset();
    
  return 0;
}