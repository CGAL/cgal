#include <fstream>

// CGAL headers
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Timer.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point_2;
typedef K::FT FT;

typedef CGAL::Delaunay_triangulation_2<K> Dt;

int main()
{
  CGAL::Timer timer;
  
  std::vector<Point_2> pts;
  
  std::ifstream f("points.cin");
  
  // read radius of the enclosing sphere
  FT r;
  f >> r;
  
  Point_2 p;
  while(f >> p) {
    pts.push_back(p);
  }
  f.close();
  
  Dt dt = Dt();
  std::vector<Point_2>::iterator ip = pts.begin();
  
  /*//
  std::vector<Point_2> pts2;
  std::ifstream f2("points2.cin");
  f2 >> r;
  
  while(f2 >> p) {
    pts2.push_back(p);
  }
  f2.close();
  
  dt.insert(pts2.begin(), pts2.end());
  //*/
  
  timer.start();
  
  while(ip != pts.end()) {
    dt.insert(*ip++);
  }
  
  timer.stop();
  
  std::cout << "R^2" << std::endl;
  std::cout << "Radius: " << r << std::endl;
  std::cout << "Number of points: " << dt.number_of_vertices() << std::endl;
  std::cout << "Time: " << timer.time() << std::endl;
  
  timer.reset();
    
  return 0;
}
