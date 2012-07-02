#include <fstream>

// CGAL headers
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Delaunay_hyperbolic_triangulation_2.h>

#include <CGAL/Triangulation_hyperbolic_traits_2.h>

#include <CGAL/Timer.h>


typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Triangulation_hyperbolic_traits_2<K> Gt;

typedef K::Point_2 Point_2;
typedef K::FT FT;

typedef CGAL::Delaunay_hyperbolic_triangulation_2<Gt> HDt;

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
  
  // Radius of the Poincare Disk
  FT r = 100;
  
  timer.start();
  
  HDt hdt = HDt(Gt(r));
  
  hdt.insert(pts.begin(), pts.end());
  
  timer.stop();
  
  std::cout << "H^2" << std::endl;
  std::cout << "Number of points: " << hdt.number_of_vertices() << std::endl;
  std::cout << "Time: " << timer.time() << std::endl;
  
  return 0;
}