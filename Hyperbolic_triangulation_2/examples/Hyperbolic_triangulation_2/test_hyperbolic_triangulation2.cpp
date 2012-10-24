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
typedef HDt::Vertex_handle Vertex_handle;
typedef HDt::Face_handle Face_handle;

int main()
{
  CGAL::Timer timer;
  
  std::vector<Point_2> pts;
  std::vector<Point_2>::iterator ip;
  
  std::ifstream f("points.cin");
  
  // read radius of the enclosing sphere
  FT r;
  f >> r;
  
  Point_2 p;
  while(f >> p) {
    pts.push_back(p);
  }
  f.close();
  
  /*
  //
  std::vector<Point_2> pts2;
  std::ifstream f2("points2.cin");
  f2 >> r;
  
  while(f2 >> p) {
    pts2.push_back(p);
  }
  f2.close();
  
  HDt hdt = HDt(Gt(r));
  hdt.insert(pts2.begin(), pts2.end());
  //*/
  
  HDt hdt = HDt(Gt(r));
  ip = pts.begin();
  
  timer.start();
  
  while(ip != pts.end()) {
    hdt.insert(*ip++);
  }
  
  timer.stop();
  
  std::cout << "H^2" << std::endl;
  std::cout << "Radius: " << r << std::endl;
  std::cout << "Number of points: " << hdt.number_of_vertices() << std::endl;
  std::cout << "Time: " << timer.time() << std::endl;
  
  return 0;
}
