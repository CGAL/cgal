

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Timer.h>
#include <iostream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel  K;
typedef CGAL::Delaunay_triangulation_3<K> DT;
typedef K::Point_3 Point_3;
typedef CGAL::Timer Timer;

int main()
{

  std::vector<Point_3> points;
  Point_3 p;

  while(std::cin >> p){
    points.push_back(p);
  } 
  DT dt;
  dt.insert(points.begin(), points.end());
  
  Timer timer;
  timer.start();
  int N = 0;
  for(int i = 0; i < 5; i++){
    for(DT::Vertex_iterator vit = dt.vertices_begin(); vit!= dt.vertices_end(); ++vit){
      std::vector<DT::Edge> E;
      E.reserve(64);
      dt.incident_edges(vit,std::back_inserter(E));
      N += E.size();
    }
  }
  timer.stop();
  
  std::cerr << N << std::endl << timer.time() << " sec" << std::endl;
  return 0;
}
