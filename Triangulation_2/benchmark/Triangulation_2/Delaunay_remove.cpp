#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Timer.h>
#include <CGAL/point_generators_2.h>

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel  K;
typedef CGAL::Delaunay_triangulation_2<K> Delaunay;
typedef Delaunay::Finite_vertices_iterator FVI;
typedef Delaunay::Vertex_handle Vertex_handle;
typedef K::Point_2                                     Point;
typedef CGAL::Creator_uniform_2<double,Point>  Creator;


int main(int argc, char **argv)
{
  int n=1000000;
  int rep=100;
  if (argc>=2)
    n=atoi(argv[1]);
  if (argc>=3)
    rep=atoi(argv[2]);
  std::vector<Point> points;
  points.reserve(n);
  CGAL::Random_points_in_disc_2<Point,Creator> g(1);
  CGAL::copy_n( g, n, std::back_inserter(points));
  Delaunay original;
  original.insert(points.begin(),points.end());

  double res=0;
  for (int r=0;r<rep;++r){
    Delaunay delaunay=original;
    std::vector<Vertex_handle> vertices;
    for(FVI fvi = delaunay.finite_vertices_begin(); fvi != delaunay.finite_vertices_end();++fvi){
      vertices.push_back(fvi);
    }
    CGAL::Timer t;
    t.start();
    for (int k=0; k<vertices.size(); ++k)
      delaunay.remove(vertices[k]);
    t.stop();
    res+=t.time();
    if (delaunay.number_of_vertices()!=0){
      std::cerr << "ERROR"<< std::endl;
      return 1;
    }
  }

  std::cout << res/rep << std::endl;

  return 0;
}
