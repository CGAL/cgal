#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <list>
#include <iostream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point_2;
typedef CGAL::Delaunay_triangulation_2<K> Delaunay_triangulation_2;


int main()
{
  std::vector<Point_2> input = { Point_2(0, 0), Point_2(1,1), Point_2(2,0), Point_2(2,2), Point_2(1,2), Point_2(0,2) };

  Delaunay_triangulation_2 dt(input.begin(), input.end());

  std::list<Point_2> result;
  Delaunay_triangulation_2::Vertex_circulator vc = dt.incident_vertices(dt.infinite_vertex()), done(vc);
  do{
    std ::cout << vc->point() << std::endl;
    // push_front in order to obtain the counterclockwise sequence
    result.push_front(vc->point());
    ++vc;
  }while(vc != done);

  return 0;
}
