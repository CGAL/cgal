// file: examples/Triangulation_2/Voronoi.C

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>

#include <fstream>

struct K : CGAL::Exact_predicates_inexact_constructions_kernel {};

typedef CGAL::Delaunay_triangulation_2<K>  Triangulation;
typedef Triangulation::Edge_iterator  Edge_iterator;
typedef Triangulation::Point          Point;

int main( )
{
  std::ifstream in("data/voronoi.cin");
  std::istream_iterator<Point> begin(in);
  std::istream_iterator<Point> end;
  Triangulation T;
  T.insert(begin, end);

  int ns = 0;
  int nr = 0;  
  Edge_iterator eit =T.edges_begin();
  for ( ; eit !=T.edges_end(); ++eit) {
    CGAL::Object o = T.dual(eit);
    K::Segment_2 s;
    K::Ray_2     r;
    if (CGAL::assign(s,o)) {++ns;}
    if (CGAL::assign(r,o)) {++nr;}
  }
  std::cout << "The voronoi diagram as " << ns << "finite edges " 
	    << " and " << nr << " rays" << std::endl;
  return 0;
}
 
