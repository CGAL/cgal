// file : example/Triangulation_2/Voronoi.C
#include <CGAL/basic.h>
#include <fstream>

#include <CGAL/Cartesian.h>
#include <CGAL/Delaunay_triangulation_2.h>

typedef double coord_type;
typedef CGAL::Cartesian<coord_type>  Gt;
typedef CGAL::Delaunay_triangulation_2<Gt>  Triangulation;
typedef Triangulation::Edge_iterator  Edge_iterator;

int main( )
{
  std::ifstream in("data/voronoi.cin");
  std::istream_iterator<Gt::Point_2> begin(in);
  std::istream_iterator<Gt::Point_2> end;
  Triangulation T;
  T.insert(begin, end);

  int ns = 0;
  int nr = 0;  
  Edge_iterator eit =T.edges_begin();
  for ( ; eit !=T.edges_end(); ++eit) {
    CGAL::Object o = T.dual(eit);
    Gt::Segment_2 s;
    Gt::Ray_2     r;
    if (CGAL::assign(s,o)) {++ns;}
    if (CGAL::assign(r,o)) {++nr;}
  }
  std::cout << "The voronoi diagram as " << ns << "finite edges " 
	    << " and " << nr << " rays" << std::endl;
  return 0;
}
 
