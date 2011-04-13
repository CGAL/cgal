//  file : example/Triangulation_2/triangulation_prog1.C
#include <CGAL/basic.h>
#include <fstream>

#include <CGAL/Cartesian.h>
#include <CGAL/Triangulation_2.h>

typedef CGAL::Cartesian<double> Gt;
typedef CGAL::Triangulation_2<Gt> Triangulation;
typedef Triangulation::Vertex_circulator Vertex_circulator;
typedef Gt::Point_2   Point;

int main() {
  std::ifstream in("data/triangulation_prog1.cin");
  std::istream_iterator<Point> begin(in);
  std::istream_iterator<Point> end;

  Triangulation t;
  t.insert(begin, end);
    
  Vertex_circulator vc = t.incident_vertices(t.infinite_vertex()),
    done(vc);
  if (vc != 0) {
    do { std::cout << vc->point() << std::endl;
    }while(++vc != done);
  }
  return 0;
}
