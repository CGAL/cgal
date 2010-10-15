#include <fstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_2.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef CGAL::Triangulation_2<K>         Triangulation;
typedef Triangulation::Vertex_circulator Vertex_circulator;
typedef Triangulation::Point             Point;

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
